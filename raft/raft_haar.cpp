#include "raft_haar.h"

#include <mutex>
#include <condition_variable>
#include <vector>
#include <thread>
#include <cmath>

class semaphore {
   std::mutex m_;
   std::condition_variable cv_;
   unsigned i_;

public:

   semaphore( unsigned i ) : i_( i ) {}
   void decrement ( void )
   {
      std::unique_lock< std::mutex > lock( m_ );
      --i_;

      cv_.notify_all();

      while ( i_ )
         cv_.wait( lock );
   }
};

// Haar transform, single level, single thread.
// To be called by each thread with start being the
// thread id and skip being the number of threads.
void raft_haar_sl_st( raft_matrix data, int start, int skip, semaphore & sem, raft_matrix tmp )
{
   // Columns are processed in pairs.
   skip *= 2; start *= 2;

   // Distance between components:
   int I_del = data.lines / 2;
   int J_del = data.columns / 2;
   // Counters:
   int i, I, j, J;

   // Compute coeficients:
   for ( I = 0, i = 0; i < data.lines; ++I, i += 2 )
      for( J = 0, j = start; j < data.columns; J += 2, j += skip )
      {
         raft_matrix_element( tmp, I,         J     ) = 0.5 *(   raft_matrix_element( data, i,     j     ) 
                                                               + raft_matrix_element( data, i + 1, j     )
                                                               + raft_matrix_element( data, i,     j + 1 )
                                                               + raft_matrix_element( data, i + 1, j + 1 )
                                                             );
         raft_matrix_element( tmp, I + I_del, J     ) = 0.5 *(   raft_matrix_element( data, i,     j     )
                                                               + raft_matrix_element( data, i + 1, j     )
                                                               - raft_matrix_element( data, i,     j + 1 )
                                                               - raft_matrix_element( data, i + 1, j + 1 )
                                                             );
         raft_matrix_element( tmp, I,         J + 1 ) = 0.5 *(   raft_matrix_element( data, i,     j     )
                                                               - raft_matrix_element( data, i + 1, j     )
                                                               + raft_matrix_element( data, i,     j + 1 )
                                                               - raft_matrix_element( data, i + 1, j + 1 )
                                                             );
         raft_matrix_element( tmp, I + I_del, J + 1 ) = 0.5 *(   raft_matrix_element( data, i,     j     )
                                                               - raft_matrix_element( data, i + 1, j     )
                                                               - raft_matrix_element( data, i,     j + 1 )
                                                               + raft_matrix_element( data, i + 1, j + 1 )
                                                             );
      }

   // Waits for computations involving original matrix end:
   sem.decrement();

   // Copy coeficients:
   skip /= 2; start /= 2;
   for ( i = 0; i < data.lines; ++i )
      for( J = 0, j = start; j < data.columns / 2; J += 2, j += skip )
      {
         raft_matrix_element( data, i, j         ) = raft_matrix_element( tmp, i, J );
         raft_matrix_element( data, i, j + J_del ) = raft_matrix_element( tmp, i, J + 1 );
      }
}

// Haar transform, single level, multi-thread.
void raft_haar_sl( raft_matrix data, int nthreads )
{
   semaphore sem( nthreads );

   std::vector< std::thread > threads;
   threads.reserve( nthreads );

   raft_matrix storage = raft_matrix_create( data.lines, data.columns );
   int cur_column = 0;

   for ( int cur_thread = 0; cur_thread < nthreads; ++cur_thread )
   {
      raft_matrix tmp;
      tmp.p_data = storage.p_data + cur_column * storage.line_stride;
      tmp.line_stride = storage.line_stride;
      tmp.column_stride = storage.column_stride;
      cur_column += 2 * ( ( data.columns / ( 2 * nthreads ) ) + ( cur_thread < data.columns % ( 2 * nthreads ) ) );
      threads.push_back( std::thread( raft_haar_sl_st, 
                                      data,
                                      cur_thread,
                                      nthreads,
                                      std::ref( sem ),
                                      tmp
                                    )
                       );
   }
   for ( int cur_thread = 0; cur_thread < nthreads; ++cur_thread )
      threads[ cur_thread ].join();

   raft_matrix_destroy( &storage );
}

// Haar transform, multi-level, multi-thread.
void raft_haar( raft_matrix data, int nlevels, int nthreads )
{
   for ( int i = 0; i < nlevels; ++i )
   {
      raft_haar_sl( data, nthreads );
      data.lines /= 2;
      data.columns /= 2;
   }
}

// Inverse Haar transform, single level, single thread.
// To be called by each thread with start being the
// thread id and skip being the number of threads.
void raft_ihaar_sl_st( raft_matrix data, int start, int skip, semaphore & sem, raft_matrix tmp )
{
   // Distance between components:
   int i_del = data.lines / 2;
   int j_del = data.columns / 2;
   // Counters:
   int i, I, j, J;

   // Compute coeficients:
   for ( I = 0, i = 0; i < data.lines / 2; I += 2, ++i )
      for( J = 0, j = start; j < data.columns / 2; J += 2, j += skip )
      {
         raft_matrix_element( tmp, I,     J     ) = 0.5 *(   raft_matrix_element( data, i,         j         ) 
                                                           + raft_matrix_element( data, i + i_del, j         )
                                                           + raft_matrix_element( data, i,         j + j_del )
                                                           + raft_matrix_element( data, i + i_del, j + j_del )
                                                         );
         raft_matrix_element( tmp, I + 1, J     ) = 0.5 *(   raft_matrix_element( data, i,         j         ) 
                                                           + raft_matrix_element( data, i + i_del, j         )
                                                           - raft_matrix_element( data, i,         j + j_del )
                                                           - raft_matrix_element( data, i + i_del, j + j_del )
                                                         );
         raft_matrix_element( tmp, I,     J + 1 ) = 0.5 *(   raft_matrix_element( data, i,         j         ) 
                                                           - raft_matrix_element( data, i + i_del, j         )
                                                           + raft_matrix_element( data, i,         j + j_del )
                                                           - raft_matrix_element( data, i + i_del, j + j_del )
                                                         );
         raft_matrix_element( tmp, I + 1, J + 1 ) = 0.5 *(   raft_matrix_element( data, i,         j         ) 
                                                           - raft_matrix_element( data, i + i_del, j         )
                                                           - raft_matrix_element( data, i,         j + j_del )
                                                           + raft_matrix_element( data, i + i_del, j + j_del )
                                                         );
      }

   // Waits for computations involving original matrix end:
   sem.decrement();

   // Copy coeficients:
   skip *= 2; start *= 2;
   for ( i = 0; i < data.lines; ++i )
      for( J = 0, j = start; j < data.columns; J += 2, j += skip )
      {
         raft_matrix_element( data, i, j     ) = raft_matrix_element( tmp, i, J     );
         raft_matrix_element( data, i, j + 1 ) = raft_matrix_element( tmp, i, J + 1 );
      }
}

// Inverse Haar transform, single level, multi-thread.
void raft_ihaar_sl( raft_matrix data, int nthreads )
{
   semaphore sem( nthreads );

   std::vector< std::thread > threads;
   threads.reserve( nthreads );

   raft_matrix storage = raft_matrix_create( data.lines, data.columns );
   int cur_column = 0;

   for ( int cur_thread = 0; cur_thread < nthreads; ++cur_thread )
   {
      raft_matrix tmp;
      tmp.p_data = storage.p_data + cur_column * storage.line_stride;
      tmp.line_stride = storage.line_stride;
      tmp.column_stride = storage.column_stride;
      cur_column += 2 * ( ( data.columns / ( 2 * nthreads ) ) + ( cur_thread < data.columns % ( 2 * nthreads ) ) );
      threads.push_back( std::thread( raft_ihaar_sl_st,
                                      data,
                                      cur_thread,
                                      nthreads,
                                      std::ref( sem ),
                                      tmp
                                    )
                       );
   }
   for ( int cur_thread = 0; cur_thread < nthreads; ++cur_thread )
      threads[ cur_thread ].join();

   raft_matrix_destroy( &storage );
}

// Inverse Haar transform, multi-level, multi-thread.
void raft_ihaar( raft_matrix data, int nlevels, int nthreads )
{
   data.lines /= std::pow( 2.0, nlevels - 1 );
   data.columns /= std::pow( 2.0, nlevels - 1 );
   for ( int i = 0; i < nlevels; ++i )
   {
      raft_ihaar_sl( data, nthreads );
      data.lines *= 2;
      data.columns *= 2;
   }
}
