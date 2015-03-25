#include "raft_translate_linear.h"

#include <thread>    // for std::thread;
#include <vector>    // for std::vector (will hold the threads).
#include <cmath>     // for std::floor

// Declare blas copy:
extern "C" {
   void dcopy_( int const * n, double const * px, int const * incx, double * py, int const * incy );
}

void raft_translate_linear( double delta, raft_vector x )
{
   // Temporary space for interpolation:
   raft_vector tmp = raft_vector_create( x.size );

   // Failed allocation test:
   if ( !tmp.size )
      return;

   // Translate each element:
   for ( int i = 0; i < x.size; ++i  )
   {
      // Interpolate:
      double index = i - delta;
      int idx = std::floor( index );
      double alpha = ( index - idx );
      double left = ( ( idx >= 0 ) && ( idx < x.size ) ) ?
                     raft_vector_element( x, idx ) :
                     0.0;
      ++idx;
      double right = ( ( idx >= 0 ) && ( idx < x.size ) ) ?
                     raft_vector_element( x, idx ) :
                     0.0;
      raft_vector_element( tmp, i ) = alpha * right + ( 1.0 - alpha ) * left;
   }

   // Copy translated signal:
   dcopy_( &(x.size), tmp.p_data, &(tmp.stride), x.p_data, &(x.stride) );

   // Frees Temporary space:
   raft_vector_destroy( &tmp );
}

// Translates many signals, single-threaded:
void raft_translate_many_linear_st( raft_vector    deltas,
                                    raft_matrix    mat,
                                    int            start,
                                    int            inc
                                    )
{
   for ( int i = start; i < mat.lines; i += inc )
      raft_translate_linear( raft_vector_element( deltas, i + start ),
                             raft_matrix_get_line( mat, i + start )
                           );
}

void raft_translate_many_linear( raft_vector deltas, 
                                 raft_matrix mat, 
                                 int nthreads
                               )
{
   // Number of translations to be performed:
   int size = ( mat.lines < deltas.size ) ? mat.lines : deltas.size;
   if ( !size )
      return;
   // Make sure we do not have too many threads:
   nthreads = ( nthreads <= size ) ? nthreads : size;

   // Create working threads:
   std::vector< std::thread > threads;
   threads.reserve( nthreads );
   for ( int cur_thread = 0; cur_thread < nthreads; ++cur_thread )
   {
      threads.push_back( std::thread( raft_translate_many_linear_st,
                                      deltas,
                                      mat,
                                      cur_thread,
                                      nthreads
                                    )
                       );
   }

   // Wait for threads to finish the job:
   for ( auto& thread : threads )
      if ( thread.joinable() )
         thread.join();
}
