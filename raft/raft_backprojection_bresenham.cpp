#include "raft_backprojection_bresenham.h"

#include "raft_ddaiterator.h"
#include "raft_matrix.h"
#include <vector>
#include <thread>

extern "C" {

// BLAS daxpy:
void daxpy_( int const *, double const *, double const *, int const *, double *, int const * );
// BLAS dcopy:
void dcopy_( int const *, double const *, int const *, double *, int const * );

}

void raft_backprojection_bresenham_summer( std::vector< raft_image > const & storage,
                                           raft_image image,
                                           int start_column, int end_column,
                                           std::vector< raft_matrix > const & scale_storage,
                                           raft_matrix scale
                                         )
{
   // One:
   const double alpha = 1.0;

   // Counters:
   int i, j;
   int cur_thread;

   // For each column of partial backprojection:
   for ( j = start_column; j < end_column; ++j )
   {
      // For each working thread in the previous step (partial
      // backprojection):
      for ( cur_thread = 0; cur_thread < ( (int) storage.size() ); ++cur_thread )
      {
         // Sum the column to the full backprojection:
         daxpy_( &( image.data.lines ),
                 &alpha,
                 &( raft_image_sample( storage[ cur_thread ], 0, j ) ),
                 &( storage[ cur_thread ].data.column_stride ),
                 &( raft_image_sample( image, 0, j ) ),
                 &( image.data.column_stride )
               );
         // Sum the column to the full backprojection:
         daxpy_( &( scale.lines ),
                 &alpha,
                 &( raft_matrix_element( scale_storage[ cur_thread ], 0, j ) ),
                 &( scale_storage[ cur_thread ].column_stride ),
                 &( raft_matrix_element( scale, 0, j ) ),
                 &( scale.column_stride )
               );
      }

      // Scale the full backprojection:
      for ( i = 0; i < image.data.lines; ++i )
         raft_image_sample( image, i, j ) /=  raft_matrix_element( scale, i, j );
   }
}

void raft_backprojection_bresenham_worker( raft_image image,
                                           raft_image sino,
                                           int start_angle, int end_angle,
                                           raft_matrix scale
                                         )
{
   // Counters:
   int i, j;

   // Zero out output storage:
   const double dzero = 0.0;
   const int izero = 0;
   for ( j = 0; j < image.data.columns; ++j )
   {
      dcopy_( &( image.data.lines ),
              &dzero,
              &izero,
              &( raft_image_sample( image, 0, j ) ),
              &( image.data.column_stride )
            );
      dcopy_( &( scale.lines ),
              &dzero,
              &izero,
              &( raft_matrix_element( scale, 0, j ) ),
              &( scale.column_stride )
            );
   }

   // Iterator:
   raft_ddaiterator it;
   raft_ddaiterator_set_image( &it, image );

   // Sampling factors:
   double delta_theta = raft_image_hsamplingdistance( sino );
   double delta_t = raft_image_vsamplingdistance( sino );

   // Run through views:
   for ( j = start_angle; j < end_angle; ++j )
   {
      // Set view in iterator:
      raft_ddaiterator_set_theta( &it, sino.tl_x + ( j * delta_theta ) );
      // Run through rays in view:
      for ( i = 0; i < sino.data.lines; ++i )
      {
         // Set ray in iterator:
         raft_ddaiterator_set_t( &it, sino.tl_y + ( i * delta_t ) );
         // Trace the ray!
         while ( !raft_ddaiterator_end( &it ) )
         {
            raft_image_sample( image, it.i, it.j ) += raft_image_sample( sino, i, j ) * it.sigma;
            raft_matrix_element( scale, it.i, it.j ) += it.sigma;
            raft_ddaiterator_update( &it );
         }
      }
   }
}

void raft_backprojection_bresenham( raft_image sino,
                                    raft_image image,
                                    int nthreads
                                  )
{
   // Actual number of threads:
   int actual_nthreads;
   // Make sure we do not have too many or too little threads:
   actual_nthreads = ( nthreads <= sino.data.columns ) ? nthreads : sino.data.columns;
   actual_nthreads = ( nthreads > 0 ) ? nthreads : 1;

   // Storage for working threads:
   std::vector< raft_image > storage( actual_nthreads );
   int cur_thread;
   for ( cur_thread = 0; cur_thread < ( (int) storage.size() ); ++cur_thread )
   {
      storage[ cur_thread ] = raft_image_create_withviewport( image.data.lines,
                                                              image.data.columns,
                                                              image.tl_x, image.tl_y,
                                                              image.br_x, image.br_y
                                                            );
      // Allocation failure, leave without doing anything:
      // TODO: think of warning mechanism.
      if ( !( storage[ cur_thread ].data.p_data ) )
      {
         for ( cur_thread = cur_thread - 1; cur_thread >= 0; --cur_thread )
            raft_image_destroy( &( storage[ cur_thread ] ) );
         return;
      }
   }
   // Storage for working threads (scale factors):
   std::vector< raft_matrix > scale_storage( actual_nthreads );
   for ( cur_thread = 0; cur_thread < ( (int) scale_storage.size() ); ++cur_thread )
   {
      scale_storage[ cur_thread ] = raft_matrix_create( image.data.lines, image.data.columns );
      // Allocation failure, leave without doing anything:
      // TODO: think of warning mechanism.
      if ( !( scale_storage[ cur_thread ].p_data ) )
      {
         for ( cur_thread = cur_thread - 1; cur_thread >= 0; --cur_thread )
            raft_matrix_destroy( &( scale_storage[ cur_thread ] ) );
         for ( cur_thread = 0; cur_thread < ( (int) scale_storage.size() ); ++cur_thread )
            raft_image_destroy( &( storage[ cur_thread ] ) );
         return;
      }
   }
   // Accumulated scale factors:
   raft_matrix scale = raft_matrix_create( image.data.lines, image.data.columns );
   int i, j;
   for ( i = 0; i < image.data.lines; ++i )
      for ( j = 0; j < image.data.columns; ++j )
         raft_matrix_element( scale, i, j ) = 0.0;

   // Base number of angles per thread:
   int base_nangles = sino.data.columns / actual_nthreads;
   // Remainder, i.e., number of threads with an extra angle:
   int remainder_nangles = sino.data.columns % actual_nthreads;
   // Current starting view for working thread:
   int cur_starting_angle = 0;
   // Number of angles to be processed by current thread:
   int cur_nangles;

   // Create working threads:
   std::vector< std::thread > threads;
   threads.reserve( actual_nthreads );
   for ( cur_thread = 0; cur_thread < actual_nthreads; ++cur_thread )
   {
      cur_nangles = base_nangles + ( cur_thread < remainder_nangles );

      threads.push_back( std::thread( raft_backprojection_bresenham_worker,
                                      storage[ cur_thread ],
                                      sino,
                                      cur_starting_angle,
                                      cur_starting_angle + cur_nangles,
                                      scale_storage[ cur_thread ]
                                    )
                       );

      cur_starting_angle += cur_nangles;
   }

   // Wait for threads to finish the job:
   for ( auto& thread : threads )
      if ( thread.joinable() )
         thread.join();

   // Accumulate the partial backprojection results:

   // New number of working threads:
   actual_nthreads = ( nthreads <= image.data.columns ) ? nthreads : image.data.columns;

   // Base number of columns per thread:
   int base_ncolumns = image.data.columns / actual_nthreads;
   // Remainder, i.e., number of threads with an extra column:
   int remainder_ncolumns = image.data.columns % actual_nthreads;
   // Current starting column for working thread:
   int cur_starting_column = 0;
   // Number of columns to be processed by current thread:
   int cur_ncolumns;

   threads.resize( 0 );
   threads.reserve( actual_nthreads );
   for ( cur_thread = 0; cur_thread < actual_nthreads; ++cur_thread )
   {
      cur_ncolumns = base_ncolumns + ( cur_thread < remainder_ncolumns );

      threads.push_back( std::thread( raft_backprojection_bresenham_summer,
                                      storage,
                                      image,
                                      cur_starting_column,
                                      cur_starting_column + cur_ncolumns,
                                      scale_storage,
                                      scale
                                    )
                       );

      cur_starting_column += cur_ncolumns;
   }

   // Wait for threads to finish the job:
   for ( auto& thread : threads )
      if ( thread.joinable() )
         thread.join();

   // Release resources:
   for ( cur_thread = 0; cur_thread < ( (int) storage.size() ); ++cur_thread )
   {
      raft_image_destroy( &( storage[ cur_thread ] ) );
      raft_matrix_destroy( &( scale_storage[ cur_thread ] ) );
   }

   raft_matrix_destroy( &scale );
}
