#include "raft_radon_bresenham.h"

#include "raft_ddaiterator.h"
#include <vector>
#include <thread>

extern "C" {

// BLAS daxpy:
void daxpy_( int const *, double const *, double const *, int const *, double *, int const * );
// BLAS dcopy:
void dcopy_( int const *, double const *, int const *, double *, int const * );

}

void raft_radon_bresenham_worker( raft_image image,
                                  raft_image sino,
                                  int start_angle, int end_angle
                                )
{
   // Counters:
   int i, j;

   // Zero out output storage:
   const double dzero = 0.0;
   const int izero = 0;
   for ( j = start_angle; j < end_angle; ++j )
      dcopy_( &( sino.data.lines ),
              &dzero,
              &izero,
              &( raft_image_sample( sino, 0, j ) ),
              &( sino.data.column_stride )
            );

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
            raft_image_sample( sino, i, j ) += raft_image_sample( image, it.i, it.j ) * it.sigma;
            raft_ddaiterator_update( &it );
         }
      }
   }
}

void raft_radon_bresenham( raft_image image,
                           raft_image sino,
                           int nthreads
                          )
{
   // Actual number of threads:
   int actual_nthreads;
   // Make sure we do not have too many or too little threads:
   actual_nthreads = ( nthreads <= sino.data.columns ) ? nthreads : sino.data.columns;
   actual_nthreads = ( nthreads > 0 ) ? nthreads : 1;

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
   for ( int cur_thread = 0; cur_thread < actual_nthreads; ++cur_thread )
   {
      cur_nangles = base_nangles + ( cur_thread < remainder_nangles );

      threads.push_back( std::thread( raft_radon_bresenham_worker,
                                      image,
                                      sino,
                                      cur_starting_angle,
                                      cur_starting_angle + cur_nangles
                                    )
                       );

      cur_starting_angle += cur_nangles;
   }

   // Wait for threads to finish the job:
   for ( auto& thread : threads )
      if ( thread.joinable() )
         thread.join();
}
