#include "raft_backprojection_georgebressler.h"

#include <cmath>
#include <vector>
// #include "raft_image_rotate_sep3sinc.h"
#include "raft_rotate_sep3linear.h"

inline static void upsample ( raft_matrix & data, int max )
{
   if ( data.lines >= max )
      return;

   double cur_data_lines = data.lines;

   data.lines *= 2;
   data.lines = ( data.lines > max ) ? max : data.lines;

   double factor = data.lines / cur_data_lines;

   for ( int j = 0; j < data.columns; ++j )
      for ( int i = data.lines - 1; i >= 0; --i  )
         raft_matrix_element( data, i, j ) = 0.5 * ( raft_matrix_element( data, int( i / factor ), j ) +
                                                     raft_matrix_element( data, int( ( i + 1 ) / factor ), j )
                                                   );
}

// Sums image b over image a:
inline static void sum( raft_image a, raft_image b );
// Backprojects one single projection, overwriting previous content:
inline static void single_backprojection( raft_image image, // Where to backproject
                                          raft_image sino, // Data to backproject
                                          int i, // Which projection in sinogram
                                          double theta // Which angle to use
                                        );

void raft_backprojection_georgebressler( raft_image sino,
                                         raft_image image,
                                         int nthreads
                                       )
{
   // Number of hierarchy levels:
   int n_levels = std::ceil( std::log( sino.data.columns ) / std::log( 3.0 ) );
   // Initial sampling rate:
   int original_sampling = image.data.lines / std::pow( 2.0, n_levels - 2 );
   original_sampling = ( original_sampling > image.data.lines ) ? image.data.lines : original_sampling;
   // LIFO queue storage:
   std::vector< double > lifo_storage( ( 2 * n_levels + 1 ) * image.data.lines * image.data.columns );
   // LIFO container:
   std::vector< raft_image > lifo_data( 2 * n_levels + 1 );
   for ( unsigned i = 0; i < lifo_data.size(); ++i ) // Prepare images
   {
      lifo_data[ i ].data.p_data = &(lifo_storage[0]) + i * image.data.lines * image.data.columns;
      lifo_data[ i ].data.lines = original_sampling;
      lifo_data[ i ].data.columns = image.data.columns;
      lifo_data[ i ].data.column_stride = 1;
      lifo_data[ i ].data.line_stride = lifo_data[ i ].data.columns;
      lifo_data[ i ].tl_x = image.tl_x;
      lifo_data[ i ].tl_y = image.tl_y;
      lifo_data[ i ].br_x = image.br_x;
      lifo_data[ i ].br_y = image.br_y;
   }
   // Smallest index in the sum of LIFO item:
   std::vector< int > lifo_least_projection( 2 * n_levels + 1 );

   // Next hierarchy level:
   int next_level;
   // Current backprojection being processed:
   int curr_projection = 0;
   // LIFO queue front:
   int lifo_front = 0;

   int last = std::round( std::pow( 3, n_levels ) );

   double delta_theta = raft_image_hsamplingdistance( sino );
   while ( curr_projection < last )
   {
      if ( curr_projection < sino.data.columns )
         single_backprojection( lifo_data[ lifo_front ],
                                sino,
                                curr_projection,
                                0.0
                              );

      lifo_least_projection[ lifo_front ] = curr_projection;

      ++lifo_front;
      ++curr_projection;

      next_level = 1;
      int next_skip = 3;
      while ( curr_projection % next_skip == 0 )
      {
         if ( lifo_least_projection[ lifo_front - 3 ] < sino.data.columns )
         {
            upsample( lifo_data[ lifo_front - 3 ].data, image.data.lines );
         }

         if ( lifo_least_projection[ lifo_front - 2 ] < sino.data.columns )
         {
//             raft_image_rotate_sep3sinc( lifo_data[ lifo_front - 2 ],
//                                         std::pow( 3, next_level - 1 ) * delta_theta,
//                                         nthreads
//                                       );
            raft_rotate_sep3linear( lifo_data[ lifo_front - 2 ],
                                    std::pow( 3, next_level - 1 ) * delta_theta,
                                    nthreads
                                  );
            upsample( lifo_data[ lifo_front - 2 ].data, image.data.lines );
            sum( lifo_data[ lifo_front - 3 ], lifo_data[ lifo_front - 2 ] );
            lifo_data[ lifo_front - 2 ].data.lines = original_sampling;
         }
         if ( lifo_least_projection[ lifo_front - 1 ] < sino.data.columns )
         {
//             raft_image_rotate_sep3sinc( lifo_data[ lifo_front - 1 ],
//                                         2.0 * std::pow( 3, next_level - 1 ) * delta_theta,
//                                         nthreads
//                                       );
            raft_rotate_sep3linear( lifo_data[ lifo_front - 1 ],
                                    2.0 * std::pow( 3, next_level - 1 ) * delta_theta,
                                    nthreads
                                  );
            upsample( lifo_data[ lifo_front - 1 ].data, image.data.lines );
            sum( lifo_data[ lifo_front - 3 ], lifo_data[ lifo_front - 1 ] );
            lifo_data[ lifo_front - 1 ].data.lines = original_sampling;
         }

         lifo_front -= 2;
         ++next_level;
         next_skip *= 3;
      }
   }

   // Normalization and sum:
   double del = 0.5 / sino.data.columns;
   for ( int i = 0; i < image.data.lines; ++i )
      for ( int j = 0; j < image.data.columns; ++j )
      {
         raft_image_sample( image, i, j ) += raft_image_sample( lifo_data[ 0 ], i, j ) * del;
      }
}

static void single_backprojection( raft_image image,
                                   raft_image sino,
                                   int k,
                                   double theta
                                 )
{
   // Image sampling intervals:
   double delta_x = raft_image_hsamplingdistance( image );
   double delta_y = raft_image_vsamplingdistance( image );
   // Projection sampling interval:
   double delta_t = raft_image_vsamplingdistance( sino );

   // For all samples in image:
   for ( int i = 0; i < image.data.lines; ++i )
      for ( int j = 0; j < image.data.columns; ++j )
      {
         // Compute current sample coordinates:
         double x = image.tl_x + j * delta_x;
         double y = image.tl_y + i * delta_y;

         // Required position in projection:
         double t = x * cos( theta ) + y * sin( theta );

         // Index of sample in projection data
         // holding the desired value:
         double index = ( t - sino.tl_y ) / delta_t;

         // Interpolate:
         int idx = index;
         double alpha = ( index - idx );
         double left = ( ( idx >= 0 ) && ( idx < sino.data.lines ) ) ?
                       raft_image_sample( sino, idx, k ) :
                       0.0;
         double right = ( ( idx + 1 >= 0 ) && ( idx + 1 < sino.data.lines ) ) ?
                        raft_image_sample( sino, idx + 1, k ) :
                        0.0;
         raft_image_sample( image, i, j ) = alpha * right + ( 1.0 - alpha ) * left;
      }
}

static void sum( raft_image a, raft_image b )
{
   for ( int i = 0; i < a.data.lines; ++i )
      for ( int j = 0; j < a.data.columns; ++j )
         raft_image_sample( a, i, j ) += raft_image_sample( b, i, j );
}
