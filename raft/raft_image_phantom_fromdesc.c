#include "raft_image_phantom_fromdesc.h"
#include <math.h>

void raft_image_phantom_fromdesc( raft_image image, raft_matrix desc )
{
   // Is there anything to do?
   if ( raft_image_is_empty( image ) )
      return;

   // Not enough data to describe phantom, empty image:
   if ( desc.columns < 6 )
   {
      raft_image_destroy( &image );
      return;
   }

   // Sampling distances:
   double x_sampling_dist = ( image.br_x - image.tl_x ) / ( image.data.columns - 1 );
   double y_sampling_dist = ( image.br_y - image.tl_y ) / ( image.data.lines - 1 );

   // Loop over ellipses:
   int cur_ellipse = 0;
   double x, y, x0, y0, theta, alpha, val, a, b, c, s, temp;
   for ( ; cur_ellipse < desc.lines; ++cur_ellipse )
   {
      alpha = raft_matrix_element( desc, cur_ellipse, 0 );
      a     = raft_matrix_element( desc, cur_ellipse, 1 );
      b     = raft_matrix_element( desc, cur_ellipse, 2 );
      x0    = raft_matrix_element( desc, cur_ellipse, 3 );
      y0    = raft_matrix_element( desc, cur_ellipse, 4 );
      theta = raft_matrix_element( desc, cur_ellipse, 5 );
      c     = cos( theta );
      s     = sin( theta );

      // Loop over samples:
      int i=0, j;
      for ( ; i < image.data.lines; ++i )
         for ( j = 0; j < image.data.columns; ++j )
         {
            x = image.tl_x + j * x_sampling_dist - x0;
            y = image.tl_y + i * y_sampling_dist - y0;

            val = ( y * s ) + ( x * c );
            val /= a;
            val *= val;

            temp = ( y * c ) - ( x * s );
            temp /= b;
            temp *= temp;

            val += temp;
            val = sqrt( val );

            // This sample is inside current ellipse:
            if ( val <= 1 )
               raft_matrix_element( image.data, i, j ) += alpha;
         }
   }
}
