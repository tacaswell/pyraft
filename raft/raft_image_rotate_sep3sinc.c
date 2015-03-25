#include "raft_image_rotate_sep3sinc.h"

#include <math.h>
#include "raft_image_shear_sinc.h"

void raft_image_rotate_sep3sinc( raft_image image, double theta, int nthreads )
{
   static const double DPI = 8.0 * atan( 1.0 );
   static const double HPI = 2.0 * atan( 1.0 );
   int i;

   // Ignores 360 degrees rotations:
   int n = theta / DPI;
   theta -= n * DPI;
   n = ( theta < 0.0 ) ? -1.0 : 1.0;

   // Decompose into at most 90 degrees rotations:
   while ( fabs( theta ) > HPI )
   {
      raft_image_rotate_sep3sinc( image, n * HPI, nthreads );
      theta -= n * HPI;
   }

   // Horizontal shear factor:
   double tan_half_theta = -tan( theta / 2.0 );

   // First horizontal shear:
   raft_image_shear_sinc( image, tan_half_theta, nthreads );

   // Vertical shear:
   raft_image_shear_sinc( raft_image_transpose( image ), sin( theta ), nthreads );

   // Second horizontal shear:
   raft_image_shear_sinc( image, tan_half_theta, nthreads );
}
