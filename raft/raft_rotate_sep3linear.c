#include "raft_rotate_sep3linear.h"

#include <math.h>
#include "raft_shear_linear.h"

void raft_rotate_sep3linear( raft_image image, double theta, int nthreads )
{
   static const double DPI = 8.0 * atan( 1.0 );
   static const double QPI = atan( 1.0 );
   int i;

   // Ignores 360 degrees rotations:
   int n = theta / DPI;
   theta -= n * DPI;
   n = ( theta < 0.0 ) ? -1.0 : 1.0;

   // Decompose into at most 45 degrees rotations:
   while ( fabs( theta ) > QPI )
   {
      raft_rotate_sep3linear( image, n * QPI, nthreads );
      theta -= n * QPI;
   }

   // Horizontal shear factor:
   double tan_half_theta = -tan( theta / 2.0 );

   // First horizontal shear:
   raft_shear_linear( image, tan_half_theta, nthreads );

   // Vertical shear:
   raft_shear_linear( raft_image_transpose( image ), sin( theta ), nthreads );

   // Second horizontal shear:
   raft_shear_linear( image, tan_half_theta, nthreads );
}
