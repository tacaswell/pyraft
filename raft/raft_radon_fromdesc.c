#include "raft_radon_fromdesc.h"
#include <math.h>

void raft_radon_fromdesc( raft_image sino,
                          raft_matrix desc
                        )
{
   // Sampling distances:
   double theta_sampling_dist = raft_image_hsamplingdistance( sino );
   double t_sampling_dist = raft_image_vsamplingdistance( sino );

   // Auxiliary varables:
   double m,
          mab,
          delt,
          tp,
          c,
          s,
          cp,
          sp,
          cur_theta,
          ell_alpha,
          ell_a,
          ell_b,
          ell_theta,
          ell_x,
          ell_y;

   int cur_ellipse = 0, i, j;
   for ( ; cur_ellipse < desc.lines; ++cur_ellipse )
   {
      ell_alpha = raft_matrix_element( desc, cur_ellipse, 0 );
      ell_a = raft_matrix_element( desc, cur_ellipse, 1 );
      ell_b = raft_matrix_element( desc, cur_ellipse, 2 );
      ell_theta = raft_matrix_element( desc, cur_ellipse, 5 );
      ell_x = raft_matrix_element( desc, cur_ellipse, 3 );
      ell_y = raft_matrix_element( desc, cur_ellipse, 4 );

      for ( j = 0; j < sino.data.columns; ++j )
      {
         cur_theta = sino.tl_x + j * theta_sampling_dist;
         c = cos( cur_theta );
         s = sin( cur_theta );

         cur_theta -= ell_theta;
         cp = cos( cur_theta );
         sp = sin( cur_theta );

         m  = 1.0 / sqrt( ( ( cp * cp ) / ( ell_b * ell_b ) ) + ( ( sp * sp ) / ( ell_a * ell_a ) ) );
         mab = m / ( ell_a * ell_b );
         m *= 2.0; m *= ell_alpha;

         delt = mab * ( ( ell_x * c ) + ( ell_y * s ) - sino.tl_y );
         mab *= t_sampling_dist;

         for ( i = 0; i < sino.data.lines; ++i )
         {
            tp = mab * i - delt;
            tp *= tp;
            if ( tp <= 1.0 )
            {
               raft_image_sample( sino, i, j ) += m * sqrt( 1.0 - tp );
            }
         }
      }
   }
}
