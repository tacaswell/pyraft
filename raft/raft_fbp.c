#include "raft_fbp.h"
#include "raft_filter.h"
#include "raft_backprojection.h"

void raft_fbp( double            cutoff,
               raft_image        sino,
               raft_image        image,
               raft_radon_method method,
               int               nthreads
             )
{
   raft_filter1d_ramp( cutoff, sino, nthreads );

   switch ( method )
   {
      case RAFT_BRESENHAM :
         raft_backprojection_bresenham( sino, image, nthreads );
         break;

      case RAFT_SLANTSTACK :
         raft_backprojection_slantstack( sino, image, nthreads );
         break;
   }
}
