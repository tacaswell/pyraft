#ifndef RAFT_IMAGE_POLAR_H
#define RAFT_IMAGE_POLAR_H

#include "raft_image.h"

#ifdef __cplusplus
extern "C" {
#endif

  void 
  raft_image_c2p(raft_image cartesian,
		 raft_image polar,
		 int nthreads);
  
  void 
  raft_image_p2c(raft_image polar,
		 raft_image cartesian,
		 int nthreads);

  
  void 
  raft_image_c2lp(raft_image cartesian,
		  raft_image logpolar,
		  int nthreads);

  void 
  raft_image_lp2c(raft_image logpolar,
		  raft_image cartesian,
		  int nthreads);
  
  void 
  raft_image_s2p(raft_image sinogram,
		 raft_image polar,
		 int nthreads);

  void 
  raft_image_s2lp(raft_image sinogram,
		  raft_image logpolar,
		  int nthreads);

  
#ifdef __cplusplus
}
#endif

#endif // #ifndef RAFT_IMAGE_POLAR_H
