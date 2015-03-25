#ifndef RAFT_RADON_SLANTSTACK_WITHMESH_H
#define RAFT_RADON_SLANTSTACK_WITHMESH_H

#include "raft_image.h"
#include "raft_vector.h"

#ifdef __cplusplus
extern "C" {
#endif
  
  
  void raft_radon_slantstack_withmesh(raft_image radon,
				      raft_image phantom,
				      raft_vector t,
				      raft_vector theta,
				      int nthreads);
  
  

#ifdef __cplusplus
} // extern "C" {
#endif

#endif // #ifndef RAFT_RADON_SLANTSTACK_VIEW_H

