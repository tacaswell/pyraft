#ifndef RAFT_RADON_SLANTSTACK_VIEW_H
#define RAFT_RADON_SLANTSTACK_VIEW_H

#include "raft_image.h"
#include "raft_vector.h"

#ifdef __cplusplus
extern "C" {
#endif

  void raft_radon_slantstack_view(raft_vector view,
				  raft_image phantom,
				  raft_vector t,
				  double theta);
  
  

#ifdef __cplusplus
} // extern "C" {
#endif

#endif // #ifndef RAFT_RADON_SLANTSTACK_VIEW_H

