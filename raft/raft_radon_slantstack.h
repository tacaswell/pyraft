#ifndef RAFT_RADON_SLANTSTACK_H
#define RAFT_RADON_SLANTSTACK_H

#include "raft_image.h"

#ifdef __cplusplus
extern "C" {
#endif

void raft_radon_slantstack(raft_image phantom,
			   raft_image radon,
			   int nthreads);

  

#ifdef __cplusplus
} // extern "C" {
#endif

#endif // #ifndef RAFT_RADON_SLANTSTACK_H

