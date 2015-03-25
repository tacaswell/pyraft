#ifndef RAFT_IMAGE_BACKPROJECTION_SLANTSTACK_H
#define RAFT_IMAGE_BACKPROJECTION_SLANTSTACK_H

#include "raft_image.h"

#ifdef __cplusplus
extern "C" {
#endif

void 
raft_backprojection_slantstack(raft_image sinogram,
			       raft_image backprojection,
			       int nthreads);

#ifdef __cplusplus
}
#endif

#endif // #ifndef RAFT_IMAGE_BACKPROJECTION_SLANTSTACK_H
