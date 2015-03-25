#ifndef RAFT_RADON_BRESENHAM_H
#define RAFT_RADON_BRESENHAM_H

#include "raft_image.h"

#ifdef __cplusplus
extern "C" {
#endif

void raft_radon_bresenham( raft_image image,
                           raft_image sino,
                           int nthreads
                         );

#ifdef __cplusplus
}
#endif

#endif // #ifndef RAFT_IMAGE_BACKPROJECTION_BRESENHAM_H
