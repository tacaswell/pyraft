#ifndef RAFT_KERNEL_LP_H
#define RAFT_KERNEL_LP_H

#include "raft_image.h"
#ifdef __cplusplus
extern "C" {
#endif

raft_image raft_backprojection_kernel_lp(raft_image source, double acc);

raft_image raft_projection_kernel_lp(raft_image source, double acc);

#ifdef __cplusplus
}
#endif

#endif //#ifndef RAFT_IMAGE_BACKPROJECTION_ANDERSSON_H

