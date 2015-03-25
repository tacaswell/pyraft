#ifndef RAFT_FBP_H
#define RAFT_FBP_H

#include "raft_image.h"
#include "raft_enums.h"

#ifdef __cplusplus
extern "C" {
#endif

void raft_fbp( double            cutoff,
               raft_image        sino,
               raft_image        image,
               raft_radon_method method,
               int               nthreads
             );

#ifdef __cplusplus
} // extern "C"
#endif

#endif // #ifndef RAFT_FBP_H
