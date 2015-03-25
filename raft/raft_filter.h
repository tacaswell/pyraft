#ifndef RAFT_FILTER1D_H
#define RAFT_FILTER1D_H

#include "raft_image.h"
#include "raft_matrix.h"
#include "raft_vector.h"

#ifdef __cplusplus
extern "C" {
#endif

void raft_filter1d( raft_vector filter,
                    raft_matrix data,
                    int nthreads
                  );

void raft_filter1d_ramp( double cutoff,
                         raft_image sino,
                         int nthreads
                       );

#ifdef __cplusplus
} // extern "C"
#endif

#endif // #ifndef RAFT_FILTER1D_H
