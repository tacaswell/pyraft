#ifndef RAFT_IMAGE_BACKPROJECTION_ANDERSSON_H
#define RAFT_IMAGE_BACKPROJECTION_ANDERSSON_H

#include "raft_image.h"
#ifdef __cplusplus
extern "C" {
#endif

void raft_backprojection_andersson(raft_image, raft_image, int);

void raft_backprojection_andersson_sectoral(raft_image, raft_image, 
		double, double, int);

void raft_backprojection_andersson_sectoral2(raft_image sino, raft_image res_, 
		double t0, double t1, int nthreads);
void full_sectoral(raft_image sino, raft_image res_, int nthreads);
#ifdef __cplusplus
}
#endif

#endif //#ifndef RAFT_IMAGE_BACKPROJECTION_ANDERSSON_H

