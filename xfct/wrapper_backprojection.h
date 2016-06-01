#ifndef RAFT_WRAPPER_BACK_H
#define RAFT_WRAPPER_BACK_H

#ifdef __cplusplus
extern "C" {
#endif

#include "../raft/raft_image.h"

void 
oldraft_backprojection(raft_image back, raft_image sino);

void 
oldraft_backprojection_xfct(raft_image back,
			raft_image sino,
			raft_image trans,
			raft_image fluor);

void 
oldraft_fbp360(	raft_image fbp, 
		raft_image sino);

#ifdef __cplusplus
} //extern "C" {
#endif

#endif // #ifndef RAFT_WRAPPER_BACK_H
