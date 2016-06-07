#ifndef RAFT_WRAPPER_EM_XFCT_H
#define RAFT_WRAPPER_EM_XFCT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "../raft/raft_image.h"

void oldraft_em_xfct(	raft_image rec,
			raft_image sino,
			raft_image trans,
			raft_image fluor,
			int niter);	

void oldraft_em360(raft_image rec,
		   raft_image sino, 
		   int niter);	

#ifdef __cplusplus
} //extern "C" {
#endif

#endif // #ifndef RAFT_WRAPPER_EM_XFCT_H
