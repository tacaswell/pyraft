#ifndef RAFT_WRAPPER_AKT_XFCT_H
#define RAFT_WRAPPER_AKT_XFCT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "../raft/raft_image.h"

void oldraft_akt_xfct(	raft_image rec,
			raft_image sino,
			raft_image trans,
			raft_image fluor,
			int niter,
			int wtype);	


#ifdef __cplusplus
} //extern "C" {
#endif

#endif // #ifndef RAFT_WRAPPER_AKT_XFCT_H
