#ifndef RAFT_WRAPPER_RADON_H
#define RAFT_WRAPPER_RADON_H

#ifdef __cplusplus
extern "C" {
#endif

#include "../raft/raft_image.h"

void 
oldraft_radon(raft_image sino, raft_image phantom);

void 
oldraft_radon_xfct(raft_image sino, 
	        raft_image dens,
	        raft_image trans, 
		raft_image fluor);	


#ifdef __cplusplus
} //extern "C" {
#endif

#endif // #ifndef RAFT_WRAPPER_RADON_H
