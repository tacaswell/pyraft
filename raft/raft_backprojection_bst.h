#ifndef RAFT_IMAGE_BACKPROJECTION_BST
#define RAFT_IMAGE_BACKPROJECTION_BST

#include "raft_image.h"
#include <fftw3.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef struct{
	double padding_coeff;
	raft_image polarsino; // Holds sinogram in polar coordinates;
	raft_image zero_padded_sino; // Holds zero-padded sinogram;
	raft_image fftp_re, fftp_im; // Holds Fourier image in polar coordinates;
	raft_image fftc_re, fftc_im; // Holds Fourier image in Cartesian coordinates;
	raft_image int1, cut;
}raft_bst;


raft_bst raft_bst_plan_create( raft_image sino, double padding_coeff );

void raft_bst_plan_destroy( raft_bst *plan );

// void raft_bst_plan_set_corners( rafr_bst * plan );

void raft_backprojection_bst (raft_image sino,	// Sinogram
			raft_image res, 	// Result
			raft_bst plan,
			int threads		// Number of threads
			);
  

#ifdef __cplusplus
}
#endif

#endif //#ifndef RAFT_IMAGE_BACKPROJECTION_BST

