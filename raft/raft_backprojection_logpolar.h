#ifndef RAFT_IMAGE_BACKPROJECTION_ANDERSSON_H
#define RAFT_IMAGE_BACKPROJECTION_ANDERSSON_H

#include "raft_image.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct{
	double padding_coeff;
	double r0;
	raft_image polar_sino; // Holds sinogram in polar coordinates;
	raft_image logpolar_sino; // Holds sinogram in logpolar coordinates;
	raft_image kernel; // Holds kernel;
// 	raft_image fft_kernel_re, fft_kernel_im; // Holds FT(kernel)
	raft_image logpolar_res;
} raft_plan_logpolar;

raft_plan_logpolar raft_plan_logpolar_create(raft_image sino,
		raft_image res,
		double padding_coeff);

void raft_plan_logpolar_destroy(raft_plan_logpolar *plan);

void raft_kernel_lp_create(raft_image kernel, double r0);

void raft_backprojection_logpolar(raft_image sino,  // Sinogram
		raft_image res,			     // Image for storing the result
		raft_plan_logpolar plan,
		int nthreads = 8		     // Numbner of threads 		
		);

#ifdef __cplusplus
}
#endif

#endif //#ifndef RAFT_IMAGE_BACKPROJECTION_ANDERSSON_H

