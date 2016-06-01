#ifndef RAFT_IMAGE_FUNC_2 
#define RAFT_IMAGE_FUNC_2

#include "raft_image.h"
#include <fftw3.h>

#ifdef __cplusplus
extern "C" {
#endif

int iDivUp_bst(int a, int b);

int iAlignUp_bst(int a, int b);

int snapTransformSize_bst(int dataSize);


#ifdef __cplusplus
}
#endif


void sino2sp_bst(raft_image source, raft_image res);

void zero_padding_bst(raft_image source, raft_image res);

void sp2c_miqueles_bst(raft_image source_r, raft_image source_i,
		raft_image res_r, raft_image res_i, int nthreads);



raft_image sp2lp(raft_image source, raft_image res, double r0, int nthreads);

raft_image lp2c(raft_image source, raft_image res, int nthreads);

void convolution_2d(raft_image source, raft_image kernel, raft_image res, int nthreads); 










void ifftw_2d_bst(raft_image xr, raft_image xi, int threads);

void fft_shift_2d_bst(raft_image x);

void ifftw_2d_C2R_bst(raft_image source_r, raft_image source_i, raft_image res, int nthreads);

void bl_interpolate_mt_bst(raft_image source, raft_image res, int nthreads);

void symmetric_copy_bst(raft_image source, raft_image res);






double mean_dc(raft_matrix data);

void subtract_mean_dc(raft_image img, double meandc);

void add_mean_dc(raft_image img, double meandc);

double raft_matrix_maximum_value(raft_matrix data, int &i_max, int &j_max);

double raft_matrix_minimum_value(raft_matrix data, int &i_min, int &j_min);

#endif
