#ifndef RAFT_IMAGE_FUNC 
#define RAFT_IMAGE_FUNC

#include "raft_image.h"

double const pi = 3.1415926535897932384626433832795;

double mean(raft_matrix data);

double mean_dc(raft_matrix data);

void to_pos(raft_matrix &x);

double raft_matrix_maximum_value(raft_matrix data, int &i_max, int &j_max); 

double raft_matrix_minimum_value(raft_matrix data, int &i_min, int &j_min);

void raft_image_normalize(raft_image &source); 

double desc_l2(raft_image img1, raft_image img2);

raft_image sino2sp(raft_image source);

raft_image get_sector(raft_image source, double t1, double t2);

raft_image filter_sector(raft_image source, double t1, double t2);

raft_image rotate(raft_image source, double t);

raft_image cut(raft_image source, double x0, double x1, 
		double y0, double y1);

raft_image raft_image_flip_x(raft_image source);

raft_image raft_image_flip_y(raft_image &source);

raft_image zero_padding(raft_image source, int Nl, int Nc);

raft_image zero_padding_on_s(raft_image source, int Ns);

raft_image zero_padding_on_s2(raft_image source, int Ns);

raft_image zero_padding_sino(raft_image source, int Ns);

raft_image bl_interpolate(raft_image source, int Nx, int Ny);

raft_image sp2c(raft_image source, int Nx, int Ny);

raft_image c2sp(raft_image source, int Nt, int Ns);

raft_image sp2lp(raft_image source, double Nr, double r0);

void lp2sp(raft_image source, raft_image &res, int Ns);

raft_image lp2c(raft_image source, int Nx, int Ny);

raft_image c2lp(raft_image source, int Nt, int Nr, double r0);

raft_image raft_kernel_lp_create(raft_image source, double acc);

// FFT functions
void fft_shift_2d(raft_matrix &x);

void convolution_2d(raft_matrix x, raft_matrix k, raft_matrix &res, int threads);

void convolution_2d_2(raft_matrix x, raft_matrix k, raft_matrix &res, int threads);

void deconvolution_2d(raft_matrix x, raft_matrix k, raft_matrix &res, double a, int threads);

void ifftw_2d(raft_matrix &xr, raft_matrix &xi, int threads);

void fftw_2d(raft_matrix &xr, raft_matrix &xi, int threads);

raft_image raft_kernel_lp_create(raft_image source, double acc);

raft_image raft_straight_kernel_lp_create(raft_image source, double acc);




///// Multithread functions
raft_image bl_interpolate_mt(raft_image source, int Nx, int Ny, int nthreads = 8);
raft_image lp2c_mt(raft_image source, int Nx, int Ny, int nthreads = 8);
raft_image sp2lp_mt(raft_image source, double Nr, double r0, int nthreads = 8);
void sp2c_miqueles(raft_image source, raft_image source_i,
		raft_image &res_r, raft_image &res_i, 
		int Nx, int Ny, int nthreads=8);










raft_image mo_sino(raft_image source, double dx, double dy, int Ns);


#endif
