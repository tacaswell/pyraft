#ifndef RAFT_IMAGE_FUNC 
#define RAFT_IMAGE_FUNC

#include "raft_image.h"
#include <fftw3.h>
#include <algorithm>


int iDivUp(int a, int b);

int iAlignUp(int a, int b);

int snapTransformSize(int dataSize);

double norm(double *data, int size);

void raft_image_set_corners(raft_image source, raft_image &res);

raft_image sino2sp(raft_image source);

raft_image sp2sino(raft_image source);

raft_image raft_image_rescale(raft_image source, double scale);

raft_image get_sector(raft_image source, double t1, double t2);

raft_image filter_sector(raft_image source, double t1, double t2);

raft_image filter_sector(raft_image source, double t0, double t1, double beta_nat);

raft_image pad_sector(raft_image sector, int Nt, double t0, double t1);

raft_image pad_sector(raft_image sector);

raft_image rotate(raft_image source, double t);

raft_image cut(raft_image source, double x0, double x1, 
		double y0, double y1);

raft_image raft_image_flip_x(raft_image source);

raft_image raft_image_flip_y(raft_image &source);

raft_image zero_padding(raft_image source, int Nl, int Nc);

raft_image zero_padding_on_s(raft_image source, int Ns);

void zero_padding_on_s(raft_image source, raft_image res_, int Ns);

raft_image zero_padding_on_s2(raft_image source, int Ns);

raft_image zero_padding_sino(raft_image source, int Ns);

// FFT functions
void fft_shift_2d(raft_matrix &x);

void fft_shift_2d(double *x, int Nx, int Ny);

double *convolution_2d_C2R(double *x, double *k, int Nx, int Ny, int threads );
/////////////
void convolution_2d(raft_matrix x, raft_matrix k, raft_matrix &res, int threads);

void convolution_2d_2(raft_matrix x, raft_matrix k, raft_matrix &res, int threads);

void deconvolution_2d(raft_matrix x, raft_matrix k, raft_matrix &res, double a, int threads);

void ifftw_2d(raft_matrix &xr, raft_matrix &xi, int threads);

void fftw_2d(raft_matrix &xr, raft_matrix &xi, int threads);

raft_matrix ifftw_2d_c2r(raft_matrix &xr, raft_matrix &xi, int threads); 

fftw_complex *fftw_2d_R2C(double *x, int Nx, int Ny, int threads);

double *ifftw_2d_C2R(fftw_complex *sp, int Nx, int Ny, int nthreads);
double *semi_convolution_2d_C2R(double *x, 
		fftw_complex *spectrum_kernel, int Nx, int Ny, int threads);


raft_image mo_sino(raft_image source, double dx, double dy, int Ns);


/////// NEW ////////////
raft_image mo_dec(raft_image src, double dx, double dy);

#endif
