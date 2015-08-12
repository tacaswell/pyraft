#include "raft_backprojection_andersson.h" 
#include <fftw3.h>
#include <algorithm>
#include <iostream>
#include <string.h>
#include "raft_image_functions.h"

////////////////////// MAIN FUNCTION //////////////////////////////
void raft_backprojection_andersson(raft_image sino, raft_image res_, int nthreads)
{
	raft_image res;
	int Nx = res_.data.lines;
	int Ny = res_.data.columns;
	// Just flipping the negative s to positive to have clear semi-polar coordinates
	res = sino2sp(sino);
	// Since after this transformation the size will differ, we renew it
	int Nt = res.data.columns;
	int Ns = Nt; 
	// Calculating the rho_0 on the base of result image step (assuming that its size is Ns X Nx)
	double r0 = 2*log(1.0/(Ns));
	// Transform to log-polar system
	res = sp2lp_mt(res, Ns, r0);
	// Calculating the kernel with sizes of sinogram in log-polar
	raft_image kernel = raft_kernel_lp_create(res, 0.005);
	// Convolve sino and kernel
	convolution_2d(res.data, kernel.data, res.data, nthreads);
	// Transform the coordinates from log-polar to Cartesian
	res = lp2c_mt(res, Nx, Ny);
	memcpy((double*)res_.data.p_data, (double*)res.data.p_data, sizeof(double)*Nx*Ny); 
}
