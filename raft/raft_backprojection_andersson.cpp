#include "raft_backprojection_andersson.h" 
#include <fftw3.h>
#include <algorithm>
#include <iostream>
#include <string.h>
#include "raft_image_functions.h"

////////////////////// MAIN FUNCTION //////////////////////////////
void raft_backprojection_andersson_sectoral(raft_image sino, raft_image res_, 
		double t0, double t1, int nthreads)
{
	raft_image res;
	int Nx = res_.data.columns;
	int Ny = res_.data.lines;
	double beta = (t1 - t0)/2;
	double sinb =  sin( beta );
 	double ar = sinb/(1 + sinb);
	double k = 1 - ar;
	double r0 = log( 1-2*ar );
	double arg = (t1 + t0)/2;

	double dy = 0.0; 
	double dx = 1 - ar;
	std::cout<<"beta="<<beta<<"; sin(beta)="<<sinb<<std::endl;;
	std::cout<<"ar="<<ar<<"; r0="<<r0<<std::endl;;
	std::cout<<"dx="<<dx<<"; dy="<<dy<<std::endl;

	// Just flipping the negative s to positive to have clear semi-polar coordinates
	res = sino2sp(sino);
	res.tl_y = ar;
	res = rotate(res, -arg);
	// Since after this transformation the size will differ, we renew it
	int Nt = res.data.columns;
	int Ns = Nt; 
	res = mo_sino(res, dx, dy, Ns);
	//show_image(res, "", "", "", "mo_sino.png");
	// Transform to log-polar system
	res = sp2lp_mt(res, Ns, r0);
	// Calculating the kernel with sizes of sinogram in log-polar
	raft_image kernel = raft_kernel_lp_create(res, 0.001);
	res = filter_sector(res, -1.5*beta,1.5*beta);
	//show_image(res, "", "Sector", "", NULL);
	//show_image(kernel, "", "", "", NULL);
// 	// Convolve sino and kernel
	convolution_2d(res.data, kernel.data, res.data, nthreads);
// 	res = filter_sector(res, -beta, beta);
	res = rotate(res, arg + pi/2);
	//show_image(res, "", "Res_LP", "", NULL);
// // 	// Transform the coordinates from log-polar to Cartesian
	res = lp2c_mt(res, Nx, Ny);
	//show_image(res, "", "", "", "pbp.res");
	memcpy((double*)res_.data.p_data, (double*)res.data.p_data, sizeof(double)*Nx*Ny); 
}

void raft_backprojection_andersson_sectoral2(raft_image sino, raft_image res_, 
		double t0, double t1, int nthreads)
{
	raft_image res;
	double beta = (t1 - t0)/2;
	int Nx = res_.data.columns;
	int Ny = res_.data.lines;
	// Just flipping the negative s to positive to have clear semi-polar coordinates
	res = sino2sp(sino);
	// Since after this transformation the size will differ, we renew it
	int Nt = res.data.columns;
	int Ns = Nt; 
	// Transform to log-polar system
	double r0 = 2*log(1.0/std::max(Nx, Ny));
	res = sp2lp_mt(res, Ns, r0);
	double arg = (t1 + t0)/2;
	res = rotate(res, -arg);
	// Calculating the kernel with sizes of sinogram in log-polar
	raft_image kernel = raft_kernel_lp_create(res, 0.01);
	res = filter_sector(res, -beta-pi/3, beta+pi/3);
	kernel =  filter_sector(kernel, -beta-pi/3, beta+pi/3);
// 	show_image(res, "", "Sector", "", NULL);
// 	show_image(kernel, "", "", "", NULL);
// 	// Convolve sino and kernel
	convolution_2d(res.data, kernel.data, res.data, nthreads);
	res = rotate(res, arg);
// 	res = filter_sector(res, -beta, beta);
// 	show_image(res, "", "Res_LP", "", NULL);
// // 	// Transform the coordinates from log-polar to Cartesian
	res = lp2c_mt(res, Nx, Ny);
	raft_image_normalize(res);
// 	show_image(res, "", "", "", "pbp.res");
	memcpy((double*)res_.data.p_data, (double*)res.data.p_data, sizeof(double)*Nx*Ny); 
}

void full_sectoral(raft_image sino, raft_image res_, int nthreads)
{
	int nu = 2;
	double beta = pi/nu;
	raft_image res = res_, tmp_res = res_;
	int Nx = res_.data.columns;
	int Ny = res_.data.lines;
	double t0 = 0;
	for(int k = 0; k<nu; k++) {
		raft_backprojection_andersson_sectoral2(sino, tmp_res, t0, t0+beta, nthreads);
		
		for(int j=0; j<Ny; j++) {
			for(int i=0; i<Nx; i++) {
				double a = raft_matrix_element(res.data, i, j) + 
					raft_matrix_element(tmp_res.data, i, j);
				raft_matrix_element(res.data, i, j) = a;
					
			}
		}
		//QString huy;
		//huy.setNum(k);
		//show_image(tmp_res, "", huy, "", NULL);
		t0 += beta;
	}
	//show_image(res, "FULL", "FULL", "", NULL);
}

void raft_backprojection_andersson(raft_image sino, raft_image res_, int nthreads)
{
	//show_image(sino, "", "", "", "shepp_logan_sino.png", 1);
	double meandc = mean_dc(sino.data);
	raft_image res;
	int Nx = res_.data.lines;
	int Ny = res_.data.columns;
	// Just flipping the negative s to positive to have clear semi-polar coordinates
	res = sino2sp(sino);
// 	show_image(res, "", "", "", "shepp_logan_sino_sp.png", 1);
	// Since after this transformation the size will differ, we renew it
	int Nt = res.data.columns;
	int Ns = Nt; 
	// Calculating the rho_0 on the base of result image step (assuming that its size is Ns X Nx)
	double r0 = 2*log(1.0/std::max(Nx, Ny));
	// Transform to log-polar system
	res = sp2lp_mt(res, Ns, r0);
// 	show_image(res, "", "", "", "shepp_logan_sino_lp.png", 1);
	// Calculating the kernel with sizes of sinogram in log-polar
	raft_image kernel = raft_kernel_lp_create(res, 0.005);
	// Convolve sino and kernel
	convolution_2d(res.data, kernel.data, res.data, nthreads);
// 	raft_image_normalize(res);
// 	// Transform the coordinates from log-polar to Cartesian
// 	show_image(res, "", "", "", "bp_andersson_lp.png", 1);
	res = lp2c_mt(res, Nx, Ny);
	int r = 2;
	for(int i=0; i<Nx; i++) {
		for(int j=0; j<Ny; j++) {
			if((i-Nx/2)*(i-Nx/2) + (j-Ny/2)*(j-Ny/2)< r*r ) {
				raft_matrix_element(res.data, i,j) = 
					raft_matrix_element(res.data, Nx/2-r, Ny/2-r);
			}
		}
	}
	raft_image_normalize(res);
	memcpy((double*)res_.data.p_data, (double*)res.data.p_data, sizeof(double)*Nx*Ny); 
}
