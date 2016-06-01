#ifndef RAFT_FUNCTIONS_MT_H
#define RAFT_FUNCTIONS_MT_H


#include "raft_image.h"
#include <fftw3.h>

double const pi = 3.1415926535897932384626433832795;
void c2sp_mt(	raft_image source, 
		raft_image &res, 
		int Nt, int Ns, int nthreads);

void sp2c_mt(	raft_image source, 
		raft_image &res, 
		int Nt, int Ns, int nthreads);

void c2lp_mt(raft_image source, 
		raft_image &res, 
		int Nt, int Nr, double r0, int nthreads);

void sp2c_miqueles(	raft_image source_r, 
			raft_image source_i,
			raft_image &res_r, 
			raft_image &res_i, 
			int Nx, int Ny, int nthreads);

void sp2c_miqueles2(	raft_image source_r, 
			raft_image source_i,
			raft_image &res_r, 
			raft_image &res_i, 
			int nthreads);

raft_image lp2c_mt(	raft_image source, 
			int Nx, 
			int Ny, 
			int nthreads	);

raft_image sp2lp_mt(	raft_image source, 
			double Nr, 
			double r0, 
			int nthreads	);

raft_image lp2sp_mt(	raft_image source, 
			double Ns, 
			int nthreads	);

raft_image bl_interpolate_mt(	raft_image source, 
				int Nx, 
				int Ny, 
				int nthreads	);


raft_vector interpolate_1d_mt(raft_vector f, int N, int nthreads);


fftw_complex *sp2c_special(fftw_complex *data, int Ns, int Nt, int nthreads);





#endif

