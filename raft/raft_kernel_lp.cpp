#include "raft_kernel_lp.h"
#include <algorithm>
#include <iostream>

raft_image raft_backprojection_kernel_lp(raft_image source, double acc)
{
	int Nt = source.data.columns;
	int Nr = source.data.lines;
	raft_image res = source;
	res.data = raft_matrix_create(Nr, Nt);

	double t0 = std::min(source.tl_x, source.br_x);
	double t1 = std::max(source.tl_x, source.br_x);
	double dt = (t1-t0)/Nt;

	double r0 = std::min(source.tl_y, source.br_y);
	double r1 = std::max(source.tl_y, source.br_y);
	double dr = (r1-r0)/Nr;

	double t=t0;
	r0 = -std::max(source.tl_y, source.br_y);
	double r  = r0;

// 	std::cout<<"Creating kernel: t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
// 	std::cout<<"r0="<<r0<<"; r1="<<r1<<"; dr="<<dr<<std::endl;

	for(int j=0; j<Nt; j++) {
		for(int i=0; i<Nr; i++) {
			double arg = fabs(exp(r)*cos(t) - 1);
			if(arg< acc) raft_matrix_element(res.data, i, j) = 1/acc;
			r+=dr;
		}
		r=r0;
		t+=dt;
	}
	return res;
}

raft_image raft_projection_kernel_lp(raft_image source, double acc)
{
	int Nt = source.data.columns;
	int Nr = source.data.lines;
	raft_image res = source;
	res.data = raft_matrix_create(Nr, Nt);

	double t0 = std::min(source.tl_x, source.br_x);
	double t1 = std::max(source.tl_x, source.br_x);
	double dt = (t1-t0)/Nt;

	double r0 = std::min(source.tl_y, source.br_y);
	double r1 = std::max(source.tl_y, source.br_y);
	double dr = (r1-r0)/Nr;

	double t  = t0;
	double r  = r0;

	std::cout<<"Creating straight kernel: t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
	std::cout<<"r0="<<r0<<"; r1="<<r1<<"; dr="<<dr<<std::endl;
	std::cout<<"Nr="<<Nr<<"; Nt="<<Nt<<std::endl;
	for(int j=0; j<Nt; j++) {
		for(int i=0; i<Nr; i++) {
			double arg = fabs(cos(t) - exp(r));
			if(arg< acc) raft_matrix_element(res.data, i, j) = 1/acc;
			r+=dr;
		}
		r=r0;
		t+=dt;
	}
	return res;
}

