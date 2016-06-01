#include "raft_backprojection_logpolar.h" 
#include <fftw3.h>
#include <algorithm>
#include <iostream>
#include <string.h>
#include "raft_image_functions2.h"

double const pi = 3.1415926535897932384626433832795;

inline int N_rho(double r0, double ds)
{
	std::cout<<"r0="<<r0<<"; ds="<<ds<<std::endl;
	int res = (r0/log(1-ds));
	res = snapTransformSize_bst(res-1);
	std::cout<<"Calculated N_rho="<<res<<std::endl;
	return res;
}

raft_plan_logpolar raft_plan_logpolar_create(raft_image sino,
		raft_image res,
		double padding_coeff)
{
	int nrays, nviews;
	raft_plan_logpolar result;
	
	nrays  = sino.data.lines;
	nviews = sino.data.columns; 

	// Sinogram in polar coordinates
	result.polar_sino = raft_image_create( nrays/2, 2*nviews );
	if ( raft_image_is_empty( result.polar_sino ) )
	{
		result.polar_sino.tl_x = result.polar_sino.tl_y = 
			result.polar_sino.br_x = result.polar_sino.br_y = 0.0;
		return result;
	}

	result.polar_sino.tl_x = -pi;
	result.polar_sino.br_x = pi;
	result.polar_sino.tl_y = 1;
	result.polar_sino.br_y = 0;

	// Sinogram in logpolar coordinates
	double ds = result.polar_sino.tl_y/result.polar_sino.data.lines;
	double dx = 2.0/res.data.lines;
	double dy = 2.0/res.data.columns;
	double dec_st = std::min(dx,dy);
	double r0 = padding_coeff*log(dec_st);
	result.r0 = r0;

	int Nr = N_rho(r0, ds);

	result.logpolar_sino = raft_image_create( Nr, 2*nviews );
	if ( raft_image_is_empty( result.logpolar_sino ) )
	{
		result.logpolar_sino.tl_x = result.logpolar_sino.tl_y = 
			result.logpolar_sino.br_x = result.logpolar_sino.br_y = 0.0;
		return result;
	}

	result.logpolar_sino.tl_x = -pi;
	result.logpolar_sino.br_x = pi;
	result.logpolar_sino.tl_y = 0;
	result.logpolar_sino.br_y = r0;
	
	// Kernel (create it with its FT (?) here)
	result.kernel = raft_image_create(result.logpolar_sino.data.lines,
			result.logpolar_sino.data.columns);

	double acc = 0.005;
	double t0 = -pi;
	double t1 = pi;
	double dt = (t1-t0)/result.logpolar_sino.data.columns;

	double r1 = 0;
	double dr = (r1-r0)/Nr;

	double t=t0;
	r0 = -std::max(result.logpolar_sino.tl_y, result.logpolar_sino.br_y);
	double r  = r0;

	for(int j=0; j<result.kernel.data.columns; j++) {
		for(int i=0; i<Nr; i++) {
			double arg = fabs(exp(r)*cos(t) - 1);
			if(arg< acc) raft_matrix_element(result.kernel.data, i, j) = 1/acc;
			r+=dr;
		}
		r=r0;
		t+=dt;
	}

	result.kernel.tl_x = -pi;
	result.kernel.br_x = pi;
	result.kernel.tl_y = 0;
	result.kernel.br_y = result.r0;
	// FIXME Could we create fft_plan here and calculate the kernel FT here? 
	// Or just use the function fftw_C2R like the following code does?

// 	result.fft_kernel_re = raft_image_create(result.logpolar_sino.data.lines,
// 			logpolar_sino.data.columns);


	// Result in Logpolar coordinates
	result.logpolar_res = raft_image_create(result.logpolar_sino.data.lines,
			result.logpolar_sino.data.columns);
	
	result.logpolar_res.tl_x = -pi;
	result.logpolar_res.br_x = pi;
	result.logpolar_res.tl_y = 0;
	result.logpolar_res.br_y = result.r0;
	return result;
}

void raft_plan_logpolar_destroy(raft_plan_logpolar *plan)
{
	/* destroy sinogram @ polar coordinates */

	raft_image_destroy( &( plan->polar_sino ) );
	raft_image_destroy( &( plan->logpolar_sino ) );
	raft_image_destroy( &( plan->kernel ) );
	raft_image_destroy( &( plan->logpolar_res ) );
	
	/* */
}

void raft_kernel_lp_create(raft_image kernel, double r0)
{
	int Nr = kernel.data.lines;
	int Nt = kernel.data.columns;
	std::cout<<"kernel size: <"<<Nr<<","<<Nt<<">"<<std::endl;
	std::cout<<"r0="<<r0<<std::endl;
// 	std::cout<<
	double acc = 0.05;
	double t0 = -pi;
	double t1 = pi;
	double dt = (t1-t0)/Nt;

	double r1 = -r0;
	double dr = -r0/(double)Nr;

	double t=t0;
	double r  = 0;

	for(int j=0; j<Nt; j++) {
		for(int i=0; i<Nr; i++) {
			double arg = fabs(exp(r)*cos(t) - 1);
			if(arg< acc) kernel.data.p_data[j*Nr+i] = 1.0/acc;
				
// 				raft_matrix_element(kernel.data, i, j) = 1.0/acc;
			r+=dr;
		}
		r=0;
		t+=dt;
	}
	
	std::cout<<"kernel size: <"<<Nr<<","<<Nt<<">"<<std::endl;
// 	kernel.tl_x = -pi;
// 	kernel.br_x = pi;
// 	kernel.tl_y = 0;
// 	kernel.br_y = r0;
}

void raft_plan_logpolar_set_korners(raft_plan_logpolar *plan)
{
	plan->polar_sino.tl_x = -pi;
	plan->polar_sino.br_x = pi;
	plan->polar_sino.tl_y = 1;
	plan->polar_sino.br_y = 0;


	plan->logpolar_sino.tl_x = -pi;
	plan->logpolar_sino.br_x = pi;
	plan->logpolar_sino.tl_y = 0;
	plan->logpolar_sino.br_y = plan->r0;

	plan->kernel.tl_x = -pi;
	plan->kernel.br_x = pi;
	plan->kernel.tl_y = 0;
	plan->kernel.br_y = plan->r0;

	plan->logpolar_res.tl_x = -pi;
	plan->logpolar_res.br_x = pi;
	plan->logpolar_res.tl_y = 0;
	plan->logpolar_res.br_y = plan->r0;

}


void raft_backprojection_logpolar(raft_image sino, 
		raft_image res, 
		raft_plan_logpolar plan, 
		int nthreads)
{
	raft_plan_logpolar_set_korners( &plan );
	sino2sp_bst(sino, plan.polar_sino);
	sp2lp(plan.polar_sino, plan.logpolar_sino, plan.r0, nthreads);
	convolution_2d(plan.logpolar_sino, plan.kernel, plan.logpolar_res, nthreads);
	lp2c(plan.logpolar_res, res, nthreads);
	
}

