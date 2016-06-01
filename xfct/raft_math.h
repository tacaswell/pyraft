#ifndef _RAFT_MATH_H_
#define _RAFT_MATH_H_

#include "raft_param.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <math.h>


typedef struct{
  
  int N;
  double step;
  gsl_fft_complex_wavetable *gslwavetable;
  gsl_fft_complex_workspace *gslworkspace;
  gsl_vector *real, *imag;

}fft_t;


typedef struct{

  fft_t *fourier;
  
  gsl_vector *kernel;
  gsl_vector *fftREsignal, *fftIMsignal;
  gsl_vector *fftREkernel, *fftIMkernel;
  gsl_vector *fftREproduct, *fftIMproduct;
  gsl_vector *zero, *aux;

}hilb_t;


int
intRandomGen(int n);


void
maxVector(gsl_vector *x);


void 
eval_chebyshev2(int k, 
		double x, 
		double *u);


void 
eval_array_chebyshev2(int k, 
		      gsl_vector *x,
		      gsl_vector *U);


void 
svd_solve(gsl_matrix *A, 
	  gsl_vector *b,
	  int n,
	  gsl_matrix *V,
	  gsl_vector *S,
	  gsl_vector *x);


void
hilbert_workspace_alloc(hilb_t *workspace,
			int N, 
			int filter,
			double step);


void 
hilbert_workspace_free(hilb_t *workspace);


void
hilbert_transform(gsl_vector *g,
		  gsl_vector *h,
		  hilb_t *workspace);  





#endif

