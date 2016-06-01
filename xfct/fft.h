#ifndef _FFT_H_
#define _FFT_H_

#include "raft_math.h"

void
fft_workspace_alloc(fft_t *workspace, 
		    int N,
		    double step);


void 
fft_workspace_free(fft_t *workspace);


void 
ifft_complex(gsl_vector *rex, 
	     gsl_vector *imx,
	     gsl_vector *reX,
	     gsl_vector *imX,
	     fft_t *workspace);


void 
fft_complex(gsl_vector *reX,
	    gsl_vector *imX,
	    gsl_vector *rex,
	    gsl_vector *imx,
	    fft_t *workspace);


void
fft2_complex(gsl_vector *reQ,
	     gsl_vector *imQ,
	     gsl_vector *req,
	     gsl_vector *imq,
	     fft_t *workspace);


void
ifft2_complex(gsl_vector *req,
	      gsl_vector *imq,
	      gsl_vector *reQ,
	      gsl_vector *imQ,
	      fft_t *workspace);


void 
fft_shift(gsl_vector *x, int N);


#endif
