#ifndef _RAFT_BACKP_H_
#define _RAFT_BACKP_H_

#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_trig.h>

#include "raft_scan.h"
#include "raft_errno.h"
#include "raft_projection.h"
#include "raft_math.h"
#include "raft_param.h"

/*######################################################
  Title: Backprojection
  ####################################################*/

/*+====================================================+
  
  TYPEDEF: raft_backp_t

  Purpose:

  RAFT backprojection structure data.   
  +====================================================+
*/

typedef struct{

  int nviews;  
  int filter;
  raft_proj_t projwork;

  struct{    
    fft_t *fftwork;
    gsl_vector *q;
    gsl_vector *signal, *impulse;
    gsl_vector *fftRE1, *fftRE2, *fftIM1, *fftIM2;
    gsl_vector *fftREimpulse, *fftIMimpulse, *zero;
    gsl_vector *fftREimpulseTretiak, *fftIMimpulseTretiak;
  }fourier;
  
  struct{    
    gsl_vector *vector;
    gsl_vector **matrix;    
  }integrand;

  struct{    
    gsl_vector *radon, *expradon;
    int defined;       
  }novikov;

  struct{
    gsl_vector *radon, *radonF, *expradon, *expradonF;
    gsl_vector *Phant, *RadonPhant, *kernel;
    gsl_vector *aux;
    int defined;    
  }invxfct;

  struct{
    gsl_vector *modproj;
    gsl_vector *hilbert, *coshilbert, *sinhilbert;
    gsl_vector *aux1, *aux2, *aux3, *aux4, *haux1, *haux2;    
    gsl_vector *data, *kernel, *att;
    hilb_t *hilbwork;  
  }inversion;

  struct{
    gsl_vector *next, *rad, *direction;
  }neumann;

  struct{
    fft_t *fft2work;
    gsl_vector *zero, *impulse, *filter;
    gsl_vector *fftREimpulse, *fftIMimpulse;
    gsl_vector *fftREback, *fftIMback;
  }fob;
  
}raft_backp_t;


void 
raft_backp(raft_scan_t *data, 
	   gsl_vector *p,
	   gsl_vector *b);


void 
raft_backp_partial(raft_scan_t *data, 
		   gsl_vector *p,
		   gsl_vector *b,
		   gsl_vector *T,
		   int sT);

int 
raft_backp_fob(raft_scan_t *data,
	       raft_backp_t *workspace,
	       gsl_vector *b,
	       gsl_vector *f);


int 
raft_backp_fbp(raft_scan_t *data, 	       
	       raft_backp_t *workspace,
	       gsl_vector *p,
	       gsl_vector *f);


void 
raft_backp_attenuated_ert(raft_scan_t *data, 
			  raft_backp_t *workspace,
			  double *att,
			  gsl_vector *p,
			  gsl_vector *b);


void
raft_backp_attenuated_pet(raft_scan_t *data, 
			  raft_backp_t *workspace,
			  gsl_vector *p,
			  gsl_vector *att,
			  gsl_vector *b);


void 
raft_backp_attenuated_spect(raft_scan_t *data, 
			    raft_backp_t *workspace,
			    gsl_vector *p,
			    gsl_vector *att,
			    gsl_vector *b);


void
raft_backp_attenuated_xfct(raft_scan_t *data, 
			   raft_backp_t *workspace,
			   gsl_vector *p,
			   gsl_vector *attT,
			   gsl_vector *attF,
			   gsl_vector *b);


void 
raft_backp_generalized(raft_scan_t *data, 
		       raft_backp_t *workspace,
		       gsl_vector *p,
		       gsl_vector **weight,				
		       gsl_vector *b);


int 
raft_backp_novikov(raft_scan_t *data, 
		   raft_backp_t *workspace,
		   gsl_vector *att,
		   gsl_vector *p,
		   gsl_vector *f);


int
raft_backp_tretiakMetz(raft_scan_t *data, 
		       raft_backp_t *workspace,
		       double *att,
		       gsl_vector *p,
		       gsl_vector *f);


int 
raft_backp_invxfct(raft_scan_t *data, 
		   raft_backp_t *workspace,
		   gsl_vector *attT,
		   gsl_vector *attF,
		   gsl_vector *p,
		   gsl_vector *f);


int 
raft_backp_invxfct_partial(raft_scan_t *data, 
			   raft_backp_t *workspace,
			   gsl_vector *attT,
			   gsl_vector *attF,
			   gsl_vector *p,
			   gsl_vector *f);


int
raft_backp_invxfct_neumann(raft_scan_t *data, 
			   raft_backp_t *workspace,
			   gsl_vector *attT,
			   gsl_vector *attF,
			   gsl_vector *p,
			   gsl_vector *f,
			   int n);


int
raft_backp_workspace_alloc(raft_scan_t *data,
			   int filter,
			   raft_backp_t *workspace);


void
raft_backp_workspace_free(raft_backp_t *workspace);


#endif
