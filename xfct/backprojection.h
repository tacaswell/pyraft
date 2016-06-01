#ifndef _FBACKPROJECTION_H_
#define _FBACKPROJECTION_H_

#include "raft_param.h"
#include "fft.h"
#include "interp.h"

typedef struct{

  raft_scan_t *data;
  gsl_vector *p;
  gsl_vector *b;
  int nthread;
  int rowIndex[2];
  int colIndex[2];  
  
}parBackLoop_t;

typedef struct{

  raft_scan_t *data;
  raft_backp_t *workspace;
  gsl_vector *p;
  gsl_vector *b;
  int nthread;
  int colIndex[2];  
  
}parBackLoop_xfct_t;

void 
*backp_pixelLoop_xfct(void *t);


void 
matrix_blockIndex(int n, int T, int M, 
		  int *fC, int *lC, int *fR, int *lR);


void 
*backp_pixelLoop_Mblock(void *t);

void 
*backp_pixelLoop_Lblock(void *t);


double 
eval_anglesum(raft_scan_t *data,
	      gsl_vector *projection,
	      int j,
	      int k,
	      double sum[4]);


double 
eval_anglesum_partial(raft_scan_t *data,
		      gsl_vector *projection,
		      int j,
		      int k,
		      gsl_vector *T,
		      int sT);


double
eval_anglesum_ert(raft_scan_t *data,
		  gsl_vector *projection,
		  double *att,
		  int j,
		  int k);


double
eval_anglesum_generalized(raft_scan_t *data,
			  gsl_vector **weight,
			  gsl_vector *projection,
			  int j,
			  int k);


double 
eval_anglesum_inversion(raft_scan_t *data,
			gsl_vector *mpr,
			gsl_vector **weight,
			int i,
			int k);


double 
eval_pixel_backprojection(raft_scan_t *data,
			  gsl_vector *projection,
			  int j);


int
filtered_projection(raft_scan_t *data, 
		    raft_backp_t *workspace,
		    gsl_vector *p,
		    gsl_vector *q);



void 
modified_projection_inversion(raft_scan_t *data, 
			      raft_backp_t *backp,
			      gsl_vector *p, 
			      gsl_vector *d,
			      gsl_vector *ep2);



void
define_novikov_workspace(raft_scan_t *data, 
			 raft_backp_t *backp,
			 gsl_vector *att, 
			 gsl_vector *p);



void
define_invxfctPartial_workspace(raft_scan_t *data, 
				raft_backp_t *backp,
				gsl_vector *attT,
				gsl_vector *attF,
				gsl_vector *p,
				double *apert);


void
define_invxfct_workspace(raft_scan_t *data, 
			 raft_backp_t *backp,
			 gsl_vector *attT,
			 gsl_vector *attF,
			 gsl_vector *p);



#endif



