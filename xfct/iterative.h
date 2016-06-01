#ifndef _ITERATIVE_H_
#define _ITERATIVE_H_

#include "raft_param.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


double
eval_ct_loglikelihood(raft_scan_t *data,
		      raft_iterative_t *workspace,
		      gsl_vector *f,
		      gsl_vector *ph);


void 
em_ct_iteration_knownAtt(raft_scan_t *data,
			 raft_iterative_t *workspace,
			 gsl_vector *p,
			 gsl_vector *n,
			 gsl_vector *b);


void 
em_pet_iteration_knownAtt(raft_scan_t *data,
			  raft_iterative_t *workspace,
			  gsl_vector *p,
			  gsl_vector *n,
			  gsl_vector *b,
			  gsl_vector *att);


void 
em_spect_iteration_knownAtt(raft_scan_t *data,
			    raft_iterative_t *workspace,
			    gsl_vector *p,
			    gsl_vector *n,
			    gsl_vector *b,
			    gsl_vector *att);


void
em_fluor_iteration_knownAtt(raft_scan_t *data,
			    raft_iterative_t *workspace,
			    gsl_vector *p,
			    gsl_vector *n,
			    gsl_vector *b,
			    gsl_vector *attT,
			    gsl_vector *attF);


void
em_ct_iteration_grad(raft_scan_t *data,
		     raft_iterative_t *workspace,
		     gsl_vector *p,
		     gsl_vector *n,
		     gsl_vector *b);


void
em_ct_iteration_convex(raft_scan_t *data,
		       raft_iterative_t *workspace,
		       gsl_vector *p,
		       gsl_vector *n,
		       gsl_vector *b);


void 
em_ct_iteration_standard(raft_scan_t *data,
			 raft_iterative_t *workspace,
			 gsl_vector *p,
			 gsl_vector *n,
			 gsl_vector *ph);



void 
art_ct_iteration(raft_scan_t *data,	   
		 raft_iterative_t *workspace,
		 gsl_vector *n,
		 gsl_vector *b,
		 int k);


void 
art_pet_iteration_knownAtt(raft_scan_t *data,
			   raft_iterative_t *workspace,
			   gsl_vector *n,
			   gsl_vector *b,
			   gsl_vector *att,
			   int k);


void
art_spect_iteration_knownAtt(raft_scan_t *data,	   
			     raft_iterative_t *workspace,
			     gsl_vector *n,
			     gsl_vector *b,
			     gsl_vector *att,
			     int k);


void
art_xfct_iteration_knownAtt(raft_scan_t *data,	   
			    raft_iterative_t *workspace,
			    gsl_vector *n,
			    gsl_vector *b,
			    gsl_vector *attT,
			    gsl_vector *attF,
			    int k);


void
kun_set_workspace(raft_scan_t *data,
		  raft_iterative_t *workspace,
		  gsl_vector **weight,
		  gsl_vector *p,
		  int acctype);


void
kun_set_workspace_precond(raft_scan_t *data,
			  raft_iterative_t *workspace,
			  gsl_vector **weight,
			  gsl_vector *precond,
			  gsl_vector *p,
			  int acctype);


void
kun_acceleration_standard(raft_scan_t *data,
			  raft_iterative_t *workspace,
			  gsl_vector **weight,
			  gsl_vector *prev,
			  gsl_vector *direction,
			  double *alpha);

void
kun_acceleration_gsjb(raft_scan_t *data,
		      raft_iterative_t *workspace,
		      gsl_vector **weight,
		      gsl_vector *prev,
		      gsl_vector *direction,
		      double par[2]);



#endif
