#ifndef _RAFT_WEIGHT_H_
#define _RAFT_WEIGHT_H_

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>

#include "raft_scan.h"
#include "raft_phantom.h"
#include "raft_projection.h"
#include "raft_errno.h"
#include "raft_cbf.h" 

void
raft_weight_pet(raft_scan_t *data,
		raft_proj_t *proj,
		gsl_vector *att);


void
raft_weight_spect(raft_scan_t *data,
		  raft_proj_t *proj,
		  gsl_vector *att);


void
raft_weight_xfct(raft_scan_t *data,
		 raft_proj_t *proj,
		 gsl_vector *attT,
		 gsl_vector *attF);


void 
raft_weight_xfct_partial(raft_scan_t *data,
			 raft_proj_t *proj,
			 gsl_vector *attT,
			 gsl_vector *attF);


void 
raft_weight_spect_view(raft_scan_t *data,
		       raft_proj_t *proj,
		       gsl_vector *att,
		       int j);


void 
raft_weight_xfct_view(raft_scan_t *data,
		      raft_proj_t *proj,
		      gsl_vector *attT,
		      gsl_vector *attF,
		      int j);


gsl_vector *
raft_weight_pet_get(raft_proj_t *proj);


void
raft_weight_spect_get(raft_proj_t *proj,
		      gsl_vector **weight);

void
raft_weight_xfct_get(raft_proj_t *proj,
		     gsl_vector **weight);



gsl_vector *
raft_weight_spect_get_view(raft_proj_t *proj,
			   int view);


gsl_vector *
raft_weight_xfct_get_view(raft_proj_t *proj,
			  int view);


void
raft_weight_pet_cbf(raft_scan_t *data,
		    raft_proj_t *proj,
		    raft_cbf_t *cbf,
		    gsl_vector *att);


void
raft_weight_spect_cbf(raft_scan_t *data,
		      raft_proj_t *proj,
		      raft_cbf_t *cbf,
		      gsl_vector *att);


void 
raft_weight_xfct_cbf(raft_scan_t *data,
		     raft_proj_t *proj,
		     raft_cbf_t *cbf,
		     gsl_vector *attT,
		     gsl_vector *attF,
		     int id);

void 
raft_weight_average(raft_scan_t *data,
		    gsl_vector **weight,
		    gsl_vector *a);


double 
raft_weight_boundKw(raft_scan_t *data,
		    gsl_vector **weight);


double 
raft_weight_boundKuw(raft_scan_t *data,
		     gsl_vector **weight);



#endif
