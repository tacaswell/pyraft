#ifndef _CBF_H_
#define _CBF_H_

#include "interp.h"
#include "raft_param.h"
#include "raft_phantom.h"
#include "raft_scan.h"
#include "raft_cbf.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

void
eval_rayintegral_dbt_cbf(raft_scan_t *data, 
			 gsl_vector *att,
			 int j,
			 int xi,
			 int yi,
			 raft_cbf_t *cbf,
			 double *current,
			 double *forward,
			 double *backward);

#endif

