#ifndef _RIDGE_H_
#define _RIDGE_H_

#include "raft_ridge.h"
#include "raft_param.h"
#include "raft_math.h"
#include "interp.h"

#include <gsl/gsl_linalg.h>
#include <math.h>


double 
eval_ridgesum(raft_scan_t *data,
	      raft_ridge_t *workspace,
	      int j,
	      int k);


double 
eval_intprojc(raft_scan_t *data,
	      gsl_vector *p,
	      gsl_vector *U,
	      gsl_vector *s0);


#endif
