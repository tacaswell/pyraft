#ifndef _RAFT_LIKELIHOOD_H_
#define _RAFT_LIKELIHOOD_H_

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include "raft_scan.h"
#include "raft_errno.h"
#include "raft_projection.h"


double
raft_loglikelihood_pet(raft_scan_t *data,
		       raft_proj_t *workspace,
		       gsl_vector *p,
		       gsl_vector *act,
		       gsl_vector *att);


double
raft_loglikelihood_spect(raft_scan_t *data,
			 raft_proj_t *workspace,
			 gsl_vector *p,
			 gsl_vector *act,
			 gsl_vector *att);


double 
raft_loglikelihood_xfct(raft_scan_t *data,
			raft_proj_t *workspace,
			gsl_vector *p,
			gsl_vector *act,
			gsl_vector *attT,
			gsl_vector *attF);

#endif
