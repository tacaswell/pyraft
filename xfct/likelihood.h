#ifndef _LIKELIHOOD_H_
#define _LIKELIHOOD_H_

#include "raft_param.h"
#include "raft_projection.h"
#include "raft_iterative.h"
#include <gsl/gsl_linalg.h>

double 
eval_loglikelihood_emission(raft_proj_t *proj,
			    gsl_vector *radon,
			    gsl_vector *logradon,
			    gsl_vector *p);


double
eval_loglikelihood_transmission(raft_scan_t *data,
				raft_iterative_t *workspace,
				gsl_vector *f,
				gsl_vector *ph);



#endif
