#ifndef _INTERP_H_
#define _INTERP_H_

#include <gsl/gsl_vector.h>
#include <stdlib.h>
#include <math.h>

int 
interp_search(gsl_vector *v, 
	      double z);



int 
interp_nearest(gsl_vector *x,
	       double z,
	       int linear);



double
interp_linear(gsl_vector *f,
	      gsl_vector *x,
	      int index,
	      double z);


double 
interp_linear_reverse(gsl_vector *f,
		      gsl_vector *x,
		      int index,
		      double z);



#endif
