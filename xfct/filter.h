#ifndef _FILTER_H_
#define _FILTER_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <math.h>

double
filter_lowpass(double w, 
	       double wc);


double 
filter_ramlak(double w, 
	      double wc);



double 
filter_shepplogan(double w, 
		  double wc);



double 
filter_cosine(double w, 
	      double wc);


double
filter_lowpass_2D(double wx, 
		  double wy,
		  double wc);


double 
filter_ramlak_2D(double wx, 
		 double wy,
		 double wc);



double 
filter_shepplogan_2D(double wx,
		     double wy,
		     double wc);


double 
filter_cosine_2D(double wx, 
		 double wy,
		 double wc);


double 
filter_impulse_ramlak(double t,
		      double wc);


double 
filter_impulse_cosine(double t,
		      double wc);



double
filter_impulse_shepp(double t,
		     double wc);



void
filter(int type,
       int N,
       double dt,
       double wc,
       gsl_vector *filter);


void
filter_2D(int type,
	  int N,
	  double dx,
	  double wc,
	  gsl_vector *filter);


#endif
