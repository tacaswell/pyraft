#include "filter.h"
#include "raft_param.h"
#include "raft_backprojection.h"

/*######################################################
  Private

  $Id: filter.c,v 1.18 2008-11-24 16:13:08 miqueles Exp $ - Last Update  
  ####################################################*/

/*+====================================================+
  
  FUNCTION filter_lowpass
  
  Lowpass filter.
  
  Input 
  
  w - filter argument 
  wc - cutoff frequency

  Return

  Filter evaluated at omega.
    
  +====================================================+
*/

double filter_lowpass(double w, 
		      double wc)
{
  if(fabs(w) < wc || fabs(w-wc)<ZERO)
    return 1;
  else
    return 0;
} 

/*+====================================================+
  
  FUNCTION filter_ramlak
  
  Ram-Lak filter
  
  Input 
  
  w - filter argument 
  wc - cutoff frequency

  Return

  Filter evaluated at omega.  

  +====================================================+
*/

double filter_ramlak(double w, 
		     double wc)
{
  return fabs(w) * filter_lowpass(w, wc);  
}



/*+====================================================+
  
  FUNCTION filter_shepplogan
  
  Shepp Logan filter.
  
  Input 
  
  w - filter argument 
  wc - cutoff frequency

  Return

  Filter evaluated at omega.

  +====================================================+
*/

double filter_shepplogan(double w,
			 double wc)
{
  if(fabs(w)<ZERO)
    return 0.0;
  else{
    double u, sinc;

    u = w/(2*wc);
    u = MPI*u;
    
    sinc = sin(u)/u;

    return fabs(w) * sinc * filter_lowpass(w,wc);
  }
}


/*+====================================================+
  
  FUNCTION filter_cosine
  
  Cosine filter.
  
  Input 
  
  w - filter argument 
  wc - cutoff frequency

  Return

  Filter evaluated at omega.
  
  +====================================================+
*/

double filter_cosine(double w,
		     double wc)
{
  double u, cosu;

  u = MPI*w/(2*wc);
  cosu = cos(u);

  return fabs(w) * cosu * filter_lowpass(w,wc);
}


/*+====================================================+
  
  FUNCTION filter
  
  Compute an array with filter values
  
  Input 
  
  type - filter type 
  dt - step size
  N - number os samples in the signal 
  wc - cutoff frequency

  Output

  filter - filter array

  Remark

  The Ram-Lak filter is used as the default one.

  +====================================================+
*/

void filter(int type,  
	    int N, 
	    double dt,
	    double wc, 
	    gsl_vector *filter)
{
  int i,k,half;
  double f,w,wn,dw;
  double (*filterfun)(double w, double wc);

  wn = 1/(2*dt);  
  dw = (2*wn)/N;  
  half = floor(N/2); 

  switch(type)
    {
    case RAFT_RAMLAK:
      filterfun = filter_ramlak;
      break;
    case RAFT_SHEPP:
      filterfun = filter_shepplogan;
      break;
    case RAFT_COSINE:
      filterfun = filter_cosine;
      break;
    default:
      filterfun = filter_ramlak;
    }

  for(k=0; k<N; k++)
    {
      i = k - SIGN(k-half)*half;
      w = -wn + dw*i; 

      f = filterfun(w, wc);
      
      gsl_vector_set(filter, k, f);
    }
}

/*+====================================================+
  
  FUNCTION filter_lowpass_2D
  
  Two-dimensional lowpass filter.
  
  Input 
  
  wx - filter argument (x direction)  
  wy - filter argument (y direction)
  wc - cutoff frequency

  Return

  Evaluated filter.
  
  +====================================================+
*/

double filter_lowpass_2D(double wx,
			 double wy,
			 double wc)
{
  if(fabs(wx)<=wc && fabs(wy)<=wc)
    return 1;
  else
    return 0;
} 

/*+====================================================+
  
  FUNCTION filter_ramlak_2D
  
  Two-dimensional Ram-Lak filter
  
  Input 
  
  wx - filter argument (x direction)  
  wy - filter argument (y direction)
  wc - cutoff frequency

  Return

  Evaluated filter.
  
  +====================================================+
*/

double filter_ramlak_2D(double wx,
			double wy, 
			double wc)
{
  double norm;

  norm = gsl_hypot(wx,wy);

  return norm * filter_lowpass_2D(wx,wy,wc);  
}


/*+====================================================+
  
  FUNCTION filter_shepplogan_2D
  
  Two-dimensional Shepp-Logan filter.
  
  Input 
  
  wx - filter argument (x direction)  
  wy - filter argument (y direction) 
  wc - cutoff frequency

  Return

  Evaluated filter.
        
  +====================================================+
*/

double filter_shepplogan_2D(double wx,
			    double wy,
			    double wc)
{
  double norm;

  norm = gsl_hypot(wx,wy);

  if(norm<ZERO)
    return 0.0;
  else{
    double ux, uy, sincx, sincy;

    ux = wx/(2*wc);
    ux = MPI*ux;
    
    uy = wy/(2*wc);
    uy = MPI*uy;

    if(fabs(ux)<ZERO)
      sincx=1;
    else
      sincx = sin(ux)/ux;

    if(fabs(uy)<ZERO)
      sincy=1;
    else
      sincy = sin(uy)/uy;
    
    return norm*sincx*sincy*filter_lowpass_2D(wx,wy,wc);
  }
}


/*+====================================================+
  
  FUNCTION filter_cosine_2D
  
  Two-dimensional Cosine filter.
  
  Input 
  
  wx - filter argument (x direction)  
  wy - filter argument (y direction) 
  wc - cutoff frequency

  Return

  Evaluated filter.
    
  +====================================================+
*/

double filter_cosine_2D(double wx,
			double wy,
			double wc)
{
  double norm, ux, uy, cosx, cosy;

  norm = gsl_hypot(wx,wy);
  
  ux = MPI*wx/(2*wc);
  uy = MPI*wy/(2*wc);
  cosx = cos(ux);
  cosy = cos(uy);
  
  return norm*cosx*cosy*filter_lowpass_2D(wx,wy,wc);
}


/*+====================================================+
  
  FUNCTION filter_2D
  
  Compute an array with two-dimensional filter values
  
  Input 
  
  type - filter type
  dx - step size (equal to dy)
  N - matrix size 
  wc - cutoff frequency

  Output

  filter - filter array

  Remark

  The Ram-Lak filter is used as the default one.

  +====================================================+
*/

void filter_2D(int type,  
	       int N, 
	       double dx,
	       double wc, 
	       gsl_vector *filter)
{
  int j,k,half,ij,ik;
  double f,wx,wy,wn,dw;
  double (*filterfun)(double wx, double wy, double wc);

  wn = 1/(2*dx);
  dw = 2*wn/N;  
  half = N/2; 

  switch(type)
    {
    case RAFT_RAMLAK:
      filterfun = filter_ramlak_2D;
      break;
    case RAFT_SHEPP:
      filterfun = filter_shepplogan_2D;
      break;
    case RAFT_COSINE:
      filterfun = filter_cosine_2D;
      break;
    default:
      filterfun = filter_ramlak_2D;
    }

  for(k=0; k<N; k++)
    {
      ik = k - SIGN(k-half)*half;
      wy = -wn + dw*ik;
	        
      for(j=0; j<N; j++)
	{
	  ij = j - SIGN(j-half)*half;
	  wx = -wn + dw * ij;
	  
	  f = filterfun(wx, wy, wc);
	  
	  gsl_vector_set(filter, j + k*N ,f);
	}
    }
}


/*+====================================================+
  
  FUNCTION filter_impulse_ramlak
  
  Shepp-Logan impulse response.
  
  Input 
  
  t  - domain 
  wc - cutoff frequency

  Output 
  
  Impulse response.
  
  +====================================================+
*/

double filter_impulse_ramlak(double t,
			     double wc)
{
  double T, T2, imp;
  
  T  = MPI * wc * t;
  T2 = 2*T;
        
  if(fabs(t) < ZERO)
    imp = 1;
  else
    imp = sin(T2)/T - SQR((sin(T)/T));
  
  return SQR(wc) * imp;
}

/*+====================================================+
  
  FUNCTION  filter_impulse_cosine
  
  Cosine impulse response.
  
  Input  
  
  t  - domain 
  wc - cutoff frequency

  Return 
  
  Impulse response.        
  
  +====================================================+
*/

double filter_impulse_cosine(double t,
			     double wc)
{
  double dt, imp;
  
  dt = 1/(2*wc);

  imp = (filter_impulse_ramlak(t - dt/2, wc) + 
	 filter_impulse_ramlak(t + dt/2, wc));

  return imp/2;
}

/*+====================================================+
  
  FUNCTION  filter_impulse_shepp
  
  Shepp-Logan impulse response.
  
  Input  
  
  t  - domain 
  wc - cutoff frequency

  Return 
  
  Impulse response.
  
  [from http://iria.pku.edu.cn/~jiangm/courses/IRIA/IRIA.html]
  
  +====================================================+
*/

double filter_impulse_shepp(double t,
			    double wc)
{
  double imp, bw, p, a, d, b;

  bw = 5*wc; /* new bandwidth */
  
  p = MPI/2;
  
  d = 1/(2*bw);

  if( fabs(fabs(t) - p) < ZERO)
    imp = 1/MPI;
  else
    {
      b = SQR(p) - SQR(bw*t);

      a = p - (bw*t)*sin(bw*t); 

      imp = a/b;
    }
  
  a = SQR(bw)/(SQR(MPI)*MPI);
      
  return a*imp;  
}


