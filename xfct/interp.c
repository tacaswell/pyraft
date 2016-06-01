#include "interp.h"
#include "raft_param.h"
#include "raft_scan.h"

/*######################################################
  Private
  
  $Id: interp.c,v 1.11 2008-11-07 21:13:49 miqueles Exp $ - Last Update
  ####################################################*/

/*+====================================================+
  
  FUNCTION interp_search
  
  Find index of a given value within a
  one-dimensional array.

  Input 

  v - array
  z - value  
    
  Return

  Return -1 whenever the index do not exist or out 
  of domain. Otherwise, the requested index.

  +====================================================+
*/

int interp_search(gsl_vector *v, double z)
{
  int i, N, ind;
  double vi, vii;

  N = v->size;

  vi = gsl_vector_get(v,0);
  vii = gsl_vector_get(v,N-1);

  ind = -1;

  if(z < vi || z > vii)
    return -1;
  else
    {
      for(i=0; i < N-1; i++)
	{
	  vi = gsl_vector_get(v,i);
	  vii = gsl_vector_get(v,i+1);
	  
	  if(z>= vi && z<= vii)
	    {
	      ind = i;
	      break;
	    }
	}
    }
  
  if(SIGN(ind)<0)
    return -1;
  else
    return ind;
}


/*+====================================================+
  
  FUNCTION interp_nearest
  
  Find nearest integer of a given value, within a
  one-dimensional array.

  Input 

  x - array
  z - value  
  linear - identify a linear array
  
  Return

  Return -1 whenever the nearest value is negative or out 
  of domain. Otherwise, the requested nearest value.
  
  +====================================================+
*/

int interp_nearest(gsl_vector *x,
		   double z,
		   int linear)
  
{
  if(linear == RAFT_REG)
    {
      int nearest, N;
      double y, dx, x0;

      N = x->size;
      
      x0 = gsl_vector_get(x,0);
      dx = gsl_vector_get(x,1) - x0;
      
      y = (z-x0)/dx;
      
      nearest = (int) floor(y);
      
      if(nearest<0 || nearest>=N)
	return -1;
      else
	return nearest;
    }
  else
    {
      return interp_search(x,z);
    }
}


/*+====================================================+
  
  FUNCTION interp_linear
  
  Linear interpolation at a given index 
  of a vector.

  Input 

  f - function 
  x - domain
  index - nearest index
  z - value  
  
  Return

  Interpolated value.
  
  +====================================================+
*/

double interp_linear(gsl_vector *f,
		     gsl_vector *x,
		     int index,
		     double z)
{
  int size,k;
  
  k = index;
  size = f->size;

  if(k+1==size)
    return gsl_vector_get(f, k);
  else
    {
      double fk, xk, df, dx, L;
      
      fk = gsl_vector_get(f, k);
      df = gsl_vector_get(f, k+1) - fk;
      
      xk = gsl_vector_get(x, k);
      dx = gsl_vector_get(x, k+1) - xk;
      
      L = (df/dx)*(z - xk) + fk;    
      
      return L;
    }
}

/*+====================================================+
  
  FUNCTION interp_linear_reverse
  
  Reverse linear interpolation at a given index 
  of a vector.

  Input 

  f - function 
  x - domain
  index - nearest index
  z - value  
  
  Return

  Interpolated value.
  
  +====================================================+
*/

double interp_linear_reverse(gsl_vector *f,
			     gsl_vector *x,
			     int index,
			     double z)
{
  int size,k;
  
  k = index;
  size = f->size;

  if(k+1==size)
    return gsl_vector_get(f, k);
  else
    {
      double fk, xk, df, dx, L;
      
      fk = gsl_vector_get(f, size-1-k);
      df = gsl_vector_get(f, size-1-k-1) - fk;
      
      xk = gsl_vector_get(x, k);
      dx = gsl_vector_get(x, k+1) - xk;
      
      L = (df/dx)*(z - xk) + fk;    
      
      return L;
    }
}
