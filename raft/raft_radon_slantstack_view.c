#include "raft_radon_slantstack_view.h"
#include "raft_vector.h"
#include <math.h>

#define MAX( x, y ) ( ( ( x ) > ( y ) ) ? ( x ) : ( y ) )
#define SIGN( x ) ( ( ( x ) > 0.0 ) ? 1.0 : ( ( ( x ) < 0.0 ) ? -1.0 : 0.0 ) )
#define PI 3.1415926535897932384626433832795

double eval_integral_view(raft_image phantom,
			  double t,
			  double theta,
			  int nrays);

/*========================================================*/

/*!
 * \brief View of Radon transform: slant stack
 * \param phantom raft image
 * \param nthreads number of threads
 */

void raft_radon_slantstack_view(raft_vector view,
				raft_image phantom,
				raft_vector t,
				double theta)
{
  int k;
  double sum;
  
  for( k = 0; k < t.size; k++)
    {
      sum = eval_integral_view(phantom,
			       raft_vector_element(t,k),
			       theta,
			       t.size);

      raft_vector_element(view, k) = sum;
    }
}

double eval_integral_view(raft_image phantom,
			  double t,
			  double theta,
			  int nrays)
{
  int i, size, m;
  double sum, cost, sint, tant;
  double x, xmin, xmax, dx;
  double y, ymin, ymax, dy;
  double fsint, fcost, tol, value, Delta;
   
  //
  tol  = 1.0/sqrt(2);
  size = raft_matrix_nlines(phantom.data);
  
  // image corners
  xmin = phantom.tl_x;
  xmax = phantom.br_x;
  
  ymin = phantom.br_y;
  ymax = phantom.tl_y;

  //sampling distance for phantom
  dx = (xmax - xmin)/(size - 1);
  dy = (ymax - ymin)/(size - 1);
  
  //angle 
  cost  = cos(theta);
  sint  = sin(theta);
  tant  = tan(theta);  
  fsint = fabs(sint);
  fcost = fabs(cost);
  
  //ray @ angle
  
  if(fcost < tol) 
    {
      sum   = 0;
      Delta = (dy/fsint);

      for(i=0; i < size; i++)
	{
	  x = (double) (xmin + i * dx);
	  y = (double) (t/sint - x/tant);
	  
	  m = (int) floor((y-ymin)/dy);

	  if(m>-1 && m < size-1)	  
	    {
	      sum += (raft_matrix_element(phantom.data,m+1,i)-raft_matrix_element(phantom.data,m,i))*(y-(ymin+m*dy))/dy + raft_matrix_element(phantom.data,m,i);
	    }
	}
      value = sum*Delta;
    } 
  else
    {
      sum   = 0;
      Delta = (dx/fcost);

      for(i=0; i < size; i++)
	{
	  y = (double) (xmin + i * dy); 
	  x = (double) (t/cost - y*tant);
	  
          m = (int) floor ( (x-xmin)/dx);
	  
	  if(m>-1 && m < size-1)	  
	    {
	      sum += (raft_matrix_element(phantom.data,i,m+1)-raft_matrix_element(phantom.data,i,m))*(x-(xmin+m*dx))/dx + raft_matrix_element(phantom.data,i,m);
	    }	  
	}
      value = sum*Delta;
    }
  
  return value;
}


