#include "raft_image_phantom_gear.h"
#include <math.h>

#define MAX( x, y ) ( ( ( x ) > ( y ) ) ? ( x ) : ( y ) )
#define MPI 3.14159265358979

/*=====================================================*/

/*!
 * \brief Gear phantom: micro-ct
 * \param image raft image
 */

void raft_image_phantom_gear(raft_image image)
{
  int i,j,k, N;
  double max, delta, dim;
  double x, y, th0, th, dth=45*MPI/180;
  double h, a, b;
  double B=1.5, R=1.9, T=1, t=0.25, r=0.4, q=0.6;
  
  N = raft_matrix_nlines(image.data);
  
  //square size!
  //max = sqrt(2)/2; 
  max = 3;
  
  
  //sampling distance!
  delta = 2*max/(N-1); 
  
  //scaling factor for phantom!
  //dim = 1/image.tl_x;
  dim = 1;

  for (i=0; i< N; i++) 
    {
      x = -max + i*delta;

      for (j=0; j < N; j++)
	{
	  y = -max + j *delta;
	  
	  /* 0th circle */
	  
	  h = hypot(x,y);
	  
	  if ( h*h <= B*B ){
	    //matrix[i + j*N] = 2;
	    raft_matrix_element(image.data, j, i) = 2;
	  }
	  
	  /* 1st circle */
	  
	  h = hypot(x,y);

	  if ( h*h <= r*r ){
	    //matrix[i + j*N] = 0;
	    raft_matrix_element(image.data, j, i) = 0;
	  }

	  /* 2nd circles */
	  
	  th0 = 22*MPI/180;
	  for (k=0; k<8; k++)
	    {
	      th = th0 + k*dth;

	      a = T*cos(th);
	      b = T*sin(th);

	      h = hypot( x-a, y-b );
	      
	      if ( h*h <= t*t )
		{
		  //matrix[ i + j*N] = 0.8 * fabs(th);
		  raft_matrix_element(image.data, j, i) = 0.8 * fabs(th);
		}
	    }
	  
	  /* 3rd circles */

	  th0 = 0;
	  for (k=0; k<8; k++)
	    {
	      th = th0 + k*dth;

	      a = R*cos(th);
	      b = R*sin(th);

	      h = hypot( x-a, y-b );
	      
	      if ( h*h <= q*q )
		{
		  //matrix[ i + j*N] = 0;	      
		  raft_matrix_element(image.data, j, i) = 0;
		}
	    }

	  /* sample holder */
	 
	  /* 
	  if ( fabs(x)<0.2 && fabs(y-2.3) < 1  )
	    {
	      //matrix[ i + j*N] = 0.5;
	      raft_matrix_element(image.data, j, i) = 0.5;
	    }
          */
	  	  
	}
    }
}  
