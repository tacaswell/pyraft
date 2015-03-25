#include "raft_image_phantom_golosio.h"
#include <math.h>

#define MAX( x, y ) ( ( ( x ) > ( y ) ) ? ( x ) : ( y ) )

/*=====================================================*/

/*!
 * \brief Golosio's phantom: fluorescence density
 * \param image raft image
 */

void raft_image_phantom_golosio_I(raft_image image)
{
  int i,j, size;
  double max, r, x, y, att, radius, X, Y, trans, delta, dim;
  
  //desc for phantom
  att    = 1.5;
  radius = 0.16;
  trans  = 0.22;
  
  size  = raft_matrix_nlines(image.data);

  //square size!
  max = sqrt(2)/2; 
  
  //sampling distance!
  delta = 2*max/(size-1); 
  
  //scaling factor for phantom!
  dim = 1/image.tl_x;

  for(i=0; i < size; i++)
    {
      x = -max + i*delta;
      x = x/dim;
      
      for(j=0; j < size; j++)
	{
	  y = -max + j*delta;
	  y = y/dim;  
	  
	  /* 1st structure (up,left) */
	  
	  X = x+trans;
	  Y = y+trans;

	  if( MAX(fabs(X),fabs(Y)) < radius)
	    raft_matrix_element(image.data, i, j) = att;
	  
	  /* 2nd structure (down,right) */
	  
	  Y = y-trans;
	  X = x-trans;

	  r = sqrt(X*X + Y*Y);

	  if( MAX(fabs(X),fabs(Y)) < radius)
	    raft_matrix_element(image.data, i, j) = att;

	  if(fabs(r) < 0.7 * radius)
	    raft_matrix_element(image.data, i, j) = 0.7 * att;
	  
	  /* 3rd structure (down,left) */
	  
	  X = x-trans;
	  Y = y+trans;	 

	  r = sqrt(X*X + Y*Y);
	  
	  if( r < radius)
	    raft_matrix_element(image.data, i, j) = att;

	  if( fabs(X) + fabs(Y) < 0.9 * radius)
	    raft_matrix_element(image.data, i, j) = 0;
	  
	  /* 4th structure (up,right) */

	  X = x+trans;
	  Y = y-trans;

	  r = sqrt(X*X + Y*Y);
      
	  if( r < radius)
	    raft_matrix_element(image.data, i, j) = att;
	}
    }
}

/*=====================================================*/

/*!
 * \brief Golosio's phantom: transmission attenuation
 * \param image raft image
 */

void raft_image_phantom_golosio_II(raft_image image)
{
  int i,j, size;
  double max, r, x, y, att, att2, radius, X, Y, trans, square, delta, dim;
  
  //desc for phantom!
  att    = 0.2;
  att2   = 10*att;
  radius = 0.16;
  square = 0.45;
  trans  = 0.22;
  
  size  = raft_matrix_nlines(image.data);

  //square size!
  max = sqrt(2)/2; 
  
  //sampling distance!
  delta = 2*max/(size-1); 
  
  //scaling factor for phantom!
  dim = 1/image.tl_x;
  
  for(i=0; i < size; i++)
    {
      x = -max + i*delta;
      x = x/dim;
      
      for(j=0; j < size; j++)
	{
	  y = -max + j*delta;
	  y = y/dim;  
	  
	  /* outer square */

	  if( (MAX( fabs(x), fabs(y)) < square) &&
	      (MAX( fabs(x), fabs(y)) > square*(1-0.03)))
	    raft_matrix_element(image.data, i, j) = att;
	  
	  if( MAX( fabs(x), fabs(y)) < square*(1-0.03) )
	    raft_matrix_element(image.data, i, j) = att;

	  /* 1st structure (up,left) */
	  
	  X = x+trans;
	  Y = y+trans;

	  if( MAX(fabs(X),fabs(Y)) < radius)
	    raft_matrix_element(image.data, i, j) = att2;
	  
	  /* 2nd structure (down,right) */
	  
	  Y = y-trans;
	  X = x-trans;

	  r = sqrt(X*X + Y*Y);

	  if( MAX(fabs(X),fabs(Y)) < radius)
	    raft_matrix_element(image.data, i, j) = att2;
	  
	  /* 3rd structure (down,left) */
	  
	  X = x-trans;
	  Y = y+trans;	 

	  r = sqrt(X*X + Y*Y);
	  
	  if( r < radius)
	    raft_matrix_element(image.data, i, j) =  att2;

	  if( fabs(X) + fabs(Y) < 0.9 * radius)
	    raft_matrix_element(image.data, i, j) = att;

	  /* 4th structure (up,right) */

	  X = x+trans;
	  Y = y-trans;

	  r = sqrt(X*X + Y*Y);
      
	  if( r < radius)
	    raft_matrix_element(image.data, i, j) = att2;
	}
    }
}

/*=====================================================*/

/*!
 * \brief Golosio's phantom: fluorescence attenuation
 * \param image raft image
 */

void raft_image_phantom_golosio_III(raft_image image)
{
  int i,j, size;
  double max, r, x, y, att, att2, radius, X, Y, trans, square, delta, dim;
  
  //desc for phantom!
  att    = 0.5;
  att2   = 2*att;
  radius = 0.16;
  square = 0.45;
  trans  = 0.22;
  
  size  = raft_matrix_nlines(image.data);
    
  //square size!
  max = sqrt(2)/2; 
  
  //sampling distance!
  delta = 2*max/(size-1); 
  
  //scaling factor for phantom!
  dim = 1/image.tl_x;
  
  for(i=0; i < size; i++)
    {
      x = -max + i*delta;
      x = x/dim;
      
      for(j=0; j < size; j++)
	{
	  y = -max + j*delta;
	  y = y/dim;  
	  
	  /* outer square */

	  if( (MAX( fabs(x), fabs(y)) < square) &&
	      (MAX( fabs(x), fabs(y)) > square*(1-0.03)))
	    raft_matrix_element(image.data, i, j) = att;
	  
	  if( MAX( fabs(x), fabs(y)) < square*(1-0.03) )
	    raft_matrix_element(image.data, i, j) = att;

	  /* 1st structure (up,left) */
	  
	  X = x+trans;
	  Y = y+trans;

	  if( MAX(fabs(X),fabs(Y)) < radius)
	    raft_matrix_element(image.data, i, j) = att2;
	  
	  /* 2nd structure (down,right) */
	  
	  Y = y-trans;
	  X = x-trans;

	  r = sqrt(X*X + Y*Y);

	  if( MAX(fabs(X),fabs(Y)) < radius)
	    raft_matrix_element(image.data, i, j) = att2;

	  if(fabs(r) < 0.7 * radius)
	    raft_matrix_element(image.data, i, j) = 1.1 * att2;
	  
	  /* 3rd structure (down,left) */
	  
	  X = x-trans;
	  Y = y+trans;	 

	  r = sqrt(X*X + Y*Y);
	  
	  if( r < radius)
	    raft_matrix_element(image.data, i, j) = att2;

	  if( fabs(X) + fabs(Y) < 0.9 * radius)
	    raft_matrix_element(image.data, i, j) = att;
	  
	  /* 4th structure (up,right) */

	  X = x+trans;
	  Y = y-trans;

	  r = sqrt(X*X + Y*Y);
      
	  if( r < radius)
	    raft_matrix_element(image.data, i, j) = att2;
	}
    }

}
