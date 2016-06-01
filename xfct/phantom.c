#include "raft_phantom.h"
#include "raft_param.h"
#include <gsl/gsl_math.h>

#define EPS 10e-08
#define DEFINE_BOUNDARY 0

/*######################################################
  Title: Phantom
  
  Header - <raft/raft_phantom.h>
  Type - <raft_phantom_t>
  $Id: phantom.c,v 1.20 2009-12-06 15:13:16 miqueles Exp $ - Last Update  
  ####################################################*/

/*######################################################
  Public
  ####################################################*/

/*######################################################
  Section: Allocation
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_phantom_alloc

  Allocattes a phantom structure for RAFT simulation. 

  Input: 

  Empty.

  Return:

  Phantom structure. See <raft_phantom_t>.
  
  +====================================================+
*/

raft_phantom_t *raft_phantom_alloc(void)
{
  return (raft_phantom_t *)malloc(sizeof(raft_phantom_t));  
}

/*+====================================================+
  
  FUNCTION: raft_phantom_free

  Frees the phantom structure. 

  Input: 

  Phantom structure. See <raft_phantom_t>.

  +====================================================+
*/

void raft_phantom_free(raft_phantom_t *phantom)
{
  free(phantom);
}

/*+====================================================+
  
  FUNCTION: raft_phantom_alloc_data

  Allocattes a phantom structure data for RAFT simulation. 

  Input: 
  
  phantom - RAFT phantom data type. See <raft_phantom_t>.
  size - number of pixels vertical/horizontal direction
  min - lower bound for horizontal axis
  max - upper bound for vertical axis
  etype - energy type (RAFT_MONO or RAFT_POLY)

  Return:

  RAFT_SUCCESS - Successfull operation
  RAFT_ENOMEM - Not enough memory for allocattion
  RAFT_EDOM - Domain error, wrong input data.
   
  +====================================================+
*/

int raft_phantom_alloc_data(raft_phantom_t *phantom, 
			    int size,
			    double min, 
			    double max)
{
  int i;
  double z;
  
  gsl_set_error_handler_off();
  
  if(min >= max || size<=0)
    return RAFT_EDOM;
  
  phantom->x = gsl_vector_alloc(size);
  if(!phantom->x)
    return RAFT_ENOMEM;
  
  phantom->y = gsl_vector_alloc(size);
  if(!phantom->y)
    return RAFT_ENOMEM;
  
  phantom->vector = gsl_vector_calloc(size*size);
  if(!phantom->vector)
    return RAFT_ENOMEM;
  
  phantom->min = min;
  phantom->max = max;
  phantom->size = size;
  
  /* x = linspace(a,b,n) & z = linspace(c,d,n) */
  
  for(i=0; i< size-1; i++)
    {
      z = min + i*(max - min)/(size-1);
      
      gsl_vector_set(phantom->x, i, z);
      gsl_vector_set(phantom->y, i, z);
    }
  
  gsl_vector_set(phantom->x, size-1, max);
  gsl_vector_set(phantom->y, size-1, max);

  phantom->etype = RAFT_MONO;
  phantom->nsubs = 0;
  
  return RAFT_SUCCESS;
}


/*+====================================================+
  
  FUNCTION: raft_phantom_alloc_data_basis

  Allocattes a phantom basis for polychromatic experiments. 

  Input: 
  
  phantom - RAFT phantom data type. See <raft_phantom_t>.
  nsubs - number of basis substances

  Return:

  RAFT_SUCCESS - Successfull operation
  RAFT_ENOMEM - Not enough memory for allocattion
  RAFT_EDOM - Domain error, wrong input data.
  
  +====================================================+
*/

int raft_phantom_alloc_data_basis(raft_phantom_t *phantom, 
				  int nsubs)
{
  int i, pixels;
  
  if(nsubs<=0)
    return RAFT_EDOM;

  pixels = (phantom->size)*(phantom->size);
  
  phantom->etype = RAFT_POLY;
  phantom->nsubs = nsubs;

  phantom->basis = (basis_t *)malloc(sizeof(basis_t)*nsubs);
  if(!phantom->basis)
    return RAFT_ENOMEM;
  
  for(i=0; i < nsubs; i++)
    {
      phantom->basis[i].vector = gsl_vector_alloc(pixels);
    }
    
  return RAFT_SUCCESS;
}


/*+====================================================+
  
  FUNCTION: raft_phantom_free_data

  Frees the phantom structure data. 

  Input: 

  RAFT phantom structure. See <raft_phantom_t>.
  
  +====================================================+
*/
 
void raft_phantom_free_data(raft_phantom_t *phantom)
{
  gsl_vector_free(phantom->x);
  gsl_vector_free(phantom->y);
  gsl_vector_free(phantom->vector);      

  if(phantom->etype == RAFT_POLY)
    {
      int i;
      
      for(i=0; i < phantom->nsubs; i++)
	gsl_vector_free(phantom->basis[i].vector);
      
      free(phantom->basis);
    }
}

/*######################################################
  Section: Getting phantom parameters
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_phantom_get_size

  Get the phantom size. 

  Input: 
  
  P - RAFT phantom data type. See <raft_phantom_t>.

  Return:

  Phantom size (i.e, number of pixels in both directions)
  
  +====================================================+
*/

int raft_phantom_get_size(raft_phantom_t *P)
{
  return P->size;
}

/*+====================================================+
  
  FUNCTION: raft_phantom_get_min

  Get the axis minimum coordinate from a phantom. 

  Input: 
  
  P - RAFT phantom data type. See <raft_phantom_t>.

  Return:

  Minimum coordinate.   
  
  +====================================================+
*/

double raft_phantom_get_min(raft_phantom_t *P)
{
  return P->min;
}

/*+====================================================+
  
  FUNCTION: raft_phantom_get_max

  Get the axis maximum coordinate from a phantom. 

  Input: 
  
  P - RAFT phantom data type. See <raft_phantom_t>.

  Return:

  Maximum coordinate.
  
  +====================================================+
*/

double raft_phantom_get_max(raft_phantom_t *P)
{
  return P->max;
}

/*+====================================================+
  
  FUNCTION: raft_phantom_get_x

  Get the ith element of the x-axis from a phantom. 

  Input: 
  
  P - RAFT phantom data type. See <raft_phantom_t>.
  i - index

  Return:

  The ith position of the x-axis.
  
  +====================================================+
*/


double raft_phantom_get_x(raft_phantom_t *P,
			  int i)
{
  return gsl_vector_get(P->x,i);
}

/*+====================================================+
  
  FUNCTION: raft_phantom_get_y

  Get the ith element of the y-axis from a phantom. 

  Input: 
  
  P - RAFT phantom data type. See <raft_phantom_t>.
  i - index

  Return:

  The ith position of the y-axis.
  
  +====================================================+
*/

double raft_phantom_get_y(raft_phantom_t *P,
			 int i)
{
  return gsl_vector_get(P->y, i);
}

/*+====================================================+
  
  FUNCTION: raft_phantom_get_step

  Get the step in the x,y axis. 

  Input: 
  
  P - RAFT phantom data type. See <raft_phantom_t>.
  
  Return:

  The step in x,y axis.
  
  +====================================================+
*/

double raft_phantom_get_step(raft_phantom_t *phantom)
{
  return raft_phantom_get_x(phantom,1) - 
    raft_phantom_get_x(phantom,0);
}


/*+====================================================+
  
  FUNCTION: raft_phantom_get

  Get a position of the phantom matrix. 

  Input: 
  
  P - phantom vector
  size - matrix size
  i - row index 
  j - column index
  
  Return:

  Matrix entry.
  
  +====================================================+
*/

double raft_phantom_get(gsl_vector *P,
			int i,
			int j)
{
  int size;

  size = sqrt(P->size);

  return gsl_vector_get(P, j + i*size);
}

/*+====================================================+
  
  FUNCTION: raft_phantom_get_vector_basis

  Get the kth phantom vector. 

  Input: 
  
  P - RAFT phantom data. See <raft_phantom_t>.
  k - basis index

  Return:

  Phantom vector or a null pointer if there is no 
  previously allocated basis.
  
  +====================================================+
*/

gsl_vector *
raft_phantom_get_vector_basis(raft_phantom_t *phantom,
			      int k)
{
  if(phantom->etype == RAFT_POLY)
    {
      return phantom->basis[k].vector;
    }
  else
    {
      return NULL;
    }
}


/*+====================================================+
  
  FUNCTION: raft_phantom_get_vector

  Get the phantom vector. 

  Input: 
  
  P - RAFT phantom data type
    
  Return:

  Phantom vector.
  
  +====================================================+
*/

gsl_vector *raft_phantom_get_vector(raft_phantom_t *P)
{
  return P->vector;
}


/*######################################################
  Section: Setting phantom parameters
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_phantom_set

  Set an entry of the phantom matrix. 

  Input: 
  
  P - phantom vector
  i - row index 
  j - column index
  z - matrix entry
  
  +====================================================+
*/

void raft_phantom_set(gsl_vector *P, 
		      int i,
		      int j, 
		      double z)
{
  int size;

  size = sqrt(P->size);

  gsl_vector_set(P, j + i*size, z);
}

/*+====================================================+
  
  FUNCTION: raft_phantom_set_basis

  Set an entry of the kth phantom matrix. 

  Input: 
  
  phantom - RAFT phantom data type. See <raft_phantom_t>.
  k - basis index
  i - row index
  j - column index
  z - matrix entry
    
  Return:

  RAFT_EDOM - Domain error.
  RAFT_SUCCESS - Successfull operation.
  
  +====================================================+
*/

int raft_phantom_set_basis(raft_phantom_t *phantom, 
			   int k, 
			   int i,
			   int j, 
			   double z)
{
  if(k < 0 || k > phantom->nsubs)
    return RAFT_EDOM;
  else
    {
      raft_phantom_set(phantom->basis[k].vector, i, j, z);
      
      return RAFT_SUCCESS;
    }
}

/*+====================================================+
  
  FUNCTION: raft_phantom_set_default_basis
  
  Define a default basis for a monochromatic phantom.
  
  The basis is computed using a singular value decomposition
  of the energy matrix.
   
  Input: 

  phantom - phantom data . See <raft_phantom_t>..
  data- scan data . See <raft_scan_t>..
  f - phantom vector.
  
  +====================================================+
*/

void raft_phantom_set_default_basis(raft_phantom_t *phantom,
				       raft_scan_t *data,
				       gsl_vector *f)
{
  int j, k, i, size, nrays, nviews, nergs, nsubs;
  double fun;
  
  gsl_vector *q, *aux;
    
  size = raft_scan_get_size(data);
  nrays = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);
  nergs = raft_scan_get_nenergies(data);
  nsubs = raft_scan_get_nsubs(data);
  
  q = gsl_vector_alloc(nsubs);
  aux = gsl_vector_alloc(nsubs);
    
  gsl_vector_set_all(data->energy->ones, 1.0);

  gsl_matrix_memcpy(data->energy->U, data->energy->G);
  
  gsl_linalg_SV_decomp_jacobi(data->energy->U, 
			      data->energy->V, 
			      data->energy->S);

  data->energy->decomposed = 1;

  /*--------------------*/
  /* q = pinv(G) * ones */
  
  gsl_blas_dgemv(CblasTrans, 1.0, data->energy->U, data->energy->ones, 0.0, q);
	  
  for(k=0; k<nsubs; k++)
    {
      double sk, div;
      
      sk = gsl_vector_get(data->energy->S,k);
      
      if(fabs(sk) < 1e-06)
	gsl_vector_set(q,k,0.0);
      else
	{
	  div = gsl_vector_get(q,k)/sk;
	  gsl_vector_set(q, k, div);
	}
    }
  
  gsl_blas_dgemv(CblasNoTrans, 1.0, data->energy->V, q, 0.0, aux);
  
  gsl_vector_memcpy(q, aux);

  /*-----------*/
  /* Set basis */
  
  for(i=0; i<size; i++)
    {
      for(j=0; j<size; j++)
	{
	  fun = raft_phantom_get(f, i, j);
	  
	  for(k=0; k<nsubs; k++)
	    {
	      double qq;
	      
	      qq = gsl_vector_get(q, k) * fun;

	      raft_phantom_set(phantom->basis[k].vector, i, j, qq);
	    }
	}
    }
}

/*######################################################
  Section: Special phantoms
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_phantom_define_ghost
  
  Define a phantom for RAFT simulation, with a ghost 
  shape. 

  Input: 

  P - RAFT phantom data type. See <raft_phantom_t>.
  
  +====================================================+
*/

void raft_phantom_define_ghost(raft_phantom_t *P)
{
  int i,j, size;
  double aux1, aux2, S1, S2, S3, S4, S5, x, y, scale;
  double aa[2] = { 0.0, -0.45}; 
  double bb[2] = {-0.5, -0.8};
  double cc[2] = { 0.5, -0.8};
  
  size = raft_phantom_get_size(P);
  scale = raft_phantom_get_max(P);

  for(i=0; i < size; i++)
    {
      x = raft_phantom_get_x(P,i);
      x = x/scale;
      
      for(j=0; j < size; j++)
	{
	  y = raft_phantom_get_y(P, j);
	  y = y/scale;
      
	  aux1 = ((aa[0]-bb[0])/(aa[1]-bb[1]))*(x-bb[0])+bb[1];
	  aux2 = ((aa[0]-cc[0])/(aa[1]-cc[1]))*(x-cc[0])+cc[1];
	  
	  if((y<=aux1) && (y<=aux2) && (y>=bb[1]))   
	    raft_phantom_set(P->vector, size-1-j, size-1-i, WHITE);
	  
	  S1 = (x+0.045-0.4)*(x+0.045-0.4) + y*y;
	  S2 = (x-0.1-0.4)*(x-0.1-0.4) + (y-0.2)*(y-0.2);
	  S3 = (x-0.1-0.4)*(x-0.1-0.4) + (y+0.2)*(y+0.2);
	  S4 = (x+0.125-0.4)*(x+0.125-0.4) + y*y;
	  S5 = (x+0.345-0.4)*(x+0.345-0.4) + 0.125*y*y;
	  
	  if((S1<=0.25) && (S2>=0.0125) && (S3>=0.0125) &&
	     (S4>=0.00325) && (S5>=0.00325))
	    raft_phantom_set(P->vector, size-1-i, size-1-j, WHITE);
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_phantom_define_ring

  Define a phantom for RAFT simulation, with a ring shape. 

  Input: 

  P - RAFT phantom data type. See <raft_phantom_t>. 
  
  +====================================================+
*/

void raft_phantom_define_ring(raft_phantom_t *P)
{
  int i,j, size;
  double x, y, aux, scale;
    
  size = raft_phantom_get_size(P);
  scale = raft_phantom_get_max(P);

  for(i=0; i < size; i++)
    {
      x = raft_phantom_get_x(P, i);
      x = x/scale;
      
      for(j=0; j < size; j++)
	{
	  y = raft_phantom_get_y(P,j);
	  y = y/scale;
	  
	  aux = x*x + y*y;
	  
	  if ((aux <=0.1) && (aux >= 0.03))
	    raft_phantom_set(P->vector, i, j, WHITE);
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_phantom_define_elipring

  Define a phantom for RAFT simulation, with an elliptic 
  ring shape. 

  Input: 

  P - RAFT phantom data type. See <raft_phantom_t>. 
   
  +====================================================+
*/

void raft_phantom_define_elipring(raft_phantom_t *P)
{
  int i,j, size;
  double x, y, aux, scale;
  
  size = raft_phantom_get_size(P);
  scale = raft_phantom_get_max(P);

  for(i=0;i < size; i++)
    {
      x = raft_phantom_get_x(P,i);
      x = x/scale;
      
      for(j=0;j < size; j++)
	{
	  y = raft_phantom_get_y(P,j);
	  y = y/scale;

	  aux = 0.125*x*x + y*y;
	  
	  if ((aux <=0.05) && (aux >= 0.03))
	    raft_phantom_set(P->vector, i, j, WHITE);	
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_phantom_define_ball

  Define a phantom for RAFT simulation, with a ball shape. 

  Input: 

  P - RAFT phantom data type. See <raft_phantom_t>. 
  
  +====================================================+
*/

void raft_phantom_define_ball(raft_phantom_t *P)
{
  int i,j, size;
  double x, y, scale;

  size = raft_phantom_get_size(P);
  scale = raft_phantom_get_max(P);
  
  for(i=0; i < size; i++)
    {
      x = raft_phantom_get_x(P,i);
      x = x/scale;
            
      for(j=0; j < size; j++)
	{
	  y = raft_phantom_get_y(P, j);
	  y = y/scale;
	  
	  if((x-0.31)*(x-0.31) + (y-0.31)*(y-0.31) <= 0.004)
	    raft_phantom_set(P->vector,i,j,WHITE);	  

	  if((x-0.41)*(x-0.41) + (y-0.15)*(y-0.15) <= 0.004)
	    raft_phantom_set(P->vector,i,j,WHITE);

	  if((x-0.51)*(x-0.51) + (y-0.31)*(y-0.31) <= 0.004)
	    raft_phantom_set(P->vector,i,j,WHITE);
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_phantom_define_rectang
  
  Define a phantom for RAFT simulation, with a rectangle
  shape. 

  Input: 

  P - RAFT phantom data type. See <raft_phantom_t>. 
  
  +====================================================+
*/

void raft_phantom_define_rectang(raft_phantom_t *P)
{
  int i,j, size;
  double x, y, scale;
  
  size = raft_phantom_get_size(P);
  scale = raft_phantom_get_max(P);
  
  for(i=0; i < size; i++)
    {
      x = raft_phantom_get_x(P, i);
      x = x/scale;
      
      for(j=0; j < size; j++)
	{
	  y = raft_phantom_get_y(P, j);
	  y = y/scale;

	  if(fabs(x)<=0.5 && fabs(y)<=0.05)
	    raft_phantom_set(P->vector, i, j, WHITE);
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_phantom_define_square

  Define a phantom for RAFT simulation, with a square shape. 

  Input: 

  P - RAFT phantom data type. See <raft_phantom_t>. 
  
  +====================================================+
*/

void raft_phantom_define_square(raft_phantom_t *P)
{
  int i,j, size;
  double x, y, scale;
  
  size = raft_phantom_get_size(P);
  scale = raft_phantom_get_max(P);

  for(i=0; i < size; i++)
    {
      x = raft_phantom_get_x(P, i);
      x = x/scale;
      
      for(j=0; j < size; j++)
	{
	  y = raft_phantom_get_y(P, j);
	  y = y/scale;
	  
	  if(fabs(x)<=0.10 && fabs(y)<=0.10)
	    raft_phantom_set(P->vector, i, j, WHITE);
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_phantom_define_shiftsquare

  Define a phantom for RAFT simulation, with a shifted 
  square shape. 

  Input: 

  P - RAFT phantom data type. See <raft_phantom_t>. 
  
  +====================================================+
*/

void raft_phantom_define_shiftsquare(raft_phantom_t *P)
{
  int i,j,size;
  double x, y, scale;
  
  size = raft_phantom_get_size(P);
  scale = raft_phantom_get_max(P);

  for(i=0; i < size; i++)
    {
      x = raft_phantom_get_x(P, i);
      x = x/scale;

      for(j=0; j < size; j++)
	{
	  y = raft_phantom_get_y(P, j);
	  y = y/scale;

	  if(fabs(x-0.31)<=0.09 && fabs(y+0.31)<=0.09)
	    raft_phantom_set(P->vector, i, j, WHITE);
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_phantom_define_star

  Define a phantom for RAFT simulation, with a star shape. 

  Input: 

  P - RAFT phantom data type. See <raft_phantom_t>. 
  
  +====================================================+
*/

void raft_phantom_define_star(raft_phantom_t *P)
{
  int i,j, size;
  double S1, S2, S3, S4, S5, x, y, scale;
  double a1=0.5, a2=0.3, b1=2, b2=0.6, cc=0.1;
  
  size = raft_phantom_get_size(P);
  scale = raft_phantom_get_max(P);

  for(i=0; i < size; i++)
    {
      x = raft_phantom_get_x(P,i);
      x = x/scale;

      for(j=0; j < size; j++)
	{
	  y = raft_phantom_get_y(P, j);
	  y = y/scale;
	  
	  S1 = y - (-a1*x + a2);
	  S2 = y - ( a1*x - a2);
	  S3 = y - (-b1*x - b2);
	  S4 = y - ( b1*x + b2);
	  S5 = x - cc;
	
	  if((S1<=0) && (S2>=0) && (S3>=0) && (S4<=0) && (S5<=0))
	    raft_phantom_set(P->vector, size-1-i, size-1-j, WHITE);
	  	  
	  if((S1<=0) && (S2>=0) && (S5>=0))
	    raft_phantom_set(P->vector, size-1-i, size-1-j, WHITE);
	      
	  if((S1<=0) && (S3>=0) && (S4>=0))
	    raft_phantom_set(P->vector, size-1-i, size-1-j, WHITE);
	  
	  if((S2>=0) && (S4<=0) && (S3<=0))
	    raft_phantom_set(P->vector, size-1-i, size-1-j, WHITE);
	  	  
	  if((S4<=0) && (S5<=0) && (S1>=0))
	    raft_phantom_set(P->vector, size-1-i, size-1-j, WHITE);
	    	  
	  if((S5<=0) && (S3>=0) && (S2<=0))
	    raft_phantom_set(P->vector, size-1-i, size-1-j, WHITE);
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_phantom_define_sheeplogan

  Define the Sheep-Logan phantom for RAFT simulation. 
  
  Input: 

  P - RAFT phantom data type. See <raft_phantom_t>. 
  
  +====================================================+
*/
 
void raft_phantom_define_sheeplogan(raft_phantom_t *P)
{
  int i, j, k, size;
  double fun, u1,u2, aa , bb, cosa, sina, x, y, z, scale;

  if(DEFINE_BOUNDARY)
    {
      double table_head[10][6] = {
	{0,0,0,0,0,0},
	{0.0, -1.84e-002, 6.624e-001, 8.74e-001, 0.0, 1},
	{0,0,0,0,0,0},
	{0,0,0,0,0,0},
	{0,0,0,0,0,0},
	{0,0,0,0,0,0},
	{0,0,0,0,0,0},
	{0,0,0,0,0,0},
	{0,0,0,0,0,0},
	{0,0,0,0,0,0}
      };

      struct{
	double c1, c2, A, B, theta, density;
      }data[10];   
      
      size = raft_phantom_get_size(P);
      scale = raft_phantom_get_max(P);
      
      for(k=0;k<10;k++){
	
	data[k].c1 = table_head[k][0];
	data[k].c2 = table_head[k][1];
	data[k].A  = table_head[k][2];
	data[k].B  = table_head[k][3];
	data[k].theta = table_head[k][4];
	data[k].density = table_head[k][5];
	
	aa = data[k].A * data[k].A;
	bb = data[k].B * data[k].B;
	
	for(j=size-1; j>=0; j--)
	  {
	    y = raft_phantom_get_y(P, j);
	    
	    y=y/scale; 
	    
	    for(i=0; i < size; i++)
	      {
		double U1, U2; 
		
		x = raft_phantom_get_x(P, i);
		
		x=x/scale;
		
		u1 = x + data[k].c1;
		u2 = y + data[k].c2;
		
		cosa = cos(data[k].theta);
		sina = sin(data[k].theta);
		U1 =  cosa*u1 + sina*u2;
		U2 = -sina*u1 + cosa*u2;
		
		fun = U1*U1/aa + U2*U2/bb;
		
		if(fun<=1)
		  {
		    z  = data[k].density;		
		    z += raft_phantom_get(P->vector, j, size-1-i); 
		    raft_phantom_set(P->vector, j, size-1-i, z);
		  }	      
	      }
	  }
      }
      
    }
  else
    {
      double table_head[10][6] = {
	{0.0,  0.0, 6.9e-001, 9.2e-001, 0.0, 1.0},
	{0.0, -1.84e-002, 6.624e-001, 8.74e-001, 0.0,-9.80e-001},
	{2.2e-001, 0.0, 1.1e-001, 3.1e-001,-3.1415927e-001,-8.0e-001}, 
	{-2.2e-001, 0.0, 1.6e-001, 4.1e-001, 3.1415927e-001,-8.0e-001},
	{0.0,  3.5e-001, 2.1e-001, 2.5e-001, 0.0, 1.0e-002},
	{0.0,  1.0e-001, 4.6e-002, 4.6e-002, 0.0, 2.0e-001},   
	{0.0, -1.0e-001, 4.6e-002, 4.6e-002, 0.0, 1.0e-002},
	{-8.0e-002,-6.05e-001, 4.60e-002, 2.3e-002, 0.0, 1.0e-002},
	{0.00, -6.06e-001, 2.3e-002, 2.3e-002, 0.0, 1.0e-002},
	{6.00e-002, -6.05e-001, 2.30e-002, 4.60e-002, 0.0, 1.0e-002}
      };

       struct{
	 double c1, c2, A, B, theta, density;
       }data[10];   
       
       size = raft_phantom_get_size(P);
       scale = raft_phantom_get_max(P);
       
       for(k=0;k<10;k++){
	 
	 data[k].c1 = table_head[k][0];
	 data[k].c2 = table_head[k][1];
	 data[k].A  = table_head[k][2];
	 data[k].B  = table_head[k][3];
	 data[k].theta = table_head[k][4];
	 data[k].density = table_head[k][5];
	 
	 aa = data[k].A * data[k].A;
	 bb = data[k].B * data[k].B;
	 
	 for(j=size-1; j>=0; j--)
	   {
	     y = raft_phantom_get_y(P, j);
	     
	     y=y/scale; 
	     
	     for(i=0; i < size; i++)
	       {
		 double U1, U2; 
		 
		 x = raft_phantom_get_x(P, i);
		 
		 x=x/scale;
		 
		 u1 = x + data[k].c1;
		 u2 = y + data[k].c2;
		 
		 cosa = cos(data[k].theta);
		 sina = sin(data[k].theta);
		 U1 =  cosa*u1 + sina*u2;
		 U2 = -sina*u1 + cosa*u2;
		 
		 fun = U1*U1/aa + U2*U2/bb;
		 
		 if(fun<=1)
		   {
		     z  = 20 + 30 * data[k].density;		
		     z += raft_phantom_get(P->vector, j, size-1-i); 
		     raft_phantom_set(P->vector, j, size-1-i, z/2);
		   }	      
	       }
	   }
       }
       
    }
}

/*+====================================================+
  
  FUNCTION: raft_phantom_define_xfct_density_miqdep

  Define a density phantom for tomographic fluorescence 
  simulations. Miqueles & De Pierro (2009)
  
  Input: 

  P - RAFT phantom data type. See <raft_phantom_t>.  
  
  +====================================================+
*/

void raft_phantom_define_xfct_density_miqdep(raft_phantom_t *P)
{
  int i,j, size, leafs;
  double scale, theta, r, f, x, y, att, radius2;
  double apert1, apert2, trans, X, Y;
  
  leafs = 2;
  apert1 = 0.14;
  apert2 = 0.16;
  trans = 0.22;
  radius2 = 0.16;
  
  size  = raft_phantom_get_size(P);
  scale = raft_phantom_get_max(P);
  
  att = 10e-06;

  for(i=0; i < size; i++)
    {
      x = raft_phantom_get_x(P, i);
      x = x/scale;
      
      for(j=0; j < size; j++)
	{
	  y = raft_phantom_get_y(P, j);
	  y = y/scale;  
	  
	  /* Boundary for another XFCT phantom density */
	  
	  if(DEFINE_BOUNDARY)
	    {
	      r = gsl_hypot(x,y);
	      if( (r<SQR(.82)) && (r>SQR(.80)))
		raft_phantom_set(P->vector, i, j, 10e-05);
	      
	      if( r<SQR(.80) )
		raft_phantom_set(P->vector, i, j, 1.5*att);
	    }
	  
	  /* 1st rose */
	  
	  X = x+trans;
	  Y = y+trans;

	  if(fabs(x) < EPS)
	    theta = 0;
	  else
	    theta = atan(Y/X);
	  
	  r = gsl_hypot(X,Y);
	  
	  f = apert1 * cos(leafs * theta);
	  
	  if( MAX(fabs(X),fabs(Y)) < radius2)
	    raft_phantom_set(P->vector, i, j, 4.8*att);
	  
	  if(fabs(r) < fabs(f))
	    raft_phantom_set(P->vector, i, j, 3.0*att);
	  
	  /* 2nd rose */
	  
	  Y = y-trans;
	  X = x-trans;

	  if(fabs(x) < EPS)
	    theta = 0;
	  else
	    theta = atan(Y/X);

	  r = gsl_hypot(Y,X);
	  
	  f = apert1 * cos((leafs+1) * theta);
	  
	  if( MAX(fabs(X),fabs(Y)) < radius2)
	    raft_phantom_set(P->vector, i, j, 4.8*att);

	  if(fabs(r) < fabs(f))
	    raft_phantom_set(P->vector, i, j, 3.0*att);

	  /* 3rd rose */

	  X = x-trans;
	  Y = y+trans;

	  if(fabs(x) < EPS)
	    theta = 0;
	  else
	    theta = atan(Y/X);

	  r = gsl_hypot(Y,X);
	  
	  f = apert1 * cos(leafs * (theta + MPI/4));
	  
	  if( r < radius2)
	    raft_phantom_set(P->vector, i, j, 4.8*att);

	  if(fabs(r) < fabs(f))
	    raft_phantom_set(P->vector, i, j, 10*att);
	  
	  /* 4th rose */

	  X = x+trans;
	  Y = y-trans;

	  if(fabs(x) < EPS)
	    theta = 0;
	  else
	    theta = atan(Y/X);

	  r = gsl_hypot(Y,X);
	  
	  f = apert1 * cos((leafs+1) * (theta + MPI/4));
	  
	  if( r < radius2)
	    raft_phantom_set(P->vector, i, j, 4.8*att);

	  if(fabs(r) < fabs(f))
	    raft_phantom_set(P->vector, i, j, 10*att);
	  
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_phantom_define_xfct_Tattenuation_miqdep

  Define the transmission attenuation phantom for 
  tomographic fluorescence simulations. Miqueles & De Pierro
  (2009).
  
  Input: 

  P - RAFT phantom data type. See <raft_phantom_t>.  
  
  +====================================================+
*/

void raft_phantom_define_xfct_Tattenuation_miqdep(raft_phantom_t *P)
{
  int i,j, size, leafs;
  double scale, theta, r, f, x, y, att, radius, radius2;
  double apert, trans,X,Y;
  
  leafs = 2;
  apert = 0.14;
  radius = 0.5;
  trans = 0.22;
  radius2 = 0.16;
  
  size  = raft_phantom_get_size(P);
  scale = raft_phantom_get_max(P);
  
  att = 0.15;

  for(i=0; i < size; i++)
    {
      x = raft_phantom_get_x(P, i);
      x = x/scale;
      
      for(j=0; j < size; j++)
	{
	  y = raft_phantom_get_y(P, j);
	  y = y/scale;

	  if(!DEFINE_BOUNDARY)
	    {
	      /* outer square */

	      if( (MAX( fabs(x), fabs(y)) < radius) &&
		  (MAX( fabs(x), fabs(y)) > radius*(1-0.03)))
		raft_phantom_set(P->vector, i, j, 4.8*att);

	      if( MAX( fabs(x), fabs(y)) < radius*(1-0.03) )
		raft_phantom_set(P->vector, i, j, 1.5*att);
	      
	      /* 1st rose */
	      
	      X = x+trans;
	      Y = y+trans;
	      
	      if(fabs(x) < EPS)
		theta = 0;
	      else
		theta = atan(Y/X);
	      
	      r = gsl_hypot(X,Y);
	      
	      f = apert * cos(leafs * theta);
	      
	      if( MAX(fabs(X),fabs(Y)) < radius2)
		raft_phantom_set(P->vector, i, j, 4*att);
	      
	      if(fabs(r) < fabs(f))
		raft_phantom_set(P->vector, i, j, 0);
	      
	      /* 2nd rose */
	      
	      Y = y-trans;
	      X = x-trans;
	      
	      if(fabs(x) < EPS)
		theta = 0;
	      else
		theta = atan(Y/X);
	      
	      r = gsl_hypot(Y,X);
	      
	      f = apert * cos((leafs+1) * theta);
	      
	      if( MAX(fabs(X),fabs(Y)) < radius2)
		raft_phantom_set(P->vector, i, j, 4*att);
	      
	      if(fabs(r) < fabs(f))
		raft_phantom_set(P->vector, i, j, 0);
	      
	      /* 3rd rose */
	      
	      X = x-trans;
	      Y = y+trans;
	      
	      if(fabs(x) < EPS)
		theta = 0;
	      else
		theta = atan(Y/X);
	      
	      r = gsl_hypot(Y,X);
	      
	      f = apert * cos(leafs * (theta + MPI/4));
	      
	      if( r < radius2)
		raft_phantom_set(P->vector, i, j, 4*att);
	      
	      if(fabs(r) < fabs(f))
		raft_phantom_set(P->vector, i, j, 0);
	      
	      /* 4th rose */
	      
	      X = x+trans;
	      Y = y-trans;
	      
	      if(fabs(x) < EPS)
		theta = 0;
	      else
		theta = atan(Y/X);
	      
	      r = gsl_hypot(Y,X);
	      
	      f = apert * cos((leafs+1) * (theta + MPI/4));
	      
	      if( r < radius2)
		raft_phantom_set(P->vector, i, j, 4*att);
	      
	      if(fabs(r) < fabs(f))
		raft_phantom_set(P->vector, i, j, 0);

	    }
	  else
	    {
	      r = gsl_hypot(x,y);
	      
	      if( r<SQR(.80) )
		raft_phantom_set(P->vector, i, j, 1);	      
	    }
	  

	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_phantom_define_xfct_attenuation_miqdep

  Define the fluorescence attenuation phantom for 
  tomographic fluorescence simulations. Miqueles &
  De Pierro (2009).
  
  Input: 

  P - RAFT phantom data type. See <raft_phantom_t>.  
  
  +====================================================+
*/

void raft_phantom_define_xfct_attenuation_miqdep(raft_phantom_t *P)
{
  int i,j, size, leafs;
  double scale, theta, r, f, x, y, att, radius, radius2;
  double apert, trans,X,Y;
  
  leafs = 2;
  apert = 0.14;
  radius = 0.5;
  trans = 0.22;
  radius2 = 0.16;
  
  size  = raft_phantom_get_size(P);
  scale = raft_phantom_get_max(P);
  
  att = 0.3;

  for(i=0; i < size; i++)
    {
      x = raft_phantom_get_x(P, i);
      x = x/scale;
      
      for(j=0; j < size; j++)
	{
	  y = raft_phantom_get_y(P, j);
	  y = y/scale;
	  
	  if(!DEFINE_BOUNDARY)
	    {
	      /* outer square */
	      
	      if( (MAX( fabs(x), fabs(y)) < radius) &&
		  (MAX( fabs(x), fabs(y)) > radius*(1-0.03)))
		raft_phantom_set(P->vector, i, j, 4.8*att);

	      if( MAX( fabs(x), fabs(y)) < radius*(1-0.03) )
		raft_phantom_set(P->vector, i, j, 1.5*att);
	      
	      /* 1st rose */
	      
	      X = x+trans;
	      Y = y+trans;
	      
	      if(fabs(x) < EPS)
		theta = 0;
	      else
		theta = atan(Y/X);
	      
	      r = gsl_hypot(X,Y);
	      
	      f = apert * cos(leafs * theta);
	      
	      if( MAX(fabs(X),fabs(Y)) < radius2)
		raft_phantom_set(P->vector, i, j, 4*att);
	      
	      if(fabs(r) < fabs(f))
		raft_phantom_set(P->vector, i, j, 3.0*att);
	      
	      /* 2nd rose */
	      
	      Y = y-trans;
	      X = x-trans;
	      
	      if(fabs(x) < EPS)
		theta = 0;
	      else
		theta = atan(Y/X);
	      
	      r = gsl_hypot(Y,X);
	      
	      f = apert * cos((leafs+1) * theta);
	      
	      if( MAX(fabs(X),fabs(Y)) < radius2)
		raft_phantom_set(P->vector, i, j, 4*att);
	      
	      if(fabs(r) < fabs(f))
		raft_phantom_set(P->vector, i, j, 3.0*att);
	      
	      /* 3rd rose */
	      
	      X = x-trans;
	      Y = y+trans;
	      
	      if(fabs(x) < EPS)
		theta = 0;
	      else
		theta = atan(Y/X);
	      
	      r = gsl_hypot(Y,X);
	      
	      f = apert * cos(leafs * (theta + MPI/4));
	      
	      if( r < radius2)
		raft_phantom_set(P->vector, i, j, 4*att);
	      
	      if(fabs(r) < fabs(f))
		raft_phantom_set(P->vector, i, j, att/2);
	      
	      /* 4th rose */
	      
	      X = x+trans;
	      Y = y-trans;
	      
	      if(fabs(x) < EPS)
		theta = 0;
	      else
		theta = atan(Y/X);
	      
	      r = gsl_hypot(Y,X);
	      
	      f = apert * cos((leafs+1) * (theta + MPI/4));
	      
	      if( r < radius2)
		raft_phantom_set(P->vector, i, j, 4*att);
	      
	      if(fabs(r) < fabs(f))
		raft_phantom_set(P->vector, i, j, att/2);
	    }
	  else
	    {
	      r = gsl_hypot(x,y);
	      
	      if( r<SQR(.80) )
		raft_phantom_set(P->vector, i, j, 1);
	    }
	}
    }
}


/*+====================================================+
  
  FUNCTION: raft_phantom_define_xfct_density_golosio

  Define a density phantom for tomographic fluorescence 
  simulations. Golosio et al (2003), density distribution
  for Calcium.
  
  Input: 

  P - RAFT phantom data type. See <raft_phantom_t>.  
  
  +====================================================+
*/

void raft_phantom_define_xfct_density_golosio(raft_phantom_t *P)
{
  int i,j, size;
  double scale, r, x, y, att, radius, X, Y, trans;
  
  att    = 1.5;
  radius = 0.16;
  trans  = 0.22;
  
  size  = raft_phantom_get_size(P);
  scale = raft_phantom_get_max(P);
  
  for(i=0; i < size; i++)
    {
      x = raft_phantom_get_x(P, i);
      x = x/scale;
      
      for(j=0; j < size; j++)
	{
	  y = raft_phantom_get_y(P, j);
	  y = y/scale;  
	  
	  /* 1st structure (up,left) */
	  
	  X = x+trans;
	  Y = y+trans;

	  if( MAX(fabs(X),fabs(Y)) < radius)
	    raft_phantom_set(P->vector, i, j, att);
	  
	  /* 2nd structure (down,right) */
	  
	  Y = y-trans;
	  X = x-trans;

	  r = gsl_hypot(X,Y);

	  if( MAX(fabs(X),fabs(Y)) < radius)
	    raft_phantom_set(P->vector, i, j, att);

	  if(fabs(r) < 0.7 * radius)
	    raft_phantom_set(P->vector, i, j, 0.7 * att);
	  
	  /* 3rd structure (down,left) */
	  
	  X = x-trans;
	  Y = y+trans;	 

	  r = gsl_hypot(Y,X);
	  
	  if( r < radius)
	    raft_phantom_set(P->vector, i, j, att);

	  if( fabs(X) + fabs(Y) < 0.9 * radius)
	    raft_phantom_set(P->vector, i, j, 0);
	  
	  /* 4th structure (up,right) */

	  X = x+trans;
	  Y = y-trans;

	  r = gsl_hypot(Y,X);
      
	  if( r < radius)
	    raft_phantom_set(P->vector, i, j, att);
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_phantom_define_xfct_Tattenuation_golosio

  Define the transmission attenuation phantom for 
  tomographic fluorescence simulations. Golosio et al
  (2003), attenuation for Calcium.
  
  Input: 

  P - RAFT phantom data type. See <raft_phantom_t>.  
  
  +====================================================+
*/

void raft_phantom_define_xfct_Tattenuation_golosio(raft_phantom_t *P)
{
  int i,j, size;
  double scale, r, x, y, att, att2, radius, X, Y, trans, square;
  
  att    = 0.2;
  att2   = 10*att;
  radius = 0.16;
  square = 0.45;
  trans  = 0.22;
  
  size  = raft_phantom_get_size(P);
  scale = raft_phantom_get_max(P);
  
  for(i=0; i < size; i++)
    {
      x = raft_phantom_get_x(P, i);
      x = x/scale;
      
      for(j=0; j < size; j++)
	{
	  y = raft_phantom_get_y(P, j);
	  y = y/scale;  
	  
	  /* outer square */

	  if( (MAX( fabs(x), fabs(y)) < square) &&
	      (MAX( fabs(x), fabs(y)) > square*(1-0.03)))
	    raft_phantom_set(P->vector, i, j, att);
	  
	  if( MAX( fabs(x), fabs(y)) < square*(1-0.03) )
	    raft_phantom_set(P->vector, i, j, att);

	  /* 1st structure (up,left) */
	  
	  X = x+trans;
	  Y = y+trans;

	  if( MAX(fabs(X),fabs(Y)) < radius)
	    raft_phantom_set(P->vector, i, j, att2);
	  
	  /* 2nd structure (down,right) */
	  
	  Y = y-trans;
	  X = x-trans;

	  r = gsl_hypot(X,Y);

	  if( MAX(fabs(X),fabs(Y)) < radius)
	    raft_phantom_set(P->vector, i, j, att2);
	  
	  /* 3rd structure (down,left) */
	  
	  X = x-trans;
	  Y = y+trans;	 

	  r = gsl_hypot(Y,X);
	  
	  if( r < radius)
	    raft_phantom_set(P->vector, i, j, att2);

	  if( fabs(X) + fabs(Y) < 0.9 * radius)
	    raft_phantom_set(P->vector, i, j, att);

	  /* 4th structure (up,right) */

	  X = x+trans;
	  Y = y-trans;

	  r = gsl_hypot(Y,X);
      
	  if( r < radius)
	    raft_phantom_set(P->vector, i, j, att2);
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_phantom_define_xfct_attenuation_golosio

  Define the fluorescence attenuation phantom for 
  tomographic fluorescence simulations. Golosio et al
  (2003), attenuation for Calcium.
  
  Input: 

  P - RAFT phantom data type. See <raft_phantom_t>.  
  
  +====================================================+
*/

void raft_phantom_define_xfct_attenuation_golosio(raft_phantom_t *P)
{
  int i,j, size;
  double scale, r, x, y, att, att2, radius, X, Y, trans, square;
  
  att    = 0.5;
  att2   = 2*att;
  radius = 0.16;
  square = 0.45;
  trans  = 0.22;
  
  size  = raft_phantom_get_size(P);
  scale = raft_phantom_get_max(P);
  
  for(i=0; i < size; i++)
    {
      x = raft_phantom_get_x(P, i);
      x = x/scale;
      
      for(j=0; j < size; j++)
	{
	  y = raft_phantom_get_y(P, j);
	  y = y/scale;  
	  
	  /* outer square */

	  if( (MAX( fabs(x), fabs(y)) < square) &&
	      (MAX( fabs(x), fabs(y)) > square*(1-0.03)))
	    raft_phantom_set(P->vector, i, j, att);
	  
	  if( MAX( fabs(x), fabs(y)) < square*(1-0.03) )
	    raft_phantom_set(P->vector, i, j, att);

	  /* 1st structure (up,left) */
	  
	  X = x+trans;
	  Y = y+trans;

	  if( MAX(fabs(X),fabs(Y)) < radius)
	    raft_phantom_set(P->vector, i, j, att2);
	  
	  /* 2nd structure (down,right) */
	  
	  Y = y-trans;
	  X = x-trans;

	  r = gsl_hypot(X,Y);

	  if( MAX(fabs(X),fabs(Y)) < radius)
	    raft_phantom_set(P->vector, i, j, att2);

	  if(fabs(r) < 0.7 * radius)
	    raft_phantom_set(P->vector, i, j, 1.1 * att2);
	  
	  /* 3rd structure (down,left) */
	  
	  X = x-trans;
	  Y = y+trans;	 

	  r = gsl_hypot(Y,X);
	  
	  if( r < radius)
	    raft_phantom_set(P->vector, i, j, att2);

	  if( fabs(X) + fabs(Y) < 0.9 * radius)
	    raft_phantom_set(P->vector, i, j, att);
	  
	  /* 4th structure (up,right) */

	  X = x+trans;
	  Y = y-trans;

	  r = gsl_hypot(Y,X);
      
	  if( r < radius)
	    raft_phantom_set(P->vector, i, j, att2);
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_phantom_define_spect_attenuation

  Define the attenuation phantom for SPECT simulations.
  
  It is similar to the one proposed by F.Natterer at "Inversion of
  the attenuated Radon transform" (Inverse Problems)
  
  Input: 

  P - RAFT phantom data type. See <raft_phantom_t>.  
  
  +====================================================+
*/

void raft_phantom_define_spect_attenuation(raft_phantom_t *P)
{
  int i,j, size;
  double scale, r, x, y, att;
    
  size  = raft_phantom_get_size(P);
  scale = raft_phantom_get_max(P);
  
  att = 0.1;

  for(i=0; i < size; i++)
    {
      x = raft_phantom_get_x(P, i);
      x = x/scale;
      
      for(j=0; j < size; j++)
	{
	  y = raft_phantom_get_y(P, j);
	  y = y/scale;
		 
	  r = gsl_hypot(x,y);
	  
	  if(fabs(r) < 0.7)
	    raft_phantom_set(P->vector, i, j, att);
	  
	  r = gsl_hypot(x,y-0.3);

	  if(fabs(r) < 0.2)
	    raft_phantom_set(P->vector, i, j, 0);
	  
	  r = gsl_hypot(x,2*(y+0.3));

	  if(fabs(r) < 0.35)
	    raft_phantom_set(P->vector, i, j, 0);
	}
    }
}


/*+====================================================+
  
  FUNCTION: raft_phantom_define_spect_activity

  Define the activity phantom for SPECT simulations.

  It is similar to the one proposed by F.Natterer at 
  "Inversion of the attenuated Radon transform" (Inverse Problems)
  
  Input: 

  P - RAFT phantom data type. See <raft_phantom_t>.  
  
  +====================================================+
*/

void raft_phantom_define_spect_activity(raft_phantom_t *P)
{
  int i,j, size;
  double scale, r, x, y, act;
    
  size  = raft_phantom_get_size(P);
  scale = raft_phantom_get_max(P);
  
  act = 100;

  for(i=0; i < size; i++)
    {
      x = raft_phantom_get_x(P, i);
      x = x/scale;
      
      for(j=0; j < size; j++)
	{
	  y = raft_phantom_get_y(P, j);
	  y = y/scale;
		 
	  r = gsl_hypot(x,y);
	  
	  if(fabs(r) < 0.7)
	    raft_phantom_set(P->vector, i, j, act);
	  
	  r = gsl_hypot(x,y-0.3);

	  if(fabs(r) < 0.2)
	    raft_phantom_set(P->vector, i, j, 0);
	  
	  r = gsl_hypot(x,2*(y+0.3));

	  if(fabs(r) < 0.35)
	    raft_phantom_set(P->vector, i, j, 0);

	  /* activity */

	  r = gsl_hypot(x-0.3,y-0.1);

	  if(fabs(r) < 0.1)
	    raft_phantom_set(P->vector, i, j, 10*act);
	}
    }
}
