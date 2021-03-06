#include "raft_radon_slantstack.h"
#include "raft_vector.h"
#include <pthread.h>
#include <math.h>

#define MAX_NTHREADS 20000

#define MIN( a, b ) ( ( ( a ) < ( b ) ) ? ( a ) : ( b ) )
#define MAX( x, y ) ( ( ( x ) > ( y ) ) ? ( x ) : ( y ) )
#define SIGN( x ) ( ( ( x ) > 0.0 ) ? 1.0 : ( ( ( x ) < 0.0 ) ? -1.0 : 0.0 ) )
#define PI 3.1415926535897932384626433832795

typedef struct{

  raft_image phantom;
  raft_image radon;
  int size, nviews, nrays;
  int nthread;
  int rowIndex[2];
  int colIndex[2];  
  double sinogram_tl[2];
  double sinogram_br[2];
  
}param_t;

void *slantstack_loop(void *t);

double eval_rayintegral(raft_image phantom,
			double sinogram_tl[2],
			double sinogram_br[2],
			int nviews,
			int nrays,
			int j,
			int k);

/*========================================================*/

/*!
 * \brief Radon transform: slant stack
 * \param phantom raft image
 * \param nthreads number of threads
 */

void raft_radon_slantstack(raft_image phantom,
			   raft_image radon,
			   int nthreads)
{
  pthread_t thread[MAX_NTHREADS];
  pthread_attr_t attr;
  int size, nviews, nrays, e, n, rc;    
  param_t param[MAX_NTHREADS];  
  void *status;
  
  // Initialize and set thread detached attribute
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  size   = raft_matrix_nlines(phantom.data);
  nviews = raft_matrix_ncolumns(radon.data);
  nrays  = raft_matrix_nlines(radon.data);

  e = (int) floor((nviews*nrays)/nthreads);
      
  for(n = 0; n < nthreads + 1; n++) 
    {
      param[n].phantom = phantom;
      param[n].radon   = radon;
      param[n].size    = size; 
      param[n].nviews  = nviews;
      param[n].nrays   = nrays;
      param[n].nthread = n;      
      param[n].sinogram_tl[0] = radon.tl_x;
      param[n].sinogram_tl[1] = radon.tl_y;
      param[n].sinogram_br[0] = radon.br_x;
      param[n].sinogram_br[1] = radon.br_y;      
      
      param[n].colIndex[0] = e * n;
      param[n].colIndex[1] = (n+1) * e;
      
      rc = pthread_create(&thread[n], &attr, slantstack_loop, (void *)&param[n]);
    }
  
  // Free attribute and wait for the other threads
  pthread_attr_destroy(&attr);
  for(n = 0; n < nthreads + 1; n++) 
    {
      rc = pthread_join(thread[n], &status);
    }  
}

/*=========== PRIVATE FUNCTIONS ========================*/

void *slantstack_loop(void *t)
{
  int w, j, k, ndata;
  param_t *param;
  double sum;
  
  param = (param_t *)t;
  
  ndata = (param->nrays)*(param->nviews);

  for(w = param->colIndex[0]; w < MIN( ndata, param->colIndex[1]); w++)
    {
      j = w/param->nrays;
      k = w%param->nrays;

      sum = eval_rayintegral(param->phantom, 
			     param->sinogram_tl,
			     param->sinogram_br,
			     param->nviews, 
			     param->nrays, j, k);
      
      raft_matrix_element(param->radon.data, k, j) = sum;
    }  
  
  pthread_exit(NULL);
}

/*=======================================================*/

double eval_rayintegral(raft_image phantom,
			double sinogram_tl[2],
			double sinogram_br[2],
			int nviews,
			int nrays,
			int j,
			int k)
{
  int i, size, m;
  double sum, cost, sint, tant;
  double x, xmin, xmax, dx;
  double y, ymin, ymax, dy;
  double theta, th0, thf, dth;
  double t, tmin, tmax, dt; 
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
  th0 = sinogram_tl[0];
  thf = sinogram_br[0];
  dth = (thf-th0)/(nviews);

  theta = th0 + j * dth;

  cost  = cos(theta);
  sint  = sin(theta);
  tant  = tan(theta);  
  fsint = fabs(sint);
  fcost = fabs(cost);
  
  //ray @ angle
  
  tmin = sinogram_br[1];
  tmax = sinogram_tl[1];
  dt = (tmax - tmin)/(nrays - 1);

  t = tmin + k * dt;
   
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


