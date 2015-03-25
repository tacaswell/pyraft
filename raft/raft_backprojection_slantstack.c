#include "raft_backprojection_slantstack.h"

#include "raft_image.h"
#include <pthread.h>
#include <math.h>

#define MAX_NTHREADS 20000
#define MIN( a, b ) ( ( ( a ) < ( b ) ) ? ( a ) : ( b ) )
#define MAX( x, y ) ( ( ( x ) > ( y ) ) ? ( x ) : ( y ) )
#define SIGN( x ) ( ( ( x ) > 0.0 ) ? 1.0 : ( ( ( x ) < 0.0 ) ? -1.0 : 0.0 ) )
#define PI 3.1415926535897932384626433832795

typedef struct{

  raft_image sinogram;
  raft_image backprojection;
  int size;
  int nthread;
  int rowIndex[2];
  int colIndex[2];  
  double phantom_tl[2];
  double phantom_br[2];
  
}param_t;

void *back_slantstack_loop(void *t);

double eval_anglesum_slantstack(raft_image sinogram,
				double phantom_tl[2],
				double phantom_br[2],
				int size,
				int j,
				int k);

/*==============================================================*/

/*!
 * \brief Backprojection transform: slant stack
 * \param sinogram raft image
 * \param nthreads number of threads 
 */

void raft_backprojection_slantstack(raft_image sinogram,
				    raft_image backprojection,
				    int nthreads)
{
  pthread_t thread[MAX_NTHREADS];
  pthread_attr_t attr;
  int size, e, n, rc;    
  param_t param[MAX_NTHREADS];  
  void *status;
  
  // Initialize and set thread detached attribute
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  size = raft_matrix_nlines(backprojection.data);
  
  e = (int) floor((size*size)/nthreads);;
      
  for(n = 0; n < nthreads+1; n++) 
    {
      param[n].size           = size;
      param[n].sinogram       = sinogram;
      param[n].backprojection = backprojection; 
      param[n].nthread        = n;      
      
      param[n].colIndex[0]   = e * n;
      param[n].colIndex[1]   = (n+1) * e;
      param[n].phantom_tl[0] = backprojection.tl_x;
      param[n].phantom_tl[1] = backprojection.tl_y;
      param[n].phantom_br[0] = backprojection.br_x;
      param[n].phantom_br[1] = backprojection.br_y;
      
      rc = pthread_create(&thread[n], 
			  &attr, 
			  back_slantstack_loop, 
			  (void *)&param[n]);
    }

  // Free attribute and wait for the other threads
  pthread_attr_destroy(&attr);
  for(n = 0; n < nthreads + 1; n++) 
    {
      rc = pthread_join(thread[n], &status);
    }  
}

/*==============================================================*/

void *back_slantstack_loop(void *t)
{
  int w, j, k, size, npixels;
  param_t *param;
  double sum;
  
  param = (param_t *)t;
  
  size = param->size;
  npixels = size*size;
  
  for(w = param->colIndex[0]; w < MIN(param->colIndex[1], npixels) ; w++)
    {
      j = w/size;
      k = w%size;
      
      sum = eval_anglesum_slantstack(param->sinogram, 
				     param->phantom_tl,
				     param->phantom_br,
				     size, 
				     j, k);	  
      
      raft_matrix_element(param->backprojection.data, j, k) = sum;
    }  
  
  pthread_exit(NULL);
}

/*=========== PRIVATE FUNCTIONS ========================*/

double eval_anglesum_slantstack(raft_image sinogram,
				double phantom_tl[2],
				double phantom_br[2],
				int size,
				int j,
				int k)
{
  int i, m, n, nrays, nviews;
  double cost, sint, sum, s;
  double th0, thf, th, dth, dt;
  double t, tmin, tmax;
  double x, xmin, xmax, dx;
  double y, ymin, ymax, dy;
  
  nrays  = raft_matrix_nlines(sinogram.data);
  nviews = raft_matrix_ncolumns(sinogram.data);

  tmin = sinogram.br_y;
  tmax = sinogram.tl_y;
  dt = (tmax-tmin)/(nrays-1);

  th0 = sinogram.tl_x;
  thf = sinogram.br_x;
  
  dth = (thf-th0)/(nviews);

  xmin = phantom_tl[0];
  xmax = phantom_br[0];
  dx   = (xmax - xmin)/(size - 1);

  ymin = phantom_br[1];
  ymax = phantom_tl[1];
  dy   = (ymax - ymin)/(size - 1);
  
  x = xmin + k * dx;

  y = ymin + j * dy;
    
  // stacking 

  sum = 0;

  for(i=0; i < nviews; i++)
    {
      th = th0 + i*dth;
	  
      cost  = cos(th);
      sint  = sin(th);
      
      t =  x*cost + y*sint;
      
      n = (int) floor((t-tmin)/dt);
      
      if(n > -1 && n<nrays-1)
	{
	  sum += (raft_matrix_element(sinogram.data,n+1,i)-raft_matrix_element(sinogram.data,n,i))*(t-(tmin+n*dt))/dt + raft_matrix_element(sinogram.data,n,i);
	}
    }
  
  return sum*dth;
}

