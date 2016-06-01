#include "raft_backprojection.h"
#include "raft_projection.h"
#include "raft_weight.h"
#include "raft_math.h"

#include "backprojection.h"
#include "projection.h"
#include "filter.h"
#include <pthread.h>

#define BACKP_LBLOCK 0
#define BACKP_MBLOCK 1
#define BACKP_TTYPE  BACKP_LBLOCK
#define BACKP_THREADS 7

/*######################################################
  Title: Backprojection

  Header - <raft/raft_backprojection.h>
  Type - <raft_backp_t>
  $Id: backprojection.c,v 1.45 2011-03-02 19:23:19 miqueles Exp $ - Last Update  
  ####################################################*/

/*######################################################
                         Private
  ####################################################*/

/*+====================================================+
  
  FUNCTION  eval_anglesum

  Compute the angle sum for the backprojection formula.
  
  Input  

  data - scan data . See <raft_scan_t>
  projection - projection vector. 
  j - pixel index (direction y).
  k - pixel index (direction x).   
    
  Return 

  Value for the reconstructed pixel.
    
  +====================================================+
*/

double eval_anglesum(raft_scan_t *data,
		     gsl_vector *projection,
		     int j,
		     int k,
		     double sum[4])
{
  if(data->viewsid == RAFT_REG)
    {
      int i, row[4], nrays, size;
      double x, y, da, cost, sint, t[4];
      gsl_vector_view P;
      
      nrays = raft_scan_get_nrays(data);
      size  = raft_scan_get_size(data);

      da = raft_scan_get_angle(data,1)-raft_scan_get_angle(data,0);
      
      if(data->scan != RAFT_CT)
	da = da/2;
      
      x  = raft_scan_get_x(data,k);
      y  = raft_scan_get_y(data,size-1-j);
      
      sum[0]  = 0; 
      sum[1]  = 0; 
      sum[2]  = 0; 
      sum[3]  = 0; 
      
      for(i=0; i < raft_scan_get_nviews(data); i++)
	{
	  /* ith row of the projection matrix */
	  P = gsl_vector_subvector(projection, i*nrays, nrays);  
	  
	  cost  = gsl_vector_get(data->costheta, i);
	  sint  = gsl_vector_get(data->sintheta, i);
	  	  
	  t[0] =  x*cost + y*sint;
	  t[1] = -x*cost - y*sint;
	  t[2] =  x*cost - y*sint;
	  t[3] = -x*cost + y*sint;
	  
	  row[0] = interp_nearest(data->t, t[0], data->raysid);
	  row[1] = interp_nearest(data->t, t[1], data->raysid);
	  row[2] = interp_nearest(data->t, t[2], data->raysid);
	  row[3] = interp_nearest(data->t, t[3], data->raysid);
	  
	  if(SIGN(row[0])>0)
	    {
	      sum[0] += interp_linear(&P.vector, data->t, row[0], t[0]);
	      sum[1] += interp_linear(&P.vector, data->t, row[1], t[1]);
	      sum[2] += interp_linear(&P.vector, data->t, row[2], t[2]);
	      sum[3] += interp_linear(&P.vector, data->t, row[3], t[3]);	      
	    }
	}
      
      //return MAX(0,(sum[0]*da));
      return sum[0]*da;
    }
  else
    {
      int i, row[4], nrays, size;
      double x, y, da, cost, sint, t[4];
      gsl_vector_view P;
      
      nrays = raft_scan_get_nrays(data);
      size  = raft_scan_get_size(data);

      x = raft_scan_get_x(data,k);
      y = raft_scan_get_y(data,size-1-j);
      
      sum[0]  = 0; 
      sum[1]  = 0; 
      sum[2]  = 0; 
      sum[3]  = 0;        
      
      for(i=0; i < raft_scan_get_nviews(data) - 1; i++)
	{
	  /* ith row of the projection matrix */
	  P = gsl_vector_subvector(projection, i*nrays, nrays);  
	  
	  da = raft_scan_get_angle(data,i+1) - raft_scan_get_angle(data,i);
	  
	  cost  = gsl_vector_get(data->costheta, i);
	  sint  = gsl_vector_get(data->sintheta, i);

	  t[0] =  x*cost + y*sint;
	  t[1] = -x*cost - y*sint;
	  t[2] =  x*cost - y*sint;
	  t[3] = -x*cost + y*sint;
	  
	  row[0] = interp_nearest(data->t, t[0], data->raysid);
	  row[1] = interp_nearest(data->t, t[1], data->raysid);
	  row[2] = interp_nearest(data->t, t[2], data->raysid);
	  row[3] = interp_nearest(data->t, t[3], data->raysid);
	  
	  if(SIGN(row[0])>0)	
	    {
	      sum[0] += interp_linear(&P.vector, data->t, row[0], t[0]) * da;
	      sum[1] += interp_linear(&P.vector, data->t, row[1], t[1]) * da;
	      sum[2] += interp_linear(&P.vector, data->t, row[2], t[2]) * da;
	      sum[3] += interp_linear(&P.vector, data->t, row[3], t[3]) * da;
	    }
	}
      
      if(data->scan != RAFT_CT)
	{
	  sum[0] = sum[0]/2;
	  sum[1] = sum[1]/2;
	  sum[2] = sum[2]/2;
	  sum[3] = sum[3]/2;
	}
      
      /* return MAX(0,sum); */
      return sum[0];
    }
}


/*+====================================================+
  
  FUNCTION  eval_anglesum_partial

  Compute the angle sum for the partial backprojection 
  formula.
  
  Input  

  data - scan data . See <raft_scan_t>.
  projection - projection vector. 
  T - vector with angle indexes.
  sT - size of vector T.
  j - pixel index (direction y).
  k - pixel index (direction x).   
    
  Return 

  Value for the reconstructed pixel.
    
  +====================================================+
*/

double eval_anglesum_partial(raft_scan_t *data,
			     gsl_vector *projection,
			     int j,
			     int k,
			     gsl_vector *T,
			     int sT)
{
  int l, i, row, nrays, size;
  double x, y, s, cost, sint, t;
  gsl_vector_view P;
  
  nrays = raft_scan_get_nrays(data);
  size  = raft_scan_get_size(data);
  
  x = raft_scan_get_x(data,k);
  y = raft_scan_get_y(data,size-1-j);
  s = 0; 
  
  for(l=0; l < sT; l++)
    {
      i = (int) gsl_vector_get(T,l);

      /* ith row of the projection matrix */
      P = gsl_vector_subvector(projection, i*nrays, nrays); 
      
      cost  = gsl_vector_get(data->costheta, i);
      sint  = gsl_vector_get(data->sintheta, i);
      
      t = x*cost + y*sint;
      
      row = interp_nearest(data->t, t, data->raysid);
      
      if(SIGN(row)>0)	
	s += interp_linear(&P.vector, data->t, row, t);
    }
 
  s = s/sT;
 
  if(data->scan != RAFT_CT)
    s = s/2;
  
  /* return MAX(0,s); */
  return s;
}


/*+====================================================+
  
  FUNCTION  eval_anglesum_ert

  Compute the angle sum for the ERT backprojection formula.
  
  ERT stands for the Exponential Radon Transform.

  Input  

  data - scan data . See <raft_scan_t>.
  projection - projection vector. 
  j - pixel index (direction y).
  k - pixel index (direction x).   
    
  Return 

  Value for the reconstructed pixel.
    
  +====================================================+
*/

double eval_anglesum_ert(raft_scan_t *data,
			 gsl_vector *projection,
			 double *att,
			 int j,
			 int k)
{
  if(data->viewsid == RAFT_REG)
    {
      int i, row, nrays, size;
      double x, y, sum, s, da, cost, sint, t, interp, prod;
      gsl_vector_view P;
      
      nrays = raft_scan_get_nrays(data);
      size  = raft_scan_get_size(data);

      da   = raft_scan_get_angle(data,1)-raft_scan_get_angle(data,0);
      x    = raft_scan_get_x(data,k);
      y    = raft_scan_get_y(data,size-1-j);      
      sum  = 0; 
      
      for(i=0; i < raft_scan_get_nviews(data); i++)
	{
	  /* ith row of the projection matrix */
	  P = gsl_vector_subvector(projection, i*nrays, nrays);  
	  
	  cost  = gsl_vector_get(data->costheta, i);
	  sint  = gsl_vector_get(data->sintheta, i);
	  
	  t =  x*cost + y*sint;
	  s = -x*sint + y*cost;
	  
	  row = interp_nearest(data->t, t, data->raysid);
	  
	  if(SIGN(row)>0)	
	    {
	      interp = interp_linear(&P.vector, data->t, row, t);
	      prod = s*(*att);	      
	      sum +=  interp * exp(prod);	      
	    }
	}
      
      return sum*da; 
    }
  else
    {
      int i, row, nrays, size;
      double x, y, s, sum, da, cost, sint, t, interp, prod;
      gsl_vector_view P;
      
      nrays = raft_scan_get_nrays(data);
      size  = raft_scan_get_size(data);

      x   = raft_scan_get_x(data,k);
      y   = raft_scan_get_y(data,size-1-j);
      sum = 0; 
      
      for(i=0; i < raft_scan_get_nviews(data) - 1; i++)
	{
	  /* ith row of the projection matrix */
	  P = gsl_vector_subvector(projection, i*nrays, nrays);  
	  
	  da = raft_scan_get_angle(data,i+1) - raft_scan_get_angle(data,i);

	  cost  = gsl_vector_get(data->costheta, i);
	  sint  = gsl_vector_get(data->sintheta, i);

	  t =  x*cost + y*sint;
	  s = -x*sint + y*cost;
	  
	  row = interp_nearest(data->t, t, data->raysid);
	  
	  if(SIGN(row)>0)	
	    {
	      interp = interp_linear(&P.vector, data->t, row, t);
	      prod = s*(*att);	      
	      sum += interp * exp(prod) * da;	      
	    }
	}
      
      return sum; 
    }
}

/*+====================================================+
  
  FUNCTION  eval_anglesum_generalized

  Compute the angle sum for the generalized 
  attenuated backprojection formula.
  
  Input  

  data - scan data . See <raft_scan_t>.
  weight - weight function
  projection - projection vector 
  j - pixel index (direction y)
  k - pixel index (direction x)   
    
  Return 

  Pixel backprojection.
  
  +====================================================+
*/

double eval_anglesum_generalized(raft_scan_t *data,
				 gsl_vector **weight,
				 gsl_vector *projection,
				 int j,
				 int k)
{
  if(data->viewsid == RAFT_REG)
    {
      int i, row, nrays, size, Pixel, Ray;
      double x, y, s, da, cost, sint, t, w;
      gsl_vector_view P;
      
      nrays = raft_scan_get_nrays(data);
      size  = raft_scan_get_size(data);

      da = raft_scan_get_angle(data,1)-raft_scan_get_angle(data,0);
      x = raft_scan_get_x(data,k);
      y = raft_scan_get_y(data,size-1-j);
      s = 0; 

      Pixel = k + j*size;      
      
      for(i=0; i < raft_scan_get_nviews(data); i++)
	{
	  /* ith row of the projection matrix */
	  P = gsl_vector_subvector(projection, i*nrays, nrays);  
	  
	  cost  = gsl_vector_get(data->costheta, i);
	  sint  = gsl_vector_get(data->sintheta, i);
	  	  
	  t = x*cost + y*sint;
	  
	  row = interp_nearest(data->t, t, data->raysid);
	  
	  Ray = row + i*nrays;

	  w = raft_projection_radon_dbt_get(data, weight, Ray, Pixel);
	  
	  if(SIGN(row)>0)	
	    s += w * interp_linear(&P.vector, data->t, row, t);
	}
      
      return MAX(0,(s*da));
      //return s*da;
    }
  else
    {
      int i, row, nrays, size, Pixel, Ray;
      double x, y, s, da, cost, sint, t, w;
      gsl_vector_view P;
      
      nrays = raft_scan_get_nrays(data);
      size  = raft_scan_get_size(data);

      x = raft_scan_get_x(data,k);
      y = raft_scan_get_y(data,size-1-j);
      s = 0; 
      
      Pixel = k + j*size;
      
      for(i=0; i < raft_scan_get_nviews(data) - 1; i++)
	{
	  /* ith row of the projection matrix */
	  P = gsl_vector_subvector(projection, i*nrays, nrays);  
	  
	  da = raft_scan_get_angle(data,i+1) - raft_scan_get_angle(data,i);

	  cost  = gsl_vector_get(data->costheta, i);
	  sint  = gsl_vector_get(data->sintheta, i);

	  t = x*cost + y*sint;
	  
	  row = interp_nearest(data->t, t, data->raysid);
	
	  Ray   = row + i*nrays;
	  
	  w = raft_projection_radon_dbt_get(data, weight, Ray, Pixel);
	  
	  if(SIGN(row)>0)	
	    s += w * interp_linear(&P.vector, data->t, row, t) * da;
	}
      
      /*return MAX(0,s);*/
      return s;
    }
}


/*+====================================================+
  
  FUNCTION  eval_anglesum_inversion

  Compute the angle sum for inversion.
  
  Input  

  data - scan data . See <raft_scan_t>.
  mpr - modified projections
  weight - weight matrix
  i - row (pixel y)
  k - column (pixel x)
    
  Return 

  Value for reconstructed pixel.
  
  +====================================================+
*/

double eval_anglesum_inversion(raft_scan_t *data,
			       gsl_vector *mpr,
			       gsl_vector **weight,
			       int i,
			       int k)
{
  double Z2[2], Z1[2] = {0, -0.5}; 
  /* double Z2[2] = {-0.5, 0}; */
  double absZ1, absZ2, cost1, sint1, cost2, sint2, PC, PS;
  
  Z2[0] =  Z1[1];
  Z2[1] = -Z1[0];
  
  absZ1 = gsl_hypot(Z1[0],Z1[1]);
  absZ2 = gsl_hypot(Z2[0],Z2[1]);
  
  cost1 = -Z1[1]/absZ1;
  sint1 =  Z1[0]/absZ1;
  
  cost2 = -Z2[0]/absZ2;
  sint2 = -Z2[1]/absZ2;
  
  if(data->viewsid == RAFT_REG)
    {
      int T, j, nrays, nviews, size;
      double dm, dwght[3], dt, dx, da, ida, idx, idt, ipi, m;
      double cost, sint, t, g, x, y, sum;
      gsl_vector_view P;
      
      nrays  = raft_scan_get_nrays(data);
      nviews = raft_scan_get_nviews(data);
      size   = raft_scan_get_size(data);
      
      dx  = raft_scan_get_x(data,1) - raft_scan_get_x(data,0);
      dt  = raft_scan_get_ray(data,1)-raft_scan_get_ray(data,0);
      da  = (2*MPI)/nviews;
      idx = 1/dx;
      idt = 1/dt;
      ida = 1/da;
      ipi = 1/(2*MPI);
      
      y = raft_scan_get_y(data, size-1-i);
      x = raft_scan_get_x(data, k);
      
      sum = 0;
      
      for(j=0; j < nviews; j++)
	{
	  cost  = gsl_vector_get(data->costheta, j);
	  sint  = gsl_vector_get(data->sintheta, j);
	  
	  PC = cost1 * cost - sint1*sint;
	  PS = sint2 * cost + cost2*sint; 

	  t = x*cost + y*sint;
	  T = interp_nearest(data->t, t, data->raysid);
	  
	  if(SIGN(T) > 0)
	    {
	      /* jth row of the projection matrix */
	      P = gsl_vector_subvector(mpr, j*nrays, nrays);	         
	      	      
	      /* diff modified projection */
	      
	      m = interp_linear(&P.vector, data->t, T, t);
	      
	      dm = -m;
	      	      
	      if((T+1)<nrays)
		{
		  dm += interp_linear(&P.vector, data->t, T+1, t + dt);
		}
		
	      dm *= idt;
	      
	      /* diff inverse weight matrix */
	      
	      g = gsl_vector_get(weight[j], k + i*size);   	      

	      dwght[0] = -g;
	      dwght[1] = -g;
	      
	      if((k+1)<size && (i+1)<size)
		{
		  dwght[0] += gsl_vector_get(weight[j], (k+1) + i * size);
		  dwght[1] += gsl_vector_get(weight[j], k + (i+1) * size);
		}
	      else
		{
		  dwght[0] = 0;
		  dwght[1] = 0;
		}

	      dwght[0] *= idx;
	      dwght[1] *= idx;
	      
	      /* sum */
	      
	      g = 1/g;
	      
	      /* cost, sint */

	      sum += (dm * cost * g - (m*g) * (dwght[0]*g)) * absZ1* PC +
		     (dm * sint * g - (m*g) * (dwght[1]*g)) * absZ2* PS;
	    }
	}
      
      return MAX(0,(sum*ipi*da));
    }
  else
    {
      int T, j, nrays, nviews, size;
      double dm, dwght[3], dt, dx, da, ida, idt, idx, ipi, m;
      double cost, sint, t, g, x, y, sum, theta;
      gsl_vector_view P;
      
      nrays  = raft_scan_get_nrays(data);
      nviews = raft_scan_get_nviews(data);
      size   = raft_scan_get_size(data);
      
      dx  = raft_scan_get_x(data,1) - raft_scan_get_x(data,0);  
      dt  = raft_scan_get_ray(data,1)-raft_scan_get_ray(data,0);
      idx = 1/dx;
      idt = 1/dt;
      ipi = 1/(2*MPI);
      
      y = raft_scan_get_y(data, size-1-i);
      x = raft_scan_get_x(data, k);
      
      sum = 0;
      
      for(j=0; j < nviews-1; j++)
	{
	  theta = raft_scan_get_angle(data,j);
	  cost  = gsl_vector_get(data->costheta, j);
	  sint  = gsl_vector_get(data->sintheta, j);
	  
	  PC = cost1 * cost - sint1*sint;
	  PS = sint2 * cost + cost2*sint; 

	  da  = raft_scan_get_angle(data,j+1) - theta;
	  ida = 1/da;

	  t = x*cost + y*sint;
	  T = interp_nearest(data->t, t, data->raysid);
	  
	  if(SIGN(T) > 0)
	    {
	      /* jth row of the projection matrix */
	      P = gsl_vector_subvector(mpr, j*nrays, nrays);
	      
	      /* diff modified projection */
	      
	      m = interp_linear(&P.vector, data->t, T, t);

	      dm = -m;
	      
	      if((T+1)<nrays)
		{
		  dm += interp_linear(&P.vector, data->t, T+1, t+dt);
		}
		
	      dm *= idt;
	            
	      /* diff inverse weight matrix */
	      
	      g = gsl_vector_get(weight[j], k + i*size);   	      

	      dwght[0] = -g;
	      dwght[1] = -g;
	      
	      if((k+1)<size && (i+1)<size)
		{
		  dwght[0] += gsl_vector_get(weight[j], k+1 + i * size);
		  dwght[1] += gsl_vector_get(weight[j], k + (i+1) * size);
		}
	      else
		{
		  dwght[0] = 0;
		  dwght[1] = 0;
		}
	      
	      dwght[0] *= idx;
	      dwght[1] *= idx;     
	      
	      /* sum */
	      
	      g = 1/g;
	      
	      sum += ((dm * cost * g - (m*g) * (dwght[0]*g)) * absZ1* PC +
		      (dm * sint * g - (m*g) * (dwght[1]*g)) * absZ2* PS) * da;	     
	    }
	}
      
      return MAX(0,(sum*ipi));
    }
}

/*+====================================================+
  
  FUNCTION  eval_pixel_backprojection

  Evaluate the backprojection for a given pixel. 
  
  Input  

  scan - scan data . See <raft_scan_t>.
  projection - projection vector
  j - pixel number
   
  Return 

  Evaluated backprojection formula.
  
  +====================================================+
*/

double eval_pixel_backprojection(raft_scan_t *data,
				 gsl_vector *projection,
				 int j)
{
  double s[4];
  int x, y, size;
  div_t D;  

  size = raft_scan_get_size(data); 

  D = div(j, raft_scan_get_size(data));
  x = D.rem;
  y = D.quot;
  
  return eval_anglesum(data, projection, y, x, s);
}

/*+====================================================+
  
  FUNCTION  filtered_projection

  Compute a filter of the projection matrix.

  Compute the convolution of projection matrix with a 
  transfer function.
  
  Input  

  data - scan data . See <raft_scan_t>.
  workspace - backprojection workspace. See <raft_backp_t>.
  ftype - filter type
  p - projection vector

  Output 

  q - filtered projection vector

  Return 

  RAFT_EDOM - Domain error. Nonequispaced rays
  RAFT_SUCCESS - Successfull operation.

  _Remark_ 
  
  It only works for equispaced rays.
  
  +====================================================+
*/

int filtered_projection(raft_scan_t *data, 
			raft_backp_t *workspace,
			gsl_vector *p,
			gsl_vector *q)
{
  if(data->raysid == RAFT_REG)
    {
      int j, nrays, nviews;
      double dt;
      gsl_vector_view v, w, V, W;
      
      nrays = raft_scan_get_nrays(data);
      nviews = raft_scan_get_nviews(data);
      dt = raft_scan_get_raystep(data);
      
      for(j=0; j < nviews; j++)
	{
	  v = gsl_vector_subvector(p, j*nrays, nrays);
	  w = gsl_vector_subvector(q, j*nrays, nrays);

	  gsl_vector_set_all(workspace->fourier.signal, 0.0);
	  V = gsl_vector_subvector(workspace->fourier.signal, 0, nrays);
	  gsl_blas_dcopy(&v.vector, &V.vector);	  
	  
	  fft_complex(workspace->fourier.fftRE1, 
		      workspace->fourier.fftIM1,
		      workspace->fourier.signal, 
		      workspace->fourier.zero,
		      workspace->fourier.fftwork);      
	  
	  gsl_vector_mul(workspace->fourier.fftRE1, workspace->fourier.fftREimpulse);
	  gsl_vector_mul(workspace->fourier.fftIM1, workspace->fourier.fftREimpulse);
	  
	  gsl_vector_set_all(workspace->fourier.signal, 0.0);
	  ifft_complex(workspace->fourier.impulse,
		       workspace->fourier.fftIM2, 
		       workspace->fourier.fftRE1, 
		       workspace->fourier.fftIM1,
		       workspace->fourier.fftwork);     	  
	  
	  W = gsl_vector_subvector(workspace->fourier.impulse, 0, nrays); 
	  gsl_blas_dcopy(&W.vector, &w.vector);	  
	    
	  /*fft_shift(&w.vector, nrays);      
	  gsl_vector_scale(&w.vector, dt);*/
	}
      
      return RAFT_SUCCESS;
    }
  else
    return RAFT_EDOM;
}

/*+====================================================+
  
  FUNCTION  modified_projection_inversion

  Compute modified projections for inversion.
    
  Input  

  data - scan data . See <raft_scan_t>.
  backp - backprojection workspace. See <raft_backp_t>.
  p - projection vector (sinogram/inversion).  
  d - projection vector (sinogram/data).
  ep2 - function exp(p/2) 

  Return 

  RAFT_EDOM - Domain error. Nonequispaced rays
  RAFT_SUCCESS - Successfull operation.

  _Remark_ 
  
  It only works for equispaced rays.
  
  +====================================================+
*/

void modified_projection_inversion(raft_scan_t *data, 
				   raft_backp_t *backp,
				   gsl_vector *p, 
				   gsl_vector *d,
				   gsl_vector *ep2)
{
  int i, j, nviews, nrays;
  double Hc, Hs;
  gsl_vector_view D, P, EP, M;

  nviews = raft_scan_get_nviews(data);
  nrays  = raft_scan_get_nrays(data);
  
  for(j = 0; j < nviews; j++)
    {
      D  = gsl_vector_subvector(d,   j*nrays, nrays);      
      P  = gsl_vector_subvector(p,   j*nrays, nrays);
      EP = gsl_vector_subvector(ep2, j*nrays, nrays);
      M  = gsl_vector_subvector(backp->inversion.modproj, j*nrays, nrays);
      
      gsl_vector_set_all(backp->inversion.aux1, 0.0);
      gsl_vector_set_all(backp->inversion.aux2, 0.0);
      gsl_vector_set_all(backp->inversion.aux3, 0.0);
      gsl_vector_set_all(backp->inversion.aux4, 0.0);
      gsl_vector_set_all(&M.vector, 0.0);
            
      hilbert_transform(&P.vector, 
			backp->inversion.hilbert, 
			backp->inversion.hilbwork);
      
      for(i=0; i < nrays; i++)
	{
	  Hc = cos(gsl_vector_get(backp->inversion.hilbert, i)/2);
	  Hs = sin(gsl_vector_get(backp->inversion.hilbert, i)/2);
	  
	  gsl_vector_set(backp->inversion.coshilbert, i, Hc);
	  gsl_vector_set(backp->inversion.sinhilbert, i, Hs);
	} 
      
      gsl_vector_add(backp->inversion.aux1, backp->inversion.coshilbert);
      gsl_vector_mul(backp->inversion.aux1, &EP.vector);
      gsl_vector_mul(backp->inversion.aux1, &D.vector);
      
      gsl_vector_add(backp->inversion.aux2, backp->inversion.sinhilbert);
      gsl_vector_mul(backp->inversion.aux2, &EP.vector);
      gsl_vector_mul(backp->inversion.aux2, &D.vector);
      
      hilbert_transform(backp->inversion.aux1, 
			backp->inversion.haux1, 
			backp->inversion.hilbwork);  
      
      hilbert_transform(backp->inversion.aux2, 
			backp->inversion.haux2, 
			backp->inversion.hilbwork);  
      
      gsl_vector_add(backp->inversion.aux3, backp->inversion.coshilbert);
      gsl_vector_mul(backp->inversion.aux3, backp->inversion.haux1);
      
      gsl_vector_add(backp->inversion.aux4, backp->inversion.sinhilbert);
      gsl_vector_mul(backp->inversion.aux4, backp->inversion.haux2);
      
      gsl_vector_add(&M.vector, backp->inversion.aux3);
      gsl_vector_add(&M.vector, backp->inversion.aux4);	  
      gsl_vector_div(&M.vector, &EP.vector); 
    }    
}

/*+====================================================+
  
  FUNCTION  define_novikov_workspace

  Define the workspace for Novikov's inversion.
    
  Input  

  data - scan data . See <raft_scan_t>.
  backp - backprojection workspace. See <raft_backp_t>.
  att - attenuation matrix
  p - SPECT data.  
  
  +====================================================+
*/

void define_novikov_workspace(raft_scan_t *data, 
			      raft_backp_t *backp,
			      gsl_vector *att, 
			      gsl_vector *p)
{
  if(!backp->novikov.defined)
    {
      /* weight function */
  
      raft_weight_spect(data, &backp->projwork, att);  
      
      /* novikov's sinogram */ 
  
      raft_projection_radon(data, &backp->projwork,
			    att, backp->novikov.radon);
      
      gsl_blas_dcopy(backp->projwork.likelihood.exprad2, 
		     backp->novikov.expradon);
      
      backp->novikov.defined = 1; 
    }
  
  /* modified projections */
  
  modified_projection_inversion(data, 
				backp,
				backp->novikov.radon, 
				p,
				backp->novikov.expradon);
}

/*+====================================================+
  
  FUNCTION  define_invxfctPartial_workspace

  Define the workspace for partial XFCT inversion.
    
  Input  

  data - scan data . See <raft_scan_t>.
  backp - backprojection workspace. See <raft_backp_t>.
  attT - transmission attenuation matrix
  attF - fluorescence attenuation matrix
  p - XFCT data.  
  
  Output

  apert - angle section 
  
  +====================================================+
*/

void define_invxfctPartial_workspace(raft_scan_t *data, 
				     raft_backp_t *backp,
				     gsl_vector *attT,
				     gsl_vector *attF,
				     gsl_vector *p,
				     double *apert)
{
  if(!backp->invxfct.defined)
    {
      div_t D;
      int N, J, j, k, size, nviews, nrays;
      double kernel, radon, value[2];
      gsl_vector_view R,E;
      
      size   = raft_scan_get_size(data);
      nviews = raft_scan_get_nviews(data);
      nrays  = raft_scan_get_nrays(data);
      
      N = (int)(nviews/2);           
      
      /* weight */
      
      raft_weight_xfct_partial(data, &backp->projwork, attT, attF);      
      
      *apert = backp->projwork.xfct.apert;

      /* inversion sinogram */
      
      for(j = 0; j < nviews; j++)
	{
	  D = div(j + N, nviews);
	  J = D.rem; 
	  
	  R = gsl_vector_subvector(backp->invxfct.radon,    j*nrays, nrays); 
	  E = gsl_vector_subvector(backp->invxfct.expradon, j*nrays, nrays);
	  
	  gsl_blas_dscal(0, backp->invxfct.Phant);
	  gsl_blas_daxpy(-1.0, backp->projwork.xfct.LogWeightP[j], backp->invxfct.Phant);
	  gsl_blas_daxpy(-1.0, backp->projwork.xfct.LogWeightP[J], backp->invxfct.Phant);
	  
	  for(k=0; k<nrays; k++)
	    {
	      kernel = gsl_vector_get(backp->invxfct.kernel, k);
	      
	      radon = eval_rayintegral_(data, backp->invxfct.Phant, j, k, 
				       data->memory.a, &data->memory.norma, &data->memory.dotp);
	      
	      gsl_vector_set(backp->invxfct.RadonPhant, k + j*nrays, kernel*radon);

	      value[0] = kernel * radon;
	      value[1] = exp(value[0]/2);
	      
	      gsl_vector_set(&R.vector, k, value[0]);
	      gsl_vector_set(&E.vector, k, value[1]);
	    }
	}
      
      backp->invxfct.defined = 1;
    }
  
  /* modified projections */
  
  *apert = backp->projwork.xfct.apert;

  modified_projection_inversion(data, 
				backp,
				backp->invxfct.radon, 
				p,
				backp->invxfct.expradon);
    
}

/*+====================================================+
  
  FUNCTION  define_invxfct_workspace

  Define the workspace for XFCT inversion.
    
  Input  

  data - scan data . See <raft_scan_t>.
  backp - backprojection workspace. See <raft_backp_t>.
  attT - transmission attenuation matrix
  attF - fluorescence attenuation matrix
  p - XFCT data.  
  
  +====================================================+
*/

void define_invxfct_workspace(raft_scan_t *data, 
			      raft_backp_t *backp,
			      gsl_vector *attT,
			      gsl_vector *attF,
			      gsl_vector *p)
{
  if(!backp->invxfct.defined)
    {
      div_t D;
      int N, J, k, j, nviews, nrays, size;
      double value[2], apert, kernel, radon;
      gsl_vector_view R, E;
      
      size   = raft_scan_get_size(data);
      nviews = raft_scan_get_nviews(data);
      nrays  = raft_scan_get_nrays(data);
      N      = (int)(nviews/2);
      
      /* weight */
      
      raft_weight_xfct(data, &backp->projwork, attT, attF);      
      
      apert = backp->projwork.xfct.apert;
      
      /* inversion sinogram */
      
      gsl_vector_set_all(backp->invxfct.expradon, 1.0);
      
      for(j=0; j<nviews; j++)
	{
	  D = div(j + N, nviews);
	  J = D.rem;
	  
	  R = gsl_vector_subvector(backp->invxfct.radon,    j*nrays, nrays); 
	  E = gsl_vector_subvector(backp->invxfct.expradon, j*nrays, nrays);
	  
	  /*gsl_blas_dscal(0, backp->invxfct.Phant);*/
	  gsl_vector_set_all(backp->invxfct.Phant, 0);
	  gsl_blas_daxpy(-1.0, backp->projwork.xfct.LogWeightF[j], backp->invxfct.Phant);
	  gsl_blas_daxpy(-1.0, backp->projwork.xfct.LogWeightF[J], backp->invxfct.Phant);
	  
	  for(k=0; k<nrays; k++)
	    {
	      kernel = gsl_vector_get(backp->invxfct.kernel, k);
	      
	      radon = eval_rayintegral_(data, backp->invxfct.Phant, j, k, 
				       data->memory.a, &data->memory.norma, &data->memory.dotp);
	      
	      gsl_vector_set(backp->invxfct.RadonPhant, k + j*nrays, kernel*radon);

	      value[0] = kernel*radon;
	      value[1] = exp(value[0]/2);
	      
	      gsl_vector_set(&R.vector, k, value[0]);
	      gsl_vector_set(&E.vector,  k, value[1]);
	    }
	}
      
      backp->invxfct.defined = 1; 
    }
  
  /* modified projections */
  
  modified_projection_inversion(data, 
				backp,
				backp->invxfct.radon, 
				p,
				backp->invxfct.expradon);
  
}


/*######################################################
  Public
  ####################################################*/

/*######################################################
  Section: Allocation
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_backp_workspace_alloc
  
  Workspace allocation for filtered backprojection
  
  Input: 

  data - scan data
  filter - filter type

  Output:

  workspace - backprojection workspace   
  
  +====================================================+
*/

int raft_backp_workspace_alloc(raft_scan_t *data,
			       int filter,
			       raft_backp_t *workspace)  
{ 
  int nrays, ndata, npixels, nviews, size, N, j;
  double dt, wc, dx, w, wn, dw, t, fun, s;

  if(filter<0)
    return RAFT_EDOM;

  size    = raft_scan_get_size(data);
  nrays   = raft_scan_get_nrays(data);
  nviews  = raft_scan_get_nviews(data);
  ndata   = raft_scan_get_ndata(data);
  npixels = raft_scan_get_npixels(data);

  workspace->nviews = nviews;

  dx = raft_scan_get_x(data,1)-raft_scan_get_x(data,0);
  dt = raft_scan_get_raystep(data);
  wc = 1/(2*dt);
  
  /* integrand workspace */

  workspace->integrand.vector = gsl_vector_alloc(ndata);
  
  /* Fourier workspace */
  
  N = nrays;
  
  workspace->fourier.fftwork = (fft_t *)malloc(sizeof(fft_t));
  fft_workspace_alloc(workspace->fourier.fftwork, N, dt);
  
  workspace->fourier.signal  = gsl_vector_calloc(N);
  workspace->fourier.impulse = gsl_vector_calloc(N);
  workspace->fourier.zero    = gsl_vector_calloc(N);
  workspace->fourier.q       = gsl_vector_calloc(ndata);
  
  workspace->fourier.fftREimpulse = gsl_vector_calloc(N);
  workspace->fourier.fftIMimpulse = gsl_vector_calloc(N);
  workspace->fourier.fftRE1       = gsl_vector_calloc(N);
  workspace->fourier.fftRE2       = gsl_vector_calloc(N);
  workspace->fourier.fftIM1       = gsl_vector_calloc(N);
  workspace->fourier.fftIM2       = gsl_vector_calloc(N);
  
  workspace->fourier.fftREimpulseTretiak = gsl_vector_calloc(N);
  workspace->fourier.fftIMimpulseTretiak = gsl_vector_calloc(N);
  
  switch(filter)
    {
    case RAFT_COSINE:
      
      wn = wc;
      dw = (2*wn)/N;
  
      for(j=0; j < N; j++)
	{
	  w = -wn + j*dw;	  
	  fun = fabs(w) * cos(MPI*w*dt);	  
	  gsl_vector_set(workspace->fourier.fftREimpulse, j, fun);
	}      
      fft_shift(workspace->fourier.fftREimpulse, N);
      
      break;
      
    case RAFT_SHEPP:
      
      wn = wc;
      dw = (2*wn)/N;
  
      for(j=0; j < N; j++)
	{
	  w = -wn + j*dw;	  
	  fun = fabs(w) * gsl_sf_sinc(w*dt);	  
	  gsl_vector_set(workspace->fourier.fftREimpulse, j, fun);
	}       
      fft_shift(workspace->fourier.fftREimpulse, N);

      break;
      
    default:

      wn = wc;
      dw = (2*wn)/N;
      
      for(j=0; j < N; j++)
	{
	  w = -wn + j*dw;	  
	  fun = fabs(w); 	  
	  gsl_vector_set(workspace->fourier.fftREimpulse, j, fun);
	}       
      fft_shift(workspace->fourier.fftREimpulse, N);
    }

  raft_projection_workspace_alloc(data, &workspace->projwork);

  /* Inversion workspace */ 

  workspace->inversion.hilbwork = (hilb_t *)malloc(sizeof(hilb_t));
  
  hilbert_workspace_alloc(workspace->inversion.hilbwork, 
			  nrays, RAFT_RAMLAK, dt);

  workspace->inversion.modproj  = gsl_vector_alloc(ndata);  
  workspace->inversion.hilbert    = gsl_vector_alloc(nrays);
  workspace->inversion.coshilbert = gsl_vector_alloc(nrays);
  workspace->inversion.sinhilbert = gsl_vector_alloc(nrays);
  
  workspace->inversion.aux1  = gsl_vector_alloc(nrays);
  workspace->inversion.aux2  = gsl_vector_alloc(nrays);
  workspace->inversion.aux3  = gsl_vector_alloc(nrays);
  workspace->inversion.aux4  = gsl_vector_alloc(nrays);
  workspace->inversion.haux1 = gsl_vector_alloc(nrays);
  workspace->inversion.haux2 = gsl_vector_alloc(nrays);    
  workspace->inversion.kernel = gsl_vector_alloc(nrays);
  workspace->inversion.data   = gsl_vector_alloc(ndata);
  workspace->inversion.att    = gsl_vector_alloc(npixels);

  /* Novikov's workspace */

  workspace->novikov.radon    = gsl_vector_alloc(ndata);
  workspace->novikov.expradon = gsl_vector_alloc(ndata);  
  workspace->novikov.defined  = 0;

  /* XFCT inversion workspace */
  
  workspace->invxfct.RadonPhant = gsl_vector_alloc(ndata);
  workspace->invxfct.Phant      = gsl_vector_alloc(npixels);
  workspace->invxfct.kernel     = gsl_vector_alloc(nrays);  

  for(j=1; j<nrays-1; j++)
    {
      t   = gsl_vector_get(data->t, j);
      s   = sqrt(1-t*t);
      fun = 1.0/(2*s);
      gsl_vector_set(workspace->invxfct.kernel, j, fun);
      gsl_vector_set(workspace->inversion.kernel, j, s);
    }
  gsl_vector_set(workspace->invxfct.kernel, 0, 0);
  gsl_vector_set(workspace->invxfct.kernel, nrays-1, 0);
  gsl_vector_set(workspace->inversion.kernel, 0, 0);
  gsl_vector_set(workspace->inversion.kernel, nrays-1, 0);

  workspace->invxfct.radon     = gsl_vector_alloc(ndata);
  workspace->invxfct.radonF    = gsl_vector_alloc(ndata);
  workspace->invxfct.expradon  = gsl_vector_alloc(ndata);
  workspace->invxfct.expradonF = gsl_vector_alloc(ndata);
  workspace->invxfct.aux       = gsl_vector_alloc(npixels);
  workspace->invxfct.defined   = 0;

  workspace->neumann.rad       = gsl_vector_alloc(ndata);
  workspace->neumann.next      = gsl_vector_alloc(npixels);
  workspace->neumann.direction = gsl_vector_alloc(npixels);
 
  /* Matrices */
  
  workspace->integrand.matrix = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);        

  for(j=0; j < nviews; j++)
    {
      workspace->integrand.matrix[j] = gsl_vector_alloc(npixels);            
    }

  /* Filter of the backprojection */

  workspace->fob.fft2work = (fft_t *)malloc(sizeof(fft_t));
  fft_workspace_alloc(workspace->fob.fft2work, size, dx);
  
  workspace->fob.zero    = gsl_vector_calloc(npixels);
  workspace->fob.impulse = gsl_vector_alloc(npixels);
  workspace->fob.filter = gsl_vector_alloc(npixels);
  
  workspace->fob.fftREimpulse = gsl_vector_alloc(npixels);
  workspace->fob.fftIMimpulse = gsl_vector_alloc(npixels);
  workspace->fob.fftREback    = gsl_vector_alloc(npixels);
  workspace->fob.fftIMback    = gsl_vector_alloc(npixels);

  wc = 1/(2*dx);  

  switch(filter)
    {
    case RAFT_COSINE:
      
      filter_2D(RAFT_COSINE, size, dx, wc, workspace->fob.filter);
      break;

    case RAFT_SHEPP:

      filter_2D(RAFT_SHEPP, size, dx, wc, workspace->fob.filter);
      break;

    default:      
      filter_2D(RAFT_RAMLAK, size, dx, wc, workspace->fob.filter);
      
    }
  
  return RAFT_SUCCESS;
}

/*+====================================================+
  
  FUNCTION: raft_backp_workspace_free
  
  Frees backprojection workspace 
  
  Input: 

  workspace - backprojection workspace  
  
  +====================================================+
*/

void raft_backp_workspace_free(raft_backp_t *workspace)
{
  int j;

  /* integrand workspace */
  
  gsl_vector_free(workspace->integrand.vector);
    
  /* Fourier workspace */

  gsl_vector_free(workspace->fourier.signal);
  gsl_vector_free(workspace->fourier.impulse);
  gsl_vector_free(workspace->fourier.zero);
  gsl_vector_free(workspace->fourier.q);
  gsl_vector_free(workspace->fourier.fftREimpulse);
  gsl_vector_free(workspace->fourier.fftIMimpulse);
  gsl_vector_free(workspace->fourier.fftRE1);
  gsl_vector_free(workspace->fourier.fftRE2);
  gsl_vector_free(workspace->fourier.fftIM1);
  gsl_vector_free(workspace->fourier.fftIM2);

  gsl_vector_free(workspace->fourier.fftREimpulseTretiak);
  gsl_vector_free(workspace->fourier.fftIMimpulseTretiak);
  
  fft_workspace_free(workspace->fourier.fftwork);
  free(workspace->fourier.fftwork);

  raft_projection_workspace_free(&workspace->projwork);

  /* Inversion workspace */ 
  
  gsl_vector_free(workspace->inversion.modproj);
  gsl_vector_free(workspace->inversion.hilbert);
  gsl_vector_free(workspace->inversion.coshilbert);
  gsl_vector_free(workspace->inversion.sinhilbert);
  gsl_vector_free(workspace->inversion.aux1);
  gsl_vector_free(workspace->inversion.aux2);
  gsl_vector_free(workspace->inversion.aux3);
  gsl_vector_free(workspace->inversion.aux4);  
  gsl_vector_free(workspace->inversion.haux1);
  gsl_vector_free(workspace->inversion.haux2);
  gsl_vector_free(workspace->inversion.kernel);
  gsl_vector_free(workspace->inversion.data);
  gsl_vector_free(workspace->inversion.att);
  
  hilbert_workspace_free(workspace->inversion.hilbwork);
  free(workspace->inversion.hilbwork);
  
  /* Novikov's workspace */
  
  gsl_vector_free(workspace->novikov.radon);
  gsl_vector_free(workspace->novikov.expradon);  
  
  /* XFCT workspace */
  
  gsl_vector_free(workspace->invxfct.Phant);
  gsl_vector_free(workspace->invxfct.RadonPhant);
  gsl_vector_free(workspace->invxfct.kernel);
  
  gsl_vector_free(workspace->invxfct.radon);
  gsl_vector_free(workspace->invxfct.radonF); 
  gsl_vector_free(workspace->invxfct.expradon);
  gsl_vector_free(workspace->invxfct.expradonF);
  gsl_vector_free(workspace->invxfct.aux);

  gsl_vector_free(workspace->neumann.rad);
  gsl_vector_free(workspace->neumann.next);
  gsl_vector_free(workspace->neumann.direction);
  
  /* Matrices */

  for(j = 0; j < workspace->nviews; j++)
    {
      gsl_vector_free(workspace->integrand.matrix[j]);        
    }

  free(workspace->integrand.matrix);

  /* Filter of the backprojection */

  gsl_vector_free(workspace->fob.zero);
  gsl_vector_free(workspace->fob.impulse);
  gsl_vector_free(workspace->fob.filter);
  gsl_vector_free(workspace->fob.fftREimpulse);
  gsl_vector_free(workspace->fob.fftIMimpulse);
  gsl_vector_free(workspace->fob.fftREback);
  gsl_vector_free(workspace->fob.fftIMback);

  fft_workspace_free(workspace->fob.fft2work);
  free(workspace->fob.fft2work);
}

/*######################################################
  Section: Backprojection operators
  ####################################################*/

/*-----------------------------*/
/* Block Index for submatrices */

void matrix_blockIndex(int n, int T, int M, int *fC, int *lC, int *fR, int *lR)
{
  int e, i, L, C, J, K, NC, NR;

  e = (int) floor(sqrt(T)); 

  C = e + 1;
  L = (int) floor(T/C);

  J = n/C;
  K = n%C;
  
  NC = (int) floor(M/(C-1));
  NR = (int) floor(M/(L-1));
  
  if((J <= L-2) & (K <= C-2))
    {
      *fC = NC*K;
      *lC = *fC + NC - 1;
      
      *fR = NR*J;
      *lR = *fR + NR - 1;
    }
  else
    {      
      if((J > L-2) && (K > C-2))
	{
	  *fC = NC*K;
	  *lC = M - 1;
	  
	  *fR = NR*J;
	  *lR = M - 1;
	}
      
      if((J <= L-2) && (K > C-2))
	{
	  *fC = NC*K;
	  *lC = M - 1;
	  
	  *fR = NR*J;
	  *lR = *fR + NR - 1;
	}
      
      if((J > L-2) && (K <= C-2))
	{
	  *fC = NC*K;
	  *lC = *fC + NC - 1;
	  
	  *fR = NR*J;
	  *lR = M - 1;
	}
    } 
}

/*--------------------------------------------------*/
/* Backprojection with threads using: matrix blocks */

void *backp_pixelLoop_Mblock(void *t)
{
  int j, k, size;
  parBackLoop_t *param;
  double s[4], sum;
  
  param = (parBackLoop_t *)t;
  
  size = raft_scan_get_size(param->data);
  
  for(j = param->colIndex[0]; j < param->colIndex[1]; j++)
    {
      for(k = param->rowIndex[0]; k < param->rowIndex[1]; k++)
	{
	  sum = eval_anglesum(param->data, param->p, j, k, s);	  
	  
	  gsl_vector_set(param->b, k + j*size, sum);
	}
    }  
  
  pthread_exit(NULL);
}

/*------------------------------------------------*/
/* Backprojection with threads: using line blocks */

void *backp_pixelLoop_Lblock(void *t)
{
  int w, j, k, size;
  parBackLoop_t *param;
  double sum, s[4];
  
  param = (parBackLoop_t *)t;
  
  size = raft_scan_get_size(param->data);
  
  for(w = param->colIndex[0]; w < MIN(param->colIndex[1], size*size); w++)
    {
      j = w/size;
      k = w%size;
      
      sum = eval_anglesum(param->data, param->p, j, k, s);	  
      
      gsl_vector_set(param->b, k + j*size, sum);
    }  
  
  pthread_exit(NULL);
}

/*+====================================================+
  
  FUNCTION: raft_backp

  Compute backprojection matrix.

  Input: 

  data - scan data . See <raft_scan_t>.  
  p - projection vector
	       
  Output:

  b - backprojection vector

  _Remark_:
  
  It works either with equispaced or nonequispaced rays.
  The same for the views.

  _Contributor_:

  Hugo H.Slepicka @ LNLS - Parallel threads for the 
  backprojection operator (April/2013).

  +====================================================+
*/

void raft_backp(raft_scan_t *data, 
		gsl_vector *p,
		gsl_vector *b)
{
  pthread_t thread[BACKP_THREADS+1];
  pthread_attr_t attr;
  int size, n, rc, lpt, apt, fC, fR, lC, lR;    
  parBackLoop_t param[BACKP_THREADS+1];  
  void *status;
  
  /* Initialize and set thread detached attribute */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  size = raft_scan_get_size(data);
  
  if(BACKP_TTYPE==BACKP_MBLOCK)
    {
      for(n = 0; n < BACKP_THREADS + 1; n++) 
	{
	  param[n].data        = data;
	  param[n].p           = p;
	  param[n].b           = b; 
	  param[n].nthread     = n;      
	  
	  matrix_blockIndex(n, BACKP_THREADS, size, &fR, &lR, &fC, &lC);
	    
	  param[n].rowIndex[0] = fR;
	  param[n].rowIndex[1] = lR;
	  
	  param[n].colIndex[0] = fC;
	  param[n].colIndex[1] = lC;

	  rc = pthread_create(&thread[n], &attr, backp_pixelLoop_Mblock, (void *)&param[n]);
	}
    }
  else
    {
      int e = (int) floor((size*size)/BACKP_THREADS);
      
      for(n = 0; n < BACKP_THREADS+1; n++) 
	{
	  param[n].data        = data;
	  param[n].p           = p;
	  param[n].b           = b; 
	  param[n].nthread     = n;      
	  
	  param[n].colIndex[0] = e * n;
	  param[n].colIndex[1] = (n+1) * e;
	  
	  rc = pthread_create(&thread[n], &attr, backp_pixelLoop_Lblock, (void *)&param[n]);
	}
    }
      
  /* Free attribute and wait for the other threads */
  pthread_attr_destroy(&attr);
  for(n = 0; n < BACKP_THREADS + 1; n++) 
    {
      rc = pthread_join(thread[n], &status);
    }  
}

/*+====================================================+
  
  FUNCTION: raft_backp_partial

  Compute partial backprojection matrix, i.e., only for 
  a finite set of angles. 

  Input: 

  data - scan data . See <raft_scan_t>.  
  p - projection vector.
  T - vector with angle indexes.
  sT - length of vector T.
	       
  Output:

  b - backprojection vector

  _Remark_:
  
  It works either with equispaced or nonequispaced rays.
  The same for the views.

  +====================================================+
*/

void raft_backp_partial(raft_scan_t *data, 
			gsl_vector *p,
			gsl_vector *b,
			gsl_vector *T,
			int sT)
{
  int j, k, size;
  double sum;
  
  size = raft_scan_get_size(data);
  
  for(j = 0; j < size; j++)
    {
      for(k = 0; k < size; k++)
	{
	  sum = eval_anglesum_partial(data, p, j, k, T, sT);
	  
	  gsl_vector_set(b, k + j*size, sum);
	}
    }
}


/*+====================================================+
  
  FUNCTION: raft_backp_attenuated_ert

  Compute the ERT attenuated backprojection matrix.

  ERT stands for the Exponential Radon Transform.

  Input: 

  data - scan data . See <raft_scan_t>.  
  p - projection vector
  att - attenuation value
	       
  Output:

  b - backprojection vector

  _Remark_:
  
  It works either with equispaced or nonequispaced rays.
  The same for the views.
   
  +====================================================+
*/

void raft_backp_attenuated_ert(raft_scan_t *data, 
			       raft_backp_t *workspace,
			       double *att,
			       gsl_vector *p,
			       gsl_vector *b)
{
  int j, k, size;
  double sum;
  
  size = raft_scan_get_size(data);
  
  for(j = 0; j < size; j++)
    {
      for(k = 0; k < size; k++)
	{
	  sum = eval_anglesum_ert(data, p, att, j, k);
	  
	  gsl_vector_set(b, k + j*size, sum);
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_backp_attenuated_pet

  Compute the PET attenuated backprojection matrix.

  Input: 

  data - scan data . See <raft_scan_t>.  
  p - projection vector
  att - attenuation map
	       
  Output:

  b - backprojection vector

  _Remark_:
  
  It works either with equispaced or nonequispaced rays.
  The same for the views.
   
  +====================================================+
*/

void raft_backp_attenuated_pet(raft_scan_t *data, 
			       raft_backp_t *workspace,
			       gsl_vector *p,
			       gsl_vector *att,
			       gsl_vector *b)
{
  int j, k, size;
  double sum, s[4];

  size = raft_scan_get_size(data);

  if(!workspace->projwork.pet.defined)
    {
      raft_weight_pet(data, &workspace->projwork, att);
    }

  gsl_vector_memcpy(workspace->integrand.vector, 
		    workspace->projwork.pet.weight);
  
  gsl_vector_mul(workspace->integrand.vector, p);
  
  for(j = 0; j < size; j++)
    {
      for(k = 0; k < size; k++)
	{
	  sum = eval_anglesum(data, workspace->integrand.vector, j, k, s);
	  
	  gsl_vector_set(b, k + j*size, sum);
	}
    }
}


/*+====================================================+
  
  FUNCTION: raft_backp_attenuated_spect

  Compute the SPECT attenuated backprojection matrix.

  Input: 

  data - scan data . See <raft_scan_t>.  
  p - projection vector
  att - attenuation map
	       
  Output:

  b - backprojection vector

  _Remark_:
  
  It works either with equispaced or nonequispaced rays.
  The same for the views.
  
  +====================================================+
*/

void raft_backp_attenuated_spect(raft_scan_t *data, 
				 raft_backp_t *workspace,
				 gsl_vector *p,
				 gsl_vector *att,
				 gsl_vector *b)
{
  int j, k, size;
  double sum;
  
  if(!workspace->projwork.spect.defined)
    {
      raft_weight_spect(data, &workspace->projwork, att);
    }
  
  size = raft_scan_get_size(data);
  
  for(j = 0; j < size; j++)
    {
      for(k = 0; k < size; k++)
	{
	  sum = eval_anglesum_generalized(data, 
					  workspace->projwork.spect.weight,
					  p, 
					  j,
					  k);
	  
	  gsl_vector_set(b, k + j*size, sum);
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_backp_attenuated_xfct

  Compute the XFCT attenuated backprojection matrix.

  Input: 

  data - scan data . See <raft_scan_t>.  
  p - projection vector
  att - attenuation map
  attF - fluorescence attenuation map
	       
  Output:

  b - backprojection vector

  _Remark_:
  
  It works either with equispaced or nonequispaced rays.
  The same for the views.
  
  +====================================================+
*/

void *backp_pixelLoop_xfct(void *t)
{
  int w, j, k, size;
  parBackLoop_xfct_t *param;
  double sum, s;
  
  param = (parBackLoop_xfct_t *)t;
  
  size = raft_scan_get_size(param->data);
  
  for(w = param->colIndex[0]; w < MIN(param->colIndex[1], size*size); w++)
    {
      j = w/size;
      k = w%size;
      
      sum = eval_anglesum_generalized(param->data, param->workspace->projwork.xfct.weight, param->p, j, k);	  

      gsl_vector_set(param->b, k + j*size, sum);
    }  
  
  pthread_exit(NULL);
}

void raft_backp_attenuated_xfct(raft_scan_t *data, 
				raft_backp_t *workspace,
				gsl_vector *p,
				gsl_vector *attT,
				gsl_vector *attF,
				gsl_vector *b)
{
  pthread_t thread[BACKP_THREADS+1];
  pthread_attr_t attr;
  int size, n, rc;    
  parBackLoop_xfct_t param[BACKP_THREADS+1];  
  void *status;
  
  if(!workspace->projwork.xfct.defined)
    {
      raft_weight_xfct(data, &workspace->projwork, attT, attF);
    }

  /* Initialize and set thread detached attribute */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  size = raft_scan_get_size(data);
  
  int e = (int) floor((size*size)/BACKP_THREADS);;
  
  for(n = 0; n < BACKP_THREADS+1; n++) 
    {
      param[n].data        = data;
      param[n].workspace   = workspace;
      param[n].p           = p;
      param[n].b           = b; 
      param[n].nthread     = n;      
      
      param[n].colIndex[0] = e * n;
      param[n].colIndex[1] = (n+1) * e;
      
      rc = pthread_create(&thread[n], &attr, backp_pixelLoop_xfct, (void *)&param[n]);
    }
  
  /* Free attribute and wait for the other threads */
  pthread_attr_destroy(&attr);
  for(n = 0; n < BACKP_THREADS + 1; n++) 
    {
      rc = pthread_join(thread[n], &status);
    }  

  /*
  int j, k, size;
  double sum;
  
  if(!workspace->projwork.xfct.defined)
    raft_weight_xfct(data, &workspace->projwork, attT, attF);
  
  size = raft_scan_get_size(data);
      
  for(j = 0; j < size; j++)
    {
      for(k = 0; k < size; k++)
	{
	  sum = eval_anglesum_generalized(data, workspace->projwork.xfct.weight, p, j, k);
	  gsl_vector_set(b, k + j*size, sum);
	}
    }
  */
}


/*+====================================================+
  
  FUNCTION: raft_backp_generalized

  Compute the generalized backprojection matrix.

  Input: 

  data - scan data . See <raft_scan_t>.  
  p - projection vector
  weight - weight vector
	       
  Output:

  b - backprojection vector

  _Remark_:
  
  It works either with equispaced or nonequispaced rays.
  The same for the views.
   
  +====================================================+
*/

void raft_backp_generalized(raft_scan_t *data, 
			    raft_backp_t *workspace,
			    gsl_vector *p,
			    gsl_vector **weight,
			    gsl_vector *b)
{
  int j, k, size;
  double sum;
  
  size = raft_scan_get_size(data);
      
  for(j = 0; j < size; j++)
    {
      for(k = 0; k < size; k++)
	{
	  sum = eval_anglesum_generalized(data, 
					  weight,
					  p, 
					  j,
					  k);
	  
	  gsl_vector_set(b, k + j*size, sum);
	}
    }
}

/*######################################################
  Section: Inversion operators
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_backp_fbp

  Compute the (f)iltered (b)ack(p)rojection matrix.
  
  It applies the backprojection operator 
  (see <raft_backprojection>) on the filtered 
  projection matrix.
  
  Input: 

  data - scan data . See <raft_scan_t>.
  workspace - backprojection workspace 
  p - projection vector

  Ouput:

  f - filtered backprojection vector (reconstructed image)
  
  Return:

  RAFT_EDOM - Domain error. Nonequispaced rays
  RAFT_SUCCESS - Successfull operation.

  _Remark_:
  
  It only works for equispaced rays.
  
  +====================================================+
*/

int raft_backp_fbp(raft_scan_t *data, 
		   raft_backp_t *workspace,
		   gsl_vector *p,
		   gsl_vector *f)
{
  if(data->raysid == RAFT_REG)
    {
      filtered_projection(data, workspace, p, workspace->fourier.q);
      
      raft_backp(data, workspace->fourier.q, f);  

      return RAFT_SUCCESS;
    }
  else
    return RAFT_EDOM;
}

/*+====================================================+
  
  FUNCTION: raft_backp_fob

  Compute a (f)ilter (o)f the (b)ackprojection matrix.

  It applies the 2D inverse Fourier transform in the
  filtered backprojection matrix. 
  
  Input: 

  data - scan data . See <raft_scan_t>.
  workspace - backprojection workspace
  b - backprojection matrix

  Ouput:

  f - reconstructed matrix

  Return:

  RAFT_EDOM - Domain error. Nonequispaced rays
  RAFT_SUCCESS - Successfull operation.

  _Remark_:
  
  It only works for equispaced rays.
  
  +====================================================+
*/

int raft_backp_fob(raft_scan_t *data,
		   raft_backp_t *workspace,
		   gsl_vector *b,
		   gsl_vector *f)
{
  if(data->raysid == RAFT_REG)
    {
      fft2_complex(workspace->fob.fftREback, 
		   workspace->fob.fftIMback, 
		   b, 
		   workspace->fob.zero, 
		   workspace->fob.fft2work);
      
      gsl_vector_mul(workspace->fob.fftREback, workspace->fob.filter);
      gsl_vector_mul(workspace->fob.fftIMback, workspace->fob.filter);

      ifft2_complex(f, 
		    workspace->fob.zero, 
		    workspace->fob.fftREback, 
		    workspace->fob.fftIMback,
		    workspace->fob.fft2work);      
      
      return RAFT_SUCCESS;
    }
  else
    return RAFT_EDOM;
}


/*+====================================================+
  
  FUNCTION: raft_backp_novikov
  
  Novikov's inversion for SPECT.
  
  Input: 

  data - scan data . See <raft_scan_t>.
  workspace - backprojection workspace
  att - attenuation matrix
  p - projection matrix

  Ouput:

  f - SPECT reconstructed image (activity)
  
  Return:

  RAFT_EDOM - Domain error. Nonequispaced rays
  RAFT_SUCCESS - Successfull operation.

  _Remark_:
  
  It only works for equispaced rays.
       
  +====================================================+
*/

int raft_backp_novikov(raft_scan_t *data, 
		       raft_backp_t *workspace,
		       gsl_vector *att,
		       gsl_vector *p,
		       gsl_vector *f)
{
  if(data->raysid == RAFT_REG)
    {
      int i, k, size, nviews;
      double anglesum;
      
      size   = raft_scan_get_size(data);
      nviews = raft_scan_get_nviews(data);

      /* define workspace */

      define_novikov_workspace(data, workspace, att, p);
            
      /* inversion */
	
      for(i=0; i < size; i++)
	{
	  for(k=0; k<size; k++)
	    {
	      anglesum = eval_anglesum_inversion(data,
						 workspace->inversion.modproj,
						 workspace->projwork.spect.weight,
						 i, k);
	      
	      gsl_vector_set(f, k + i*size, anglesum);
	    }
	}
      
      return RAFT_SUCCESS;
    }
  else
    return RAFT_EDOM;
}

/*+====================================================+
  
  FUNCTION: raft_backp_tretiakMetz
  
  Tretiak & Metz inversion for ERT.
  
  ERT stands for the Exponential Radon Tranform

  Input: 

  data - scan data . See <raft_scan_t>.
  workspace - backprojection workspace
  att - attenuation value
  p - projection matrix

  Ouput:

  f - ERT reconstructed image (activity)
  
  Return:

  RAFT_EDOM - Domain error. Nonequispaced rays
  RAFT_SUCCESS - Successfull operation.

  _Remark_:
  
  It only works for equispaced rays.
       
  +====================================================+
*/

int raft_backp_tretiakMetz(raft_scan_t *data, 
			   raft_backp_t *workspace,
			   double *att,
			   gsl_vector *p,
			   gsl_vector *f)
{
  /*
  if(data->raysid == RAFT_REG)
    {
      int j, nrays, nviews;
      double Att, aux;
      gsl_vector_view v;
      
      Att = *att;
      
      gsl_vector_set_all(workspace->inversion.att, Att);
      
      nrays  = raft_scan_get_nrays(data);
      nviews = raft_scan_get_nviews(data);
      
      for(j=0; j<nrays; j++)
	{
	  aux = exp(-Att*gsl_vector_get(workspace->inversion.kernel,j));
	  gsl_vector_set(workspace->inversion.kernel, j, aux);
	}
      
      gsl_blas_dcopy(p, workspace->inversion.data);
      
      for(j=0; j<nviews; j++)
	{
	  v = gsl_vector_subvector(workspace->inversion.data, j*nrays, nrays);
	  gsl_vector_mul(&v.vector, workspace->inversion.kernel);
	}
      
      raft_backp_novikov(data, workspace, workspace->inversion.att, p, f);

      return RAFT_SUCCESS;
    }
  else
    return RAFT_EDOM;
  */
  
  if(data->raysid == RAFT_REG)
    { 
      int k, j, nrays, nviews;
      double dt, ww, wn, wc, dw, fun, Att;
      gsl_vector_view v, w, V, W;
      
      Att = *att;

      nrays = raft_scan_get_nrays(data);
      nviews = raft_scan_get_nviews(data);
      dt = raft_scan_get_raystep(data);
      
      wc = 1/(2*dt);
      wn = wc;
      dw = (2*wn)/nrays;
  
      for(k=0; k < nrays; k++)
	{
	  ww = -wn + k*dw;	  	  

	  if(fabs(ww) > fabs(Att) || fabs(ww)==fabs(Att))
	    fun = fabs(ww);
	  else
	    fun = 0;
	  
	  gsl_vector_set(workspace->fourier.fftREimpulseTretiak, k, fun);	  
	}
      
      fft_shift(workspace->fourier.fftREimpulseTretiak, nrays);

      for(j=0; j < nviews; j++)
	{
	  v = gsl_vector_subvector(p, j*nrays, nrays);
	  w = gsl_vector_subvector(workspace->fourier.q, j*nrays, nrays);

	  gsl_vector_set_all(workspace->fourier.signal, 0.0);
	  V = gsl_vector_subvector(workspace->fourier.signal, 0, nrays);
	  gsl_blas_dcopy(&v.vector, &V.vector);	  

	  fft_complex(workspace->fourier.fftRE1, 
		      workspace->fourier.fftIM1,
		      workspace->fourier.signal, 
		      workspace->fourier.zero,
		      workspace->fourier.fftwork);      
	  
	  gsl_vector_mul(workspace->fourier.fftRE1, workspace->fourier.fftREimpulseTretiak);
	  gsl_vector_mul(workspace->fourier.fftIM1, workspace->fourier.fftREimpulseTretiak);
	  
	  gsl_vector_set_all(workspace->fourier.signal, 0.0);

	  ifft_complex(workspace->fourier.impulse,
		       workspace->fourier.fftIM2, 
		       workspace->fourier.fftRE1, 
		       workspace->fourier.fftIM1,
		       workspace->fourier.fftwork);     	  
	  
	  W = gsl_vector_subvector(workspace->fourier.impulse, 0, nrays); 
	  gsl_blas_dcopy(&W.vector, &w.vector);	  
	}
      
      /* Att *= -1; */
      
      raft_backp_attenuated_ert(data, workspace, &Att, workspace->fourier.q, f); 
      
      return RAFT_SUCCESS;
    }
  else
    return RAFT_EDOM;
  
}

/*+====================================================+
  
  FUNCTION: raft_backp_invxfct
  
  XFCT inversion.

  input: 

  data - scan data. See <raft_scan_t>.
  workspace - backprojection workspace
  attT - transmission attenuation matrix
  attF - fluorescence attenuation matrix
  p - projection matrix

  Ouput:

  f - XFCT reconstructed image (density)
  
  Return:

  RAFT_EDOM - Domain error. Nonequispaced rays
  RAFT_SUCCESS - Successfull operation.

  _Remark_:
  
  It only works for equispaced rays.
  
  +====================================================+
*/

int raft_backp_invxfct(raft_scan_t *data, 
		       raft_backp_t *workspace,
		       gsl_vector *attT,
		       gsl_vector *attF,
		       gsl_vector *p,
		       gsl_vector *f)
{
  if(data->raysid == RAFT_REG)
    {
      int i, k, size;
      double anglesum;

      size   = raft_scan_get_size(data);
      
      /* define workspace */
      
      define_invxfct_workspace(data, workspace, attT, attF, p);
      
      /* inversion */
	
      for(i=0; i < size; i++)
	{
	  for(k=0; k<size; k++)
	    {
	      anglesum = eval_anglesum_inversion(data,
						 workspace->inversion.modproj,
						 workspace->projwork.xfct.weight,
						 i, k);

	      gsl_vector_set(f, k + i*size, anglesum);
	    }
	}
      
      return RAFT_SUCCESS;
    }
  else
    return RAFT_EDOM;
}

/*+====================================================+
  
  FUNCTION: raft_backp_invxfct_partial
  
  XFCT partial inversion.

  Only for the perpendicular fluorescence ray.
  
  Input: 

  data - scan data. See <raft_scan_t>.
  workspace - backprojection workspace
  attT - transmission attenuation matrix
  attF - fluorescence attenuation matrix
  p - projection matrix

  Ouput:

  f - XFCT reconstructed image (density)
  
  Return:

  RAFT_EDOM - Domain error. Nonequispaced rays
  RAFT_SUCCESS - Successfull operation.

  _Remark_:
  
  It only works for equispaced rays.
       
  +====================================================+
*/

int raft_backp_invxfct_partial(raft_scan_t *data, 
			       raft_backp_t *workspace,
			       gsl_vector *attT,
			       gsl_vector *attF,
			       gsl_vector *p,
			       gsl_vector *f)
{
  if(data->raysid == RAFT_REG)
    {
      int i, k, size;
      double anglesum, apert, iapert; 
      
      size   = raft_scan_get_size(data);
      
      /* define workspace */
      
      define_invxfctPartial_workspace(data, workspace, attT, attF, p, &apert);
      
      iapert = 1.0/apert;
      
      /* inversion */
	
      for(i=0; i < size; i++)
	{
	  for(k=0; k<size; k++)
	    {
	      anglesum = eval_anglesum_inversion(data,
						 workspace->inversion.modproj,
						 workspace->projwork.xfct.weightP,
						 i, k);
	      
	      anglesum *= iapert;

	      gsl_vector_set(f, k + i*size, anglesum);
	    }
	}
      
      return RAFT_SUCCESS;
    }
  else
    return RAFT_EDOM;
}


/*+====================================================+
  
  FUNCTION: raft_backp_invxfct_neumann
  
  XFCT inversion using Neumann series.

  Input: 

  data - scan data. See <raft_scan_t>.
  workspace - backprojection workspace
  attT - transmission attenuation matrix
  attF - fluorescence attenuation matrix
  p - projection matrix
  n - number of partial sums

  Ouput:

  f - XFCT reconstructed image (density)
  
  Return:

  RAFT_EDOM - Domain error. Nonequispaced rays
  RAFT_SUCCESS - Successfull operation.

  _Remark_:
  
  It only works for equispaced rays.
       
  +====================================================+
*/

int raft_backp_invxfct_neumann(raft_scan_t *data, 
			       raft_backp_t *workspace,
			       gsl_vector *attT,
			       gsl_vector *attF,
			       gsl_vector *p,
			       gsl_vector *f,
			       int n)
{
  int k, status;
  
  k = 0;
  
  gsl_vector_set_all(workspace->neumann.next, 0.0);

  do{
    
    raft_projection_radon_xfct(data, 
			       &workspace->projwork,
			       workspace->neumann.next,
			       attF,
			       attT,
			       workspace->neumann.rad);
    
    gsl_vector_sub(workspace->neumann.rad, p);
    gsl_vector_scale(workspace->neumann.rad, -1.0);

    status = raft_backp_invxfct_partial(data, 
					workspace,
					attT,
					attF,
					workspace->neumann.rad,
					workspace->neumann.direction);
    
    if(status!=RAFT_SUCCESS)
      {
	return RAFT_EDOM;
      }
    
    gsl_vector_add(workspace->neumann.next, workspace->neumann.direction);
    
    k++;
    
  }while(k < n);

  return RAFT_SUCCESS;
}
