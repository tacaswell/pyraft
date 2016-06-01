#include "raft_projection.h"
#include "raft_weight.h"
#include "projection.h"
#include "likelihood.h"
#include <gsl/gsl_math.h>
#include <pthread.h>

#define RADON_THREADS 10

/*######################################################
  Title: Projections
  
  Header - <raft/raft_projection.h>
  Type - <raft_proj_t>
  $Id: projection.c,v 1.37 2011-03-02 19:23:19 miqueles Exp $ - Last Update
  ####################################################*/

/*+====================================================+
  
  FUNCTION intersection_ray_pixel
  
  Determines whether the ray intersects a given pixel. 
  
  Input

  scan - scan data . See <raft_scan_t>.
  i - ray number
  j - pixel number 
  
  Return

  1 whenever the ray intersects the pixel and 0 otherwise.
  
  +====================================================+
*/

int intersection_ray_pixel(raft_scan_t *data,
			   int i, int j)
{
  int size, row, column, J, angle, ray;
  div_t D;
  double dt, t0, t, x, y, cost, sint;

  D = div(i, raft_scan_get_nrays(data));
  angle = D.quot;
  ray = D.rem;
  
  dt = raft_scan_get_raystep(data);
  t0 = raft_scan_get_ray(data,0);
  size = raft_scan_get_size(data);

  D = div(j, size);
  row = D.quot;
  column = D.rem;

  cost = gsl_vector_get(data->costheta, angle);
  sint = gsl_vector_get(data->sintheta, angle);
  
  x = raft_scan_get_x(data, column);
  y = raft_scan_get_y(data, size-1-row);

  t = x*cost + y*sint;
    
  /* J = interp_search(data->t, t); */
  
  J = interp_nearest(data->t, t, RAFT_REG);

  if(SIGN(J)>0 && J==ray)
    return 1;
  else
    return 0;
}


/*+====================================================+
  
  FUNCTION eval_radon_charf_pixel

  Radon transform of the characteristic 
  function with support at a given pixel. 
  
  Input

  scan - scan data . See <raft_scan_t>.
  i - ray number
  pixel - pixel number 
  
  Return

  Evaluated Radon transform.
  
  +====================================================+
*/

double eval_radon_charf_pixel(raft_scan_t *data,
			      int i,
			      int pixel)
{
  int angle, ray, row, column, size, J;
  div_t D;
  double t, x, y, dx, cost, sint, fcost, fsint, tol;
  
  size = raft_scan_get_size(data);

  D = div(i, raft_scan_get_nrays(data));
  angle = D.quot;
  ray = D.rem;
  
  tol = 1/sqrt(2);
  dx = raft_scan_get_x(data,1)-raft_scan_get_x(data,0);
  
  D = div(pixel, size);
  row = D.quot;
  column = D.rem;

  cost = gsl_vector_get(data->costheta, angle);
  sint = gsl_vector_get(data->sintheta, angle);
   
  fcost=fabs(cost);
  fsint=fabs(sint);
  
  x = raft_scan_get_x(data, column);
  y = raft_scan_get_y(data, size-1-row);
  
  t = x*cost + y*sint;

  J = interp_nearest(data->t, t, RAFT_REG);

  if(SIGN(J)>0 && J==ray)
    {
      if(fcost<tol)
	return dx/fsint;
      else
	return dx/fcost;    
    }
  else
    return 0.0;  
}

/*+====================================================+
  
  FUNCTION eval_rayintegral

  Compute a ray integral for the Radon transform. 
  
  Input

  scan - scan data. See <raft_scan_t>.
  workspace - projection workspace. See <raft_proj_t>.
  phantom - phantom vector. See <raft_phantom_t>.
  j - angle index
  k - ray index
    
  Output

  a - row of the projection matrix, corresponding to 
      ray with index (j,k).
  norma - squared norm of vector a
  dotp - dot product of vector 'a' with 'phantom'

  Return 

  Evaluated ray integral.
  
  +====================================================+
*/

double eval_rayintegral_(raft_scan_t *data,
			gsl_vector *phantom,
			int j,
			int k,
			gsl_vector *a,
			double *norma,
			double *dotp)
{
  int i, size, row, column, pixel, m;
  double sum, cost, sint, tant, x, y, t, dot;
  double fsint, fcost, tol, step, value, Delta;
  gsl_vector_view phant;  

  tol = 1/sqrt(2);
  step = raft_scan_get_x(data,1)-raft_scan_get_x(data,0);
  size = raft_scan_get_size(data);
    
  cost  = gsl_vector_get(data->costheta,j);
  sint  = gsl_vector_get(data->sintheta,j);
  tant  = gsl_vector_get(data->tantheta,j);
  fsint = fabs(sint);
  fcost = fabs(cost);
  
  t = raft_scan_get_ray(data,k);
  
  if(fcost < tol) 
    {
      m   = 0; 
      sum = 0;
      dot = 0;
      
      Delta = (step/fsint);

      for(i=0; i < size; i++)
	{
	  /* column of the phantom matrix */
	  phant = gsl_vector_subvector_with_stride(phantom, i, size, size);
	  
	  x = raft_scan_get_x(data,i);
	  y = t/sint - x/tant;
	  
	  row = interp_nearest(data->y, y, RAFT_REG);
	  
	  if(SIGN(row)>0)	  
	    {
	      sum += interp_linear_reverse(&phant.vector, data->y, row, y);

	      pixel = i + size*(size-1-row);	      
	      gsl_vector_set(a, pixel, Delta);
	      
	      m++;

	      dot += raft_phantom_get(phantom, size-1-row, i);
	    }
	}
      value = sum*Delta;
    } 
  else
    {
      m   = 0;
      sum = 0;
      dot = 0;
      
      Delta = (step/fcost);

      for(i=0; i < size; i++)
	{
	  /* row of the phantom matrix */
	  phant = gsl_vector_subvector(phantom, (size-1-i)*size, size);
	  
	  y = raft_scan_get_y(data,i);
	  x = t/cost - y*tant;
	  
          column = interp_nearest(data->x, x, RAFT_REG);
	  
	  if(SIGN(column)>0)	  
	    {
	      sum += interp_linear(&phant.vector, data->x, column, x);

	      pixel = column + size*(size-1-i);
	      gsl_vector_set(a, pixel, Delta);
	      
	      m++;

	      dot += raft_phantom_get(phantom, size-1-i, column);
	    }	  
	}
      value = sum*Delta;
    }

  *norma = m*SQR(Delta);  
  *dotp  = dot*Delta;

  return value;
}

/*+====================================================+
  
  FUNCTION eval_rayintegral_ert

  Compute a ray integral for the Exponential Radon transform (ERT). 
  
  Input

  scan - scan data. See <raft_scan_t>.
  workspace - projection workspace. See <raft_proj_t>.
  att - attenuation value.
  act - activity phantom vector. See <raft_phantom_t>.
  j - angle index
  k - ray index
    
  Output

  a - row of the projection matrix, corresponding to 
      ray with index (j,k).
  norma - squared norm of vector a
  dotp - dot product of vector 'a' with 'phantom'

  Return 

  Evaluated ray integral.
  
  +====================================================+
*/

double eval_rayintegral_ert(raft_scan_t *data,
			    double *att,
			    gsl_vector *act,
			    int j,
			    int k,
			    gsl_vector *a,
			    double *norma,
			    double *dotp)
{
  int i, size, row, column, pixel, m;
  double sum, cost, sint, tant, x, y, t, dot;
  double fsint, fcost, tol, step, value, Delta, E, s, prod;
  gsl_vector_view phant;  

  tol = 1/sqrt(2);
  step = raft_scan_get_x(data,1)-raft_scan_get_x(data,0);
  size = raft_scan_get_size(data);
    
  cost  = gsl_vector_get(data->costheta,j);
  sint  = gsl_vector_get(data->sintheta,j);
  tant  = gsl_vector_get(data->tantheta,j);
  fsint = fabs(sint);
  fcost = fabs(cost);
  
  t = raft_scan_get_ray(data,k);
  
  if(fcost < tol) 
    {
      m   = 0; 
      sum = 0;
      dot = 0;
      
      Delta = (step/fsint);

      for(i=0; i < size; i++)
	{
	  /* column of the phantom matrix */
	  phant = gsl_vector_subvector_with_stride(act, i, size, size);
	  
	  x = raft_scan_get_x(data,i);
	  y = t/sint - x/tant;
	  
	  s = (t*cost-x)/sint;	 
	  
	  row = interp_nearest(data->y, y, RAFT_REG);
	  
	  if(SIGN(row)>0)	  
	    {
	      prod = - s*(*att);	      
	      E = exp(prod); 	      
	      
	      sum += interp_linear_reverse(&phant.vector, data->y, row, y) * E;

	      pixel = i + size*(size-1-row);	      
	      gsl_vector_set(a, pixel, Delta*E);
	      
	      m++;

	      dot += raft_phantom_get(act, size-1-row, i);
	    }
	}
      value = sum*Delta;
    } 
  else
    {
      m   = 0;
      sum = 0;
      dot = 0;
      
      Delta = (step/fcost);

      for(i=0; i < size; i++)
	{
	  /* row of the phantom matrix */
	  phant = gsl_vector_subvector(act, (size-1-i)*size, size);
	  
	  y = raft_scan_get_y(data,i);
	  x = t/cost - y*tant;
	  
	  s = (y - t*sint)/cost;	  

          column = interp_nearest(data->x, x, RAFT_REG);
	  
	  if(SIGN(column)>0)	  
	    {
	      prod = - s*(*att);	      
	      E = exp(prod); 	          
	           
	      sum += interp_linear(&phant.vector, data->x, column, x) * E;
	      
	      pixel = column + size*(size-1-i);
	      gsl_vector_set(a, pixel, Delta*E);
	      
	      m++;
	      
	      dot += raft_phantom_get(act, size-1-i, column);
	    }	  
	}
      value = sum*Delta;
    }

  *norma = m*SQR(Delta);  
  *dotp  = dot*Delta;

  return value;
}

/*+====================================================+
  
  FUNCTION eval_rayintegral_pet

  Compute a ray integral for the PET Radon transform. 
  
  Input

  scan - scan data. See <raft_scan_t>.
  proj - projection workspace. See <raft_proj_t>.
  phantom - phantom vector
  j - angle index
  k - ray index
    
  Output

  a - row of the projection matrix, corresponding to 
      ray with index (j,k).
  norma - squared norm of vector a
  dotp - dot product of vector 'a' with 'phantom'

  Return 

  Evaluated ray integral.
  
  +====================================================+
*/

double eval_rayintegral_pet(raft_scan_t *data,
			    raft_proj_t *proj,
			    gsl_vector *phantom,
			    int j,
			    int k,
			    gsl_vector *a,
			    double *norma,
			    double *dotp)
{
  int i, size, row, column, pixel, nrays;
  double sum, cost, sint, tant, x, y, t, weight, dw;
  double fsint, fcost, tol, step, value, Delta;
  gsl_vector_view phant;  

  tol = 1/sqrt(2);
  step = raft_scan_get_x(data,1)-raft_scan_get_x(data,0);
  size = raft_scan_get_size(data);
  nrays= raft_scan_get_nrays(data);

  cost  = gsl_vector_get(data->costheta,j);
  sint  = gsl_vector_get(data->sintheta,j);
  tant  = gsl_vector_get(data->tantheta,j);
  fsint = fabs(sint);
  fcost = fabs(cost);
  
  t = raft_scan_get_ray(data,k);
  
  weight = gsl_vector_get(proj->pet.weight, k + nrays*j);
  
  if(fcost < tol) 
    {
      *norma = 0;
      *dotp  = 0;
      sum    = 0;
            
      Delta = (step/fsint);

      dw = (Delta*weight);

      for(i=0; i < size; i++)
	{
	  /* column of the phantom matrix */
	  phant = gsl_vector_subvector_with_stride(phantom, i, size, size);
	  
	  x = raft_scan_get_x(data,i);
	  y = t/sint - x/tant;
	  
	  row = interp_nearest(data->y, y, RAFT_REG);
	  
	  if(SIGN(row)>0)	  
	    {
	      sum += interp_linear_reverse(&phant.vector, data->y, row, y);

	      pixel = i + size*(size-1-row);	      
	      gsl_vector_set(a, pixel, dw);
	      
	      (*norma) += SQR(dw);

	      (*dotp) += dw * raft_phantom_get(phantom, size-1-row, i);
	    }
	}
      value = sum*Delta;
    } 
  else
    {
      *norma = 0;
      *dotp  = 0;
      sum    = 0;
      
      Delta = (step/fcost);
      
      dw = (Delta*weight);

      for(i=0; i < size; i++)
	{
	  /* row of the phantom matrix */
	  phant = gsl_vector_subvector(phantom, (size-1-i)*size, size);
	  
	  y = raft_scan_get_y(data,i);
	  x = t/cost - y*tant;
	  
          column = interp_nearest(data->x, x, RAFT_REG);
	  
	  if(SIGN(column)>0)	  
	    {
	      sum += interp_linear(&phant.vector, data->x, column, x);

	      pixel = column + size*(size-1-i);
	      gsl_vector_set(a, pixel, dw);
	      
	      (*norma) += SQR(dw);

	      (*dotp) += dw * raft_phantom_get(phantom, size-1-i, column);
	    }	  
	}
      value = sum*Delta;
    }

  return value;
}

/*+====================================================+
  
  FUNCTION eval_rayintegral_spect

  Compute a ray integral for the SPECT Radon transform. 
  
  Input

  scan - scan data. See <raft_scan_t>.
  proj - projection workspace. See <raft_proj_t>.
  phantom - phantom vector
  j - angle index
  k - ray index
    
  Output

  a - row of the projection matrix, corresponding to 
      ray with index (j,k).
  norma - squared norm of vector a
  dotp - dot product of vector 'a' with 'phantom'

  Return 

  Evaluated ray integral.
  
  +====================================================+
*/

double eval_rayintegral_spect(raft_scan_t *data,
			      raft_proj_t *proj,
			      gsl_vector *phantom,
			      int j,
			      int k,
			      gsl_vector *a,
			      double *norma,
			      double *dotp)
{
  int i, size, row, column, pixel, nrays;
  double sum, cost, sint, tant, x, y, t, dw;
  double fsint, fcost, tol, step, value, Delta;
  gsl_vector_view phant; 
  
  tol = 1/sqrt(2);
  step = raft_scan_get_x(data,1)-raft_scan_get_x(data,0);
  size = raft_scan_get_size(data);
  nrays= raft_scan_get_nrays(data);

  cost  = gsl_vector_get(data->costheta,j);
  sint  = gsl_vector_get(data->sintheta,j);
  tant  = gsl_vector_get(data->tantheta,j);
  fsint = fabs(sint);
  fcost = fabs(cost);
  
  t = raft_scan_get_ray(data,k);
  
  gsl_blas_dcopy(phantom, proj->spect.integrand);
  
  gsl_vector_mul(proj->spect.integrand, proj->spect.weight[j]);
  
  if(fcost < tol) 
    {
      *norma = 0;
      *dotp  = 0;
      sum    = 0;
            
      Delta = (step/fsint);

      for(i=0; i < size; i++)
	{
	  /* column of the phantom matrix */
	  phant = gsl_vector_subvector_with_stride(proj->spect.integrand, i, size, size);
	  
	  x = raft_scan_get_x(data,i);
	  y = t/sint - x/tant;
	  
	  row = interp_nearest(data->y, y, RAFT_REG);
	  
	  if(SIGN(row)>0)	  
	    {
	      sum += interp_linear_reverse(&phant.vector, data->y, row, y);
	      
	      pixel = i + size*(size-1-row);	      
	      
	      dw = Delta * raft_phantom_get(proj->spect.weight[j], size-1-row, i);
	      
	      gsl_vector_set(a, pixel, dw);
	      
	      (*norma) += SQR(dw);
	      
	      (*dotp) += dw * raft_phantom_get(phantom, size-1-row, i);
	    }	  
	}
      value = sum*Delta;
    } 
  else
    {
      *norma = 0;
      *dotp  = 0;
      sum    = 0;
      
      Delta = (step/fcost);
      
      for(i=0; i < size; i++)
	{
	  /* row of the phantom matrix */
	  phant = gsl_vector_subvector(proj->spect.integrand, (size-1-i)*size, size);
	  
	  y = raft_scan_get_y(data,i);
	  x = t/cost - y*tant;
	  
          column = interp_nearest(data->x, x, RAFT_REG);
	  
	  if(SIGN(column)>0)	  
	    {
	      sum += interp_linear(&phant.vector, data->x, column, x);

	      pixel = column + size*(size-1-i);

	      dw = Delta * raft_phantom_get(proj->spect.weight[j], size-1-i, column);
	      
	      gsl_vector_set(a, pixel, dw);
	      
	      (*norma) += SQR(dw);
	      
	      (*dotp) += dw * raft_phantom_get(phantom, size-1-i, column);
	    }	  	    
	}
      value = sum*Delta;
    }

  return value;
}

/*+====================================================+
  
  FUNCTION eval_rayintegral_xfct

  Compute a ray integral for the XFCT Radon transform. 
  
  Input

  scan - scan data. See <raft_scan_t>.
  proj - projection workspace. See <raft_proj_t>.
  phantom - phantom vector
  j - angle index
  k - ray index
    
  Output

  a - row of the projection matrix, corresponding to 
      ray with index (j,k).
  norma - squared norm of vector a
  dotp - dot product of vector 'a' with 'phantom'

  Return 

  Evaluated ray integral.
  
  +====================================================+
*/

double eval_rayintegral_xfct(raft_scan_t *data,
			     raft_proj_t *proj,
			     gsl_vector *phantom,
			     int j,
			     int k,
			     gsl_vector *a,
			     double *norma,
			     double *dotp)
{
  int i, size, row, column, pixel, nrays;
  double sum, cost, sint, tant, x, y, t, dw;
  double fsint, fcost, tol, step, value, Delta;
  gsl_vector_view phant;  

  tol = 1/sqrt(2);
  step = raft_scan_get_x(data,1)-raft_scan_get_x(data,0);
  size = raft_scan_get_size(data);
  nrays= raft_scan_get_nrays(data);

  cost  = gsl_vector_get(data->costheta,j);
  sint  = gsl_vector_get(data->sintheta,j);
  tant  = gsl_vector_get(data->tantheta,j);
  fsint = fabs(sint);
  fcost = fabs(cost);
  
  t = raft_scan_get_ray(data,k);

  gsl_blas_dcopy(phantom, proj->xfct.integrand[0]);
  
  gsl_vector_mul(proj->xfct.integrand[0], proj->xfct.weight[j]);
  
  if(fcost < tol) 
    {
      *norma = 0;
      *dotp  = 0;
      sum    = 0;
            
      Delta = (step/fsint);

      for(i=0; i < size; i++)
	{
	  /* column of the phantom matrix */
	  phant = gsl_vector_subvector_with_stride(proj->xfct.integrand[0], i, size, size);
	    
	  x = raft_scan_get_x(data,i);
	  y = t/sint - x/tant;
	  
	  row = interp_nearest(data->y, y, RAFT_REG);
	  
	  if(SIGN(row)>0)	  
	    {
	      sum += interp_linear_reverse(&phant.vector, data->y, row, y);
	      
	      pixel = i + size*(size-1-row);	      
	      
	      dw = Delta * raft_phantom_get(proj->xfct.weight[j], size-1-row, i);
	      
	      gsl_vector_set(a, pixel, dw);
	      
	      (*norma) += SQR(dw);
	      
	      (*dotp) += dw * raft_phantom_get(phantom, size-1-row, i);
	    }
	}
      value = sum*Delta;
    } 
  else
    {
      *norma = 0;
      *dotp  = 0;
      sum    = 0;
      
      Delta = (step/fcost);
      
      for(i=0; i < size; i++)
	{
	  /* row of the phantom matrix */
	  phant = gsl_vector_subvector(proj->xfct.integrand[0], (size-1-i)*size, size);
	  
	  y = raft_scan_get_y(data,i);
	  x = t/cost - y*tant;
	  
          column = interp_nearest(data->x, x, RAFT_REG);
	  
	  if(SIGN(column)>0)	  
	    {
	      sum += interp_linear(&phant.vector, data->x, column, x);

	      pixel = column + size*(size-1-i);

	      dw = Delta * raft_phantom_get(proj->xfct.weight[j], size-1-i, column);
	      
	      gsl_vector_set(a, pixel, dw);
	      
	      (*norma) += SQR(dw);
	      
	      (*dotp) += dw * raft_phantom_get(phantom, size-1-i, column);
	    }	  
	}
      value = sum*Delta;
    }

  return value;
}

double eval_rayintegral_xfct_bythread(int nthread,
				      raft_scan_t *data,
				      raft_proj_t *proj,
				      gsl_vector *phantom,
				      int j,
				      int k,
				      gsl_vector *a,
				      double *norma,
				      double *dotp)
{
  int i, size, row, column, pixel, nrays;
  double sum, cost, sint, tant, x, y, t, dw;
  double fsint, fcost, tol, step, value, Delta;
  gsl_vector_view phant;  

  tol = 1/sqrt(2);
  step = raft_scan_get_x(data,1)-raft_scan_get_x(data,0);
  size = raft_scan_get_size(data);
  nrays= raft_scan_get_nrays(data);

  cost  = gsl_vector_get(data->costheta,j);
  sint  = gsl_vector_get(data->sintheta,j);
  tant  = gsl_vector_get(data->tantheta,j);
  fsint = fabs(sint);
  fcost = fabs(cost);
  
  t = raft_scan_get_ray(data,k);

  gsl_blas_dcopy(phantom, proj->xfct.integrand[nthread]);
  
  gsl_vector_mul(proj->xfct.integrand[nthread], proj->xfct.weight[j]);
  
  if(fcost < tol) 
    {
      *norma = 0;
      *dotp  = 0;
      sum    = 0;
            
      Delta = (step/fsint);

      for(i=0; i < size; i++)
	{
	  /* column of the phantom matrix */
	  phant = gsl_vector_subvector_with_stride(proj->xfct.integrand[nthread], i, size, size);
	    
	  x = raft_scan_get_x(data,i);
	  y = t/sint - x/tant;
	  
	  row = interp_nearest(data->y, y, RAFT_REG);
	  
	  if(SIGN(row)>0)	  
	    {
	      sum += interp_linear_reverse(&phant.vector, data->y, row, y);
	      
	      pixel = i + size*(size-1-row);	      
	      
	      dw = Delta * raft_phantom_get(proj->xfct.weight[j], size-1-row, i);
	      
	      gsl_vector_set(a, pixel, dw);
	      
	      (*norma) += SQR(dw);
	      
	      (*dotp) += dw * raft_phantom_get(phantom, size-1-row, i);
	    }
	}
      value = sum*Delta;
    } 
  else
    {
      *norma = 0;
      *dotp  = 0;
      sum    = 0;
      
      Delta = (step/fcost);
      
      for(i=0; i < size; i++)
	{
	  /* row of the phantom matrix */
	  phant = gsl_vector_subvector(proj->xfct.integrand[nthread], (size-1-i)*size, size);
	  
	  y = raft_scan_get_y(data,i);
	  x = t/cost - y*tant;
	  
          column = interp_nearest(data->x, x, RAFT_REG);
	  
	  if(SIGN(column)>0)	  
	    {
	      sum += interp_linear(&phant.vector, data->x, column, x);

	      pixel = column + size*(size-1-i);

	      dw = Delta * raft_phantom_get(proj->xfct.weight[j], size-1-i, column);
	      
	      gsl_vector_set(a, pixel, dw);
	      
	      (*norma) += SQR(dw);
	      
	      (*dotp) += dw * raft_phantom_get(phantom, size-1-i, column);
	    }	  
	}
      value = sum*Delta;
    }

  return value;
}


/*+====================================================+
  
  FUNCTION eval_rayintegral_generalized

  Compute a ray integral for the generalized Radon 
  transform. 
  
  Input

  scan - scan data. See <raft_scan_t>.
  proj - projection workspace. See <raft_proj_t>.
  weight - weight matrix
  phantom - phantom vector
  j - angle index
  k - ray index
    
  Output

  a - row of the projection matrix, corresponding to 
      ray with index (j,k).
  norma - squared norm of vector a
  dotp - dot product of vector 'a' with 'phantom'

  Return 

  Evaluated ray integral.
  
  +====================================================+
*/

double eval_rayintegral_generalized(raft_scan_t *data,
				    raft_proj_t *proj,
				    gsl_vector **weight,
				    gsl_vector *phantom,
				    int j,
				    int k,
				    gsl_vector *a,
				    double *norma,
				    double *dotp)
{
  int i, size, row, column, pixel, nrays;
  double sum, cost, sint, tant, x, y, t, dw;
  double fsint, fcost, tol, step, value, Delta;
  gsl_vector_view phant;  

  tol = 1/sqrt(2);
  step = raft_scan_get_x(data,1)-raft_scan_get_x(data,0);
  size = raft_scan_get_size(data);
  nrays= raft_scan_get_nrays(data);

  cost  = gsl_vector_get(data->costheta,j);
  sint  = gsl_vector_get(data->sintheta,j);
  tant  = gsl_vector_get(data->tantheta,j);
  fsint = fabs(sint);
  fcost = fabs(cost);
  
  t = raft_scan_get_ray(data,k);
  
  gsl_blas_dcopy(phantom, proj->generalized.integrand);
  
  gsl_vector_mul(proj->generalized.integrand, weight[j]);
  
  if(fcost < tol) 
    {
      *norma = 0;
      *dotp  = 0;
      sum    = 0;
            
      Delta = (step/fsint);

      for(i=0; i < size; i++)
	{
	  /* column of the phantom matrix */
	  phant = gsl_vector_subvector_with_stride(proj->generalized.integrand, 
						   i, size, size);
	  
	  x = raft_scan_get_x(data,i);
	  y = t/sint - x/tant;
	  
	  row = interp_nearest(data->y, y, RAFT_REG);
	  
	  if(SIGN(row)>0)	  
	    {
	      sum += interp_linear_reverse(&phant.vector, data->y, row, y);
	      
	      pixel = i + size*(size-1-row);	      
	      
	      dw = Delta * raft_phantom_get(weight[j], size-1-row, i);
	      
	      gsl_vector_set(a, pixel, dw);
	      
	      (*norma) += SQR(dw);
	      
	      (*dotp) += dw * raft_phantom_get(phantom, size-1-row, i);
	    }
	}
      value = sum*Delta;
    } 
  else
    {
      *norma = 0;
      *dotp  = 0;
      sum    = 0;
      
      Delta = (step/fcost);
      
      for(i=0; i < size; i++)
	{
	  /* row of the phantom matrix */
	  phant = gsl_vector_subvector(proj->generalized.integrand, 
				       (size-1-i)*size, size);
	  
	  y = raft_scan_get_y(data,i);
	  x = t/cost - y*tant;
	  
          column = interp_nearest(data->x, x, RAFT_REG);
	  
	  if(SIGN(column)>0)	  
	    {
	      sum += interp_linear(&phant.vector, data->x, column, x);

	      pixel = column + size*(size-1-i);

	      dw = Delta * raft_phantom_get(weight[j], size-1-i, column);
	      
	      gsl_vector_set(a, pixel, dw);
	      
	      (*norma) += SQR(dw);
	      
	      (*dotp) += dw * raft_phantom_get(phantom, size-1-i, column);
	    }	  
	}
      value = sum*Delta;
    }

  return value;
}

/*+====================================================+
  
  FUNCTION  eval_rayintegral_dbt
  
  Compute the ray integral for divergent beam transform 
  
  Input  

  scan - scan data . See <raft_scan_t>.
  att - attenuation phantom vector
  j - angle index
  xi - x-coordinate index (column)
  yi - y-coordinate index (row)
    
  Return 

  Evaluated ray integral, i.e, the integral of the 
  attenuation function starting at a given point until 
  reach the detector.   
  
  +====================================================+
*/

double eval_rayintegral_dbt(raft_scan_t *data, 
			    gsl_vector *att,
			    int j,
			    int xi,
			    int yi)
{
  int i, begin, end, size, row, column;
  double sum, cost, sint, tant, x, y, X, Y;
  double fsint, fcost, tol, step, value;
  gsl_vector_view atten;  
  
  tol = 1/sqrt(2);
  step = raft_scan_get_x(data,1)-raft_scan_get_x(data,0);
  size = raft_scan_get_size(data);
  
  cost  = gsl_vector_get(data->costheta, j);
  sint  = gsl_vector_get(data->sintheta, j);
  tant  = gsl_vector_get(data->tantheta, j);
  fsint = fabs(sint);
  fcost = fabs(cost);
  
  x = raft_scan_get_x(data, xi);
  y = raft_scan_get_y(data, size-1-yi);
  
  if(fcost<tol)
    {
      if(SIGN(sint)<0)
	{
	  begin = xi;
	  end   = size;
	}
      else
	{
	  begin = 0;
	  end   = xi+1;
	}
      
      sum = 0;

      for(i=begin; i < end; i++)
	{
	  /* column of the attenuation phantom matrix*/
	  atten = gsl_vector_subvector_with_stride(att, i, size, size);
	  
	  X = raft_scan_get_x(data,i);
	  Y = y + (x-X)/tant;
	  
	  row = interp_nearest(data->y, Y, RAFT_REG);
	  
	  if(SIGN(row)>0)	  
	    sum += interp_linear_reverse(&atten.vector, data->y, row, Y);
	}
      
      value = (sum*step)/fsint;
    } 
  else
    {
      if(SIGN(cost)>0) 
	{
	  begin = 0;
	  end   = yi+1;
	}
      else
	{
	  begin = yi;
	  end   = size;
	}
      
      sum = 0;
      
      for(i=begin; i < end; i++)
	{
	  /* row of the activity phantom matrix*/
	  atten = gsl_vector_subvector(att, i*size, size);	  
	  
	  Y = raft_scan_get_y(data,size-1-i);
	  X = x + (y-Y)*tant;
	  
	  column = interp_nearest(data->x, X, RAFT_REG);
	  
	  if(SIGN(column)>0)	  
	    sum += interp_linear(&atten.vector, data->x, column, X);
	}

      value = (sum*step)/fcost;
    }
  
  return value;
}


/*+====================================================+
  
  FUNCTION  eval_ray_projection

  Evaluate the projection for a given ray. 
  
  Input  

  scan - scan data . See <raft_scan_t>.
  phantom - phantom vector
  i - ray number
   
  Output

  a - row of the projection matrix, corresponding to 
      ray with index (j,k).
  norma - squared norm of vector a
  dotp - dot product of vector 'a' with 'phantom'

  Return 

  Evaluated integral.
  
  +====================================================+
*/

double eval_ray_projection(raft_scan_t *data,
			   gsl_vector *phantom,
			   int i,
			   gsl_vector *a,
			   double *norma,
			   double *dotp)
{
  int angle, ray;
  div_t D;

  D = div(i, raft_scan_get_nrays(data));
  angle = D.quot;
  ray = D.rem;
  
  return eval_rayintegral_(data, phantom, angle, 
			  ray, a, norma, dotp);
}

/*+====================================================+
  
  FUNCTION  eval_ray_projection_ert

  Evaluate the projection for a given ray in ERT.
  
  ERT stands for the Exponential Radon Transform.
  
  Input  

  scan - scan data . See <raft_scan_t>.
  att - attenuation value
  phantom - phantom vector
  i - ray number
   
  Output

  a - row of the projection matrix, corresponding to 
      ray with index (j,k).
  norma - squared norm of vector a
  dotp - dot product of vector 'a' with 'phantom'

  Return 

  Evaluated integral.
  
  +====================================================+
*/

double eval_ray_projection_ert(raft_scan_t *data,
			       double *att,
			       gsl_vector *phantom,
			       int i,
			       gsl_vector *a,
			       double *norma,
			       double *dotp)
{
  int angle, ray;
  div_t D;

  D = div(i, raft_scan_get_nrays(data));
  angle = D.quot;
  ray = D.rem;
  
  return eval_rayintegral_ert(data, att, phantom, angle, 
			      ray, a, norma, dotp);
}


/*+====================================================+
  
  FUNCTION  eval_ray_projection_pet

  Evaluate the projection for a given ray in PET. 
  
  Input  

  scan - scan data . See <raft_scan_t>.
  workspace - projection workspace (see <raft_proj_t>)
  act - activity map
  att - attenuation map  
  i - ray number
  
  Output
 
  a - row of the projection matrix, corresponding to 
      ray with index i.
  norma - squared norm of vector a
  dotp - dot product of vector 'a' with 'phantom'

  Return 

  Evaluated integral.
  
  +====================================================+
*/

double eval_ray_projection_pet(raft_scan_t *data,
			       raft_proj_t *workspace,
			       gsl_vector *act,
			       gsl_vector *att,
			       int i,
			       gsl_vector *a,
			       double *norma,
			       double *dotp)
{
  double radon;
  int angle, ray;
  div_t D;

  D = div(i, raft_scan_get_nrays(data));
  angle = D.quot;
  ray = D.rem;
  
  if(!workspace->pet.defined)
    {
      raft_weight_pet(data, workspace, att);
    }

  radon = eval_rayintegral_pet(data, workspace, act, 
			       angle, ray, a, norma, dotp);
  
  return radon;
}

/*+====================================================+
  
  FUNCTION  eval_ray_projection_spect

  Evaluate the projection for a given ray in SPECT. 
  
  Input  

  scan - scan data . See <raft_scan_t>.
  workspace - projection workspace (see <raft_proj_t>)
  act - activity map
  att - attenuation map  
  i - ray number

  Output
 
  a - row of the projection matrix, corresponding to 
      ray with index i.
  norma - squared norm of vector a
  dotp - dot product of vector 'a' with 'phantom'
   
  Return 

  Evaluated integral.
  
  +====================================================+
*/

double eval_ray_projection_spect(raft_scan_t *data,
				 raft_proj_t *workspace,
				 gsl_vector *act,
				 gsl_vector *att,
				 int i,
				 gsl_vector *a,
				 double *norma,
				 double *dotp)
{
  int angle, ray;
  div_t D;
  double radon;
  
  D = div(i, raft_scan_get_nrays(data));
  angle = D.quot;
  ray = D.rem;
  
  if(!workspace->spect.defined)
    {
      raft_weight_spect_view(data,
			     workspace,
			     att,
			     angle);

      /* raft_weight_spect(data, workspace, att);  */
    }
  
  radon = eval_rayintegral_spect(data, workspace, act, 
				 angle, ray, a, norma, dotp);

  return radon;
}

/*+====================================================+
  
  FUNCTION  eval_ray_projection_xfct

  Evaluate the projection for a given ray in xfct. 
  
  Input  

  scan - scan data . See <raft_scan_t>.
  workspace - projection workspace (see <raft_proj_t>)
  act - activity map
  attT - transmission attenuation map  
  attF - XFCT attenuation map
  i - ray number
   
  Output
 
  a - row of the projection matrix, corresponding to 
      ray with index i.
  norma - squared norm of vector a
  dotp - dot product of vector 'a' with 'phantom'

  Return 

  Evaluated integral.   
  
  +====================================================+
*/

double eval_ray_projection_xfct(raft_scan_t *data,
				raft_proj_t *workspace,
				gsl_vector *act,
				gsl_vector *attT,
				gsl_vector *attF,
				int i,
				gsl_vector *a,
				double *norma,
				double *dotp)
{
  int angle, ray;
  div_t D;
  double radon;
  
  D = div(i, raft_scan_get_nrays(data));
  angle = D.quot;
  ray = D.rem;
  
  if(!workspace->xfct.defined)
    {
      raft_weight_xfct_view(data,
			    workspace,
			    attT,
			    attF,
			    angle);
      
      /* raft_weight_xfct(data, workspace, attT, attF); */
    }

  radon = eval_rayintegral_xfct(data, workspace, act, 
				angle, ray, a, norma, dotp);

  return radon;
}

/*+====================================================+
  
  FUNCTION  projection_workspace_memcpy

  Copy memory area from one projection workspace to 
  another. 
  
  Input  

  sworkspace - source projection workspace
  
  Output 
  
  dworkspace - destiny projection workspace  
  
  +====================================================+
*/

void projection_workspace_memcpy(raft_proj_t *dworkspace,
				 raft_proj_t *sworkspace)
{
  int j, nviews;

  nviews = sworkspace->nviews;

  dworkspace->npixels = sworkspace->npixels;
  dworkspace->ndata   = sworkspace->ndata;
  dworkspace->nrays   = sworkspace->nrays;
  dworkspace->nviews  = sworkspace->nviews;
  dworkspace->size    = sworkspace->size;
  
  dworkspace->pet.defined = sworkspace->pet.defined;
  
  gsl_blas_dcopy(sworkspace->pet.Radon,dworkspace->pet.Radon);
  gsl_blas_dcopy(sworkspace->pet.ExpRadon,dworkspace->pet.ExpRadon);
  gsl_blas_dcopy(sworkspace->pet.integrand,dworkspace->pet.integrand);
  gsl_blas_dcopy(sworkspace->pet.weight,dworkspace->pet.weight);
  
  dworkspace->spect.defined = sworkspace->spect.defined;
  
  gsl_blas_dcopy(sworkspace->spect.integrand,
		 dworkspace->spect.integrand);
  
  dworkspace->xfct.defined = sworkspace->xfct.defined;
  
  gsl_blas_dcopy(sworkspace->xfct.ones, 
		 dworkspace->xfct.ones);
  gsl_blas_dcopy(sworkspace->xfct.integrand[0], 
		 dworkspace->xfct.integrand[0]);  
  
  for(j = 0; j < nviews; j++)
    {
      gsl_blas_dcopy(sworkspace->dbt.Dbt[j], dworkspace->dbt.Dbt[j]);
      gsl_blas_dcopy(sworkspace->dbt.ExpDbt[j], dworkspace->dbt.ExpDbt[j]);

      gsl_blas_dcopy(sworkspace->spect.weight[j], dworkspace->spect.weight[j]);
      gsl_blas_dcopy(sworkspace->spect.Dbt[j], dworkspace->spect.Dbt[j]);

      gsl_blas_dcopy(sworkspace->xfct.DbtT[j], dworkspace->xfct.DbtT[j]);
      gsl_blas_dcopy(sworkspace->xfct.DbtF[j], dworkspace->xfct.DbtF[j]);
      gsl_blas_dcopy(sworkspace->xfct.ExpDbtT[j], dworkspace->xfct.ExpDbtT[j]);
      gsl_blas_dcopy(sworkspace->xfct.ExpDbtF[j], dworkspace->xfct.ExpDbtF[j]);
      gsl_blas_dcopy(sworkspace->xfct.weight[j], dworkspace->xfct.weight[j]);      
    }
}

/******************************************************/
/*                      Public                        */
/******************************************************/

/*######################################################
  Section: Allocation
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_projection_workspace_alloc
  
  Workspace allocation for projection procedures.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
   
  Output:
  
  workspace - projection workspace
  
  +====================================================+
*/

void
raft_projection_workspace_alloc(raft_scan_t *data,
				raft_proj_t *workspace)
{
  int j;
  int ndata, npixels, nviews, nrays, size;

  ndata   = raft_scan_get_ndata(data);
  npixels = raft_scan_get_npixels(data);
  nviews  = raft_scan_get_nviews(data);
  nrays   = raft_scan_get_nrays(data);
  size    = raft_scan_get_size(data);
  
  workspace->likelihood.lograd = gsl_vector_alloc(ndata);
  workspace->likelihood.invrad = gsl_vector_alloc(ndata);
  workspace->likelihood.exprad = gsl_vector_alloc(ndata);
  workspace->likelihood.exprad2= gsl_vector_alloc(ndata);
  workspace->likelihood.radon  = gsl_vector_alloc(ndata);
  workspace->likelihood.ones   = gsl_vector_alloc(ndata);
  
  workspace->likelihood.ratio = 1;
  workspace->likelihood.eps   = ZERO;
  gsl_vector_set_all(workspace->likelihood.ones, 1.0);
  
  workspace->generalized.integrand = gsl_vector_alloc(npixels);

  workspace->dbt.Dbt    = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  workspace->dbt.ExpDbt = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  
  workspace->xfct.weight  = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  workspace->xfct.weightP = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  workspace->xfct.weightF = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  workspace->xfct.LogWeightF = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  workspace->xfct.LogWeightP = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);

  workspace->xfct.DbtT    = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  workspace->xfct.DbtF    = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  workspace->xfct.ExpDbtT = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  workspace->xfct.ExpDbtF = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews); 
  
  workspace->spect.weight = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  workspace->spect.Dbt    = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);

  for(j = 0; j < nviews; j++)
    {
      workspace->dbt.Dbt[j]    = gsl_vector_alloc(npixels);
      workspace->dbt.ExpDbt[j] = gsl_vector_alloc(npixels);
      
      workspace->xfct.weight[j]  = gsl_vector_alloc(npixels);
      workspace->xfct.weightP[j] = gsl_vector_alloc(npixels);
      workspace->xfct.weightF[j] = gsl_vector_alloc(npixels);
      workspace->xfct.LogWeightF[j] = gsl_vector_alloc(npixels);
      workspace->xfct.LogWeightP[j] = gsl_vector_alloc(npixels);
      
      workspace->xfct.DbtT[j]    = gsl_vector_alloc(npixels);
      workspace->xfct.DbtF[j]    = gsl_vector_alloc(npixels);
      workspace->xfct.ExpDbtT[j] = gsl_vector_alloc(npixels);
      workspace->xfct.ExpDbtF[j] = gsl_vector_alloc(npixels);

      workspace->spect.weight[j] = gsl_vector_alloc(npixels);
      workspace->spect.Dbt[j]    = gsl_vector_alloc(npixels);
    }
  
  workspace->pet.weight    = gsl_vector_alloc(ndata);
  workspace->pet.integrand = gsl_vector_alloc(npixels);
  workspace->pet.Radon     = gsl_vector_alloc(ndata);
  workspace->pet.ExpRadon  = gsl_vector_alloc(ndata);
  workspace->pet.defined   = 0;
  
  workspace->xfct.ones  = gsl_vector_alloc(npixels);
  for (j=0; j< MAX_NTHREADS; j++)
    workspace->xfct.integrand[j] = gsl_vector_alloc(npixels);
  
  workspace->xfct.defined     = 0;
  workspace->xfct.definedDbtT = 0;

  workspace->spect.integrand         = gsl_vector_alloc(npixels);
  workspace->spect.defined = 0;

  workspace->npixels = npixels;
  workspace->ndata   = ndata;
  workspace->size    = size;
  workspace->nviews  = nviews;
  workspace->nrays   = nrays;
}

/*+====================================================+
  
  FUNCTION: raft_projection_workspace_free
  
  Frees projection workspace.
  
  Input: 
  
  workspace - projection workspace
    
  +====================================================+
*/

void
raft_projection_workspace_free(raft_proj_t *workspace)
{
  int j, nviews;

  nviews = workspace->nviews;
  
  gsl_vector_free(workspace->likelihood.lograd);
  gsl_vector_free(workspace->likelihood.invrad);
  gsl_vector_free(workspace->likelihood.exprad);
  gsl_vector_free(workspace->likelihood.exprad2);
  gsl_vector_free(workspace->likelihood.radon);
  gsl_vector_free(workspace->likelihood.ones);

  gsl_vector_free(workspace->generalized.integrand);

  gsl_vector_free(workspace->pet.integrand);
  gsl_vector_free(workspace->pet.weight);
  gsl_vector_free(workspace->pet.Radon);
  gsl_vector_free(workspace->pet.ExpRadon);
  
  gsl_vector_free(workspace->xfct.ones);

  for(j=0; j < MAX_NTHREADS; j++)
    gsl_vector_free(workspace->xfct.integrand[j]);
    
  gsl_vector_free(workspace->spect.integrand);
  
  for(j = 0; j < nviews; j++)
    {
      gsl_vector_free(workspace->dbt.Dbt[j]);
      gsl_vector_free(workspace->dbt.ExpDbt[j]);

      gsl_vector_free(workspace->xfct.weight[j]);
      gsl_vector_free(workspace->xfct.weightP[j]);
      gsl_vector_free(workspace->xfct.weightF[j]);
      gsl_vector_free(workspace->xfct.LogWeightF[j]);
      gsl_vector_free(workspace->xfct.LogWeightP[j]);
      
      gsl_vector_free(workspace->xfct.DbtT[j]);
      gsl_vector_free(workspace->xfct.DbtF[j]);
      gsl_vector_free(workspace->xfct.ExpDbtT[j]);
      gsl_vector_free(workspace->xfct.ExpDbtF[j]);
      
      gsl_vector_free(workspace->spect.Dbt[j]);
      gsl_vector_free(workspace->spect.weight[j]);
    }

  free(workspace->dbt.Dbt);
  free(workspace->dbt.ExpDbt);

  free(workspace->xfct.weight);
  free(workspace->xfct.weightP);
  free(workspace->xfct.weightF);
  free(workspace->xfct.LogWeightF);
  free(workspace->xfct.LogWeightP);
  
  free(workspace->xfct.DbtT);
  free(workspace->xfct.DbtF);
  free(workspace->xfct.ExpDbtT);
  free(workspace->xfct.ExpDbtF);

  free(workspace->spect.Dbt);
  free(workspace->spect.weight);
}

/*######################################################
  Section:  Restriction at Rays
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_projection_phantom_atray

  Restriction of a given phantom at a straight line 
  (a ray). 

  Let "f" be the input phantom, and "(t,a)" a pair in
  polar coordinates defining a ray, referenced by an index
  "j". This function returns the one-dimensional function 
  "g(s) = f(tv + sv')" with "v=(cos a,sin a)" and 
  "v'=(-sin a,cos a)".
  
  Input: 
  
  data - scan data. See <raft_scan_t>.
  workspace - projection workspace. See <raft_proj_t>.
  phantom - phantom vector
  j - ray index
  
  Output:
  
  rest - one-dimensional function (previously alloc'd with 
         sufficient dimension)

  rsize - length for vector rest
  
  +====================================================+
*/

void raft_projection_phantom_atray(raft_scan_t *data,
				   gsl_vector *phantom,
				   int j,
				   gsl_vector *rest,
				   int rsize)
{
  div_t D;
  int i, ray, angle, size, row, column;
  double cost, sint, tant, x, y, t;
  double fcost, tol, value;
  gsl_vector_view phant;  

  tol = 1/sqrt(2); 
  size = raft_scan_get_size(data);
  
  D = div(j, raft_scan_get_nrays(data));
  angle = D.quot;
  ray   = D.rem;
  
  cost  = gsl_vector_get(data->costheta,angle);
  sint  = gsl_vector_get(data->sintheta,angle);
  tant  = gsl_vector_get(data->tantheta,angle); 
  fcost = fabs(cost);
  
  t = raft_scan_get_ray(data,ray);

  if(fcost < tol) 
    {
      rsize = 0; 
      
      for(i=0; i < size; i++)
	{
	  /* column of the phantom matrix */
	  phant = gsl_vector_subvector_with_stride(phantom, i, size, size);
	  
	  x = raft_scan_get_x(data,i);
	  y = t/sint - x/tant;
	  
	  row = interp_nearest(data->y, y, RAFT_REG);
	  
	  if(SIGN(row)>0)	  
	    {
	      value = interp_linear_reverse(&phant.vector, data->y, row, y);
	      
	      gsl_vector_set(rest, rsize, value);

	      rsize++;
	    }
	}
    }
  else
    {
      rsize = 0;
      
      for(i=0; i < size; i++)
	{
	  /* row of the phantom matrix */
	  phant = gsl_vector_subvector(phantom, (size-1-i)*size, size);
	  
	  y = raft_scan_get_y(data,i);
	  x = t/cost - y*tant;
	  
          column = interp_nearest(data->x, x, RAFT_REG);
	  
	  if(SIGN(column)>0)	  
	    {
	      value = interp_linear(&phant.vector, data->x, column, x);

	      gsl_vector_set(rest, rsize, value);
	      
	      rsize++;      
	    }	  
	}     
    }
}

/*######################################################
  Section: Computing Radon transforms
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_projection_radon
  
  Compute the Radon Transform.
  
  Input: 
  
  data - scan data. See <raft_scan_t>.
  workspace - projection workspace. See <raft_proj_t>.
  phantom - phantom vector
    
  Output:
  
  radon - radon transform
  
  +====================================================+
*/

void *radon_loop(void *t)
{
  int w, j, k, ndata, ray;
  parRad_t *param;
  double r, eps, ratio, logr, invr, expr, expr2;
     
  param = (parRad_t *)t;
  
  ndata = param->nviews * param->nrays;

  eps   = param->workspace->likelihood.eps;
  ratio = param->workspace->likelihood.ratio;
  
  for(w = param->colIndex[0]; w < MIN(param->colIndex[1], ndata); w++)
    {
      j = w/param->nrays;
      k = w%param->nrays;

      r = eval_rayintegral_(param->data, 
			   param->phantom, 
			   j, k,
			   param->data->memory.a, 
			   &param->data->memory.norma,
			   &param->data->memory.dotp);
	  
      if(fabs(r) > eps*ratio){
	logr = log(r);
	invr = 1/r;
      }
      else{
	logr = 0;
	invr = 1;
      }
      
      expr  = exp(r);
      expr2 = exp(r/2);
      
      ray = k + j*param->nrays;
	  
      gsl_vector_set(param->workspace->likelihood.lograd, ray, logr);
      gsl_vector_set(param->workspace->likelihood.invrad, ray, invr);
      gsl_vector_set(param->workspace->likelihood.exprad, ray, expr);
      gsl_vector_set(param->workspace->likelihood.exprad2, ray, expr2);
      gsl_vector_set(param->radon, ray, r);
    }  
  
  pthread_exit(NULL);
}


void raft_projection_radon(raft_scan_t *data, 
			   raft_proj_t *workspace,
			   gsl_vector *phantom,
			   gsl_vector *radon)
{
  pthread_t thread[RADON_THREADS+1];
  pthread_attr_t attr;
  int size, nviews, nrays, e, n, rc;    
  parRad_t param[RADON_THREADS+1];  
  void *status;
  
  // Initialize and set thread detached attribute
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  size   = raft_scan_get_size(data);
  nviews = raft_scan_get_nviews(data);
  nrays  = raft_scan_get_nrays(data);

  e = (int) floor((nviews*nrays)/RADON_THREADS);
      
  for(n = 0; n < RADON_THREADS+1; n++) 
    {
      param[n].data = data;
      param[n].workspace = workspace;
      param[n].phantom = phantom;
      param[n].radon   = radon;
      param[n].size    = size; 
      param[n].nviews  = nviews;
      param[n].nrays   = nrays;
      param[n].nthread = n;      
      
      param[n].colIndex[0] = e * n;
      param[n].colIndex[1] = (n+1) * e;
      
      rc = pthread_create(&thread[n], &attr, radon_loop, (void *)&param[n]);
    }
  
  // Free attribute and wait for the other threads
  pthread_attr_destroy(&attr);
  for(n = 0; n < RADON_THREADS+1; n++) 
    {
      rc = pthread_join(thread[n], &status);
    }  
}

/*+====================================================+
  
  FUNCTION: raft_projection_radon_ert
  
  Compute the Exponential Radon Transform (ERT).
  
  Input: 
  
  data - scan data. See <raft_scan_t>.
  workspace - projection workspace. See <raft_proj_t>.
  act - activity vector (phantom)
  att - attenuation value
    
  Output:
  
  radon - radon transform
  
  +====================================================+
*/

void raft_projection_radon_ert(raft_scan_t *data, 
			       raft_proj_t *workspace,
			       double *att,
			       gsl_vector *act,
			       gsl_vector *radon)
{
  int j, k, nrays, nviews, ray;
  double r, logr, invr, eps, ratio, expr, expr2;
  
  nrays  = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);
  eps    = workspace->likelihood.eps;
  ratio = workspace->likelihood.ratio;
  
  for(j = 0; j < nviews; j++)
    {
      for(k = 0; k < nrays; k++)
	{
	  r = eval_rayintegral_ert(data, 
				   att,
				   act, 
				   j, k,
				   data->memory.a, 
				   &data->memory.norma,
				   &data->memory.dotp);
	  
	  if(fabs(r) > eps*ratio){
	    logr = log(r);
	    invr = 1/r;
	  }
	  else{
	    logr = 0;
	    invr = 1;
	  }

	  expr  = exp(r);
	  expr2 = exp(r/2);

	  ray = k + j*nrays;
	  
	  gsl_vector_set(workspace->likelihood.lograd, ray, logr);
	  gsl_vector_set(workspace->likelihood.invrad, ray, invr);
	  gsl_vector_set(workspace->likelihood.exprad, ray, expr);
	  gsl_vector_set(workspace->likelihood.exprad2, ray, expr2);
	  gsl_vector_set(radon, ray, r);
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_projection_radon_pet
  
  Compute the PET Radon transform.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  workspace - projection workspace 
  act - activity phantom vector
  att - attenuation phantom vector
    
  Output:
  
  radon - attenuated radon transform
  
  +====================================================+
*/

void raft_projection_radon_pet(raft_scan_t *data, 
			       raft_proj_t *workspace,
			       gsl_vector *act,
			       gsl_vector *att,
			       gsl_vector *radon)
{
  int j, k, nrays, nviews, ray;
  double r, logr, invr, eps, ratio, expr, expr2;
  
  if(!workspace->pet.defined)
    {
      raft_weight_pet(data, workspace, att);
    }

  nrays = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);
  eps    = workspace->likelihood.eps;
  ratio = workspace->likelihood.ratio;
  
  for(j = 0; j < nviews; j++)
    {
      for(k = 0; k < nrays; k++)
	{
	  r = eval_rayintegral_pet(data, 
				   workspace,
				   act, 
				   j, k,
				   data->memory.a, 
				   &data->memory.norma,
				   &data->memory.dotp);
	  
	  if(fabs(r) > eps*ratio){
	    logr = log(r);
	    invr = 1/r;
	  }
	  else{
	    logr = 0;
	    invr = 1;
	  }
	  
	  expr  = exp(r);
	  expr2 = exp(r/2);
	  
	  ray = k + j*nrays;

	  gsl_vector_set(workspace->likelihood.lograd, ray, logr);
	  gsl_vector_set(workspace->likelihood.invrad, ray, invr);
	  gsl_vector_set(workspace->likelihood.exprad, ray, expr);
	  gsl_vector_set(workspace->likelihood.exprad2, ray, expr2);
	  gsl_vector_set(radon, ray, r);
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_projection_radon_spect
  
  Compute the SPECT Radon transform.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  workspace - projection workspace 
  act - activity phantom vector
  att - attenuation phantom vector
    
  Output:
  
  radon - attenuated radon transform
  
  +====================================================+
*/

void raft_projection_radon_spect(raft_scan_t *data, 
				 raft_proj_t *workspace,
				 gsl_vector *act,
				 gsl_vector *att,
				 gsl_vector *radon)
{
  int j, k, nrays, nviews, ray;
  double r, logr, invr, eps, ratio, expr, expr2;
  
  if(!workspace->spect.defined)
    {
      raft_weight_spect(data, workspace, att);
    }
  
  nrays  = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);
  eps    = workspace->likelihood.eps;
  ratio  = workspace->likelihood.ratio;

  for(j = 0; j < nviews; j++)
    {
      for(k = 0; k < nrays; k++)
	{
	  r = eval_rayintegral_spect(data, 
				     workspace,
				     act, 
				     j, k,
				     data->memory.a, 
				     &data->memory.norma,
				     &data->memory.dotp);
	  
	  if(fabs(r) > eps*ratio){
	    logr = log(r);
	    invr = 1/r;
	  }
	  else{
	    logr = 0;
	    invr = 1;
	  }
	  
	  expr  = exp(r);
	  expr2 = exp(r/2);
	  
	  ray = k + j*nrays;

	  gsl_vector_set(workspace->likelihood.lograd, ray, logr);
	  gsl_vector_set(workspace->likelihood.invrad, ray, invr);
	  gsl_vector_set(workspace->likelihood.exprad, ray, expr);
	  gsl_vector_set(workspace->likelihood.exprad2, ray, expr2);
	  gsl_vector_set(radon, ray, r);
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_projection_radon_xfct
  
  Compute the fluorescent Radon transform.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  workspace - projection workspace 
  act - activity phantom vector
  attF - XFCT attenuation phantom vector
  attT - transmission attenuation phantom vector
    
  Output:
  
  radon - attenuated radon transform
   
  +====================================================+
*/

void *radonXfct_loop(void *t)
{
  int w, j, k, ndata, ray;
  parRad_t *param;
  double r, eps, ratio, logr, invr, expr, expr2;
     
  param = (parRad_t *)t;
  
  ndata = param->nviews * param->nrays;

  eps   = param->workspace->likelihood.eps;
  ratio = param->workspace->likelihood.ratio;
 
  for(w = param->colIndex[0]; w < MIN(param->colIndex[1], ndata); w++)
    {
      j = w/param->nrays;
      k = w%param->nrays;      
      
      r = eval_rayintegral_xfct_bythread(param->nthread,
					 param->data, 
					 param->workspace,
					 param->phantom, 
					 j, k,
					 param->data->memory.a, 
					 &param->data->memory.norma,
					 &param->data->memory.dotp);
      
  
      if(fabs(r) > eps*ratio){
	logr = log(r);
	invr = 1/r;
      }
      else{
	logr = 0;
	invr = 1;
      }
	  
      expr  = exp(r);
      expr2 = exp(r/2);
      
      ray = k + j*param->nrays;

      gsl_vector_set(param->workspace->likelihood.lograd, ray, logr);
      gsl_vector_set(param->workspace->likelihood.invrad, ray, invr);
      gsl_vector_set(param->workspace->likelihood.exprad, ray, expr);
      gsl_vector_set(param->workspace->likelihood.exprad2, ray, expr2);

      gsl_vector_set(param->radon, ray, r);
    }  
  
  pthread_exit(NULL);
}

void raft_projection_radon_xfct(raft_scan_t *data, 
				raft_proj_t *workspace,
				gsl_vector *act,
				gsl_vector *attF,
				gsl_vector *attT,
				gsl_vector *radon)
{
  pthread_t thread[RADON_THREADS+1];
  pthread_attr_t attr;
  int size, nviews, nrays, e, n, rc;    
  parRad_t param[RADON_THREADS+1];  
  void *status;

  if(!workspace->xfct.defined)
    {
      raft_weight_xfct(data, workspace, attT, attF);
    }
   
  // Initialize and set thread detached attribute
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  size   = raft_scan_get_size(data);
  nviews = raft_scan_get_nviews(data);
  nrays  = raft_scan_get_nrays(data);

  e = (int) floor((nviews*nrays)/RADON_THREADS);
      
  for(n = 0; n < RADON_THREADS+1; n++) 
    {
      param[n].data = data;
      param[n].workspace = workspace;
      param[n].phantom = act;
      param[n].radon   = radon;
      param[n].size    = size; 
      param[n].nviews  = nviews;
      param[n].nrays   = nrays;
      param[n].nthread = n;      
      
      param[n].colIndex[0] = e * n;
      param[n].colIndex[1] = (n+1) * e;
      
      rc = pthread_create(&thread[n], &attr, radonXfct_loop, (void *)&param[n]);
    }
  
  // Free attribute and wait for the other threads
  pthread_attr_destroy(&attr);
  for(n = 0; n < RADON_THREADS+1; n++) 
    {
      rc = pthread_join(thread[n], &status);
    }  

  /*
  int j, k, nrays, nviews, ray;
  double r, logr, invr, eps, ratio, expr, expr2;
  
  if(!workspace->xfct.defined)
    {
      raft_weight_xfct(data, workspace, attT, attF);
    }
  
  nrays = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);
  eps    = workspace->likelihood.eps;
  ratio  = workspace->likelihood.ratio;

  for(j = 0; j < nviews; j++)
    {
      for(k = 0; k < nrays; k++)
	{
	  r = eval_rayintegral_xfct(data, 
				    workspace,
				    act, 
				    j, k,
				    data->memory.a, 
				    &data->memory.norma,
				    &data->memory.dotp);
	  
	  if(fabs(r) > eps*ratio){
	    logr = log(r);
	    invr = 1/r;
	  }
	  else{
	    logr = 0;
	    invr = 1;
	  }
	  
	  expr  = exp(r);
	  expr2 = exp(r/2);
	  
	  ray = k + j*nrays;

	  gsl_vector_set(workspace->likelihood.lograd, ray, logr);
	  gsl_vector_set(workspace->likelihood.invrad, ray, invr);
	  gsl_vector_set(workspace->likelihood.exprad, ray, expr);
	  gsl_vector_set(workspace->likelihood.exprad2, ray, expr2);
	  gsl_vector_set(radon, ray, r);
	}
    }
  */
}

/*+====================================================+
  
  FUNCTION: raft_projection_radon_generalized
  
  Compute the generalized Radon transform.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  workspace - projection workspace 
  f - phantom vector
  weight - weight matrix
      
  Output:
  
  radon - generalized radon transform
  
  +====================================================+
*/

void raft_projection_radon_generalized(raft_scan_t *data, 
				       raft_proj_t *workspace,
				       gsl_vector *f,
				       gsl_vector **weight,
				       gsl_vector *radon)
{
  int j, k, nrays, nviews, ray;
  double r, logr, invr, eps, ratio, expr, expr2;
  
  nrays  = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);

  eps    = workspace->likelihood.eps;
  ratio  = workspace->likelihood.ratio;

  for(j = 0; j < nviews; j++)
    {
      for(k = 0; k < nrays; k++)
	{
	  r = eval_rayintegral_generalized(data, 
					   workspace,
					   weight,
					   f, 
					   j, k,
					   data->memory.a, 
					   &data->memory.norma,
					   &data->memory.dotp);
	  
	  if(fabs(r) > eps*ratio){
	    logr = log(r);
	    invr = 1/r;
	  }
	  else{
	    logr = 0;
	    invr = 0;
	  }
	  
	  expr  = exp(r);
	  expr2 = exp(r/2);

	  ray = k + j*nrays;
	  
	  gsl_vector_set(workspace->likelihood.lograd, ray, logr);
	  gsl_vector_set(workspace->likelihood.invrad, ray, invr);
	  gsl_vector_set(workspace->likelihood.exprad, ray, expr);
	  gsl_vector_set(workspace->likelihood.exprad2, ray, expr2);
	  gsl_vector_set(radon, ray, r);
	}
    }
}

/*######################################################
  Section: Computing divergent beam transform
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_projection_radon_dbt
  
  Compute the divergent beam transform.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  workspace - projection workspace
  phantom - attenuation phantom vector
    
  Output:
  
  dbt - divergent beam transform
       
  _Remark_:

   the output is a matrix with dimension nviews x (nrays)^2, 
  each row represent a projection of the divergent transform, 
  for a fixed angle.
  
  +====================================================+
*/

void *radonDbt_loop(void *t)
{
  int w, j, i, k, ndata, I;
  parDbt_t *param;
  double r, er;
     
  param = (parDbt_t *)t;
  
  ndata = param->size * param->size;
  
  for(w = param->colIndex[0]; w < MIN(param->colIndex[1], ndata); w++)
    {
      i = w/param->size;
      k = w%param->size;

      I = k + i*param->size;

      for(j=0; j < param->nviews; j++)
	{
	  r  = eval_rayintegral_dbt(param->data, param->phantom, j, k, i);
	  er = exp(r);
	  
	  gsl_vector_set(param->dbt[j], I, r);
	  gsl_vector_set(param->workspace->dbt.ExpDbt[j], I, er);
	}
    }  
  
  pthread_exit(NULL);
}


void raft_projection_radon_dbt(raft_scan_t *data, 
			       raft_proj_t *workspace,
			       gsl_vector *phantom,
			       gsl_vector **dbt)
{
  pthread_t thread[RADON_THREADS+1];
  pthread_attr_t attr;
  int size, nviews, nrays, e, n, rc;    
  parDbt_t param[RADON_THREADS+1];  
  void *status;
  
  // Initialize and set thread detached attribute
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  size   = raft_scan_get_size(data);
  nviews = raft_scan_get_nviews(data);
  nrays  = raft_scan_get_nrays(data);

  e = (int) floor((size*size)/(RADON_THREADS));
      
  for(n = 0; n < RADON_THREADS+1; n++) 
    {
      param[n].data = data;
      param[n].workspace = workspace;
      param[n].phantom = phantom;
      param[n].dbt     = dbt;
      param[n].size    = size; 
      param[n].nviews  = nviews;
      param[n].nrays   = nrays;
      param[n].nthread = n;      
      
      param[n].colIndex[0] = e * n;
      param[n].colIndex[1] = (n+1) * e;
      
      rc = pthread_create(&thread[n], &attr, radonDbt_loop, (void *)&param[n]);
    }
  
  // Free attribute and wait for the other threads
  pthread_attr_destroy(&attr);
  for(n = 0; n < RADON_THREADS+1; n++) 
    {
      rc = pthread_join(thread[n], &status);
    }  
}

/*+====================================================+
  
  FUNCTION: raft_projection_radon_dbt_get
  
  Get a matrix entry from the divergent beam transform.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  dbt - divergent beam transform
  i - ray index
  j - pixel index 
  
  +====================================================+
*/

double
raft_projection_radon_dbt_get(raft_scan_t *data,
			      gsl_vector **dbt,
			      int i,
			      int j)
{
  int angle, row, column;
  div_t D;
  
  D = div(i, raft_scan_get_nrays(data));
  angle = D.quot;
   
  D = div(j, raft_scan_get_size(data));
  row = D.quot;
  column = D.rem;

  return raft_phantom_get(dbt[angle], row, column);
}


/*+====================================================+
  
  FUNCTION: raft_projection_workspace_get_exp_dbt
  
  Get exponential of the divergent beam transform.
  
  Input: 
  
  workspace - projection workspace
  
  Output:

  exp - exponential of divergent beam transform

  +====================================================+
*/

void 
raft_projection_workspace_get_exp_dbt(raft_proj_t *workspace,
				      gsl_vector **exp)
{
  int j;

  for(j=0; j < workspace->nviews; j++)
    {
      gsl_blas_dcopy(workspace->dbt.ExpDbt[j],exp[j]);
    }
}


/*######################################################
  Section: Computing Radon transform views
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_projection_radon_view
  
  Compute the Radon transform for a fixed angle.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  phantom - phantom vector
  j - view index
      
  Output:
  
  view - projection vector
  
  +====================================================+
*/

void raft_projection_radon_view(raft_scan_t *data, 
				gsl_vector *phantom,
				int j,
				gsl_vector *view)
{
  int k, nrays, nviews;
  double r;
  
  nrays  = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);
  
  for(k = 0; k < nrays; k++)
    {
      r = eval_rayintegral_(data, phantom, j, k,
			   data->memory.a, 
			   &data->memory.norma,
			   &data->memory.dotp);
      
      gsl_vector_set(view, k, r);
    }
}


/*+====================================================+
  
  FUNCTION: raft_projection_radon_pet_view
  
  Compute the PET Radon transform for a fixed angle.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  workspace - projection workspace 
  act - activity phantom vector
  att - attenuation phantom vector
  j - view index
    
  Output:
  
  view - projection vector
  
  +====================================================+
*/

void raft_projection_radon_pet_view(raft_scan_t *data, 
				    raft_proj_t *workspace,
				    gsl_vector *act,
				    gsl_vector *att,
				    int j,
				    gsl_vector *view)
{
  gsl_vector_view W;
  int nrays;

  nrays = raft_scan_get_nrays(data);

  if(!workspace->pet.defined)
    {
      raft_weight_pet(data, workspace, att);
    }
  
  raft_projection_radon_view(data, act, j, view);
  
  W = gsl_vector_subvector(workspace->pet.weight, j*nrays, nrays);

  gsl_vector_mul(view, &W.vector);
}

/*+====================================================+
  
  FUNCTION: raft_projection_radon_spect_view
  
  Compute the SPECT Radon transform for a fixed angle.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  workspace - projection workspace 
  act - activity phantom vector
  att - attenuation phantom vector
  j - view index
    
  Output:
  
  view - projection vector    
  
  +====================================================+
*/

void raft_projection_radon_spect_view(raft_scan_t *data, 
				      raft_proj_t *workspace,
				      gsl_vector *act,
				      gsl_vector *att,
				      int j,
				      gsl_vector *view)
{
  int nrays, nviews;
  
  nrays   = raft_scan_get_nrays(data);
  nviews  = raft_scan_get_nviews(data);
  
  if(!workspace->spect.defined)
    {
      raft_weight_spect(data, workspace, att);
    }
  
  gsl_blas_dcopy(act, workspace->spect.integrand);
  
  gsl_vector_mul(workspace->spect.integrand, workspace->spect.weight[j]);
  
  raft_projection_radon_view(data, workspace->spect.integrand, j, view);  
}


/*+====================================================+
  
  FUNCTION: raft_projection_radon_xfct_view
  
  Compute the fluorescent Radon transform for a fixed angle.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  workspace - projection workspace 
  act - activity phantom vector
  attF - XFCT attenuation phantom vector
  attT - transmission attenuation phantom vector
  j - view index
    
  Output:
  
  view - projection vector
  
  +====================================================+
*/

void raft_projection_radon_xfct_view(raft_scan_t *data, 
				     raft_proj_t *workspace,
				     gsl_vector *act,
				     gsl_vector *attF,
				     gsl_vector *attT,
				     int j,
				     gsl_vector *view)
{
  int nrays, nviews; 

  nrays = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);
  
  if(!workspace->xfct.defined)
    {
      raft_weight_xfct(data, workspace, attT, attF);
    }
   
  gsl_blas_dcopy(workspace->xfct.weight[j], workspace->xfct.integrand[0]);
  
  gsl_vector_mul(workspace->xfct.integrand[0], act);
  
  raft_projection_radon_view(data, workspace->xfct.integrand[0], j, view); 
}


/*+====================================================+
  
  FUNCTION: raft_projection_radon_dbt_view
  
  Compute the divergent beam transform for a fixed angle.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  workspace - projection workspace
  phantom - attenuation phantom vector
  j - angle index
    
  Output:
  
  dbt - divergent beam transform

  Return: 

  Updated view index j.
  
  +====================================================+
*/

int raft_projection_radon_dbt_view(raft_scan_t *data, 
			           raft_proj_t *workspace,
				   int j,
				   gsl_vector *phantom,
				   gsl_vector *dbt)
{
  int I, i, k, size, nviews, nrays, view;
  double r, er; 
  div_t D;
  
  size = raft_scan_get_size(data);
  nrays  = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);

  if( j >= nviews )
    {
      D = div(j + nviews, nviews);
      view = D.rem;
    }
  else
    {
      view = j;
    }

  for(i = 0; i < size; i++)
    {
      for(k = 0; k < size; k++)
	{
	  I = k + i*size;
	  
	  r  = eval_rayintegral_dbt(data, phantom, view, k, i);
	  er = exp(r);
	  
	  gsl_vector_set(dbt, I, r);
	  gsl_vector_set(workspace->dbt.ExpDbt[view], I, er);
	}
    }

  return view;
}

/*######################################################
  Section: Computing ray sums
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_projection_radon_ray
  
  Compute the Radon transform for a fixed ray.
  
  Input: 
  
  data - scan data. See <raft_scan_t>.
  phantom - phantom vector
  i - ray index. 
      
  Output:
  
  rsum - ray sum
  row - vector such that <row,phantom> = rsum
  dotp - dot product <row,phantom>
  snorm - squared norm of vector row
  
  Return:

  RAFT_SUCCESS - sucessfull operation
  RAFT_EDOM - input domain error (wrong ray index j)

  _Remark_:

  Index i varies from 0 to D-1, with D=R*V, V the number
  of views and R the number of rays per view. We use the
  convention i = a(i) + R*b(i), where a=a(i) define the ray 
  index (per view) and b=b(i) the view index.  
  Also, due to interpolation on the ray sum, 'dotp' is an 
  approximation to 'rsum'.

  +====================================================+
*/

int
raft_projection_radon_ray(raft_scan_t *data,
			  gsl_vector *phantom,
			  int i,
			  double *rsum,
			  double *dotp,
			  double *snorm,
			  gsl_vector *row)
{
  if(i<0 || i>raft_scan_get_ndata(data))
    {
      return RAFT_EDOM;
    }
  else
    {
      gsl_vector_set_zero(row);

      *rsum = eval_ray_projection(data,
				  phantom, 
				  i,
				  row,
				  snorm, 
				  dotp);
      
      return RAFT_SUCCESS;
    }
}

/*+====================================================+
  
  FUNCTION: raft_projection_radon_ray_pet
  
  Compute the PET Radon transform for a fixed ray.
  
  Input: 
  
  data - scan data. See <raft_scan_t>.
  phantom - phantom vector
  proj - projection workspace. See <raft_proj_t>.
  att - attenuation map.
  i - ray index. 
      
  Output:
  
  rsum - ray sum
  row - vector such that <row,phantom> = rsum
  dotp - dot product <row,phantom>
  snorm - squared norm of vector row
  
  Return:

  RAFT_SUCCESS - sucessfull operation
  RAFT_EDOM - input domain error (wrong ray index j)

  _Remark_:

  Index i varies from 0 to D-1, with D=R*V, V the number
  of views and R the number of rays per view. We use the
  convention i = a(i) + R*b(i), where a=a(i) define the ray 
  index (per view) and b=b(i) the view index.  
  Also, due to interpolation on the ray sum,
  'dotp' is an approximation to 'rsum'.

  +====================================================+
*/

int
raft_projection_radon_ray_pet(raft_scan_t *data,
			      raft_proj_t *proj,
			      gsl_vector *phantom,
			      gsl_vector *att,
			      int i,
			      double *rsum,
			      double *dotp,
			      double *snorm,
			      gsl_vector *row)
{
  if(i<0 || i>raft_scan_get_ndata(data))
    {
      return RAFT_EDOM;
    }
  else
    {
      gsl_vector_set_zero(row);
      
      *rsum = eval_ray_projection_pet(data,
				      proj,
				      phantom,
				      att,
				      i,
				      row,
				      snorm,
				      dotp);
      return RAFT_SUCCESS;
    }
}

/*+====================================================+
  
  FUNCTION: raft_projection_radon_ray_spect
  
  Compute the SPECT Radon transform for a fixed ray.
  
  Input: 
  
  data - scan data. See <raft_scan_t>.
  phantom - phantom vector
  proj - projection workspace. See <raft_proj_t>.
  att - attenuation map.
  i - ray index. 
      
  Output:
  
  rsum - ray sum
  row - vector such that <row,phantom> = rsum
  dotp - dot product <row,phantom>
  snorm - squared norm of vector row
  
  Return:

  RAFT_SUCCESS - sucessfull operation
  RAFT_EDOM - input domain error (wrong ray index j)

  _Remark_:

  Index i varies from 0 to D-1, with D=R*V, V the number
  of views and R the number of rays per view. We use the
  convention i = a(i) + R*b(i), where a=a(i) define the ray 
  index (per view) and b=b(i) the view index.  
  Also, due to interpolation on the ray sum,
  'dotp' is an approximation to 'rsum'.
  
  +====================================================+
*/

int
raft_projection_radon_ray_spect(raft_scan_t *data,
				raft_proj_t *proj,
				gsl_vector *phantom,
				gsl_vector *att,
				int i,
				double *rsum,
				double *dotp,
				double *snorm,
				gsl_vector *row)
{
  if(i<0 || i>raft_scan_get_ndata(data))
    {
      return RAFT_EDOM;
    }
  else
    {
      gsl_vector_set_zero(row);

      *rsum = eval_ray_projection_spect(data,
					proj,
					phantom,
					att,
					i,
					row,
					snorm,
					dotp);
      return RAFT_SUCCESS;
    }
}

/*+====================================================+
  
  FUNCTION: raft_projection_radon_ray_xfct
  
  Compute the XFCT Radon transform for a fixed ray.
  
  Input: 
  
  data - scan data. See <raft_scan_t>.
  phantom - phantom vector
  proj - projection workspace. See <raft_proj_t>.
  attT - transmission attenuation map.
  attF - fluorescence attenuation map.
  i - ray index. 
      
  Output:
  
  rsum - ray sum
  row - vector such that <row,phantom> = rsum
  dotp - dot product <row,phantom>
  snorm - squared norm of vector row
  
  Return:

  RAFT_SUCCESS - sucessfull operation
  RAFT_EDOM - input domain error (wrong ray index j)

  _Remark_:

  Index i varies from 0 to D-1, with D=R*V, V the number
  of views and R the number of rays per view. We use the
  convention i = a(i) + R*b(i), where a=a(i) define the ray 
  index (per view) and b=b(i) the view index.  
  Also, due to interpolation on the ray sum,
  'dotp' is an approximation to 'rsum'.

  +====================================================+
*/

int
raft_projection_radon_ray_xfct(raft_scan_t *data,
			       raft_proj_t *proj,
			       gsl_vector *phantom,
			       gsl_vector *attT,
			       gsl_vector *attF,
			       int i,
			       double *rsum,
			       double *dotp,
			       double *snorm,
			       gsl_vector *row)
{
  if(i<0 || i>raft_scan_get_ndata(data))
    {
      return RAFT_EDOM;
    }
  else
    {
      gsl_vector_set_zero(row);

      *rsum = eval_ray_projection_xfct(data,
				       proj,
				       phantom,
				       attT,
				       attF,
				       i,
				       row,
				       snorm,
				       dotp);
      return RAFT_SUCCESS;
    }
}


/*+====================================================+
  
  FUNCTION: raft_projection_radon_ray_dbt
  
  Compute the DBT passing through a given point for
  a given view.
    
  Input: 
  
  data - scan data . See <raft_scan_t>.
  phantom - phantom matrix.  
  j - view index 
  x - x-coordinate index
  y - y-coordinate index
  
  +====================================================+
*/

double
raft_projection_radon_ray_dbt(raft_scan_t *data,
			      gsl_vector *phantom,
			      int j,
			      int x,
			      int y)
{
  return eval_rayintegral_dbt(data, phantom, j, x, y);  
}

/*######################################################
  Section: Computing transmission data
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_projection_monodata
  
  Compute the monochromatic transmission data.

  Input: 
  
  data - scan data . See <raft_scan_t>.
  phantom - phantom vector
    
  Output:
  
  ephotons - expected number of photons
  radon - radon transform

  +====================================================+
*/

void raft_projection_monodata(raft_scan_t *data, 
			      gsl_vector *phantom,
			      gsl_vector *ephotons,
			      gsl_vector *radon)
{
  int j, k, nrays, nviews;
  double r, b, em, tr, ef, d;
  
  nrays = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);

  em = raft_scan_get_photons_emitted(data);
  tr = raft_scan_get_photons_transmittance(data);
  ef = raft_scan_get_photons_efficiency(data);
  b  = em*tr*ef; 
  
  for(j = 0; j < nviews; j++)
    {
      for(k = 0; k < nrays; k++)
	{
	  r = eval_rayintegral_(data, phantom, j, k,
			       data->memory.a, 
			       &data->memory.norma,
			       &data->memory.dotp);	  
	  
	  d = b*exp(-r);

	  gsl_vector_set(ephotons, k + j*nrays, d);
	  gsl_vector_set(radon, k + j*nrays, r);
	}
    }
}


/*+====================================================+
  
  FUNCTION: raft_projection_polydata
  
  Compute the polychromatic transmission data.

  Input: 
  
  data - scan data. See <raft_scan_t>.
  phantom - phantom data. See <raft_phantom_t>.
  workspace - projection workspace. See <raft_proj_t>.
    
  Output:
  
  ephotons - expected number of photons    
  
  +====================================================+
*/

void raft_projection_polydata(raft_scan_t *data, 
			      raft_proj_t *workspace,
			      raft_phantom_t *phantom,
			      gsl_vector *ephotons)
{
  int j, k, i, m, size, nrays, nviews, nergs, nsubs;
  double sum, d, r, g, s, ini;
  
  size = raft_scan_get_size(data);
  nrays = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);
  nergs = raft_scan_get_nenergies(data);
  nsubs = raft_scan_get_nsubs(data);
  
  for(m=0; m < nsubs; m++)
    {
      raft_projection_radon(data, 
			    workspace,
			    phantom->basis[m].vector, 
			    data->energy->basis[m].p);
    }
  
  for(j=0; j < nviews; j++)
    {
      for(i=0; i < nrays; i++)
	{
	  d = 0;
	  ini = 0;
	  
	  for(k=0; k < nergs; k++)
	    {
	      sum = 0;
	      
	      s = gsl_vector_get(data->energy->spectrum, k);
	      
	      ini += s;
	      
	      for(m=0; m < nsubs; m++)
		{
		  r = gsl_vector_get(data->energy->basis[m].p, i+j*nrays);
		  
		  g = gsl_matrix_get(data->energy->G, k, m);
		  
		  sum += (r*g);
		}

	      d += s*exp(-sum);
	    }

	  d = d/ini;
	  
	  gsl_vector_set(ephotons, i + j*nrays, d);
	}
    }
}

/*######################################################
  Section: Computing the Radon operator matrix
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_projection_pixel_intersection_length

  Compute the intersection length of a ray with a
  pixel. 
  
  Input: 

  scan - scan data . See <raft_scan_t>.
  i - ray number
  k - pixel number
  
  Return:

  Evaluated matrix.

  +====================================================+
*/

double
raft_projection_pixel_intersection_length(raft_scan_t *data, 
					  int i, 
					  int k)
{
  return eval_radon_charf_pixel(data, i, k);
}


/*+====================================================+
  
  FUNCTION: raft_projection_pixel_intersection_find
  
  Define whether a pixel belongs to a path starting at 
  a reference pixel.
  
  Input: 

  scan - scan data . See <raft_scan_t>.
  m - pixel number
  i - ray number
  j - pixel reference
  
  Return:

  It returns 1 whenever the pixel 'm' lies upon the
  ray starting at the pixel  
  
  +====================================================+
*/

int 
raft_projection_pixel_intersection_find(raft_scan_t *data, 
					int m,
					int i, 
					int j)
{
  double aij;
  
  aij = raft_projection_pixel_intersection_length(data,i,j);
  
  if(fabs(aij) < ZERO)
    return 0;
  else
    {
      div_t D;
      int size, angle, ray, row, column;
      double cost, sint, sm, sj, xj, yj, xm, ym;

      size = raft_scan_get_size(data);

      D = div(i, raft_scan_get_nrays(data));
      angle = D.quot;
      ray = D.rem;
      
      cost  = gsl_vector_get(data->costheta, angle);
      sint  = gsl_vector_get(data->sintheta, angle);
      
      D = div(j, size);
      row = D.quot;
      column = D.rem;

      xj = raft_scan_get_x(data, column);
      yj = raft_scan_get_y(data, size-1-row);
      
      sj = -xj*sint + yj*cost;

      D = div(m, size);
      row = D.quot;
      column = D.rem;

      xm = raft_scan_get_x(data, column);
      ym = raft_scan_get_y(data, size-1-row);
      
      sm = -xm*sint + ym*cost; 

      if((sm>sj) || fabs(sm-sj) < ZERO)
	return 1;
      else
	return 0;
    }
}

