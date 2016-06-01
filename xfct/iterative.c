#include "raft_scan.h"
#include "raft_iterative.h"
#include "raft_projection.h"
#include "raft_backprojection.h"
#include "raft_weight.h"

#include "projection.h"
#include "iterative.h"
#include "backprojection.h"
#include "likelihood.h"

/*######################################################
  Title: Iterative methods

  Header - <raft/raft_iterative.h>
  Type - <raft_iterative_t>
  $Id: iterative.c,v 1.28 2011-03-02 19:23:19 miqueles Exp $ - Last Update
  ####################################################*/

/*######################################################
  Private
  ####################################################*/

/*+====================================================+
  
  FUNCTION  ramla_ct_iteration

  One iteration of the RAMLA algorithm for standard CT. 
  
  Input  

  data - scan data. See <raft_scan_t>.
  workspace - workspace for iterative methods. See <raft_iterative_t>.
  b - transmission data
  p - previous estimation
  k - iteration number
    
  Output 

  n - next estimation
   
  +====================================================+
*/

void ramla_ct_iteration(raft_scan_t *data, 
			raft_iterative_t *workspace,
			gsl_vector *p,
			gsl_vector *n, 
			gsl_vector *b,
			int k)
{
  gsl_vector *I;
  int s, i, ndata, nviews, sT;
  double lambda, bb, rr;

  I = gsl_vector_alloc(workspace->ramla.sizeRow);

  ndata  = raft_scan_get_ndata(data);  
  nviews = raft_scan_get_nviews(data); 
  sT    = workspace->ramla.sizeRow;
  
  if(k>0)
    lambda = 1/((double) k);
  else
    lambda = 1;

  
  i = 0;
  
  do{

    raft_projection_radon(data, 
			  &workspace->backwork.projwork, 
			  p, 
			  workspace->em.radon);

    /* setting I = [i + T] */
    
    for(s=0; s < sT; s++)
      {
	div_t D;
	
	D = div(i + gsl_vector_get(workspace->ramla.row,s), nviews);
	
	gsl_vector_set(I, s, D.rem);
      }
    
    /**/
    
    for(s=0; s < ndata; s++)
      {
	bb = gsl_vector_get(b, s);
	rr = gsl_vector_get(workspace->em.radon, s);
	
	if(fabs(rr) < ZERO)
	  gsl_vector_set(workspace->em.data, s, 0);
	else
	  gsl_vector_set(workspace->em.data, s, bb/rr - 1);
      }
    
    raft_backp_partial(data, workspace->em.data, n, I, sT);
    
    /**/

    gsl_vector_mul(n, p);
    gsl_blas_dscal(lambda, n);
    gsl_blas_daxpy(1.0, p, n);
    
    gsl_blas_dcopy(n,p);

    i++;
        
  }while(i < nviews); 
  
  
  gsl_vector_free(I);
}

/*+====================================================+
  
  FUNCTION  art_ct_iteration

  One iteration of the ART algorithm in standard CT. 
  
  Input  

  data - scan data. See <raft_scan_t>.
  workspace - workspace for iterative methods. See <raft_iterative_t>.
  b - projection data 
  k - iteration number
    
  Output 

  n - next estimation
   
  +====================================================+
*/

void art_ct_iteration(raft_scan_t *data,	   
		      raft_iterative_t *workspace,
		      gsl_vector *n,
		      gsl_vector *b,
		      int k)
{
  double rad, norma, L, lambda, dotp, q;
  int ndata, npixels, nviews, nrays, i, m;
    
  ndata = raft_scan_get_ndata(data);  
  npixels = raft_scan_get_npixels(data);
  nviews = raft_scan_get_nviews(data);
  nrays = raft_scan_get_nrays(data);

  lambda = 1;
  
  i = 0;
  
  do{
    
    m = i;    
    
    gsl_vector_set_zero(workspace->art.row);
    
    rad = eval_ray_projection(data, n, m, workspace->art.row, &norma, &dotp);
    
    q = CVXCOMB(0.001, rad, 0.999, dotp);
    
    if(sqrt(norma) > ZERO)
      {
	L = (lambda*(gsl_vector_get(b,m) - q))/norma;
	
	gsl_blas_daxpy(L, workspace->art.row, n);
      }
    
    maxVector(n);

    i++;
    
  }while(i < ndata);   
}

/*+====================================================+
  
  FUNCTION  art_pet_iteration_knownAtt

  One iteration of the ART algorithm in PET with
  known attenuation map. 
  
  Input  

  data - scan data . See <raft_scan_t>.
  workspace - workspace for iterative methods . See <raft_iterative_t>.
  b - PET projection data
  att - attenuation map
  k - iteration number
    
  Output 

  n - next estimation

  +====================================================+
*/

void art_pet_iteration_knownAtt(raft_scan_t *data,	   
				raft_iterative_t *workspace,
				gsl_vector *n,
				gsl_vector *b,
				gsl_vector *att,
				int k)
{
  double rad, norma, L, lambda, dotp, q;
  int ndata, npixels, nviews, i, m;

  ndata = raft_scan_get_ndata(data);  
  npixels = raft_scan_get_npixels(data);
  nviews = raft_scan_get_nviews(data);

  lambda = 1;
  
  i = 0;
  
  do{
    
    m = i;

    gsl_vector_set_zero(workspace->art.row);
    
    rad = eval_ray_projection_pet(data, &workspace->backwork.projwork, n, att, m, 
				  workspace->art.row, &norma, &dotp);
    
    q = CVXCOMB(0.001, rad, 0.999, dotp);

    if(sqrt(norma) > ZERO)
      {
	L = lambda*(gsl_vector_get(b,m) - q)/norma;
	
	gsl_blas_daxpy(L, workspace->art.row, n);
      }
    
    maxVector(n);

    i++;

  }while(i<ndata);   
  
}

/*+====================================================+
  
  FUNCTION  art_spect_iteration_knownAtt

  One iteration of the ART algorithm in SPECT with
  known attenuation map. 
  
  Input  

  data - scan data . See <raft_scan_t>.
  workspace - workspace for iterative methods. See <raft_iterative_t>.
  b - SPECT projection data
  att - attenuation map
  k - iteration number
    
  Output 

  n - next estimation
  
  +====================================================+
*/

void art_spect_iteration_knownAtt(raft_scan_t *data,	   
				  raft_iterative_t *workspace,
				  gsl_vector *n,
				  gsl_vector *b,
				  gsl_vector *att,
				  int k)
{
  double rad, norma, L, lambda, dotp, q;
  int ndata, npixels, nviews, i, m;
    
  ndata = raft_scan_get_ndata(data);  
  npixels = raft_scan_get_npixels(data);
  nviews = raft_scan_get_nviews(data);

  lambda = 1;
  
  i = 0;
  
  do{
    
    m = i;
    
    gsl_vector_set_zero(workspace->art.row);
    
    rad = eval_ray_projection_spect(data, &workspace->backwork.projwork, n, att, 
				    m, workspace->art.row, &norma, &dotp);

    q = CVXCOMB(0.001, rad, 0.999, dotp);

    if(sqrt(norma) > ZERO)
      {
	L = lambda*(gsl_vector_get(b,m) - q)/norma;
	
	gsl_blas_daxpy(L, workspace->art.row, n);
      }

    maxVector(n);
    
    i++;

  }while(i<ndata); 
  
}

/*+====================================================+
  
  FUNCTION  art_xfct_iteration_knownAtt

  One iteration of the ART algorithm in fluorescence
  with known attenuation map. 
  
  Input  

  data - scan data . See <raft_scan_t>.
  workspace - workspace for iterative methods . See <raft_iterative_t>.
  b - XFCT projection data
  attT - transmission attenuation map
  attF - XFCT attenuation map
  k - iteration number
    
  Output 

  n - next estimation
  
  +====================================================+
*/

void art_xfct_iteration_knownAtt(raft_scan_t *data,	   
				 raft_iterative_t *workspace,
				 gsl_vector *n,
				 gsl_vector *b,
				 gsl_vector *attT,
				 gsl_vector *attF,
				 int k)
{
  double rad, norma, L, lambda, dotp, q;
  int ndata, npixels, nviews, i, m;
  
  ndata = raft_scan_get_ndata(data);  
  npixels = raft_scan_get_npixels(data);
  nviews = raft_scan_get_nviews(data);
  
  lambda = 2;
  
  i = 0;
  
  do{
    
    m = i;
    
    gsl_vector_set_zero(workspace->art.row);
    
    rad = eval_ray_projection_xfct(data, &workspace->backwork.projwork, n, attT, attF, 
				   m, workspace->art.row, &norma, &dotp);

    q = CVXCOMB(0.001, rad, 0.999, dotp);

    if(sqrt(norma) > ZERO)
      {
	L = lambda*(gsl_vector_get(b,m) - q)/norma;
	
	gsl_blas_daxpy(L, workspace->art.row, n);
      }
    
    maxVector(n);
    
    i++;

  }while(i<ndata);   
}

/*+====================================================+
  
  FUNCTION  em_ct_iteration_knownAtt

  One iteration for the CT-EM algorithm. 
  
  Input  

  data - scan data . See <raft_scan_t>.
  workspace - workspace for iterative methods . See <raft_iterative_t>.
  p - previous estimation
  b - emission data
  
  Output 

  n - next estimation
  
  +====================================================+
*/

void em_ct_iteration_knownAtt(raft_scan_t *data,
			       raft_iterative_t *workspace,
			       gsl_vector *p,
			       gsl_vector *n,
			       gsl_vector *b)
{
  div_t D;
  int i, angle, ray;
  double bb, rr, lh, t,th;
 
  
  raft_projection_radon(data,
            &workspace->backwork.projwork,
            p,
            workspace->em.radon);
    
  /*lh = eval_loglikelihood_emission(&workspace->backwork.projwork,
    workspace->em.radon,
    workspace->backwork.projwork.likelihood.lograd,
    b);
  */

   
  for(i=0; i < raft_scan_get_ndata(data); i++)
    {
      bb = gsl_vector_get(b, i);
      rr = gsl_vector_get(workspace->em.radon, i);

      //gsl_vector_set(workspace->em.data, i, bb/(rr+0.01));

     
      D = div(i, raft_scan_get_nrays(data));
      angle = D.quot;
      ray = D.rem;
      t = -1.0 + ray*(2.0/raft_scan_get_nrays(data));
     
      //if(fabs(rr) < 1e-15)
      if( fabs(t) > sqrt(2)/2 )
          gsl_vector_set(workspace->em.data, i, 0.0);
      else
          gsl_vector_set(workspace->em.data, i, bb/rr);
     
      //fprintf(stderr,"%f\n",workspace->em.data);

    }


  //abort;

  raft_backp(data, 
	     workspace->em.data, 
	     n);  
  
  gsl_vector_mul(n, p);
  gsl_vector_div(n, workspace->em.back);
}

/*+====================================================+
  
  FUNCTION  em_pet_iteration_knownAtt

  One iteration for the PET-EM algorithm with known
  attenuation map. 
  
  Input  

  data - scan data . See <raft_scan_t>.
  workspace - workspace for iterative methods . See <raft_iterative_t>.
  p - previous estimation
  b - emission data
  att - attenuation map
  
  Output 

  n - next estimation
  
  +====================================================+
*/

void em_pet_iteration_knownAtt(raft_scan_t *data,
			       raft_iterative_t *workspace,
			       gsl_vector *p,
			       gsl_vector *n,
			       gsl_vector *b,
			       gsl_vector *att)
{
  int i;
  double bb, rr, lh;
  
  raft_projection_radon_pet(data, 
			    &workspace->backwork.projwork, 
			    p, 
			    att,
			    workspace->em.radon);

  
  lh = eval_loglikelihood_emission(&workspace->backwork.projwork,
				   workspace->em.radon,
				   workspace->backwork.projwork.likelihood.lograd,
				   b);
  
  for(i=0; i < raft_scan_get_ndata(data); i++)
    {
      bb = gsl_vector_get(b, i);
      rr = gsl_vector_get(workspace->em.radon, i);
      
      if(fabs(rr) < ZERO)
	gsl_vector_set(workspace->em.data, i, 0);
      else
	gsl_vector_set(workspace->em.data, i, bb/rr);
    }

  raft_backp_attenuated_pet(data, 
			    &workspace->backwork,
			    workspace->em.data, 
			    att,
			    n);  
  
  gsl_vector_mul(n, p);
  gsl_vector_div(n, workspace->em.back);
}

/*+====================================================+
  
  FUNCTION  em_spect_iteration_knownAtt

  One iteration for the SPECT-EM algorithm with known
  attenuation map. 
  
  Input  

  data - scan data . See <raft_scan_t>.
  workspace - workspace for iterative methods . See <raft_iterative_t>.
  p - previous estimation
  b - SPECT data
  att - attenuation map
  
  Output 
  
  n - next estimation  
  
  +====================================================+
*/

void em_spect_iteration_knownAtt(raft_scan_t *data,
				 raft_iterative_t *workspace,
				 gsl_vector *p,
				 gsl_vector *n,
				 gsl_vector *b,
				 gsl_vector *att)
{
  int i;
  double bb, rr, lh;

  raft_projection_radon_spect(data, 
			      &workspace->backwork.projwork, 
			      p, 
			      att,
			      workspace->em.radon);
  
  lh = eval_loglikelihood_emission(&workspace->backwork.projwork,
				   workspace->em.radon,
				   workspace->backwork.projwork.likelihood.lograd,
				   b);
  
  for(i=0; i < raft_scan_get_ndata(data); i++)
    {
      bb = gsl_vector_get(b, i);
      rr = gsl_vector_get(workspace->em.radon, i);

      if(fabs(rr) < ZERO)
	gsl_vector_set(workspace->em.data, i, 0);
      else
	gsl_vector_set(workspace->em.data, i, bb/rr);
    }

  raft_backp_attenuated_spect(data, 
			      &workspace->backwork,
			      workspace->em.data, 
			      att,
			      n);  
  
  gsl_vector_mul(n, p);
  gsl_vector_div(n, workspace->em.back);
}

/*+====================================================+
  
  FUNCTION  em_fluor_iteration_knownAtt

  One iteration for the fluorescence-EM algorithm with 
  known attenuation map. 
  
  Input  

  data - scan data . See <raft_scan_t>.
  workspace - workspace for iterative methods . See <raft_iterative_t>.
  p - previous estimation
  b - XFCT data
  attT - transmission attenuation map
  attF - XFCT attenuation map
  
  Output 

  n - next estimation
  
  +====================================================+
*/

void em_fluor_iteration_knownAtt(raft_scan_t *data,
				 raft_iterative_t *workspace,
				 gsl_vector *p,
				 gsl_vector *n,
				 gsl_vector *b,
				 gsl_vector *attT,
				 gsl_vector *attF)
{
  double lh;

  workspace->backwork.projwork.likelihood.ratio = gsl_blas_dnrm2(b);

  raft_projection_radon_xfct(data, 
			     &workspace->backwork.projwork, 
			     p, 
			     attF,
			     attT,
			     workspace->em.radon);

  lh = eval_loglikelihood_emission(&workspace->backwork.projwork,
				   workspace->em.radon,
				   workspace->backwork.projwork.likelihood.lograd,
				   b);
  
  workspace->backwork.projwork.likelihood.ratio = 1;

  gsl_blas_dcopy(b, workspace->em.data);
  
  gsl_vector_mul(workspace->em.data, 
		 workspace->backwork.projwork.likelihood.invrad);

  raft_backp_attenuated_xfct(data, 
			     &workspace->backwork,
			     workspace->em.data, 
			     attT,
			     attF,
			     n);  
  
  gsl_vector_mul(n, p);
  gsl_vector_div(n, workspace->em.back);  
}

/*+====================================================+
  
  FUNCTION  em_ct_iteration_grad

  One iteration for the Gradient-type transmission EM 
  algorithm. 
  
  Input  

  data - scan data . See <raft_scan_t>.
  workspace - workspace for iterative methods . See <raft_iterative_t>.
  p - previous estimation
  ph - transmission data (number of photons)
  
  Output 

  n - next estimation

  _Remark_ 

  Gradient algorithm, by K.Lange, M.Bahn & R.Little  "A theoretical 
  study of some maximum likelihood algorithms for 
  emission and transmission tomography", in 
  IEEE, Trans. on Medical Imaging., vol 6, 1987.
  
  +====================================================+
*/

void em_ct_iteration_grad(raft_scan_t *data,
			  raft_iterative_t *workspace,
			  gsl_vector *p,
			  gsl_vector *n,
			  gsl_vector *ph)
{
  raft_projection_monodata(data, p, 
			   workspace->em.ephotons,
			   workspace->em.radon);

  raft_backp(data, workspace->em.ephotons, n);
  raft_backp(data, ph, workspace->em.back);
  
  gsl_vector_mul(n, p);
  gsl_vector_div(n, workspace->em.back);
}


/*+====================================================+
  
  FUNCTION  em_ct_iteration_convex

  One iteration for the Lange-Fessler transmission EM 
  algorithm. 
  
  Input  

  data - scan data . See <raft_scan_t>.
  workspace - workspace for iterative methods . See <raft_iterative_t>.
  p - previous estimation
  ph - transmission data (number of photons)
  
  Output 

  n - next estimation

  _Remark_ 

  Convex algorithm, by K.Lange & J.A.Fessler  "Globally convergent
  Algorithms for maximum a posteriori transmission tomography", in 
  IEEE, Trans. on Image Proc., vol 4, no 19, 1995.
  
  +====================================================+
*/

void em_ct_iteration_convex(raft_scan_t *data,
			    raft_iterative_t *workspace,
			    gsl_vector *p,
			    gsl_vector *n,
			    gsl_vector *ph)
{
  raft_projection_monodata(data, p, workspace->em.ephotons, data->q);
  
  gsl_vector_sub(workspace->em.ephotons, ph);
  raft_backp(data, workspace->em.ephotons, n);
  
  gsl_vector_add(workspace->em.ephotons, ph);
  gsl_vector_mul(workspace->em.ephotons, data->q);
  raft_backp(data, workspace->em.ephotons, workspace->em.back);
  
  gsl_vector_mul(n, p);
  gsl_vector_div(n, workspace->em.back);
  gsl_vector_add(n, p);
}

/*+====================================================+
  
  FUNCTION  em_ct_iteration_standard

  One iteration for the Lange-Carson transmission EM 
  algorithm. 
  
  Input  

  data - scan data . See <raft_scan_t>.
  workspace - workspace for iterative methods . See <raft_iterative_t>.
  p - previous estimation
  ph - transmission data (number of photons)
  
  Output 

  n - next estimation

  _Remark_ 

  Standard EM algorithm, by K.Lange & R.Carson  "EM 
  reconstruction algorithms for emission and transmission
  tomography", in Journal of Computer Assisted Tomography, 
  vol 8, no 2, 1984.
  
  +====================================================+
*/

void em_ct_iteration_standard(raft_scan_t *data,
			      raft_iterative_t *workspace,
			      gsl_vector *p,
			      gsl_vector *n,
			      gsl_vector *ph)
{
  /* fazer de novo !!! */
}


/*+====================================================+
  
  FUNCTION  kun_set_workspace

  Set workspace for Kunyansky's method.
  
  Input  

  data - scan data . See <raft_scan_t>.
  workspace - workspace for iterative methods . See <raft_iterative_t>.
  weight - weight matrix
  p - projection data (attenuated radon transform)
  acctype - acceleration type; RAFT_KUNSTD or RAFT_KUNGSJB
  
  Ouput 
  
  alpha - line search parameter 
   
  +====================================================+
*/

void
kun_set_workspace(raft_scan_t *data,
		  raft_iterative_t *workspace,
		  gsl_vector **weight,
		  gsl_vector *p,
		  int acctype)
{
  int i,j,k,size,nviews;
  double sum, w;
  
  size   = raft_scan_get_size(data);
  nviews = raft_scan_get_nviews(data);
  
  for(i=0;i<size; i++)
    {
      for(j=0; j < size; j++)
	{
	  sum = 0;
	  
	  for(k=0; k < nviews; k++)
	    {
	      sum += raft_phantom_get(weight[k],i,j);
	    }
	      
	  w = sum/nviews;
	  
	  raft_phantom_set(workspace->kun.weight, i, j, w);
	}
    }
  
  
  if(acctype == RAFT_KUNGSJB)
    {
      raft_backp_fbp(data, 	
		     &workspace->backwork,
		     p,
		     workspace->kun.chang);
    }
  
  gsl_vector_div(workspace->kun.chang, workspace->kun.weight);
  
  raft_projection_radon_generalized(data, 
				    &workspace->backwork.projwork,
				    workspace->kun.chang,
				    weight,
				    workspace->kun.rad);

  gsl_vector_sub(workspace->kun.rad, p);
  gsl_vector_scale(workspace->kun.rad, -1.0);
  
  raft_backp_fbp(data, 
		 &workspace->backwork,
		 workspace->kun.rad,
		 workspace->kun.e);
  
  
  gsl_vector_div(workspace->kun.e, workspace->kun.weight);  
  
  workspace->kun.defined = 1;
}


/*+====================================================+
  
  FUNCTION  kun_set_workspace_precond

  Set workspace for preconditioned Kunyansky's method.
  
  Input  

  data - scan data . See <raft_scan_t>.
  workspace - workspace for iterative methods . See <raft_iterative_t>.
  weight - weight matrix
  precond - preconditioning factor
  p - projection data (attenuated radon transform)  
  acctype - acceleration type; RAFT_KUNSTD or RAFT_KUNGSJB
  
  Ouput 
  
  alpha - line search parameter  
   
  +====================================================+
*/

void
kun_set_workspace_precond(raft_scan_t *data,
			  raft_iterative_t *workspace,
			  gsl_vector **weight,
			  gsl_vector *precond,
			  gsl_vector *p,
			  int acctype)
{
  gsl_vector_memcpy(workspace->kun.weight, precond);

  if(acctype == RAFT_KUNGSJB)
    {
      raft_backp_fbp(data, 
		     &workspace->backwork,
		     p,
		     workspace->kun.chang);
	
      gsl_vector_div(workspace->kun.chang, workspace->kun.weight);

      raft_projection_radon_generalized(data, 
					&workspace->backwork.projwork,
					workspace->kun.chang,
					weight,
					workspace->kun.rad);
      
      gsl_vector_sub(workspace->kun.rad, p);
      gsl_vector_scale(workspace->kun.rad, -1.0);
      
      raft_backp_fbp(data, 
		     &workspace->backwork,
		     workspace->kun.rad,
		     workspace->kun.e);
      
      gsl_vector_div(workspace->kun.e, workspace->kun.weight);
    }
  
  workspace->kun.defined = 1;
}


/*+====================================================+
  
  FUNCTION  kun_acceleration_standard

  Standard line search for Kunyansky's method.
  
  Input  

  data - scan data . See <raft_scan_t>.
  workspace - workspace for iterative methods . See <raft_iterative_t>.
  weight - weight matrix
  prev - previous estimation for Kunyansky's method
  direction - search direction
  
  Ouput 
  
  alpha - line search parameter
  
  +====================================================+
*/

void kun_acceleration_standard(raft_scan_t *data,
			       raft_iterative_t *workspace,
			       gsl_vector **weight,
			       gsl_vector *prev,
			       gsl_vector *direction,
			       double *alpha)
{
  double dotp, norm2;
  
  raft_projection_radon_generalized(data, 
				    &workspace->backwork.projwork,
				    direction,
				    weight,
				    workspace->kun.rad);

  raft_backp_fbp(data, 
		 &workspace->backwork,
		 workspace->kun.rad,
		 workspace->kun.v);
  
  gsl_vector_div(workspace->kun.v, workspace->kun.weight);
  
  gsl_blas_ddot(workspace->kun.v, workspace->kun.v, &norm2);
  gsl_blas_ddot(workspace->kun.v, direction, &dotp);

  *alpha = dotp/norm2;
}


/*+====================================================+
  
  FUNCTION  kun_acceleration_gsjb

  GSJB acceleration for Kunyansky's method.
  
  Input  

  data - scan data . See <raft_scan_t>.
  workspace - workspace for iterative methods . See <raft_iterative_t>.
  weight - weight matrix
  prev - previous estimation for Kunyansky's method
  direction - search direction
  
  Ouput 
  
  par - acceleration parameters
  
  Remark 

  GSJB stands for the method of Guy, Sales, Brami-Depaux and Joly-Cabaret
  at "Solutions for Fredholm equations through nonlinear iterative processes"
  in J.Phys. A, 17, 1403-1413 (1984).  
  
  +====================================================+
*/

void
kun_acceleration_gsjb(raft_scan_t *data,
		      raft_iterative_t *workspace,
		      gsl_vector **weight,
		      gsl_vector *prev,
		      gsl_vector *direction,
		      double par[2])
{
  
  double normv2, dotp1, dotp2, dotp3;
  
  gsl_vector_memcpy(workspace->kun.w, workspace->kun.chang);
  gsl_vector_sub(workspace->kun.w, workspace->kun.e);
  
  /* v */
  
  gsl_vector_memcpy(workspace->kun.aux, direction);
  gsl_vector_add(workspace->kun.aux, prev);
  gsl_vector_sub(workspace->kun.aux, workspace->kun.chang);
  
  raft_projection_radon_generalized(data, 
				    &workspace->backwork.projwork,
				    workspace->kun.aux,
				    weight,
				    workspace->kun.rad);

  raft_backp_fbp(data, 
		 &workspace->backwork,
		 workspace->kun.rad,
		 workspace->kun.v);
  
  gsl_vector_div(workspace->kun.v, workspace->kun.weight);
	
  /* dotp's */
  
  gsl_blas_ddot(workspace->kun.v, workspace->kun.v, &normv2);
  gsl_blas_ddot(workspace->kun.v, workspace->kun.chang, &dotp1);
  gsl_blas_ddot(workspace->kun.w, workspace->kun.chang, &dotp2);
  gsl_blas_ddot(workspace->kun.w, workspace->kun.v, &dotp3);
	
  /* h */
  
  gsl_vector_memcpy(workspace->kun.aux, workspace->kun.v);
  gsl_vector_memcpy(workspace->kun.h, workspace->kun.w);
  
  gsl_vector_scale(workspace->kun.aux, dotp3);
  gsl_vector_scale(workspace->kun.h, normv2);
  gsl_vector_sub(workspace->kun.h, workspace->kun.aux);
  
  /* g */
  
  gsl_vector_memcpy(workspace->kun.aux, workspace->kun.v);
  gsl_vector_memcpy(workspace->kun.g, workspace->kun.chang);
  
  gsl_vector_scale(workspace->kun.g, normv2);
  gsl_vector_scale(workspace->kun.aux, dotp1);
  gsl_vector_sub(workspace->kun.g, workspace->kun.aux);
  
  /* z */
  
  gsl_vector_memcpy(workspace->kun.aux, workspace->kun.v);
  gsl_vector_memcpy(workspace->kun.z, workspace->kun.w);
  
  gsl_vector_scale(workspace->kun.aux, dotp2);
  gsl_vector_scale(workspace->kun.z, dotp1);
  gsl_vector_sub(workspace->kun.z, workspace->kun.aux);
  
  /**/
  
  gsl_blas_ddot(workspace->kun.w, workspace->kun.g, &dotp1);
  gsl_blas_ddot(workspace->kun.w, workspace->kun.h, &dotp2);
  par[0] = dotp1/dotp2;
  
  gsl_blas_ddot(workspace->kun.w, workspace->kun.z, &dotp1);
  gsl_blas_ddot(workspace->kun.w, workspace->kun.h, &dotp2);
  par[1] = dotp1/dotp2;  
}


/*######################################################
  Public
  ####################################################*/

/*######################################################
  Section: Allocation
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_iterative_workspace_alloc

  Allocattes workspace for an iterative method. 
  
  Input: 

  ndata - data size
  npixels - number of pixels
  
  Output:

  workspace - method workspace. See <raft_iterative_t>.
  
  Return:

  RAFT_SUCCESS - sucessfull operation
  RAFT_ENOMEM - not enough memory.
  RAFT_EDOM - input domain error
  
  +====================================================+
*/

int raft_iterative_workspace_alloc(raft_scan_t *data,
				   raft_iterative_t *workspace)
{
  int ndata, npixels, nviews;

  ndata   = raft_scan_get_ndata(data);
  npixels = raft_scan_get_npixels(data);
  nviews  = raft_scan_get_nviews(data);

  gsl_set_error_handler_off();
  
  raft_backp_workspace_alloc(data, RAFT_RAMLAK, &workspace->backwork);
  
  workspace->next    = gsl_vector_alloc(npixels);
  workspace->nextDir = gsl_vector_alloc(ndata);
  workspace->previous = gsl_vector_alloc(npixels);
  
  workspace->em.ephotons = gsl_vector_alloc(ndata);
  workspace->em.radon = gsl_vector_alloc(ndata);
  workspace->em.data = gsl_vector_alloc(ndata);
  workspace->em.pi = gsl_vector_alloc(npixels);
  workspace->em.ones = gsl_vector_alloc(ndata);
  workspace->em.back = gsl_vector_alloc(npixels);
  
  workspace->art.canonical = gsl_vector_alloc(ndata);
  workspace->art.row = gsl_vector_alloc(npixels);
  workspace->art.aux = gsl_vector_alloc(npixels);  
    
  workspace->kun.defined = 0;
  workspace->kun.chang = gsl_vector_alloc(npixels);
  workspace->kun.direction = gsl_vector_alloc(npixels);
  workspace->kun.weight = gsl_vector_alloc(npixels);
  workspace->kun.e = gsl_vector_alloc(npixels);
  workspace->kun.v = gsl_vector_alloc(npixels);
  workspace->kun.h = gsl_vector_alloc(npixels);
  workspace->kun.g = gsl_vector_alloc(npixels);
  workspace->kun.z = gsl_vector_alloc(npixels);
  workspace->kun.w = gsl_vector_alloc(npixels);  
  workspace->kun.aux = gsl_vector_alloc(npixels);   
  workspace->kun.rad = gsl_vector_alloc(ndata);
  
  workspace->ramla.row = gsl_vector_alloc(nviews);

  return RAFT_SUCCESS;
}

/*+====================================================+
  
  FUNCTION: raft_iterative_workspace_free

  Frees the workspace for an iterative method. 
  
  Input: 

  workspace - method workspace. See <raft_iterative_t>.
  
  +====================================================+
*/

void raft_iterative_workspace_free(raft_iterative_t *workspace)
{
  raft_backp_workspace_free(&workspace->backwork);

  gsl_vector_free(workspace->next);
  gsl_vector_free(workspace->nextDir);
  gsl_vector_free(workspace->previous);
  
  gsl_vector_free(workspace->em.ephotons);
  gsl_vector_free(workspace->em.radon);
  gsl_vector_free(workspace->em.pi);
  gsl_vector_free(workspace->em.data);
  gsl_vector_free(workspace->em.back);
  gsl_vector_free(workspace->em.ones);

  gsl_vector_free(workspace->art.canonical);
  gsl_vector_free(workspace->art.row);
  gsl_vector_free(workspace->art.aux);
  
  gsl_vector_free(workspace->kun.chang);
  gsl_vector_free(workspace->kun.direction);
  gsl_vector_free(workspace->kun.weight);
  gsl_vector_free(workspace->kun.e);
  gsl_vector_free(workspace->kun.v);
  gsl_vector_free(workspace->kun.w);
  gsl_vector_free(workspace->kun.h);
  gsl_vector_free(workspace->kun.z);
  gsl_vector_free(workspace->kun.g);
  gsl_vector_free(workspace->kun.aux);  
  gsl_vector_free(workspace->kun.rad); 

  gsl_vector_free(workspace->ramla.row);
}

/*######################################################
  Section: Setting workspace
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_iterative_workspace_ramla_set

  Set angle-access for RAMLA method. 
  
  Input: 

  T - vector with angle indexes
  
  Output:
  
  workspace - method workspace
    
  +====================================================+
*/

void raft_iterative_workspace_ramla_set(int sizeT,
					gsl_vector *T,
					raft_iterative_t *workspace)
{
  int i;

  for(i=0; i < sizeT; i++)
    {
      gsl_vector_set(workspace->ramla.row, i, gsl_vector_get(T, i));
    }    

  workspace->ramla.sizeRow = sizeT;
}


/*+====================================================+
  
  FUNCTION: raft_iterative_workspace_set

  Set usual parameters for an iterative method. 
  
  Input: 

  maxiter - maximum number of iterations  
  
  Output:
  
  workspace - method workspace
    
  +====================================================+
*/

void raft_iterative_workspace_set(int maxiter,		  
				  raft_iterative_t *workspace)
{
  workspace->maxiter = maxiter;   

  gsl_vector_set_all(workspace->em.ones, 1.0);
  gsl_vector_set_all(workspace->em.pi, MPI);  
}

/*######################################################
  Section: Iterative methods
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_iterative_em_ct

  Expectation maximization algorithm for transmission
  tomography. 
  
  Input: 

  data - scan data . See <raft_scan_t>.
  r - radon transform
    
  Output:

  workspace - workspace for iterative methods . See <raft_iterative_t>.
  
  +====================================================+
*/

void raft_iterative_em_ct(raft_scan_t *data, 
			  raft_iterative_t *workspace,
			  gsl_vector *r)
{
  int n, maxiter;
  double ini;
  
  raft_backp(data, 
	     workspace->em.ones,
	     workspace->em.back);
  
  maxiter = workspace->maxiter;  
  
  ini = gsl_vector_max(r);

  gsl_vector_set_all(workspace->previous, ini);
  
  n = 0;
  
  while(n<maxiter)
    {
      em_ct_iteration_knownAtt(data, 
			       workspace,
			       workspace->previous,
			       workspace->next, 
			       r);
      
      gsl_vector_memcpy(workspace->previous, workspace->next);
      
      n++;  
    }
}

/*+====================================================+
  
  FUNCTION: raft_iterative_em_pet

  EM algorithm for PET tomography with known 
  attenuation map.  
  
  Input: 

  data - scan data . See <raft_scan_t>.
  b - emission data
  att - attenuation map
    
  Output:

  workspace - workspace for iterative methods . See <raft_iterative_t>.
  
  +====================================================+
*/

void raft_iterative_em_pet(raft_scan_t *data, 
			   raft_iterative_t *workspace,
			   gsl_vector *b,
			   gsl_vector *att)
{
  int n, maxiter;
  double ini;
  
  raft_backp_attenuated_pet(data, 
			    &workspace->backwork,
			    workspace->em.ones,
			    att,
			    workspace->em.back);
  
  maxiter = workspace->maxiter;  

  ini = gsl_vector_max(b);

  gsl_vector_set_all(workspace->previous, ini);

  n = 0;
  
  while(n<maxiter)
    {
      em_pet_iteration_knownAtt(data, 
				workspace,
				workspace->previous,
				workspace->next, 
				b,
				att);
      
      gsl_vector_memcpy(workspace->previous, workspace->next);
      
      n++;  
    }
}


/*+====================================================+
  
  FUNCTION: raft_iterative_em_spect

  EM algorithm for SPECT tomography with known 
  attenuation map. 
  
  Input: 

  data - scan data . See <raft_scan_t>.
  b - emission data
  att - attenuation map
    
  Output:

  workspace - workspace for iterative methods . See <raft_iterative_t>.
  
  +====================================================+
*/

void raft_iterative_em_spect(raft_scan_t *data, 
			     raft_iterative_t *workspace,
			     gsl_vector *b,
			     gsl_vector *att)
{
  int n, maxiter;
  double ini;
  
  raft_backp_attenuated_spect(data, 
			      &workspace->backwork,
			      workspace->em.ones,
			      att,
			      workspace->em.back);
  
  maxiter = workspace->maxiter;  
  
  ini = gsl_vector_max(b);

  gsl_vector_set_all(workspace->previous, ini);
  
  n = 0;
  
  while(n<maxiter)
    {
      em_spect_iteration_knownAtt(data, 
				  workspace,
				  workspace->previous,
				  workspace->next, 
				  b,
				  att);
      
      gsl_vector_memcpy(workspace->previous, workspace->next);
      
      n++;  
    }
}


/*+====================================================+
  
  FUNCTION: raft_iterative_em_xfct

  EM algorithm for XFCT tomography with known 
  attenuation map. 
  
  Input: 

  data - scan data . See <raft_scan_t>.
  b - emission data
  attT - transmission attenuation map
  attF - XFCT attenuation map

  Output:

  workspace - workspace for iterative methods . See <raft_iterative_t>.
  
  +====================================================+
*/

void raft_iterative_em_xfct(raft_scan_t *data, 
				    raft_iterative_t *workspace,
				    gsl_vector *b,
				    gsl_vector *attT,
				    gsl_vector *attF)
{
  int n, maxiter;
  double ini;

  raft_backp_attenuated_xfct(data, 
			     &workspace->backwork,
			     workspace->em.ones,
			     attT,
			     attF,
			     workspace->em.back);

  maxiter = workspace->maxiter;  
  
  ini = gsl_vector_max(b);

  gsl_vector_set_all(workspace->previous, ini);
 
  n = 0;
  
  do{
    
    em_fluor_iteration_knownAtt(data, 
				workspace,
				workspace->previous,
				workspace->next, 
				b,
				attT,
				attF);
    
    gsl_vector_memcpy(workspace->previous, workspace->next);
    
    n++;  
    
  }while( n < maxiter );
}

/*+====================================================+
  
  FUNCTION: raft_iterative_ramla_ct

  RAMLA algorithm for CT. 
  
  (Row-action maximum likelihood algorithm)

  Input: 

  data - scan data . See <raft_scan_t>.
  b - transmission data
  
  Output:

  workspace - workspace for iterative methods . See <raft_iterative_t>.
  
  +====================================================+
*/

void raft_iterative_ramla_ct(raft_scan_t *data, 
			     raft_iterative_t *workspace,
			     gsl_vector *b)
{
  int k, maxiter;
  double ini;

  maxiter = workspace->maxiter;  
  
  gsl_blas_ddot(b, workspace->em.ones, &ini);
  gsl_vector_set_all(workspace->previous, 1/ini);
  
  k = 0;
  
  do{
    
    ramla_ct_iteration(data, 
		       workspace,
		       workspace->previous,
		       workspace->next, 
		       b,
		       k);
    
    gsl_blas_dcopy(workspace->next, workspace->previous);

    k++;  
    
  }while( k < maxiter );
}

/*+====================================================+
  
  FUNCTION: raft_iterative_ramla_xfct

  RAMLA algorithm for XFCT tomography with known 
  attenuation map. 
  
  Input: 

  data - scan data . See <raft_scan_t>.
  b - emission data
  attT - transmission attenuation map
  attF - XFCT attenuation map

  Output:

  workspace - workspace for iterative methods . See <raft_iterative_t>.
  
  +====================================================+
*/

/*
void raft_iterative_ramla_xfct(raft_scan_t *data, 
			       raft_iterative_t *workspace,
			       gsl_vector *b,
			       gsl_vector *attT,
			       gsl_vector *attF)
{
  int n, maxiter;
  double ini;

  projection_workspace_memcpy(&workspace->projwork, 
			      &workspace->backwork.projwork);
    
  maxiter = workspace->maxiter;  
  
  gsl_blas_ddot(b, workspace->em.ones, &ini);
  gsl_vector_set_all(workspace->previous, 1/ini);
  
  n = 0;
  
  do{
    
    ramla_fluor_iteration_knownAtt(data, 
				   workspace,
				   workspace->previous,
				   workspace->next, 
				   b,
				   attT,
				   attF);
    
    gsl_vector_memcpy(workspace->previous, workspace->next);
    
    n++;  
    
  }while( n < maxiter );
}
*/

/*+====================================================+
  
  FUNCTION: raft_iterative_art_ct

  Algebraic reconstruction technique for standard CT.

  Input: 

  data - scan data. See <raft_scan_t>.
  p - projection vector (Radon transform)
  
  Output:

  workspace - workspace for iterative methods. See <raft_iterative_t>.  
  
  +====================================================+
*/

void raft_iterative_art_ct(raft_scan_t *data, 
			   raft_iterative_t *workspace,
			   gsl_vector *p) 
{
  int n, maxiter, average;
  
  average = gsl_blas_dasum(p)/raft_scan_get_ndata(data);

  gsl_vector_set_all(workspace->next, average);
  
  maxiter = workspace->maxiter;
  
  n = 0;

  do{
    
    art_ct_iteration(data, 
		     workspace,
		     workspace->next,		    
		     p,
		     n);          
    n++;  
      
  }while(n < maxiter);
}


/*+====================================================+
  
  FUNCTION: raft_iterative_art_pet

  Algebraic reconstruction technique for PET
  with known attenuation map.

  Input: 

  data - scan data. See <raft_scan_t>.
  p - projection vector (Radon transform)
  att - attenuation map
  
  Output:

  workspace - workspace for iterative methods. See <raft_iterative_t>.  
  
  +====================================================+
*/

void raft_iterative_art_pet(raft_scan_t *data, 
			    raft_iterative_t *workspace,
			    gsl_vector *p,
			    gsl_vector *att) 
{
  int n, maxiter, average;
  
  raft_weight_pet(data, &workspace->backwork.projwork, att);

  average = gsl_blas_dasum(p)/raft_scan_get_ndata(data);

  gsl_vector_set_all(workspace->next, average);
  
  maxiter = workspace->maxiter;
  
  n = 0;

  do{
    
    art_pet_iteration_knownAtt(data, 
			       workspace,
			       workspace->next,		    
			       p,
			       att,
			       n);
    
    n++;  
    
  }while(n < maxiter);
}

/*+====================================================+
  
  FUNCTION: raft_iterative_art_spect

  Algebraic reconstruction technique for SPECT
  with known attenuation map.

  Input: 

  data - scan data. See <raft_scan_t>.
  p - projection vector (Radon transform)
  att - attenuation map
  
  Output:

  workspace - workspace for iterative methods. See <raft_iterative_t>.
  
  +====================================================+
*/

void raft_iterative_art_spect(raft_scan_t *data, 
			      raft_iterative_t *workspace,
			      gsl_vector *p,
			      gsl_vector *att) 
{
  int n, maxiter, average;
    
  raft_weight_spect(data, &workspace->backwork.projwork, att);
  
  average = gsl_blas_dasum(p)/raft_scan_get_ndata(data);
  gsl_vector_set_all(workspace->next, average);
  
  maxiter = workspace->maxiter;
  
  n = 0;

  do{
    
    art_spect_iteration_knownAtt(data, 
				 workspace,
				 workspace->next,		    
				 p,
				 att,
				 n);
    
    n++;  
    
  }while(n < maxiter);
}

/*+====================================================+
  
  FUNCTION: raft_iterative_art_xfct

  Algebraic reconstruction technique for XFCT
  with known attenuation map.

  Input: 

  data - scan data . See <raft_scan_t>.
  p - projection vector (Radon transform)
  attT - transmission attenuation map
  attF - XFCT attenuation map
  
  Output:

  workspace - workspace for iterative methods . See <raft_iterative_t>.
  
  +====================================================+
*/

void raft_iterative_art_xfct(raft_scan_t *data, 
			     raft_iterative_t *workspace,
			     gsl_vector *p,
			     gsl_vector *attT,
			     gsl_vector *attF) 
{
  int n, maxiter, average;
  
  average = gsl_blas_dasum(p)/raft_scan_get_ndata(data);
  gsl_vector_set_all(workspace->next, average);
  gsl_vector_set_zero(workspace->next);
  
  raft_weight_xfct(data, &workspace->backwork.projwork, attT, attF);

  maxiter = workspace->maxiter;
  
  n = 0;

  do{
    
    art_xfct_iteration_knownAtt(data, 
				workspace,
				workspace->next,		    
				p,
				attT,
				attF,
				n);
    
    n++;  
    
  }while(n < maxiter);
}

/*+====================================================+
  
  FUNCTION: raft_iterative_kunyansky

  Kunyansky's iterative method with known attenuation 
  factor.

  Input: 

  data - scan data . See <raft_scan_t>.
  p - projection vector (attenuated Radon transform)
  weight - weight matrix
  initial - initial point
  acctype - acceleration type; RAFT_KUNSTD or RAFT_KUNGSJB
      
  Output:

  workspace - workspace for iterative methods . See <raft_iterative_t>.
  
  +====================================================+
*/

void raft_iterative_kunyansky(raft_scan_t *data,
			      raft_iterative_t *workspace,
			      gsl_vector *p,
			      gsl_vector **weight,
			      gsl_vector *initial,
			      int acctype)
{
  int k, maxiter;
  double lh, norm;

  if(!workspace->kun.defined)
    {
      kun_set_workspace(data,
			workspace,
			weight,
			p,
			acctype);
    }

  k = 0;
  
  norm = gsl_blas_dnrm2(p);

  maxiter = workspace->maxiter;

  gsl_vector_memcpy(workspace->next, initial);
  
  do{
    
    workspace->backwork.projwork.likelihood.ratio = norm;
    
    raft_projection_radon_generalized(data, 
				      &workspace->backwork.projwork,
				      workspace->next,
				      weight,
				      workspace->kun.rad);
    
    lh = eval_loglikelihood_emission(&workspace->backwork.projwork,
				     workspace->kun.rad,
				     workspace->backwork.projwork.likelihood.lograd,
				     p);
    
    workspace->backwork.projwork.likelihood.ratio = 1;
    
    gsl_vector_sub(workspace->kun.rad, p);
    gsl_vector_scale(workspace->kun.rad, -1.0);

    
    raft_backp_fbp(data, 
		   &workspace->backwork,
		   workspace->kun.rad,
		   workspace->kun.direction);
    
    gsl_vector_div(workspace->kun.direction, workspace->kun.weight);
    
    switch(acctype)
      {
      case RAFT_KUNGSJB:
	{
	  double par[2], alpha, beta, q;
	  
	  kun_acceleration_gsjb(data,
				workspace,
				weight,
				workspace->next,
				workspace->kun.direction,
				par);
	  if(k>1){
	    alpha = par[0];
	    beta  = par[1];
	    q     = alpha-beta;
	    
	    gsl_vector_add(workspace->next, workspace->kun.direction);	
	    gsl_vector_scale(workspace->next, beta);
	    
	    gsl_vector_memcpy(workspace->kun.aux, workspace->kun.chang);
	    gsl_vector_scale(workspace->kun.aux, q);
	    gsl_vector_add(workspace->next, workspace->kun.aux);
	  }
	  else
	    gsl_vector_add(workspace->next, workspace->kun.direction);
	  
	  break;
	}
      default:
	{
	  gsl_vector_add(workspace->next, workspace->kun.direction);
	}      
      }
    
    k++;
    
  }while(k < maxiter);
}

/*+====================================================+
  
  FUNCTION: raft_iterative_kunyansky_precond

  Preconditioned Kunyansky's iterative method with known attenuation 
  factor.

  Input: 

  data - scan data . See <raft_scan_t>.
  p - projection vector (attenuated Radon transform)
  weight - weight matrix
  initial - initial point 
  precond - preconditioning factor
  acctype - acceleration type; RAFT_KUNSTD or RAFT_KUNGSJB
  
  Output:

  workspace - workspace for iterative methods . See <raft_iterative_t>.

  +====================================================+
*/

void raft_iterative_kunyansky_precond(raft_scan_t *data,
				      raft_iterative_t *workspace,
				      gsl_vector *p,
				      gsl_vector **weight,
				      gsl_vector *initial,
				      gsl_vector *precond,
				      int acctype)
{
  int k, maxiter;
    
  if(!workspace->kun.defined)
    {
      kun_set_workspace_precond(data,
				workspace,
				weight,
				precond,
				p,
				acctype);
    }

  k = 0;
  
  maxiter = workspace->maxiter;

  gsl_vector_memcpy(workspace->next, initial);
  
  do{
    
    raft_projection_radon_generalized(data, 
				      &workspace->backwork.projwork,
				      workspace->next,
				      weight,
				      workspace->kun.rad);

    gsl_vector_sub(workspace->kun.rad, p);
    gsl_vector_scale(workspace->kun.rad, -1.0);
    
    raft_backp_fbp(data, 
		   &workspace->backwork,
		   workspace->kun.rad,
		   workspace->kun.direction);
  
  
    gsl_vector_div(workspace->kun.direction, workspace->kun.weight);
    
    switch(acctype)
      {
      case RAFT_KUNGSJB:
	{
	  double par[2], alpha, beta, q;
	  
	  kun_acceleration_gsjb(data,
				workspace,
				weight,
				workspace->next,
				workspace->kun.direction,
				par);
	  if(k>1){
	    alpha = par[0];
	    beta  = par[1];
	    q     = alpha-beta;
	    
	    gsl_vector_add(workspace->next, workspace->kun.direction);	
	    gsl_vector_scale(workspace->next, beta);
	    
	    gsl_vector_memcpy(workspace->kun.aux, workspace->kun.chang);
	    gsl_vector_scale(workspace->kun.aux, q);

	    gsl_vector_mul(workspace->kun.aux, workspace->next); /*new*/

	    gsl_vector_add(workspace->next, workspace->kun.aux);
	  }
	  else{
	    gsl_vector_add(workspace->next, workspace->kun.direction);
	  }
	  
	  break;
	}
      default:
	{
	  gsl_vector_add(workspace->next, workspace->kun.direction);	  
	}      
      }
    
    k++;
    
  }while(k < maxiter);
  
}


/*+====================================================+
  
  Temp function to compute the largest eigenvalue by
  the power method
  
  Input 

  data  scan data . See <raft_scan_t>.
  p  projection vector (attenuated Radon transform)
  weight  weight matrix
  initial  initial point
  acctype  acceleration type; RAFT_KUNSTD or RAFT_KUNGSJB
      
  Output

  workspace  workspace for iterative methods. See <raft_iterative_t>.
  
  +====================================================+
*/

void eig(raft_scan_t *data,
	 raft_iterative_t *workspace,
	 gsl_vector *p,
	 gsl_vector **weight,
	 gsl_vector *initial,
	 int acctype)
{
  int k, maxiter;
  double alpha, nom, den, eig;
  
  /*
    if(!workspace->kun.defined)
    {
      kun_set_workspace(data,
			workspace,
			weight,
			p,
			acctype);
    }
  */
  
  k = 0;
  
  maxiter = workspace->maxiter;

  gsl_vector_memcpy(workspace->next, initial);
  
  do{
    
    raft_projection_radon_generalized(data, 
				      &workspace->backwork.projwork,
				      workspace->next,
				      weight,
				      workspace->kun.rad);

    raft_backp_fbp(data, 
		   &workspace->backwork,
		   workspace->kun.rad,
		   workspace->kun.direction);
    
    /* gsl_vector_div(workspace->kun.direction, workspace->kun.weight); */
  
    alpha = 1.0/gsl_blas_dnrm2(workspace->kun.direction);

    gsl_blas_dcopy(workspace->kun.direction, workspace->next);
    
    gsl_blas_dscal(alpha, workspace->next);
    
    gsl_blas_ddot(workspace->next, workspace->kun.direction, &nom);
    
    gsl_blas_ddot(workspace->next, workspace->next, &den);
    
    eig = nom/den;

    fprintf(stderr,"%d %lf\n",k,eig);    
    
    k++;
    
  }while(k < maxiter);  
}
