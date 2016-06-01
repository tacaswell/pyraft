#include "raft_scan.h"
#include "ridge.h"

/*######################################################
  Title: Ridges

  Header - <raft/raft_ridge.h>
  Type - <raft_ridge_t>
  $Id: ridge.c,v 1.14 2009-04-24 13:11:18 miqueles Exp $ - Last Update
  ####################################################*/

/*######################################################
  Private
  ####################################################*/

/*+====================================================+
  
  FUNCTION  eval_intprojc

  Compute the dot product of a projection vector with a
  Chebyshev polynomial.
  
  Input  

  data - scan data . See <raft_scan_t>.
  p - projection vector 
  U - Chebyshev poynomial
  s0 - domain
        
  Return 

  The dot product using a quadrature with five points.

  +====================================================+
*/

double eval_intprojc(raft_scan_t *data,
		     gsl_vector *p,
		     gsl_vector *U,
		     gsl_vector *s0)
{
  double value;

  switch(data->raysid)
    {
    case RAFT_CHEB: 
      {
	double dot;
	int nrays;
	
	nrays = raft_scan_get_nrays(data);
	
	gsl_blas_ddot(p, U, &dot);

	value = (dot*MPI)/nrays;

	break;
      }
      
    case RAFT_NREG:
      {
	int size, index, r;
	double P, sum, t, dt;
	
	size = U->size;
	
	sum = 0;      
	P   = 1;
	
	for(r = 0; r < size-1; r++)
	  {
	    t  = gsl_vector_get(s0, r);
	    dt = gsl_vector_get(s0, r+1) - gsl_vector_get(s0,r);
	    
	    index = interp_nearest(data->t, t, data->raysid);	
	    
	    if(SIGN(index)>0) 
	      {
		P = interp_linear(p, data->t, index, t);
	      }
	    
	    sum += P * gsl_vector_get(U,r) * dt;
	  }
	
	value = sum;  

	break;
      }
    
    default: 
      {
	double dot, dt;
	
	dt = gsl_vector_get(s0,1)-gsl_vector_get(s0,0);
	
	gsl_blas_ddot(p, U, &dot);
	
	value = dt*dot;
      }
    }

  return value;
}

/*+====================================================+
  
  FUNCTION  eval_ridgesum

  Compute the ridge sum for the reconstruction formula.
  
  Input  

  data - scan data . See <raft_scan_t>.
  projection - projection vector 
  j - pixel index (direction y)
  k - pixel index (direction x)   
    
  Return 

  Reconstructed pixel value.

  +====================================================+
*/

double eval_ridgesum(raft_scan_t *data,
		     raft_ridge_t *workspace,
		     int j,
		     int k)
{
  int i, index, size;
  double x, y, sum, cost, sint, s;

  size = raft_scan_get_size(data);

  x   = raft_scan_get_x(data,j);
  y   = raft_scan_get_y(data,size-1-k);
  sum = 0; 

  for(i=0; i < raft_scan_get_nviews(data); i++)
    {
      gsl_vector_view ridge;
      
      ridge = gsl_matrix_row(workspace->ridges, i);

      cost = gsl_vector_get(data->costheta,i);
      sint = gsl_vector_get(data->sintheta,i);
      
      s = x*cost + y*sint;
      	  
      index = interp_nearest(data->t, s, data->raysid);	
      
      if(SIGN(index)>0) 
	{
          sum += interp_linear(&ridge.vector, data->t, index, s);
	}
    }
  
  return sum;  
}

/*######################################################
  Public
  ####################################################*/

/*######################################################
  Section: Allocation
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_ridge_workspace_alloc

  Allocattes workspace for a ridge method. 
  
  Input: 
  
  order - order of Chebyshev polynomial expansion 
  size - reconstruction size
  nrays - number of rays
  nviews - number of views

  Output:

  workspace - method workspace
  
  Return:

  RAFT_SUCCESS - sucessfull operation
  RAFT_ENOMEM - not enough memory.
  RAFT_EDOM - input domain error
  
  +====================================================+
*/

int raft_ridge_workspace_alloc(int order,
			       raft_scan_t *data,
			       raft_ridge_t *workspace)
{
  int k,nviews, nrays, size;
  
  if(order<0)
    return RAFT_EDOM;

  nviews = raft_scan_get_nviews(data);
  nrays  = raft_scan_get_nrays(data);
  size   = raft_scan_get_size(data);

  gsl_set_error_handler_off();
  
  workspace->term = (ridge_kth_t *)malloc(sizeof(ridge_kth_t)*order);
  if(!workspace->term)
    return RAFT_ENOMEM;

  for(k = 0; k < order; k++)
    {
      workspace->term[k].A    = gsl_matrix_alloc(nviews,nviews);
      workspace->term[k].invA = gsl_matrix_alloc(nviews,nviews);
      workspace->term[k].U    = gsl_vector_alloc(nrays);
      workspace->term[k].b    = gsl_vector_alloc(nviews);
      workspace->term[k].c    = gsl_vector_alloc(nviews);
    }
  
  workspace->ridges = gsl_matrix_alloc(nviews, nrays);
  if(!workspace->ridges)
    return RAFT_ENOMEM;
  
  workspace->eye = gsl_matrix_alloc(nviews, nviews);
  if(!workspace->eye)
    return RAFT_ENOMEM;
  
  workspace->nrays = nrays;
  workspace->nviews = nviews;
  workspace->order = order;
  
  /*---*/
  
  raft_projection_workspace_alloc(data, &workspace->generalized.projwork);

  workspace->generalized.ones     = gsl_vector_alloc(size*size);
  workspace->generalized.GenRadon = gsl_vector_alloc(nrays*nviews);
  workspace->generalized.g        = gsl_vector_alloc(nrays*nviews);
  workspace->generalized.r        = gsl_vector_alloc(nrays*nviews);
  
  gsl_vector_set_all(workspace->generalized.ones, 1.0);
  
  return RAFT_SUCCESS;
}

/*+====================================================+
  
  FUNCTION: raft_ridge_workspace_free

  Frees workspace for a ridge method. 
  
  Input: 
  
  workspace - method workspace
  
  +====================================================+
*/

void raft_ridge_workspace_free(raft_ridge_t *workspace)
{
  int k, M;

  M=workspace->order;

  gsl_matrix_free(workspace->ridges);
  gsl_matrix_free(workspace->eye);

  for(k = 0; k < M; k++)
    {
      gsl_matrix_free(workspace->term[k].A);
      gsl_matrix_free(workspace->term[k].invA);
      gsl_vector_free(workspace->term[k].b);
      gsl_vector_free(workspace->term[k].c);
      gsl_vector_free(workspace->term[k].U);
    }

  /*--*/

  raft_projection_workspace_free(&workspace->generalized.projwork);
  
  gsl_vector_free(workspace->generalized.ones);
  gsl_vector_free(workspace->generalized.GenRadon);
  gsl_vector_free(workspace->generalized.g);
  gsl_vector_free(workspace->generalized.r);
  
  free(workspace->term);
}

/*######################################################
  Section: Setting workspace
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_ridge_workspace_set

  Set usual parameters for a ridge method. 
  
  Input: 
  
  data - scan data
  workspace - method workspace
  p - radon transform (stored in a vector)
  
  Return:

  RAFT_SUCCESS - Successfull operation
  RAFT_EDOM - Domain error
  
  +====================================================+
*/

int raft_ridge_workspace_set(raft_scan_t *data,
			     raft_ridge_t *workspace,
			     gsl_vector *p)
{
  int i, j, k, order, nviews, nrays;
  double b, a, d, max, dim;
  
  gsl_vector_view P;
  gsl_vector *S;
  gsl_matrix *V;

  dim = sqrt(2)/2;
  max = raft_scan_get_max(data);
  
  if(max > dim)
    return RAFT_EDOM;
    
  nviews = workspace->nviews;
  nrays = workspace->nrays;
  order = workspace->order;
  
  gsl_matrix_set_identity(workspace->eye);
  gsl_matrix_set_identity(workspace->term[0].A);
  gsl_matrix_set_identity(workspace->term[0].invA);

  S = gsl_vector_alloc(nviews);
  V = gsl_matrix_alloc(nviews, nviews);
  
  /*-----------------------------------------*/
  /* coefficients of the Chebyshev expansion */
  
  for(k = 1; k < order; k++)
    {
      eval_array_chebyshev2(k-1, 
			    data->t, 
			    workspace->term[k-1].U);
      
      for(i=0;i<nviews; i++)
	{
	  P = gsl_vector_subvector(p, i*nrays, nrays);
	  
	  b = eval_intprojc(data, 
			    &P.vector, 
			    workspace->term[k-1].U, 
			    data->t);
	  
	  b = b/MPI;

	  gsl_vector_set(workspace->term[k].b, i, b);
	  
	  for(j=0;j<nviews; j++)
	    {
	      if(i==j)
		{
		  a = 1.0;
		}
	      else
		{
		  d = gsl_vector_get(data->theta,i) -
		    gsl_vector_get(data->theta,j) ;
		  
		  a = sin(k*d)/(k*sin(d));
		}
	      
	      gsl_matrix_set(workspace->term[k].A, i, j, a);
	    }
	}
      
      if(data->viewsid != RAFT_REG)
	{
  	  /*--------------------------------
	    Solve Ak * ck = bk numerically
	    -------------------------------*/
	  
	  gsl_matrix_memcpy(workspace->term[k].invA, 
			    workspace->term[k].A);

	  svd_solve(workspace->term[k].invA, 
		    workspace->term[k].b,
		    nviews, 
		    V,
		    S,
		    workspace->term[k].c);
	}
      else
	{
	  /*---------------------------------
	    Solve Ak * ck = bk analytically
	    --------------------------------*/
	  
	  if(k < nviews)
	    {
	      double square;
	      
	      square = ((double) k)/nviews;
	      square = SQR(square);

	      /* invA = square*I*A + 0*invA */
	      
	      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
			     square, 
			     workspace->eye, 
			     workspace->term[k].A,
			     0.0,
			     workspace->term[k].invA);
	    }
	  else
	    {
	      div_t q, test;
	      int m, l, n;
	      double f,g, fg;
	      
	      n = nviews;
	      q = div(k,n);
	      m = q.quot;
	      l = q.rem;
	      test = div(m, 2);
	      
	      /*-------------------------------------------*/
	      /* Analytical solution for equispaced angles */
	      
	      if(test.rem==0)
		{
		  f = ((double) k)/(m*n);
		  g = ((double) l)/(n*(m+1));
		  
		  fg = -f*g;

		  /* invA =  f  *I*I + 0*invA, 
		     invA = -f*g*I*A + 1*invA */
		  
		  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				 f, 
				 workspace->eye, 
				 workspace->eye,
				 0.0,
				 workspace->term[k].invA);
		  
		  if(l==0)
		    {
		      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				     fg, 
				     workspace->eye, 
				     workspace->eye,
				     1.0,
				     workspace->term[k].invA);
		    }
		  else
		    {
		      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				     fg, 
				     workspace->eye, 
				     workspace->term[l].A,
				     1.0,
				     workspace->term[k].invA);
		    }
		  
		}
	      else
		{
		  f = ((double) k)/(n*(m+1));
		  g = ((double) (n-l))/(n*m);

		  fg = f*g;

		  /* invA = f*I, invA = f*g*A + invA */
		  
		  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				 f, 
				 workspace->eye, 
				 workspace->eye,
				 0.0,
				 workspace->term[k].invA);
		  
		  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				 fg, 
				 workspace->eye, 
				 workspace->term[n-l].A,
				 1.0,
				 workspace->term[k].invA);
		}
	    }
	  
	  /* ck = invAk * bk + 0*ck */
	  
	  gsl_blas_dgemv(CblasNoTrans, 1.0, 
			 workspace->term[k].invA,
			 workspace->term[k].b, 
			 0.0, 
			 workspace->term[k].c);
	}
    }  
  
  gsl_vector_free(S);
  gsl_matrix_free(V);

  return RAFT_SUCCESS;  
}

/*+====================================================+
  
  FUNCTION: raft_generalized_ridge_workspace_set

  Set usual parameters for a generalized ridge method. 
  
  Input: 
  
  data - scan data
  workspace - method workspace
  p - generalized radon transform (stored in a vector)
  weigth - weigth for the generalized radon transform
  
  Return:

  RAFT_SUCCESS - Successfull operation
  RAFT_EDOM - Domain error
  
  +====================================================+
*/

int raft_generalized_ridge_workspace_set(raft_scan_t *data,
					 raft_ridge_t *workspace,
					 gsl_vector *p,
					 gsl_vector **weigth)
{
  int i, j, k, order, nviews, nrays;
  double b, a, d, max, dim;
  
  gsl_vector_view G;
  gsl_vector *S;
  gsl_matrix *V;

  dim = sqrt(2)/2;
  max = raft_scan_get_max(data);
  
  if(max > dim)
    return RAFT_EDOM;
    
  nviews = workspace->nviews;
  nrays = workspace->nrays;
  order = workspace->order;
  
  gsl_matrix_set_identity(workspace->eye);
  gsl_matrix_set_identity(workspace->term[0].A);
  gsl_matrix_set_identity(workspace->term[0].invA);

  S = gsl_vector_alloc(nviews);
  V = gsl_matrix_alloc(nviews, nviews);
  
  /*---------------------------------*/
  /* Generalized radon of ones(size) */
  
  raft_projection_radon_generalized(data, 
				    &workspace->generalized.projwork,
				    workspace->generalized.ones,
				    weigth,
				    workspace->generalized.GenRadon);

  raft_projection_radon(data, 
			&workspace->generalized.projwork,
			workspace->generalized.ones,
			workspace->generalized.r);
  
  
  for(j=0; j< nviews * nrays; j++)
    {
      double gr, r;

      gr = gsl_vector_get(workspace->generalized.GenRadon, j);
      r  = gsl_vector_get(workspace->generalized.r, j);
      
      if(fabs(gr)<1e-06)
	gsl_vector_set(workspace->generalized.g, j, 0);
      else
	gsl_vector_set(workspace->generalized.g, j, r/gr);
    }
  
  gsl_vector_mul(workspace->generalized.g, p);

  /*-----------------------------------------*/
  /* coefficients of the Chebyshev expansion */
  
  for(k = 1; k < order; k++)
    {
      eval_array_chebyshev2(k-1, 
			    data->t, 
			    workspace->term[k-1].U);
      
      for(i=0;i<nviews; i++)
	{
	  G  = gsl_vector_subvector(workspace->generalized.g, i*nrays, nrays);
	  
	  b = eval_intprojc(data, 
			    &G.vector, 
			    workspace->term[k-1].U, 
			    data->t);
	  
	  b = b/MPI;

	  gsl_vector_set(workspace->term[k].b, i, b);
	  
	  for(j=0;j<nviews; j++)
	    {
	      if(i==j)
		{
		  a = 1.0;
		}
	      else
		{
		  d = gsl_vector_get(data->theta,i) -
		    gsl_vector_get(data->theta,j) ;
		  
		  a = sin(k*d)/(k*sin(d));
		}
	      
	      gsl_matrix_set(workspace->term[k].A, i, j, a);
	    }
	}
  
      if(data->viewsid != RAFT_REG)
	{
  	  /*--------------------------------
	    Solve Ak * ck = bk numerically
	    -------------------------------*/
	  
	  gsl_matrix_memcpy(workspace->term[k].invA, 
			    workspace->term[k].A);

	  svd_solve(workspace->term[k].invA, 
		    workspace->term[k].b,
		    nviews, 
		    V,
		    S,
		    workspace->term[k].c);
	}
      else
	{
	  /*---------------------------------
	    Solve Ak * ck = bk analytically
	    --------------------------------*/
	  
	  if(k < nviews)
	    {
	      double square;
	      
	      square = ((double) k)/nviews;
	      square = SQR(square);

	      /* invA = square*I*A + 0*invA */
	      
	      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
			     square, 
			     workspace->eye, 
			     workspace->term[k].A,
			     0.0,
			     workspace->term[k].invA);
	    }
	  else
	    {
	      div_t q, test;
	      int m, l, n;
	      double f,g, fg;
	      
	      n = nviews;
	      q = div(k,n);
	      m = q.quot;
	      l = q.rem;
	      test = div(m, 2);
	      
	      /*-------------------------------------------*/
	      /* Analytical solution for equispaced angles */
	      
	      if(test.rem==0)
		{
		  f = ((double) k)/(m*n);
		  g = ((double) l)/(n*(m+1));
		  
		  fg = -f*g;

		  /* invA =  f  *I*I + 0*invA, 
		     invA = -f*g*I*A + 1*invA */
		  
		  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				 f, 
				 workspace->eye, 
				 workspace->eye,
				 0.0,
				 workspace->term[k].invA);
		  
		  if(l==0)
		    {
		      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				     fg, 
				     workspace->eye, 
				     workspace->eye,
				     1.0,
				     workspace->term[k].invA);
		    }
		  else
		    {
		      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				     fg, 
				     workspace->eye, 
				     workspace->term[l].A,
				     1.0,
				     workspace->term[k].invA);
		    }
		  
		}
	      else
		{
		  f = ((double) k)/(n*(m+1));
		  g = ((double) (n-l))/(n*m);
		  
		  /* invA = f*I, invA = -f*g*A + invA */
		  
		  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				 f, 
				 workspace->eye, 
				 workspace->eye,
				 0.0,
				 workspace->term[k].invA);
		  
		  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				 -f*g, 
				 workspace->eye, 
				 workspace->term[n-l].A,
				 1.0,
				 workspace->term[k].invA);
		}
	    }
	  
	  /* ck = invAk * bk + 0*ck */
	  
	  gsl_blas_dgemv(CblasNoTrans, 1.0, 
			 workspace->term[k].invA,
			 workspace->term[k].b, 
			 0.0, 
			 workspace->term[k].c);
	}
    }  
  
  gsl_vector_free(S);
  gsl_matrix_free(V);

  return RAFT_SUCCESS;
}


/*######################################################
  Section: Reconstruction
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_ridge_function

  Compute the ith ridge function.
 
  Input: 

  scan - scan data . See <raft_scan_t>.
  i - view index
  workspace - method workspace
 
  Output:
  
  h - ridge function

  +====================================================+
*/

void raft_ridge_function(raft_scan_t *data, 
			 raft_ridge_t *workspace,
			 int i,
			 gsl_vector *h)
{
  int k;
  double c;  
    
  gsl_vector_set_all(h, 0.0);  
  
  for(k = 1; k < workspace->order; k++)
    {
      c = gsl_vector_get(workspace->term[k].c, i);
      
      gsl_vector_scale(workspace->term[k-1].U, c);

      gsl_vector_add(h, workspace->term[k-1].U);

      gsl_vector_scale(workspace->term[k-1].U, 1/c);
    }
  
  /*teste: evita problemas na fronteira!!! */
  gsl_vector_set(h, 0, gsl_vector_get(h,1));
}

/*+====================================================+
  
  FUNCTION: raft_ridge_reconstruction

  Reconstruct a phantom using ridge functions.
 
  Input: 

  scan - scan data . See <raft_scan_t>.
  workspace - method workspace
  
  Output:
  
  H - reconstructed phantom

  +====================================================+
*/

void raft_ridge_reconstruction(raft_scan_t *data,
			       raft_ridge_t *workspace,
			       gsl_vector *H)
{
  int i, j, k, size, nviews;
  double sum;
  
  size = raft_scan_get_size(data);
  nviews = raft_scan_get_nviews(data);
  
  /*---------------------*/
  /* set ridge functions */

  for(i = 0; i < nviews; i++)
    {
      gsl_vector_view h;

      h = gsl_matrix_row(workspace->ridges, i);

      raft_ridge_function(data, workspace, i, &h.vector);
    }

  /*----------------*/
  /* reconstruction */

  for(j = 0; j < size; j++)
    {
      for(k = 0; k < size; k++)
	{
	  sum = eval_ridgesum(data, workspace, j, k);
	  
	  gsl_vector_set(H, j + k*size, sum);
	}
    }
}
