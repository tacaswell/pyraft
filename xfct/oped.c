#include "raft_scan.h"
#include "raft_oped.h"
#include "oped.h"
#include "math.h"
#include "ridge.h"

/*######################################################
  Title: OPED
  
  Header - <raft/raft_oped.h>
  Type - <raft_oped_t>
  $Id: oped.c,v 1.14 2008-09-24 00:47:52 miqueles Exp $ - Last Update
  ####################################################*/

/*######################################################
  Private
  ####################################################*/

/*+====================================================+
  
  FUNCTION eval_opedsum

  Compute the reconstruction formula for the oped method.
  
  Input 

  data - scan data . See <raft_scan_t>.
  workspace - method workspace 
  j - pixel index (direction y)
  i - pixel index (direction x)   
    
  Return

  Reconstructed pixel value.
  
  +====================================================+
*/

double eval_opedsum(raft_scan_t *data,
		    raft_oped_t *workspace,
		    int j,
		    int i)
{
  int k, K, index, size;
  double x, y, sum, cost, sint, t;

  size = raft_scan_get_size(data);

  x   = raft_scan_get_x(data,j);
  y   = raft_scan_get_y(data,size-1-i);
  K   = workspace->K;
  
  sum = 0;

  for(k=0; k < K; k++)
    {
      cost  = gsl_vector_get(data->costheta,k);
      sint  = gsl_vector_get(data->sintheta,k);

      t = x*cost + y*sint; 

      index = interp_nearest(data->t, t, data->raysid);

      if(SIGN(index)>0)
	{
	  sum += interp_linear(workspace->term[k].q, data->t, index, t);
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
  
  FUNCTION: raft_oped_workspace_alloc

  Allocattes workspace for an oped method. 
  
  Input: 
  
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

int raft_oped_workspace_alloc(int nrays, 
			      int nviews,
			      raft_oped_t *workspace)
{
  int k;

  if(nrays<=0 || nviews<=0)
    return RAFT_EDOM;

  gsl_set_error_handler_off();

  workspace->K = nviews;
  workspace->M = nviews;
  workspace->nrays = nrays;
  workspace->nviews = nviews;
  
  workspace->term = (oped_kth_t *)malloc(sizeof(oped_kth_t)*nviews);
  if(!workspace->term)
    return RAFT_ENOMEM;

  workspace->U      = gsl_matrix_alloc(workspace->M, nrays);
  workspace->filter = gsl_vector_alloc(workspace->M);

  for(k = 0; k < workspace->K; k++)
    {
      workspace->term[k].c = gsl_vector_alloc(workspace->M);
      workspace->term[k].q = gsl_vector_alloc(nrays);
    }
  
  return RAFT_SUCCESS;
}


/*+====================================================+
  
  FUNCTION: raft_oped_workspace_free

  Frees workspace for an oped method. 
  
  Input: 
  
  workspace - method workspace
   
  +====================================================+
*/

void raft_oped_workspace_free(raft_oped_t *workspace)
{
  int k, K;
  
  K=workspace->K;

  gsl_vector_free(workspace->filter);
  gsl_matrix_free(workspace->U);

  for(k = 0; k < K; k++)
    {
      gsl_vector_free(workspace->term[k].c);
      gsl_vector_free(workspace->term[k].q);
    }
  
  free(workspace->term);
}

/*######################################################
  Section: Setting workspace
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_oped_workspace_set

  Set usual parameters for an oped method. 
  
  Input: 

  data- scan data
  workspace - method workspace
  p - radon transform (stored in a vector)
  
  Return:

  RAFT_SUCCESS - Successfull operation
  RAFT_EDOM - Domain error

  +====================================================+
*/

int raft_oped_workspace_set(raft_scan_t *data,
			    raft_oped_t *workspace,
			    gsl_vector *p)
{
  int m, k, nrays, K, M;
  double c, max, dim;
  
  gsl_vector_view P, u;
  
  dim = sqrt(2)/2;
  max = raft_scan_get_max(data);
  
  if(max > dim)
    return RAFT_EDOM;
  
  nrays = workspace->nrays;
  
  K = workspace->K;
  M = workspace->M;
  
  /*--------------------------------*/
  /* Filter & Chebyshev polynomials */
  
  for(m=0; m < M; m++)
    {
      double f;;
      
      f = (m+1)/((double) 2*K);
      f = f * sin(MPI*(m+1)/(2*M))*2*M/(MPI*(m+1));
      
      gsl_vector_set(workspace->filter, m, f);
      
      u = gsl_matrix_row(workspace->U, m);

      eval_array_chebyshev2(m, data->t, &u.vector);
    }

  /*-------------------------------------*/
  /* Coefficients & filtered projections */

  for(k = 0; k < K; k++)
    {
      P = gsl_vector_subvector(p, k*nrays, nrays);

      gsl_vector_set_all(workspace->term[k].q, 0.0);
      
      for(m=0; m < M; m++)
	{
	  u = gsl_matrix_row(workspace->U, m); 

	  c = (2/MPI) * eval_intprojc(data, 
				       &P.vector,
				       &u.vector,
				       data->t);
	  
	  gsl_vector_set(workspace->term[k].c, m, c);

	  if(fabs(c)>1e-07)
	    {
	      c = c * gsl_vector_get(workspace->filter, m);

	      gsl_vector_scale(&u.vector, c);
	      
	      gsl_vector_add(workspace->term[k].q, &u.vector);
	      
	      gsl_vector_scale(&u.vector, 1/c);
	    }
	}

      gsl_vector_set(workspace->term[k].q, 0, 
		     gsl_vector_get(workspace->term[k].q,1));
    }

  return RAFT_SUCCESS;  
}

/*######################################################
  Section: Reconstruction
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_oped_reconstruction

  Reconstruct a phantom using an oped method.
 
  Input: 

  scan - scan data . See <raft_scan_t>.
  workspace - method workspace
  
  Output:
  
  H - reconstructed phantom
  
  +====================================================+
*/

void raft_oped_reconstruction(raft_scan_t *data,
			      raft_oped_t *workspace,
			      gsl_vector *H)
{
  int j, k, size;
  double sum;
  
  size = raft_scan_get_size(data);
      
  for(j = 0; j < size; j++)
    {
      for(k = 0; k < size; k++)
	{
	  sum = eval_opedsum(data, workspace, j, k);
	  
	  gsl_vector_set(H, j + k*size, sum);
	}
    }
}


