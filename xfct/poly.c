#include "raft_poly.h"
#include "raft_backprojection.h"
#include "poly.h"

/*######################################################
  Title: Polychromatic

  Header - <raft/raft_poly.h>
  Type - <raft_poly_t>
  $Id: poly.c,v 1.16 2008-09-24 01:37:41 miqueles Exp $ - Last Update
  ####################################################*/

/*######################################################
  Private
  ####################################################*/

/*+====================================================+
  
  FUNCTION
  
  Input

  Output
  
  +====================================================+
*/

/*######################################################
  Public
  ####################################################*/

/*######################################################
  Section: Allocation
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_poly_workspace_alloc
  
  Allocattes workspace for polychromatic algorithms. 

  Input: 

  order - order of approximation.
  data - scan data . See <raft_scan_t>..

  Output:

  workspace - polychromatic workspace . See <raft_poly_t>.

  Return:

  RAFT_EDOM - Domain error.
  RAFT_ENOMEM - Not enough memory.
  RAFT_SUCCESS - Successfull operation.
  
  +====================================================+
*/

int raft_poly_workspace_alloc(raft_scan_t *data,			      
			      raft_poly_t *workspace)
{
  int j,nergs, nsubs, nrays, nviews, ndata;

  nergs = raft_scan_get_nenergies(data);
  nsubs = raft_scan_get_nsubs(data);
  nrays = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);
  ndata = raft_scan_get_ndata(data);
  
  workspace->y = gsl_vector_alloc(nergs);
  if(!workspace->y)
    return RAFT_ENOMEM;
  
  workspace->p = gsl_vector_alloc(nsubs);
  if(!workspace->p)
    return RAFT_ENOMEM;

  workspace->base = (base_t *)malloc(sizeof(base_t)*nergs);
  for(j=0;j<nergs; j++)
    workspace->base[j].y = gsl_vector_alloc(ndata);

  workspace->beamhard.atten = gsl_vector_alloc(nergs);
  if(!workspace->beamhard.atten)
    return RAFT_ENOMEM;
  
  workspace->nergs = nergs;
  workspace->nsubs = nsubs;
  
  return RAFT_SUCCESS;  
}

/*+====================================================+
  
  FUNCTION: raft_poly_workspace_free
  
  Frees workspace for polychromatic algorithms.

  Input: 

  workspace - polychromatic workspace . See <raft_poly_t>.
 
  +====================================================+
*/

void raft_poly_workspace_free(raft_poly_t *workspace)
{
  int i;

  gsl_vector_free(workspace->y);
  gsl_vector_free(workspace->p);  
  gsl_vector_free(workspace->beamhard.atten);

  for(i=0; i<workspace->nergs; i++)
    gsl_vector_free(workspace->base[i].y);

  free(workspace->base);  
}

/*######################################################
  Section: Setting workspace
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_poly_workspace_set

  Set usual parameters for polychromatic algorithms. 
  
  Input: 

  data - scan data
  workspace - method workspace
  p - polychromatic transmission data    

  +====================================================+
*/

void raft_poly_workspace_set(raft_scan_t *data,
			     raft_poly_t *workspace,
			     gsl_vector *p)
{
  int j, nrays, nviews, nergs;
  double s, sum;
  
  nrays = raft_scan_get_nviews(data);
  nviews = raft_scan_get_nviews(data);
  nergs = raft_scan_get_nenergies(data);

  sum = 0;
  for(j=0; j<nergs; j++)
    {
      s = gsl_vector_get(data->energy->spectrum, j);
      sum += s;
    }

  workspace->ssum = sum; 
}


/*+====================================================+
  
  FUNCTION: raft_poly_workspace_set_BH

  Set beam hardening data within workspace.
  
  Input: 

  workspace - method workspace
  ieffe - index of effective energy
  order - polynomial order
  atten - linear attenuation of reference material
  
  +====================================================+
*/

void
raft_poly_workspace_set_BH(raft_poly_t *workspace,
			   int ieffe,
			   int order,
			   gsl_vector *atten)
{
  workspace->beamhard.ieffe = ieffe;
  workspace->beamhard.order = order;

  gsl_vector_memcpy(workspace->beamhard.atten, atten);
}

/*######################################################
  Section: Correction
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_poly_BH_correction

  Beam hardening correction
  
  Input:

  data - scan data
  workspace - method workspace
  ph - polychromatic transmission data
  
  Output:

  up - updated monochromatic transmission data
 
  +====================================================+
*/

void raft_poly_BH_correction(raft_scan_t *data,
			     raft_poly_t *workspace,
			     gsl_vector *ph,
			     gsl_vector *up)
{
  int L, J, O, i, j, k, nergs, nrays, ndata;
  double max, min, dm, m, p, s, xi, sum, pp, mm;
  gsl_vector *M, *a, *S;
  gsl_matrix *A, *A2, *V;
  
  nergs = raft_scan_get_nenergies(data);
  nrays = raft_scan_get_nrays(data);
  ndata = raft_scan_get_ndata(data);
  
  J = workspace->beamhard.ieffe;
  O = workspace->beamhard.order;
  L = ndata;

  O++;
  
  gsl_vector_minmax(ph, &min, &max);

  M  = gsl_vector_alloc(L);
  a  = gsl_vector_alloc(O);
  S  = gsl_vector_alloc(O);
  A  = gsl_matrix_alloc(L,O);
  A2 = gsl_matrix_alloc(L,O);
  V  = gsl_matrix_alloc(O,O);
  
  dm = (max-min)/((double)L);
  
  for(i=0; i < L; i++)
    {
      m = min + i*dm;       
      
      p = 0;
      sum = 0;
 
      for(k=0; k < nergs; k++)
	{
	  s = gsl_vector_get(data->energy->spectrum, k);
	  
	  xi = gsl_vector_get(workspace->beamhard.atten, k)/
	    gsl_vector_get(workspace->beamhard.atten, J);
	  
	  p += s * exp(-m*xi);
	  
	  sum += s;
	}
      
      p = -log(p/sum);

      for(j=0; j < O; j++)
	{
	  pp = pow(p, j);
	  gsl_matrix_set(A, i, j, pp);
	}      

      gsl_vector_set(M, i, m);
    }

  
  /*------------------------------*/
  /* Solve A*a = M: least squares */
  
  gsl_matrix_memcpy(A2, A);

  svd_solve(A2, M, O, V, S, a);
  
  /*---------------------------*/
  /* Updated transmission data */
  
  for(i=0; i < ndata; i++)
    {
      p = gsl_vector_get(ph, i);
            
      mm = 0;
      for(j=0; j < O; j++)
	mm += pow(p,j) * gsl_vector_get(a,j);
      
      gsl_vector_set(up, i, exp(-mm));
    }
  
  gsl_matrix_free(A);
  gsl_matrix_free(V);
  gsl_vector_free(S);
  gsl_vector_free(M);  
  gsl_vector_free(a);  
}


/*+====================================================+
  
  FUNCTION  raft_poly_reconstruction

  Reconstruction of the phantom attenuation basis
  
  Input  

  data - scan data
  workspace - method workspace

  Remark
  
  Apenas função de teste...

  +====================================================+
*/

void raft_poly_reconstruction(raft_scan_t *data,
			      raft_poly_t *workspace,
			      raft_phantom_t *phantom,
			      gsl_vector *ph)
{
  int i, k, j, ndata, nergs, nsubs;
  double d, y, sv, div, S, L, beta, r, xi;
  
  gsl_vector *aux;

  nergs = raft_scan_get_nenergies(data);
  nsubs = raft_scan_get_nsubs(data);
  ndata = raft_scan_get_ndata(data);
  
  /*-------------------*/
  /* SVD decomposition */
  
  if(!data->energy->decomposed)
    {
      gsl_matrix_memcpy(data->energy->U, data->energy->G);
  
      gsl_linalg_SV_decomp_jacobi(data->energy->U, 
				  data->energy->V, 
				  data->energy->S);
    }

  /*----------------*/
  /* Reconstruction */

  aux = gsl_vector_alloc(nsubs);

  gsl_blas_ddot(data->energy->spectrum, data->energy->spectrum, &beta);
  
  xi = workspace->ssum;
  
  beta = xi/beta;
  
  for(i = 0; i < ndata; i++)  
    {
      d = gsl_vector_get(ph, i);
      
      L = log(d);
      
      for(k = 0; k < nergs; k++)
	{
	  S = gsl_vector_get(data->energy->spectrum, k);
	  
	  if(fabs(S)>1e-05)
	    y = - (L + log(beta*S));
	  else
	    y = - (L - 11.512925464970229);
	  
	  gsl_vector_set(workspace->y, k, y);
	  gsl_vector_set(workspace->base[k].y, i, y);
	}
      
      /*-----------------*/
      /* p = pinv(G) * y */
      
      gsl_blas_dgemv(CblasTrans, 1.0, data->energy->U, workspace->y, 
		     0.0, workspace->p);
      
      for(j=0; j < nsubs; j++)
	{
	  sv = gsl_vector_get(data->energy->S,j);
	  
	  if(fabs(sv) < 1e-06)
	    gsl_vector_set(workspace->p,j,0.0);
	  else
	    {
	      div = gsl_vector_get(workspace->p,j)/sv;
	      gsl_vector_set(workspace->p,j,div);
	    }
	}
      
      gsl_blas_dgemv(CblasNoTrans, 1.0, 
		     data->energy->V, workspace->p, 0.0, aux);
      
      gsl_vector_memcpy(workspace->p, aux);
      
      /*----------------*/
      /* Radon matrices */
      
      for(j=0; j < nsubs; j++)
	{
	  r = gsl_vector_get(workspace->p, j);
	  
	  gsl_vector_set(data->energy->basis[j].p, i, r);
	}
    }
  
  gsl_vector_free(aux);

  /*--------------------------------------*/
  /* Reconstruction of each basis element */

  
  for(j = 0; j < nsubs; j++)
    {
      raft_backp_fbp(data, 
		     RAFT_RAMLAK, 
		     data->energy->basis[j].p, 
		     phantom->basis[j].vector);
    }
  
}
