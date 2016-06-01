#include "raft_projection.h"
#include "raft_backprojection.h"
#include "raft_scan.h"
#include "raft.h"
#include "../raft/raft_image.h"
#include "../raft/raft_matrix.h"
#include "gsl/gsl_matrix.h"

void oldraft_akt_xfct(	raft_image rec,
			raft_image sino,
			raft_image trans,
			raft_image fluor,
			int niter,
			int wtype)	
{
  int i,j, size, nrays, nviews, acctype;  
  double x, max, min;
  
  raft_scan_t data;
  raft_iterative_t workspace;
  
  gsl_vector *proj_, *trans_, *fluor_, *ini, **weight, *a;

  max =  sqrt(2.0)/2.0;
  min = -sqrt(2.0)/2.0;

  size   = raft_matrix_nlines(trans.data);
  nviews = raft_matrix_ncolumns(sino.data);
  nrays  = raft_matrix_nlines(sino.data);

  raft_scan_alloc_data(&data, nviews, nrays, size, min, max, RAFT_XFCT);
  raft_scan_set_data(&data);

  ini    = gsl_vector_alloc(size*size);
  trans_ = gsl_vector_alloc(size*size);
  fluor_ = gsl_vector_alloc(size*size);
  proj_  = gsl_vector_alloc(nrays*nviews);
  a      = gsl_vector_alloc(size*size);

  weight = (gsl_vector **)malloc(sizeof(gsl_vector*)*nviews);
  for(i=0; i < nviews; i++)
    weight[i] = gsl_vector_alloc(size*size); 

  /* 'read' input fluorescence sinogram */

  for(i=0; i < nviews; i++){
    for(j=0; j < nrays; j++){
      x = raft_matrix_element(sino.data, j, i);
      raft_scan_set_projection(&data, proj_, i, j, x);
    }
  }

  for(i=0; i < size; i++){
    for(j=0; j < size; j++){
      x = raft_matrix_element(trans.data, i, j);
      raft_phantom_set(trans_, i, j, x);

      x = raft_matrix_element(fluor.data, i, j);
      raft_phantom_set(fluor_, i, j, x);	
    }
  }
  
  /* AKT iterations: based on Radon's inversion */
  
  gsl_vector_set_all(ini, 10e-07);
  
  raft_iterative_workspace_alloc(&data, &workspace);
  raft_iterative_workspace_set(niter, &workspace);
  
  raft_weight_xfct(&data,
		   &workspace.backwork.projwork,
		   trans_,
		   fluor_);
  
  raft_weight_xfct_get(&workspace.backwork.projwork,
		       weight);

  acctype = RAFT_KUNSTD;
  
  if(wtype)
    {
      raft_iterative_kunyansky(&data,
			       &workspace,
			       proj_,
			       weight,
			       ini,
			       acctype);
    }
  else
    {
      gsl_vector_set_all(a, 1.0);

      raft_iterative_kunyansky_precond(&data,
				       &workspace,
				       proj_,
				       weight,
				       ini,
				       a,
				       acctype);
    }
  
  for(i=0;i<size; i++){
    for(j=0; j<size; j++){
      x = MAX(0,raft_phantom_get(workspace.next,i,j));
      raft_matrix_element(rec.data, i, j) = x;
    }
  }
  
  /**/
  
  for(i=0;i<nviews;i++)
    gsl_vector_free(weight[i]);
  
  free(weight);  

  gsl_vector_free(trans_);
  gsl_vector_free(fluor_);
  gsl_vector_free(proj_);
  gsl_vector_free(ini);
  gsl_vector_free(a);
  
  raft_iterative_workspace_free(&workspace);
  raft_scan_free_data(&data);
}

