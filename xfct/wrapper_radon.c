#include "raft_projection.h"
#include "raft_backprojection.h"
#include "raft_scan.h"
#include "raft.h"
#include "../raft/raft_image.h"
#include "gsl/gsl_matrix.h"

void oldraft_radon(raft_image sino, raft_image phantom)
{
  int i, j, nrays, nviews, size;
  double max, min, x;

  gsl_vector *f, *p;
  raft_scan_t data;
  raft_proj_t proj;

  size   = raft_matrix_nlines(phantom.data);
  nviews = raft_matrix_ncolumns(sino.data);
  nrays  = raft_matrix_nlines(sino.data);

  max =  sqrt(2.0)/2.0;
  min = -sqrt(2.0)/2.0;

  /* old raft stuff */ 
  
  raft_scan_alloc_data(&data, nviews, nrays, size, min, max, RAFT_XFCT);
  raft_scan_set_data(&data);

  f = gsl_vector_alloc(size*size);
  
  for(i=0; i < size; i++){
    for(j=0; j < size; j++){
      x = raft_matrix_element(phantom.data,i,j);
      raft_phantom_set(f, i, j, x);	
    }
  }

  raft_projection_workspace_alloc(&data, &proj);

  p = gsl_vector_alloc(nrays*nviews);

  raft_projection_radon(&data, &proj, f, p);

  
  for(i=0; i < nrays; i++){
    for(j=0; j < nviews; j++){
      x = raft_scan_get_projection(&data,p,j,i);
      raft_matrix_element(sino.data,i,j) = x;
    }
  }
  
  raft_scan_free_data(&data);
  raft_projection_workspace_free(&proj);

  gsl_vector_free(p);
  gsl_vector_free(f);
}


/*--------------------------------------------------------------*/

void oldraft_radon_xfct(raft_image sino, 
		        raft_image dens,
 		        raft_image trans, 
			raft_image fluor)	
{
  int i,j, nrays, nviews, size;
  double max, min, x;
  
  raft_scan_t data;
  raft_proj_t workspace;

  gsl_vector *dens_, *trans_, *fluor_, *p;

  size   = raft_matrix_nlines(trans.data);
  nviews = raft_matrix_ncolumns(sino.data);
  nrays  = raft_matrix_nlines(sino.data);

  max =  sqrt(2.0)/2.0;
  min = -sqrt(2.0)/2.0;

  dens_  = gsl_vector_alloc(size*size);
  trans_ = gsl_vector_alloc(size*size);
  fluor_ = gsl_vector_alloc(size*size);
  
  for(i=0; i < size; i++){
    for(j=0; j < size; j++){
      x = raft_matrix_element(trans.data,i,j);
      raft_phantom_set(trans_, i, j, x);	

      x = raft_matrix_element(fluor.data,i,j);
      raft_phantom_set(fluor_, i, j, x);	
 
      x = raft_matrix_element(dens.data,i,j);
      raft_phantom_set(dens_, i, j, x);	
    }
  }

  /*
  for(i=0; i < size; i++){
    for(j=0; j < size; j++){
      fprintf(stderr,"%lf ", raft_phantom_get(dens_, i, j) );	
    }
    fprintf(stderr,"\n");
  }
  */

  raft_scan_alloc_data(&data, nviews, nrays, size, min, max, RAFT_XFCT);
  raft_scan_set_data(&data);

  /*-----------------------*/
  /* fluorescence sinogram */
  
  /*
  if(arg.apert_given){
    raft_scan_set_xfctAperture(&data, arg.apert_arg);
  } 
  */ 
   
  p  = gsl_vector_alloc(nrays*nviews);

  raft_projection_workspace_alloc(&data, &workspace);
  
  raft_projection_radon_xfct(&data, &workspace, dens_, fluor_, trans_, p);
    
  for(i=0;i<nviews; i++){
    for(j=0; j<nrays; j++){
      x = raft_scan_get_projection(&data, p, i ,j);
      raft_matrix_element(sino.data,j,i) = x;
    }
  }
  
  raft_scan_free_data(&data);
  raft_projection_workspace_free(&workspace);

  gsl_vector_free(p);
  gsl_vector_free(dens_);
  gsl_vector_free(trans_);
  gsl_vector_free(fluor_); 
}
