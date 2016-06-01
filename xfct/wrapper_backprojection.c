#include "raft_projection.h"
#include "raft_backprojection.h"
#include "raft_scan.h"
#include "raft.h"
#include "../raft/raft_image.h"
#include "gsl/gsl_matrix.h"

void oldraft_backprojection(raft_image back, raft_image sino)
{
  int nrays, nviews, size, i, j;
  double x, min, max;
  
  raft_scan_t data;
  
  gsl_vector *b, *p;

  size   = raft_matrix_nlines(back.data);
  nviews = raft_matrix_ncolumns(sino.data);
  nrays  = raft_matrix_nlines(sino.data);

  max =  sqrt(2.0)/2.0;
  min = -sqrt(2.0)/2.0;

  raft_scan_alloc_data(&data, nviews, nrays, size, min, max, RAFT_XFCT);
  raft_scan_set_data(&data);
 
  b = gsl_vector_alloc(size*size);
  p = gsl_vector_alloc(nrays*nviews);
   
  
  for(i=0; i < nviews; i++){
    for(j=0; j < nrays; j++){
	x = raft_matrix_element(sino.data, j, i);
	raft_scan_set_projection(&data, p, i, j, x);	  	
    }
  }
    
  /* backprojection */
  
  raft_backp(&data, p, b);
  
  for(i=0; i < size; i++){
    for(j=0; j < size; j++){
      x = raft_phantom_get(b,i,j);
      raft_matrix_element(back.data, i, j) = x;	
    }
  }
  
  gsl_vector_free(b);
  gsl_vector_free(p);
  raft_scan_free_data(&data);  
}

/*----------------------------------------------------------*/

void oldraft_backprojection_xfct(raft_image back,
				raft_image sino,
				raft_image trans,
				raft_image fluor)
{
  int status, nrays, nviews, size, i, j;
  double x, min, max;

  raft_scan_t data;
  raft_backp_t workspace;
 
  gsl_vector *back_, *trans_, *fluor_, *proj_;

  size   = raft_matrix_nlines(back.data);
  nviews = raft_matrix_ncolumns(sino.data);
  nrays  = raft_matrix_nlines(sino.data);

  max =  sqrt(2.0)/2.0;
  min = -sqrt(2.0)/2.0;

  raft_scan_alloc_data(&data, nviews, nrays, size, min, max, RAFT_XFCT);
  raft_scan_set_data(&data);
 
  
  back_  = gsl_vector_alloc(size*size);
  proj_  = gsl_vector_alloc(nrays*nviews);
  trans_ = gsl_vector_alloc(size*size);
  fluor_ = gsl_vector_alloc(size*size);

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

  /* backprojection */

  raft_backp_workspace_alloc(&data, RAFT_RAMLAK, &workspace);
  
  raft_backp_attenuated_xfct(&data, 
			     &workspace, 
			     proj_, 
			     trans_, 
			     fluor_, 
			     back_);
  
  for(i=0; i < size; i++){
    for(j=0; j < size; j++){
      x = raft_phantom_get(back_,i,j);
      raft_matrix_element(back.data, i, j) = x;	
    }
  }
  
  
  
  gsl_vector_free(back_);
  gsl_vector_free(proj_);
  gsl_vector_free(trans_);
  gsl_vector_free(fluor_);

  raft_scan_free_data(&data);
  raft_backp_workspace_free(&workspace);
}


/*----------------------------------------------------------------------*/

void oldraft_fbp360(raft_image fbp, raft_image sino)
{
  int i,j,size, nrays, nviews, status;  
  double x, min, max;
  
  gsl_vector *p, *rec;
  
  raft_scan_t data;
  raft_backp_t workspace;
  
  size   = raft_matrix_nlines(fbp.data);
  nviews = raft_matrix_ncolumns(sino.data);
  nrays  = raft_matrix_nlines(sino.data);

  max =  sqrt(2.0)/2.0;
  min = -sqrt(2.0)/2.0;

  raft_scan_alloc_data(&data, nviews, nrays, size, min, max, RAFT_XFCT);
  raft_scan_set_data(&data);
 
  rec = gsl_vector_alloc(size*size);
  p = gsl_vector_alloc(nrays*nviews);

  for(i=0; i < nviews; i++){
    for(j=0; j < nrays; j++){
	x = raft_matrix_element(sino.data, j, i);
	raft_scan_set_projection(&data, p, i, j, x);	  	
    }
  }
  
  /*-------------------------*/
  /* filtered backprojection */
  
  raft_backp_workspace_alloc(&data, RAFT_COSINE, &workspace);
  
  status = raft_backp_fbp(&data, &workspace, p, rec);

  for(i=0; i < size; i++){
    for(j=0; j < size; j++){
      x = raft_phantom_get(rec,i,j);
      raft_matrix_element(fbp.data, i, j) = x;	
    }
  }

  
  gsl_vector_free(rec);
  gsl_vector_free(p);
  
  raft_scan_free_data(&data);
  raft_backp_workspace_free(&workspace);
}


