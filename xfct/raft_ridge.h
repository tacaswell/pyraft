#ifndef _RAFT_RIDGE_H_
#define _RAFT_RIDGE_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "raft_scan.h"
#include "raft_projection.h"

/*######################################################
  Title: Ridges
  ####################################################*/

typedef struct{

  gsl_vector *b, *U, *c;
  gsl_matrix *A, *invA;  

}ridge_kth_t;

/*+====================================================+
  
  TYPEDEF: raft_ridge_t

  Purpose:
  
  Workspace for ridge methods. 

  +====================================================+
*/

typedef struct{

  int order, nrays, nviews;
  gsl_matrix *ridges, *eye;
  ridge_kth_t *term;

  struct{
    raft_proj_t projwork;
    gsl_vector *ones;
    gsl_vector *GenRadon, *g, *r;
    
  }generalized;
  
}raft_ridge_t;   


int 
raft_ridge_workspace_alloc(int order,
			   raft_scan_t *data,
			   raft_ridge_t *workspace);


void 
raft_ridge_workspace_free(raft_ridge_t *workspace);



int
raft_ridge_workspace_set(raft_scan_t *data,
			 raft_ridge_t *workspace,
			 gsl_vector *p);


int 
raft_generalized_ridge_workspace_set(raft_scan_t *data,
				     raft_ridge_t *workspace,
				     gsl_vector *p,
				     gsl_vector **weigth);


void 
raft_ridge_function(raft_scan_t *data, 
		    raft_ridge_t *workspace,		    
		    int i,
		    gsl_vector *h);


void 
raft_ridge_reconstruction(raft_scan_t *data,
			  raft_ridge_t *workspace,
			  gsl_vector *H);



#endif
