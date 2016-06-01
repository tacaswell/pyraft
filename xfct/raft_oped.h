#ifndef _RAFT_OPED_H_
#define _RAFT_OPED_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "raft_scan.h"

/*######################################################
  Title: OPED
  ####################################################*/

typedef struct{

  gsl_vector *c;    
  gsl_vector *q;

}oped_kth_t;

/*+====================================================+
  
  TYPEDEF: raft_oped_t

  Purpose:
  
  Workspace for an OPED reconstruction method. 

  oped - orthogonal polynomial expansion on the disk
  +====================================================+
*/

typedef struct{

  int K, M, nrays, nviews;
  gsl_vector *filter;
  gsl_matrix *U;
  oped_kth_t *term;  
  
}raft_oped_t;   



int 
raft_oped_workspace_alloc(int nrays, 
			  int nviews,
			  raft_oped_t *workspace);


void 
raft_oped_workspace_free(raft_oped_t *workspace);



int
raft_oped_workspace_set(raft_scan_t *data,
			raft_oped_t *workspace,
			gsl_vector *p);


void 
raft_oped_reconstruction(raft_scan_t *data,
			 raft_oped_t *workspace,
			 gsl_vector *H);



#endif
