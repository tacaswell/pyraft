#ifndef _RAFT_ITERATIVE_H_
#define _RAFT_ITERATIVE_H_

#include "raft_scan.h"
#include "raft_phantom.h"
#include "raft_projection.h"
#include "raft_backprojection.h"
#include "raft_param.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

/*######################################################
  Title: Iterative methods
  ####################################################*/

/*+====================================================+
  
  TYPEDEF: raft_iterative_t

  Purpose:
  
  Workspace for iterative methods.   
         
  +====================================================+
*/

typedef struct{
  
  int maxiter;

  raft_proj_t projwork;
  raft_backp_t backwork;

  gsl_vector *previous, *next, *nextDir;

  struct{    
    gsl_vector *radon;
    gsl_vector *data;
    gsl_vector *pi;        
    gsl_vector *back;
    gsl_vector *ones;
    gsl_vector *ephotons;
  }em;

  struct{    
    gsl_vector *canonical;
    gsl_vector *row;
    gsl_vector *aux;  
  }art;

  struct{
    int defined;    
    gsl_vector *chang;
    gsl_vector *v, *h, *g, *z, *e, *w, *aux;
    gsl_vector *rad;
    gsl_vector *direction;
    gsl_vector *weight;
  }kun;

  struct{
    gsl_vector *row;
    int sizeRow;
  }ramla;
  
}raft_iterative_t;   


int 
raft_iterative_workspace_alloc(raft_scan_t *data,
			       raft_iterative_t *workspace);


void
raft_iterative_workspace_free(raft_iterative_t *workspace);



void 
raft_iterative_workspace_set(int maxiter,
			     raft_iterative_t *workspace);


void 
raft_iterative_workspace_ramla_set(int sizeT,
				   gsl_vector *T,
				   raft_iterative_t *workspace);



void
raft_iterative_em_xfct(raft_scan_t *data, 
		       raft_iterative_t *workspace,
		       gsl_vector *b,
		       gsl_vector *attT,
		       gsl_vector *attF);


void
raft_iterative_em_spect(raft_scan_t *data, 
			raft_iterative_t *workspace,
			gsl_vector *b,
			gsl_vector *att);


void
raft_iterative_em_pet(raft_scan_t *data, 
		      raft_iterative_t *workspace,
		      gsl_vector *b,
		      gsl_vector *att);

void 
raft_iterative_em_ct(raft_scan_t *data, 
		     raft_iterative_t *workspace,
		     gsl_vector *r);


void
raft_iterative_ramla_ct(raft_scan_t *data, 
			raft_iterative_t *workspace,
			gsl_vector *b);



void
raft_iterative_art_ct(raft_scan_t *data, 
		      raft_iterative_t *workspace,
		      gsl_vector *p);


void 
raft_iterative_art_pet(raft_scan_t *data, 
		       raft_iterative_t *workspace,
		       gsl_vector *p,
		       gsl_vector *att);


void 
raft_iterative_art_spect(raft_scan_t *data, 
			 raft_iterative_t *workspace,
			 gsl_vector *p,
			 gsl_vector *att);

void 
raft_iterative_art_xfct(raft_scan_t *data, 
			raft_iterative_t *workspace,
			gsl_vector *p,
			gsl_vector *attT,
			gsl_vector *attF);

void 
raft_iterative_kunyansky(raft_scan_t *data,
			 raft_iterative_t *workspace,
			 gsl_vector *p,
			 gsl_vector **weight,
			 gsl_vector *initial,
			 int acctype);


void
raft_iterative_kunyansky_precond(raft_scan_t *data,
				 raft_iterative_t *workspace,
				 gsl_vector *p,
				 gsl_vector **weight,
				 gsl_vector *initial,
				 gsl_vector *precond,
				 int acctype);


void 
eig(raft_scan_t *data,
    raft_iterative_t *workspace,
    gsl_vector *p,
    gsl_vector **weight,
    gsl_vector *initial,
    int acctype);

#endif
