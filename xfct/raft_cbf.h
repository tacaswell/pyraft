#ifndef _RAFT_CBF_H_
#define _RAFT_CBF_H_

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>

#include "raft_scan.h"
#include "raft_projection.h"
#include "raft_phantom.h"
#include "raft_errno.h"

/*######################################################
  Title: CBF
  ####################################################*/

/*+====================================================+
  
  TYPEDEF: raft_cbf_t
  
  Purpose:

  Workspace for cbf Radon transforms.   

  cbf - current, backward and forward
  +====================================================+
*/

typedef struct{

  int pixel, nviews;
  double step;

  gsl_vector *attplus;
  gsl_vector *attminus;
  gsl_vector *actplus;
  gsl_vector *actminus;
  gsl_vector *canonical;
  gsl_vector *kernel;
  
  struct{    
    gsl_vector **Dbt;
    gsl_vector **forward;
    gsl_vector **backward;
    gsl_vector **ExpDbt;
    gsl_vector **ExpDbtForward;
    gsl_vector **ExpDbtBackward;    
  }dbt;
  
  struct{
    
    struct{
      gsl_vector *kernelForward;
      gsl_vector *kernelBackward;
      gsl_vector *forward;
      gsl_vector *backward;
      gsl_vector *current;      
    }FixAct;
    
    struct{
      gsl_vector *kernelForward;
      gsl_vector *kernelBackward;
      gsl_vector *forward;
      gsl_vector *backward;
      gsl_vector *current;
    }FixAtt;
    
    struct{      
      
      struct{
	gsl_vector **forward;
	gsl_vector **backward;
	gsl_vector **current;
      }spect;

      struct{
	gsl_vector **forward;
	gsl_vector **backward;
	gsl_vector **current;
      }xfct;
      
      struct{
	gsl_vector *forward;
	gsl_vector *backward;
	gsl_vector *current;
      }pet;
      
    }weight;

  }projection;

}raft_cbf_t;



void
raft_cbf_workspace_alloc(raft_scan_t *data,
			 raft_cbf_t *cbf);



void
raft_cbf_workspace_free(raft_cbf_t *cbf);




void
raft_cbf_radon_pet(raft_scan_t *data, 
		   raft_proj_t *workspace,
		   raft_cbf_t *cbf,
		   gsl_vector *act,
		   gsl_vector *att);
     

void 
raft_cbf_radon_spect(raft_scan_t *data, 
		     raft_proj_t *workspace,
		     raft_cbf_t *cbf,
		     gsl_vector *act,
		     gsl_vector *att);


void 
raft_cbf_radon_xfct(raft_scan_t *data, 
		    raft_proj_t *workspace,
		    raft_cbf_t *cbf,
		    gsl_vector *act,
		    gsl_vector *attT,
		    gsl_vector *attF,
		    int id);

void 
raft_cbf_radon_dbt(raft_scan_t *data, 
		   raft_proj_t *workspace,
		   raft_cbf_t *cbf,
		   gsl_vector *phantom);


void
raft_cbf_set_pet(raft_scan_t *data,
		 raft_cbf_t *cbf,
		 int k,
		 double h,
		 gsl_vector *att,
		 gsl_vector *act);

void
raft_cbf_set_spect(raft_scan_t *data,
		   raft_cbf_t *cbf,
		   int k,
		   double h,
		   gsl_vector *att,
		   gsl_vector *act);

void
raft_cbf_set_xfct(raft_scan_t *data,
		  raft_cbf_t *cbf,
		  int k,
		  double h,
		  gsl_vector *att,
		  gsl_vector *act);

void
raft_cbf_set_dbt(raft_scan_t *data,
		 raft_cbf_t *cbf,
		 int k,
		 double h,
		 gsl_vector *att);


void 
raft_cbf_get_forward_fact(raft_cbf_t *cbf,
			  gsl_vector *rad);


void 
raft_cbf_get_backward_fact(raft_cbf_t *cbf,
			   gsl_vector *rad);


void 
raft_cbf_get_current_fact(raft_cbf_t *cbf,
			  gsl_vector *rad);


void 
raft_cbf_get_forward_fatt(raft_cbf_t *cbf,
			  gsl_vector *rad);


void 
raft_cbf_get_backward_fatt(raft_cbf_t *cbf,
			   gsl_vector *rad);


void 
raft_cbf_get_current_fatt(raft_cbf_t *cbf,
			  gsl_vector *rad);


void 
raft_cbf_get_current_dbt(raft_cbf_t *cbf,
			 gsl_vector **dbt);


void 
raft_cbf_get_forward_dbt(raft_cbf_t *cbf,
			 gsl_vector **dbt);


void
raft_cbf_get_backward_dbt(raft_cbf_t *cbf,
			  gsl_vector **dbt);



#endif
