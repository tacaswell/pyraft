#ifndef _RAFT_POLY_H_
#define _RAFT_POLY_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "raft_scan.h"
#include "raft_phantom.h"

/*######################################################
  Title: Polychromatic
  ####################################################*/

typedef struct{

  gsl_vector *y;

}base_t; 

typedef struct{
  
  int ieffe;          /*index of effective energy               */
  int order;          /*polynomial order                        */
  gsl_vector *atten;  /*linear attenuation of reference material*/
  
}BH_t;

/*+====================================================+
  
  TYPEDEF: raft_poly_t

  Purpose:
  
  Workspace for polychromatic reconstruction algorithms. 

  +====================================================+
*/

/*
  Variables:

  nsubs - number of substances
  nergs - number of energies
  y,p - constant vectors
  ssum - spectrum sum
  beamhard - beam hardening data
*/

typedef struct{

  int nsubs;
  int nergs;
  double ssum;
  gsl_vector *y, *p;
  
  BH_t beamhard;
  base_t *base;
  
}raft_poly_t;   


int 
raft_poly_workspace_alloc(raft_scan_t *data,			      
			  raft_poly_t *workspace);


void 
raft_poly_workspace_free(raft_poly_t *workspace);



void
raft_poly_workspace_set(raft_scan_t *data,
			raft_poly_t *workspace,
			gsl_vector *ph);


void 
raft_poly_workspace_set_BH(raft_poly_t *workspace,
			   int ieffe,
			   int order,
			   gsl_vector *atten);


void 
raft_poly_BH_correction(raft_scan_t *data,
			raft_poly_t *workspace,
			gsl_vector *ph,
			gsl_vector *up);


void 
raft_poly_reconstruction(raft_scan_t *data,
			 raft_poly_t *workspace,
			 raft_phantom_t *phantom,
			 gsl_vector *ph);



#endif
