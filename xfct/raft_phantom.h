#ifndef _RAFT_PHANTOM_H_
#define _RAFT_PHANTOM_H_

#include "raft_errno.h"
#include "raft_scan.h"
#include "raft_param.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <stdlib.h>
#include <math.h>

/*######################################################
  Title: Phantom
  ####################################################*/

/*######################################################
  Public
  ####################################################*/

typedef struct{

  gsl_vector *vector;

}basis_t; 

/*+====================================================+
  
  TYPEDEF: raft_phantom_t

  Purpose:

  RAFT phantom structure data. 

  +====================================================+
*/

typedef struct
{
  gsl_vector *x, *y, *vector;
  double min, max;  
  int    size, nsubs, etype;
  basis_t *basis;
  
}raft_phantom_t;


raft_phantom_t *
raft_phantom_alloc(void);


int 
raft_phantom_alloc_data(raft_phantom_t *phantom, 
		       int size,
		       double min, 
		       double max);

int 
raft_phantom_alloc_data_basis(raft_phantom_t *phantom, 
			      int nsubs);


void
raft_phantom_free(raft_phantom_t *phantom);


void
raft_phantom_free_data(raft_phantom_t *phantom);


double 
raft_phantom_get_max(raft_phantom_t *phantom);


double 
raft_phantom_get_min(raft_phantom_t *phantom);


int 
raft_phantom_get_size(raft_phantom_t *phantom);


double
raft_phantom_get_y(raft_phantom_t *phantom,
		  int i);

double 
raft_phantom_get_x(raft_phantom_t *phantom,
		  int i);

double
raft_phantom_get_step(raft_phantom_t *phantom);


void
raft_phantom_set(gsl_vector *phantom, 
		 int i,
		 int j, 
		 double z);


int 
raft_phantom_set_basis(raft_phantom_t *phantom, 
		       int k, 
		       int i,
		       int j, 
		       double z);

double 
raft_phantom_get(gsl_vector *phantom, 
		 int i,
		 int j);


gsl_vector *
raft_phantom_get_vector_basis(raft_phantom_t *phantom,
			      int k);


gsl_vector *
raft_phantom_get_vector(raft_phantom_t *P);



void 
raft_phantom_set_default_basis(raft_phantom_t *phantom,
			       raft_scan_t *data,
			       gsl_vector *f);


void 
raft_phantom_define_ghost(raft_phantom_t *p);


void 
raft_phantom_define_ring(raft_phantom_t *p);


void 
raft_phantom_define_elipring(raft_phantom_t *p);


void 
raft_phantom_define_ball(raft_phantom_t *p);


void
raft_phantom_define_rectang(raft_phantom_t *p);


void 
raft_phantom_define_square(raft_phantom_t *p);


void 
raft_phantom_define_shiftsquare(raft_phantom_t *p);


void 
raft_phantom_define_star(raft_phantom_t *p);


void 
raft_phantom_define_sheeplogan(raft_phantom_t *p);


void
raft_phantom_define_xfct_density_miqdep(raft_phantom_t *P);


void
raft_phantom_define_xfct_attenuation_miqdep(raft_phantom_t *P);


void
raft_phantom_define_xfct_Tattenuation_miqdep(raft_phantom_t *P);


void
raft_phantom_define_xfct_density_golosio(raft_phantom_t *P);


void
raft_phantom_define_xfct_attenuation_golosio(raft_phantom_t *P);


void
raft_phantom_define_xfct_Tattenuation_golosio(raft_phantom_t *P);


void
raft_phantom_define_spect_attenuation(raft_phantom_t *P);


void
raft_phantom_define_spect_activity(raft_phantom_t *P);


#endif

