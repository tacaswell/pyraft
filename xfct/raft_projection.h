#ifndef _RAFT_PROJECTION_H_
#define _RAFT_PROJECTION_H_

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_exp.h> 
#include "raft_scan.h"
#include "raft_phantom.h"
#include "raft_errno.h"

#define MAX_NTHREADS 20000

/*######################################################
  Title: Projections
  ####################################################*/

/*+====================================================+
  
  TYPEDEF: raft_proj_t

  Purpose:

  RAFT projection structure data.   
  +====================================================+
*/

typedef struct{

  int npixels;
  int ndata;
  int nrays;
  int nviews;
  int size;
  double AvAtt;

  struct{
    double eps, ratio;
    gsl_vector *invrad;
    gsl_vector *lograd;
    gsl_vector *ones;
    gsl_vector *radon;
    gsl_vector *exprad, *exprad2;
    
  }likelihood;
  
  struct{
    gsl_vector *integrand;
  }generalized;
  
  struct{    
    gsl_vector **Dbt;
    gsl_vector **ExpDbt;
  }dbt;

  struct{
    int defined;
    gsl_vector *Radon;
    gsl_vector *ExpRadon;
    gsl_vector *weight;
    gsl_vector *integrand;   
  }pet;
  
  struct{    
    int defined;        
    gsl_vector *integrand;
    gsl_vector **weight;
    gsl_vector **Dbt;
    
  }spect;
  
  struct{    
    double apert;
    int defined, definedDbtT;
    gsl_vector *ones;
    gsl_vector *integrand[MAX_NTHREADS];
    gsl_vector **weight, **weightP, **weightF, **LogWeightF, **LogWeightP;
    gsl_vector **DbtT;
    gsl_vector **DbtF;
    gsl_vector **ExpDbtT;
    gsl_vector **ExpDbtF;  
  }xfct;
  
}raft_proj_t;


void 
raft_projection_radon(raft_scan_t *data, 
		      raft_proj_t *workspace,
		      gsl_vector *phantom,
		      gsl_vector *radon);


void
raft_projection_radon_ert(raft_scan_t *data, 
			  raft_proj_t *workspace,
			  double *att,
			  gsl_vector *act,
			  gsl_vector *radon);


void
raft_projection_radon_pet(raft_scan_t *data, 
			  raft_proj_t *workspace,
			  gsl_vector *act,
			  gsl_vector *att,
			  gsl_vector *radon);

void 
raft_projection_radon_spect(raft_scan_t *data, 
			    raft_proj_t *workspace,
			    gsl_vector *act,
			    gsl_vector *att,
			    gsl_vector *radon);


void 
raft_projection_radon_xfct(raft_scan_t *data, 
			   raft_proj_t *workspace,
			   gsl_vector *act,
			   gsl_vector *attF,
			   gsl_vector *attB,
			   gsl_vector *radon);

void 
raft_projection_radon_generalized(raft_scan_t *data, 
				  raft_proj_t *workspace,
				  gsl_vector *f,
				  gsl_vector **weight,
				  gsl_vector *radon);


void
raft_projection_radon_dbt(raft_scan_t *data, 
			  raft_proj_t *workspace,
			  gsl_vector *phantom,
			  gsl_vector **dbt);



void
raft_projection_radon_view(raft_scan_t *data, 
			   gsl_vector *phantom,
			   int j,
			   gsl_vector *view);


void 
raft_projection_radon_pet_view(raft_scan_t *data, 
			       raft_proj_t *workspace,
			       gsl_vector *act,
			       gsl_vector *att,
			       int j,
			       gsl_vector *view);


void
raft_projection_radon_spect_view(raft_scan_t *data, 
				 raft_proj_t *workspace,
				 gsl_vector *act,
				 gsl_vector *att,
				 int j,
				 gsl_vector *view);


void 
raft_projection_radon_xfct_view(raft_scan_t *data, 
					raft_proj_t *workspace,
					gsl_vector *act,
					gsl_vector *attF,
					gsl_vector *attT,
					int j,
					gsl_vector *view);


int
raft_projection_radon_dbt_view(raft_scan_t *data, 
			       raft_proj_t *workspace,
			       int j,
			       gsl_vector *phantom,
			       gsl_vector *dbt);


double
raft_projection_radon_dbt_get(raft_scan_t *data,
			      gsl_vector **dbt,
			      int i,
			      int j);


void
raft_projection_monodata(raft_scan_t *data, 
			 gsl_vector *phantom,
			 gsl_vector *ephotons,
			 gsl_vector *radon);



void
raft_projection_polydata(raft_scan_t *data, 
			 raft_proj_t *workspace,
			 raft_phantom_t *phantom,
			 gsl_vector *ephotons);



double 
raft_projection_pixel_intersection_length(raft_scan_t *data, 
					  int i, 
					  int k);



int 
raft_projection_pixel_intersection_find(raft_scan_t *data, 
					int m,
					int i, 
					int j);


void
raft_projection_workspace_alloc(raft_scan_t *data,
				raft_proj_t *workspace);


void
raft_projection_workspace_free(raft_proj_t *workspace);


void
raft_projection_workspace_get_exp_dbt(raft_proj_t *workspace,
				      gsl_vector **exp);


double
raft_projection_loglikelihood_pet(raft_scan_t *data,
				  raft_proj_t *workspace,
				  gsl_vector *p,
				  gsl_vector *act,
				  gsl_vector *att);


double
raft_projection_loglikelihood_spect(raft_scan_t *data,
				    raft_proj_t *workspace,
				    gsl_vector *p,
				    gsl_vector *act,
				    gsl_vector *att);


double 
raft_projection_loglikelihood_xfct(raft_scan_t *data,
				   raft_proj_t *workspace,
				   gsl_vector *p,
				   gsl_vector *act,
				   gsl_vector *attT,
				   gsl_vector *attF);


int
raft_projection_radon_ray(raft_scan_t *data,
			  gsl_vector *phantom,
			  int i,
			  double *rsum,
			  double *dotp,
			  double *snorm,
			  gsl_vector *row);


int
raft_projection_radon_ray_pet(raft_scan_t *data,
			      raft_proj_t *proj,
			      gsl_vector *phantom,
			      gsl_vector *att,
			      int i,
			      double *rsum,
			      double *dotp,
			      double *snorm,
			      gsl_vector *row);


int
raft_projection_radon_ray_spect(raft_scan_t *data,
				raft_proj_t *proj,
				gsl_vector *phantom,
				gsl_vector *att,
				int i,
				double *rsum,
				double *dotp,
				double *snorm,
				gsl_vector *row);


int
raft_projection_radon_ray_xfct(raft_scan_t *data,
			       raft_proj_t *proj,
			       gsl_vector *phantom,
			       gsl_vector *attT,
			       gsl_vector *attF,
			       int i,
			       double *rsum,
			       double *dotp,
			       double *snorm,
			       gsl_vector *row);

double
raft_projection_radon_ray_dbt(raft_scan_t *data,
			      gsl_vector *phantom,
			      int j,
			      int x,
			      int y);


void 
raft_projection_phantom_atray(raft_scan_t *data,
			      gsl_vector *f,
			      int j,
			      gsl_vector *rest,
			      int rsize);

#endif
