#ifndef _PROJECTION_H_
#define _PROJECTION_H_

#include "interp.h"
#include "raft_param.h"
#include "raft_phantom.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

typedef struct{

  raft_scan_t *data;
  raft_proj_t *workspace;
  gsl_vector *phantom;
  gsl_vector *radon;
  gsl_vector **weight;
  int size, nviews, nrays;
  int nthread;
  int colIndex[2];  
  
}parRad_t;


typedef struct{

  raft_scan_t *data;
  raft_proj_t *workspace;
  gsl_vector *phantom;
  gsl_vector **dbt;
  int size, nviews, nrays;
  int nthread;
  int colIndex[2];  
  
}parDbt_t;


void 
*radon_loop(void *t);

void 
*radonXfct_loop(void *t);

void 
*radonDbt_loop(void *t);


int 
intersection_ray_pixel(raft_scan_t *data,
		       int i, 
		       int j);


double 
eval_radon_charf_pixel(raft_scan_t *data,
		       int i,
		       int pixel);


double 
eval_rayintegral(raft_scan_t *data,
		 gsl_vector *phantom,
		 int j,
		 int k,
		 gsl_vector *a,
		 double *norma,
		 double *dotp);


double
eval_rayintegral_ert(raft_scan_t *data,
		     double *att,
		     gsl_vector *act,
		     int j,
		     int k,
		     gsl_vector *a,
		     double *norma,
		     double *dotp);

double 
eval_rayintegral_pet(raft_scan_t *data,
		     raft_proj_t *proj,
		     gsl_vector *phantom,
		     int j,
		     int k,
		     gsl_vector *a,
		     double *norma,
		     double *dotp);


double
eval_rayintegral_spect(raft_scan_t *data,
		       raft_proj_t *proj,
		       gsl_vector *phantom,
		       int j,
		       int k,
		       gsl_vector *a,
		       double *norma,
		       double *dotp);

double
eval_rayintegral_xfct(raft_scan_t *data,
		      raft_proj_t *proj,
		      gsl_vector *phantom,
		      int j,
		      int k,
		      gsl_vector *a,
		      double *norma,
		      double *dotp);


double 
eval_rayintegral_xfct_bythread(int nthread,
			       raft_scan_t *data,
			       raft_proj_t *proj,
			       gsl_vector *phantom,
			       int j,
			       int k,
			       gsl_vector *a,
			       double *norma,
			       double *dotp);
  

double
eval_rayintegral_generalized(raft_scan_t *data,
			     raft_proj_t *proj,
			     gsl_vector **weight,
			     gsl_vector *phantom,
			     int j,
			     int k,
			     gsl_vector *a,
			     double *norma,
			     double *dotp);


double
eval_rayintegral_dbt(raft_scan_t *data, 
		     gsl_vector *att,
		     int j,
		     int xi,
		     int yi);


double 
eval_ray_projection(raft_scan_t *data,
		    gsl_vector *phantom,
		    int i,
		    gsl_vector *a,
		    double *norma,
		    double *dotp);


double 
eval_ray_projection_ert(raft_scan_t *data,
			double *att,
			gsl_vector *phantom,
			int i,
			gsl_vector *a,
			double *norma,
			double *dotp);


double
eval_ray_projection_pet(raft_scan_t *data,
			raft_proj_t *workspace,
			gsl_vector *act,
			gsl_vector *att,
			int i,
			gsl_vector *a,
			double *norma,
			double *dotp);

double
eval_ray_projection_spect(raft_scan_t *data,
			  raft_proj_t *workspace,
			  gsl_vector *act,
			  gsl_vector *att,
			  int i,
			  gsl_vector *a,
			  double *norma,
			  double *dotp);

double 
eval_ray_projection_xfct(raft_scan_t *data,
			 raft_proj_t *workspace,
			 gsl_vector *act,
			 gsl_vector *attT,
			 gsl_vector *attF,
			 int i,
			 gsl_vector *a,
			 double *norma,
			 double *dotp);


void 
projection_workspace_memcpy(raft_proj_t *dworkspace,
			    raft_proj_t *sworkspace);


#endif

