#ifndef _RAFT_SCAN_H_
#define _RAFT_SCAN_H_

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include "raft_errno.h"
#include "raft_param.h"

/*######################################################
  Title: Scan data
  ####################################################*/

/*+====================================================+
  
  TYPEDEF energy_t

  Purpose

  Energy structure data. 

  Variables

  nergs - number of energy levels
  nsubs - number of substances
  energies - energy level
  spectrum - source spectrum
  G - energy matrix (nergs x nsubs)
  {U,V,S} - SVD decomposition of matrix G
  basis - basis attenuation functions
  ones - constant vector
  decomposed - 1 if matrix was already decomposed and 0
               othwerwise
          
  +====================================================+
*/

typedef struct{

  gsl_vector *p;

}decomp_t;
  

typedef struct{

  int nsubs;
  int nergs;
  int decomposed;
  gsl_vector *energies, *spectrum, *S, *ones;
  gsl_matrix *G, *U, *V; 
  decomp_t *basis;       
  
}energy_t;


/*+====================================================+
  
  TYPEDEF: raft_scan_t

  Purpose:

  RAFT scanning structure data.   
  +====================================================+
*/

/*
  Variables

  nviews - number of views in the tomgraphy scan
  nrays  - number of rays per view
  pixels - number of pixels
  theta - grid of points for angles
  t - grid of points for the ray domain
  size - number of rows/columns of the phantom matrix
  min - lower bound for the x,y axis
  max - upper bound for the x,y axis
  x - grid of points for the x axis
  y - grid of points for the y axis   
  q - temporary projection vector
  emitted - number of emitted photons by a source
  transmittance - probability of photon scattering
  efficiency - probability of photons counting
  raysid - identify regular mesh 
  viewsid - identify regular mesh
  energy - energy data
  etype - identifies monochromatic or polychromatic scan
  xfctAperture - aperture of XFCT angle section (integer from 0 to m/4 with 
                 m the number of views)
  ktbound - computation of bounds for Kunyansky's method
*/

typedef struct{

  int nviews;
  int nrays;
  int size;
  int pixels;
  int raysid, viewsid;
  int etype;
  int scan;
  int xfctAperture;
  int ktbound;
  double min, max; 
  double emitted, transmittance, efficiency;
  
  gsl_vector *theta, *t, *x, *y, *q;
  gsl_vector *costheta, *sintheta, *tantheta;
  
  energy_t *energy;    

  struct{
    gsl_vector *a;
    double norma;
    double dotp;    
    gsl_vector *Ave, *AveMinusWeight, *MaxAveMinusWeight, *ones;
    double c, ca, ma;
  }memory;
      
}raft_scan_t;


raft_scan_t *
raft_scan_alloc(void);


int 
raft_scan_alloc_data(raft_scan_t *data,			      
		     int nviews,  
		     int nrays,
		     int size,
		     double min,
		     double max,
		     int scan);


int 
raft_scan_alloc_data_poly(raft_scan_t *data,			      
			  int nergs,  
			  int nsubs);


void 
raft_scan_free(raft_scan_t *data);


void 
raft_scan_free_data(raft_scan_t *data);


int 
raft_scan_parse(FILE *fp, 
		raft_scan_t *data);


void 
raft_scan_set_data(raft_scan_t *data);



void 
raft_scan_set_projection(raft_scan_t *data, 
			 gsl_vector *p,
			 int i, 
			 int j, 
			 double r);  


void 
raft_scan_set_xfctAperture(raft_scan_t *data, 
			   int apert);


void
raft_scan_set_energies(raft_scan_t *data,
		       gsl_vector *energies, 
		       gsl_vector *spectrum);



void 
raft_scan_set_energy(raft_scan_t *data,
		     double energy, 
		     double spectrum,
		     int n);


void 
raft_scan_set_energy_basis(raft_scan_t *data,
			   double g,
			   int i,
			   int k);


void
raft_scan_set_photon_statistics(raft_scan_t *data,
				double emitted,
				double transmittance,
				double efficiency);


int 
raft_scan_get_xfctAperture(raft_scan_t *data);


double
raft_scan_get_projection(raft_scan_t *data, 
			 gsl_vector *p,
			 int i, 
			 int j);


double 
raft_scan_get_photons_transmittance(raft_scan_t *data);


double 
raft_scan_get_photons_efficiency(raft_scan_t *data);


double 
raft_scan_get_photons_emitted(raft_scan_t *data);


int 
raft_scan_get_nenergies(raft_scan_t *data);


int 
raft_scan_get_nsubs(raft_scan_t *data);

 
int 
raft_scan_get_nrays(raft_scan_t *data);


int 
raft_scan_get_nviews(raft_scan_t *data);


int 
raft_scan_get_npixels(raft_scan_t *data);


int 
raft_scan_get_ndata(raft_scan_t *data);


int 
raft_scan_get_size(raft_scan_t *data);


double 
raft_scan_get_min(raft_scan_t *data);


double
raft_scan_get_max(raft_scan_t *data);


double
raft_scan_get_angle(raft_scan_t *data,
		    int j);


double 
raft_scan_get_ray(raft_scan_t *data,
		  int j);


double 
raft_scan_get_x(raft_scan_t *data,
		int j);


double
raft_scan_get_y(raft_scan_t *data,
		int j);


double 
raft_scan_get_raystep(raft_scan_t *data);


#endif
