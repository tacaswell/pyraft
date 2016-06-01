#include "raft_scan.h"
#include "raft_phantom.h"
#include "raft_param.h"
#include "scan.h"

/*######################################################
  Title: Scan data

  Header - <raft/raft_scan.h>
  Type - <raft_scan_t>
  $Id: scan.c,v 1.21 2010-06-24 16:13:24 miqueles Exp $ - Last Update
  ####################################################*/

/*######################################################
  Public
  ####################################################*/

/*######################################################
  Section: Allocation
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_scan_alloc

  Allocattes scan structure. 

  Input: 

  Empty.

  Return:

  Scan structure. See <raft_scan_t>.
  
  +====================================================+
*/

raft_scan_t *raft_scan_alloc(void)
{
  return (raft_scan_t *)malloc(sizeof(raft_scan_t));
}

/*+====================================================+
  
  FUNCTION: raft_scan_alloc_data

  Allocattes scan structure data. 

  Input: 

  data - scan data . See <raft_scan_t>.
  nviews - number of views
  nrays - number of rays per view
  min - lower bound for (x,y)-axis
  max - upper bound for (x,y)-axis
  scan - scan type (RAFT_CT, RAFT_SPECT, RAFT_XFCT)
  
  Return:

  RAFT_SUCCESS - successfull operation
  RAFT_ENOMEM - not enough memory
  RAFT_EDOM - domain error  
  
  +====================================================+
*/

int raft_scan_alloc_data(raft_scan_t *data,			      
			 int nviews,  
			 int nrays,			 	 
			 int size,
			 double min,
			 double max,
			 int scan) 
{
  int ndata, npixels;

  if(nviews<=0 || nrays<=0 || size<=0 || min>=max)
    return RAFT_EDOM;

  gsl_set_error_handler_off();
  
  ndata = nviews*nrays;
  npixels = size*size;
  
  data->t = gsl_vector_alloc(nrays);
  if(!data->t)
    return RAFT_ENOMEM;
  
  data->theta = gsl_vector_alloc(nviews);
  if(!data->theta)
    return RAFT_ENOMEM;
  
  data->costheta = gsl_vector_alloc(nviews);
  if(!data->costheta)
    return RAFT_ENOMEM;

  data->sintheta = gsl_vector_alloc(nviews);
  if(!data->sintheta)
    return RAFT_ENOMEM;

  data->tantheta = gsl_vector_alloc(nviews);
  if(!data->tantheta)
    return RAFT_ENOMEM;
  
  data->x = gsl_vector_alloc(size);
  if(!data->x)
    return RAFT_ENOMEM;
  
  data->y = gsl_vector_alloc(size);
  if(!data->y)
    return RAFT_ENOMEM;
 
  data->q = gsl_vector_calloc(ndata);
  if(!data->q)
    return RAFT_ENOMEM;

  data->nrays  = nrays;
  data->nviews = nviews;
  data->min = min;
  data->max = max;
  data->size = size;
  data->pixels = npixels;
  data->etype = RAFT_MONO;
  data->scan = scan;
  
  /**/
  data->xfctAperture = 5;
  data->ktbound = RAFT_YKTBND;

  /**/
  
  data->memory.a = gsl_vector_alloc(npixels);
  gsl_vector_set_all(data->memory.a, 0);

  data->memory.ones              = gsl_vector_alloc(npixels);
  data->memory.Ave               = gsl_vector_alloc(npixels);
  data->memory.AveMinusWeight    = gsl_vector_alloc(npixels);
  data->memory.MaxAveMinusWeight = gsl_vector_alloc(nviews);
  
  return RAFT_SUCCESS;
}

/*+====================================================+
  
  FUNCTION: raft_scan_alloc_data_poly

  Allocattes scan structure polychromatic data. 

  Input: 

  nergs - number of energy levels
  nsubs - number of substances

  Return:

  RAFT_SUCCESS - successfull operation
  RAFT_ENOMEM - not enough memory
  RAFT_EDOM - domain error
  
  +====================================================+
*/

int raft_scan_alloc_data_poly(raft_scan_t *data,			      
			      int nergs,  
			      int nsubs) 
{
  int k, ndata, npixels;

  npixels = raft_scan_get_npixels(data);
  ndata = raft_scan_get_ndata(data);
    
  if(nergs<=0 || nsubs<=0)
    return RAFT_EDOM;

  gsl_set_error_handler_off();

  data->energy = (energy_t *)malloc(sizeof(energy_t));
  if(!data->energy)
    return RAFT_ENOMEM;

  data->energy->energies = gsl_vector_alloc(nergs);
  if(!data->energy->energies)
    return RAFT_ENOMEM;

  data->energy->spectrum = gsl_vector_alloc(nergs);
  if(!data->energy->spectrum)
    return RAFT_ENOMEM;

  data->energy->G = gsl_matrix_alloc(nergs,nsubs);
  if(!data->energy->G)
    return RAFT_ENOMEM;
  
  data->energy->U = gsl_matrix_alloc(nergs,nsubs);
  if(!data->energy->U)
    return RAFT_ENOMEM;
  
  data->energy->V = gsl_matrix_alloc(nsubs,nsubs);
  if(!data->energy->V)
    return RAFT_ENOMEM;
  
  data->energy->S = gsl_vector_alloc(nsubs);
  if(!data->energy->S)
    return RAFT_ENOMEM;

  data->energy->ones = gsl_vector_alloc(nergs);
  if(!data->energy->ones)
    return RAFT_ENOMEM;
  
  data->energy->basis = (decomp_t *)malloc(sizeof(decomp_t)*nsubs);
  
  for(k=0; k<nsubs; k++)
    {
      data->energy->basis[k].p = gsl_vector_alloc(ndata);
    }
  
  data->energy->nsubs = nsubs;
  data->energy->nergs = nergs;
  data->etype = RAFT_POLY;

  return RAFT_SUCCESS;
}


/*+====================================================+
  
  FUNCTION: raft_scan_free

  Frees scan structure. 

  Input: 

  data - scan data . See <raft_scan_t>.  

  +====================================================+
*/

void raft_scan_free(raft_scan_t *data)
{
  free(data);
}

/*+====================================================+
  
  FUNCTION: raft_scan_free_data

  Frees scan structure data. 

  Input: 

  data - scan data . See <raft_scan_t>.
  
  +====================================================+
*/

void raft_scan_free_data(raft_scan_t *data)
{
  int k;
  
  gsl_vector_free(data->t);
  gsl_vector_free(data->theta);
  gsl_vector_free(data->costheta);
  gsl_vector_free(data->sintheta);
  gsl_vector_free(data->tantheta);
  gsl_vector_free(data->x);
  gsl_vector_free(data->y);
  gsl_vector_free(data->q);
  
  if(data->etype == RAFT_POLY)
    {
      gsl_vector_free(data->energy->energies);
      gsl_vector_free(data->energy->spectrum);
      gsl_matrix_free(data->energy->G);
      gsl_matrix_free(data->energy->U);
      gsl_matrix_free(data->energy->V);
      gsl_vector_free(data->energy->S);
      gsl_vector_free(data->energy->ones);
      
      for(k=0; k< raft_scan_get_nsubs(data); k++)
	{
	  gsl_vector_free(data->energy->basis[k].p);
	}
      
      free(data->energy->basis);
      free(data->energy);
    }

  /**/

  gsl_vector_free(data->memory.ones);
  gsl_vector_free(data->memory.a);
  gsl_vector_free(data->memory.Ave);
  gsl_vector_free(data->memory.AveMinusWeight);
  gsl_vector_free(data->memory.MaxAveMinusWeight); 
}

/*######################################################
  Section: Parse scanning description file
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_scan_parse

  Set scan data from a description file. 

  Input:

  fp - file pointer

  Output: 

  data - scan data . See <raft_scan_t>.
  
  Returns:

  RAFT_EDOM - Domain error
  RAFT_SUCCESS - Successfull operation
  RAFT_CFG - Descriprition file error
  
  +====================================================+
*/

int raft_scan_parse(FILE *fp, 
		    raft_scan_t *data)
{
  double min, max;
  int status, nrays, nviews, nergs, nsubs, size, scantype; 
  
  cfg_t *scan, *section;

  cfg_opt_t geometry_opt[] =
    {  
      CFG_INT("size", 1, CFGF_NONE),
      CFG_INT("rays",1, CFGF_NONE),
      CFG_INT("views",1, CFGF_NONE),
      CFG_STR("type",0, CFGF_NODEFAULT),
      CFG_END()
    };

  cfg_opt_t efunction_opt[] = 
    {
      CFG_INT("id",1,CFGF_NODEFAULT),
      CFG_STR("fun",0,CFGF_NODEFAULT),
      CFG_END()
    };
  
  cfg_opt_t basis_opt[] = 
    {
      CFG_SEC("efunction", efunction_opt, CFGF_MULTI),
      CFG_END()      
    };

  cfg_opt_t energy_opt[] =
    {
      CFG_INT("energies",1,CFGF_NODEFAULT),
      CFG_INT("substances",1,CFGF_NODEFAULT),
      CFG_STR("levels", 0, CFGF_NODEFAULT),
      CFG_STR("spectrum", 0, CFGF_NODEFAULT),
      CFG_SEC("basis", basis_opt, CFGF_NODEFAULT),
      CFG_END()
    };

  cfg_opt_t meshdata_opt[] =
    {
      CFG_FLOAT_LIST("array", 0, CFGF_NODEFAULT),
      CFG_STR("file", 0, CFGF_NODEFAULT ),
      CFG_STR("id", 0, CFGF_NODEFAULT),
      CFG_END()
    };

  cfg_opt_t mesh_opt[] = 
    {      
      CFG_SEC("rays", meshdata_opt, CFGF_NONE),
      CFG_SEC("views", meshdata_opt, CFGF_NONE),      
      CFG_END()
    };
  
  cfg_opt_t statistic_opt[] = 
    {
      CFG_FLOAT("emitted", 1, CFGF_NODEFAULT),      
      CFG_FLOAT("transmittance", 1, CFGF_NODEFAULT),      
      CFG_FLOAT("efficiency", 1, CFGF_NODEFAULT),
      CFG_END()
    };

  cfg_opt_t scan_opt[] = 
    {
      CFG_SEC("geometry",geometry_opt, CFGF_NONE),
      CFG_SEC("energy",energy_opt, CFGF_NODEFAULT),
      CFG_SEC("statistic",statistic_opt, CFGF_NODEFAULT),
      CFG_SEC("mesh",mesh_opt, CFGF_NODEFAULT),
      CFG_END()
    };
  
  scan = cfg_init(scan_opt, CFGF_NOCASE);

  if(cfg_parse_fp(scan, fp)==CFG_PARSE_ERROR)
    return RAFT_EDOM;
  
  /*----------*/
  /* geometry */
  
  section = cfg_getsec(scan, "geometry");
  
  status = parse_geometry(section, &nrays, &nviews, 
			  &size, &min, &max, &scantype);
  
  if(status!=RAFT_SUCCESS)
    {
      cfg_free(scan);
      return RAFT_CFG;
    }
      
  raft_scan_alloc_data(data, nviews, nrays, size, min, max, scantype);
  
  /*--------*/  
  /* energy */

  if(cfg_size(scan,"energy")!=0)
    {
      section = cfg_getsec(scan, "energy");

      nergs = cfg_getint(section,"energies");
      nsubs = cfg_getint(section,"substances");
      
      raft_scan_alloc_data_poly(data, nergs, nsubs);
	
      status = parse_energy(section, data);
      if(status!=RAFT_SUCCESS)
	{
	  cfg_free(scan);
	  return RAFT_CFG;
	}
    }
  
  /*-----------*/
  /* statistic */
  
  if(cfg_size(scan,"statistic")!=0)
    {
      section = cfg_getsec(scan, "statistic");
      
      status = parse_statistic(section, data);
      if(status!=RAFT_SUCCESS)
	{
	  cfg_free(scan);
	  return RAFT_CFG;
	}
    }
  
  /*--------------*/
  /* views & rays */

  if(cfg_size(scan,"mesh")!=0)
    {
      section = cfg_getsec(scan, "mesh");
      
      status = parse_mesh(section, data);
      if(status!=RAFT_SUCCESS)
	{
	  cfg_free(scan);
	  return RAFT_CFG;
	}
    }
  else
    raft_scan_set_data(data);
  
  /*------------*/
  
  cfg_free(scan);
  return RAFT_SUCCESS;
}

/*######################################################
  Section: Setting scanning parameters
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_scan_set_xfctAperture

  Set aperture of XFCT angle section.
  
  The input 'apert' stands for an integer varying 
  from 0 to m/4, with 'm' the number of views. 
  Set apert=0 for an angle section [0,pi] while 
  apert=m/4 for only the perpendicular ray.

  Input: 

  data - scan data. See <raft_scan_t>.
  apert - integer aperture
  
  +====================================================+
*/

void raft_scan_set_xfctAperture(raft_scan_t *data, int apert)
{
  data->xfctAperture = apert;
}

/*+====================================================+
  
  FUNCTION: raft_scan_set_data

  Set usual scan data. 

  Input: 

  data - scan data . See <raft_scan_t>.
  
  +====================================================+
*/

void raft_scan_set_data(raft_scan_t *data)  
{
  int i, size, nrays, nviews;
  double dt, dtheta, tmax, min, max, z, th, cost, sint, tant;
  
  nrays = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);	
   
  min = raft_scan_get_min(data);
  max = raft_scan_get_max(data); 
  
  tmax = (max-min)/sqrt(2);
  dt = (2*tmax)/data->nrays;  
  for(i=0; i < raft_scan_get_nrays(data); i++)    
    gsl_vector_set(data->t, i, -tmax + i*dt);

  if(data->scan == RAFT_CT)
    dtheta = MPI/nviews;
  else
    dtheta = (2*MPI)/nviews;
  
  for(i=0; i < nviews; i++)
    {
      th   = i*dtheta;
      cost = cos(th);
      sint = sin(th);
      tant = sint/cost;

      gsl_vector_set(data->theta, i, i*dtheta);
      gsl_vector_set(data->costheta, i, cost);    
      gsl_vector_set(data->sintheta, i, sint);
      gsl_vector_set(data->tantheta, i, tant);
    }
  
  size = raft_scan_get_size(data);

  for(i=0; i < size-1; i++)
    {
      z =  min + i*(max - min)/(size-1);
      gsl_vector_set(data->x, i, z);
      gsl_vector_set(data->y, i, z);      
    }
  
  gsl_vector_set(data->x, size-1, max);
  gsl_vector_set(data->y, size-1, max);
  
  data->raysid = RAFT_REG;
  data->viewsid = RAFT_REG;

  gsl_vector_set_all(data->memory.ones, 1.0);
}

/*+====================================================+
  
  FUNCTION: raft_scan_set_photon_statistics

  Photon statistics in transmission tomography. 

  Input: 

  data - scan data . See <raft_scan_t>.
  emitted - number of photons emiited by each source
  transmittance - probability of photons scattering
  efficiency - probability of photons counting
  
  +====================================================+
*/

void raft_scan_set_photon_statistics(raft_scan_t *data,
				     double emitted,
				     double transmittance,
				     double efficiency)
{
  data->emitted = emitted;
  data->transmittance = transmittance;
  data->efficiency = efficiency;
}


/*+====================================================+
  
  FUNCTION: raft_scan_set_energies

  Define energy levels and spectrum. 

  Input: 

  data - scan data . See <raft_scan_t>.
  energies - array with energy levels
  spectrum - source spectrum
  
  +====================================================+
*/

void raft_scan_set_energies(raft_scan_t *data,
			    gsl_vector *energies, 
			    gsl_vector *spectrum)
{
  if(data->etype == RAFT_POLY)
    {
      gsl_blas_dcopy(energies, data->energy->energies);  
      gsl_blas_dcopy(spectrum, data->energy->spectrum);  
    }
}

/*+====================================================+
  
  FUNCTION: raft_scan_set_energy

  Define only an energy levels with spectrum. 

  Input: 

  data - scan data . See <raft_scan_t>.
  energy - energy
  spectrum - source spectrum
  n - nth energy level
  
  +====================================================+
*/

void raft_scan_set_energy(raft_scan_t *data,
			  double energy, 
			  double spectrum,
			  int n)
{
  if(data->etype == RAFT_POLY)
    {
      gsl_vector_set(data->energy->energies, n, energy);  
      gsl_vector_set(data->energy->spectrum, n, spectrum);  
    }
}


/*+====================================================+
  
  FUNCTION: raft_scan_set_energy_basis

  Define energy basis function for a given substance. 

  Input: 

  data - scan data . See <raft_scan_t>.
  g - function value
  i - energy index
  k - substance index
  
  +====================================================+
*/

void raft_scan_set_energy_basis(raft_scan_t *data,
				double g,
				int i,
				int k)
{
  if(data->etype == RAFT_POLY)
    {
      gsl_matrix_set(data->energy->G, i, k, g);      
    }
}


/*+====================================================+
  
  FUNCTION: raft_scan_set_projection

  Set projection within scan data. 

  Input: 

  data - scan data . See <raft_scan_t>.
  p - projection vector
  i - row index (for views)
  j - column index (for rays)
  r - matrix entry  

  +====================================================+
*/

void raft_scan_set_projection(raft_scan_t *data, 
			      gsl_vector *p,
			      int i, 
			      int j, 
			      double r)
{
  gsl_vector_set(p, j + i*raft_scan_get_nrays(data), r);  
}  


/*######################################################
  Section: Getting scanning parameters
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_scan_get_xfctAperture

  Get the aperture of an XFCT angle section. 

  Input: 

  data - scan data . See <raft_scan_t>.
    
  Return:

  Integer aperture. See function 
  <raft_scan_set_xfctAperture>.
  
  +====================================================+
*/

int raft_scan_get_xfctAperture(raft_scan_t *data)
{ 
  return data->xfctAperture;
}

/*+====================================================+
  
  FUNCTION: raft_scan_get_projection

  Get projection within scan data. 

  Input: 

  data - scan data . See <raft_scan_t>.
  p - projection vector
  i - row index (view)
  j - column index (ray)
  
  Return:

  Matrix entry.
  
  +====================================================+
*/

double raft_scan_get_projection(raft_scan_t *data, 
				gsl_vector *p,
				int i, 
				int j)
{
  return gsl_vector_get(p, j + i*raft_scan_get_nrays(data));
}  


/*+====================================================+
  
  FUNCTION: raft_scan_get_photons_emitted

  Get the number of photons emitted by a source. 

  Input: 

  data - scan data . See <raft_scan_t>.
   
  Return:

  Number of emitted photons.
  
  +====================================================+
*/

double raft_scan_get_photons_emitted(raft_scan_t *data)
{
  return data->emitted;
}


/*+====================================================+
  
  FUNCTION: raft_scan_get_photons_transmittance

  Get the transmittance of a detector. 

  Input: 

  data - scan data . See <raft_scan_t>.
   
  Return:

  Probability of photons not be absorbed/scatter.
  
  +====================================================+
*/

double raft_scan_get_photons_transmittance(raft_scan_t *data)
{
  return data->transmittance;
}


/*+====================================================+
  
  FUNCTION: raft_scan_get_photons_efficiency

  Get the efficiency of a detector. 

  Input: 

  data - scan data . See <raft_scan_t>.
   
  Return:

  Probability of photon counting.
  
  +====================================================+
*/

double raft_scan_get_photons_efficiency(raft_scan_t *data)
{
  return data->efficiency;
}


/*+====================================================+
  
  FUNCTION: raft_scan_get_nenergies

  Get the number of energy levels in scan structure. 

  Input: 

  data - scan data . See <raft_scan_t>.
  
  Return:

  Number of energy levels whenever the polychromatic
  data is defined, otherwise zero is returned.
  
  +====================================================+
*/

int raft_scan_get_nenergies(raft_scan_t *data)
{
  if(data->etype == RAFT_POLY)
    return data->energy->nergs;
  else
    return 0;
}


/*+====================================================+
  
  FUNCTION: raft_scan_get_nsubs

  Get the number of substances in scan structure. 

  Input: 

  data - scan data . See <raft_scan_t>.
  
  Return:

  Number of substances whenever the polychromatic data
  is defined, otherwise zero is returned.
  
  +====================================================+
*/

int raft_scan_get_nsubs(raft_scan_t *data)
{
  if(data->etype == RAFT_POLY)
    return data->energy->nsubs;
  return 0;
}

/*+====================================================+
  
  FUNCTION: raft_scan_get_nrays

  Get the number of rays per view in scan structure. 

  Input: 

  data - scan data . See <raft_scan_t>.
  
  Return:

  Number of rays per view.

  +====================================================+
*/

int raft_scan_get_nrays(raft_scan_t *data)
{
  return data->nrays;
}

/*+====================================================+
  
  FUNCTION: raft_scan_get_nviews

  Get the number of views in scan structure. 

  Input: 

  data - scan data . See <raft_scan_t>.
  
  Return:

  Number of views.
  
  +====================================================+
*/

int raft_scan_get_nviews(raft_scan_t *data)
{
  return data->nviews;
}


/*+====================================================+
  
  FUNCTION: raft_scan_get_size

  Get the number of pixels in x,y axis. 

  Input: 

  data - scan data . See <raft_scan_t>.
  
  Return:

  Number of pixels.
  
  +====================================================+
*/

int raft_scan_get_size(raft_scan_t *data)
{
  return data->size;
}


/*+====================================================+
  
  FUNCTION: raft_scan_get_npixels

  Get the number of pixels in scan structure. 

  Input: 

  data - scan data . See <raft_scan_t>.
  
  +====================================================+
*/

int raft_scan_get_npixels(raft_scan_t *data)
{
  return data->pixels;
}

/*+====================================================+
  
  FUNCTION: raft_scan_get_ndata

  Get the number of data in scan structure. 

  Input: 

  data - scan data . See <raft_scan_t>.
  
  Return:

  Number of data, i.e, the number of views times
  the number of rays (per view).  

  +====================================================+
*/

int raft_scan_get_ndata(raft_scan_t *data)
{
  return (data->nrays) * (data->nviews);
}


/*+====================================================+
  
  FUNCTION: raft_scan_get_min

  Get the axis minimum coordinate. 

  Input: 

  data - scan data . See <raft_scan_t>.
  
  +====================================================+
*/

double raft_scan_get_min(raft_scan_t *data)
{
  return data->min;
}


/*+====================================================+
  
  FUNCTION: raft_scan_get_max

  Get the axis maximum coordinate. 

  Input: 

  data - scan data . See <raft_scan_t>.  
  
  +====================================================+
*/

double raft_scan_get_max(raft_scan_t *data)
{
  return data->max;
}

/*+====================================================+
  
  FUNCTION: raft_scan_get_angle

  Get the jth scan angle. 

  Input: 

  data - scan data . See <raft_scan_t>.
  j - index  

  +====================================================+
*/

double raft_scan_get_angle(raft_scan_t *data,
			   int j)
{
  return gsl_vector_get(data->theta, j);
}

/*+====================================================+
  
  FUNCTION: raft_scan_get_ray

  Get the jth coordinate of the ray axis. 

  Input: 

  data - scan data . See <raft_scan_t>.
  j - index
  
  +====================================================+
*/

double raft_scan_get_ray(raft_scan_t *data,
			 int j)
{
  return gsl_vector_get(data->t, j);
}

/*+====================================================+
  
  FUNCTION: raft_scan_get_x

  Get the ith element of the x-axis. 
   
  Input: 

  data - scan data . See <raft_scan_t>.
  j - index  
  
  +====================================================+
*/

double raft_scan_get_x(raft_scan_t *data,
		       int j)
{
  return gsl_vector_get(data->x, j);
}

/*+====================================================+
  
  FUNCTION: raft_scan_get_y

  Get the ith element of the y-axis. 
   
  Input: 

  data - scan data . See <raft_scan_t>.
  j - index  

  +====================================================+
*/

double raft_scan_get_y(raft_scan_t *data,
		       int j)
{
  return gsl_vector_get(data->y, j);
}


/*+====================================================+
  
  FUNCTION: raft_scan_get_raystep

  Get the ray step distance. 

  Input: 

  data - scan data . See <raft_scan_t>.    
  
  +====================================================+
*/

double raft_scan_get_raystep(raft_scan_t *data)
{
  return raft_scan_get_ray(data,1)-raft_scan_get_ray(data,0);
}

/*######################################################
  Private
  ####################################################*/

/*+====================================================+
  
  FUNCTION  parse_geometry

  Parse scan geometry. 

  Input  

  section - section from a description file
  
  Output 

  nrays - number of rays
  nviews - number of views
  size - phantom size
  min - lower bound for y,x-coordinate 
  max - upper bound for y,x-coordinate
  
  Return 

  RAFT_CFG - Description file error
  RAFT_SUCCESS - Successfull operation
  
  +====================================================+
*/

int parse_geometry(cfg_t *section,
		   int *nrays,
		   int *nviews,
		   int *size,
		   double *min,
		   double *max,
		   int *scantype)
{
  if(!cfg_size(section,"size") ||
     !cfg_size(section,"rays") ||
     !cfg_size(section,"views") || 
     !cfg_size(section,"type"))
    {
      return RAFT_CFG;
    }
  else
    {
      int defined;

      *min = XYMIN; /*cfg_getfloat(section,"min");*/
      *max = XYMAX; /*cfg_getfloat(section,"max");*/

      *nrays = cfg_getint(section,"rays");
      *nviews = cfg_getint(section,"views");
      *size = cfg_getint(section,"size");

      defined = 0;
      
      if(!defined && !strcmp("RAFT_CT",cfg_getstr(section,"type")))
	{
	  defined = 1;
	  *scantype = RAFT_CT;
	}
      
      if(!defined && !strcmp("RAFT_SPECT",cfg_getstr(section,"type")))
	{
	  defined = 1;
	  *scantype = RAFT_SPECT;
	}
      
      if(!defined && !strcmp("RAFT_XFCT",cfg_getstr(section,"type")))
	{
	  defined = 1;
	  *scantype = RAFT_XFCT;
	}
      
      if(!defined)
	return RAFT_CFG;
      else
	return RAFT_SUCCESS;
    }
}


/*+====================================================+
  
  FUNCTION  parse_energy

  Parse scan energies. 

  Input  

  section - section from a description file
  
  Output 

  data - scan data
  
  Return 

  RAFT_CFG - Description file error
  RAFT_SUCCESS - Successfull operation
  
  +====================================================+
*/

int parse_energy(cfg_t *section,
		 raft_scan_t *data)
{
  FILE *fpl, *fps, *fp;
  int i,k, sizeL, sizeS, nergs, nsubs, nbasis, ret;
  cfg_t *basis, *function;

  nergs = data->energy->nergs;
  nsubs = data->energy->nsubs;
    
  sizeL = cfg_size(section, "levels");
  sizeS = cfg_size(section, "spectrum");
  
  if(!sizeL || !sizeS)
    return RAFT_CFG;
  
  fpl = fopen(cfg_getstr(section,"levels"),"r");
  fps = fopen(cfg_getstr(section,"spectrum"),"r");
  if(fpl==NULL || fps==NULL)
    return RAFT_CFG;
  else
    {
      for(k = 0; k < nergs; k++)
	{
	  double e, s;
	  
	  ret = fscanf(fpl,"%lf",&e);
	  gsl_vector_set(data->energy->energies, k, e);
	  
	  ret = fscanf(fps,"%lf",&s);
	  gsl_vector_set(data->energy->spectrum, k, s);
	}
    }
  fclose(fpl);
  fclose(fps);

  basis = cfg_getsec(section, "basis");
  nbasis = cfg_size(basis, "efunction");
  
  if(nbasis!=nsubs)
    return RAFT_CFG;
  
  for(k = 0; k < nsubs; k++)
    {
      int id;
      double g;
      
      function = cfg_getnsec(basis,"efunction",k);

      id = cfg_getint(function,"id");
      
      fp = fopen(cfg_getstr(function,"fun"),"r");
      if(fp==NULL)
	return RAFT_CFG;
      else
	{
	  for(i=0; i<nergs; i++)
	    {
	      ret = fscanf(fp,"%lf",&g);
	      gsl_matrix_set(data->energy->G, i, id-1, g);
	    }
	}
      fclose(fp);
    }
  
  data->energy->decomposed = 0;
  
  return RAFT_SUCCESS;
}

/*+====================================================+
  
  FUNCTION  parse_statistic

  Parse scan statistics. 

  Input  

  section - section from a description file
  
  Output 

  data - scan data
  
  Return 

  RAFT_CFG - Description file error
  RAFT_SUCCESS - Successfull operation

  +====================================================+
*/

int parse_statistic(cfg_t *section,
		    raft_scan_t *data)
{
  if(!cfg_size(section,"emitted") || 
     !cfg_size(section,"transmittance") ||
     !cfg_size(section,"efficiency"))
    {
      return RAFT_CFG;
    }
  else
    {
      data->emitted = cfg_getfloat(section,"emitted");
      data->transmittance = cfg_getfloat(section,"transmittance");
      data->efficiency = cfg_getfloat(section,"efficiency");

      return RAFT_SUCCESS;
    }
}

/*+====================================================+
  
  FUNCTION  parse_mesh

  Parse scan views and rays. 

  Input  

  section - section from a description file
  
  Output 

  data - scan data
  
  Return 

  RAFT_SUCCESS - Successfull operation
  RAFT_CFG - Description file error
  
  +====================================================+
*/

int parse_mesh(cfg_t *section,
	       raft_scan_t *data)
{
  int filer, filev, defined, mesh;
  int k, nrays, nviews, size, ret;
  double t, theta, cost, sint, tant, min, max, z;
  cfg_t *subsection;
  
  size = data->size;
  min = data->min;
  max = data->max;
  
  /*------*/
  /* rays */
  
  subsection = cfg_getsec(section,"rays");
 
  nrays = cfg_size(subsection,"array");
  filer = cfg_size(subsection,"file");

  if(filer!=0 || nrays!=0)
    {
      if(!cfg_size(subsection,"id"))
	return RAFT_CFG;
      else
	{
	  defined = 0;
	  
	  if(!defined && !strcmp("RAFT_REG",cfg_getstr(subsection,"id")))
	    {
	      defined = 1;
	      mesh = RAFT_REG;
	    }
	  
	  if(!defined && !strcmp("RAFT_NREG",cfg_getstr(subsection,"id")))
	    {
	      defined = 1;
	      mesh = RAFT_NREG;
	    }
	  
	  if(!defined && !strcmp("RAFT_CHEB",cfg_getstr(subsection,"id")))
	    {
	      defined = 1;
	      mesh = RAFT_CHEB;
	    }
	  
	  if(!defined)
	    return RAFT_CFG;
	}
      
      data->raysid = mesh;
    }

  if(nrays!=0)
    {
      if(nrays > data->nrays)
	{
	  gsl_vector_free(data->t);
	  data->t = gsl_vector_alloc(nrays);
	}
      
      data->nrays = nrays;
      
      for(k = 0; k < nrays; k++)
	{
	  t = cfg_getnfloat(subsection, "array", k);
	  gsl_vector_set(data->t, k, t);
	}
    }
  else
    {
      if(filer!=0)
	{
	  FILE *fp;
	  
	  fp = fopen(cfg_getstr(subsection,"file"),"r");
	  if(fp==NULL)
	    return RAFT_CFG;
	  else
	    {
	      for(k=0; k < raft_scan_get_nrays(data); k++)    
		{
		  ret = fscanf(fp,"%lf",&t);
		  gsl_vector_set(data->t, k, t);
		}
	    }
	  fclose(fp);
	}
      else
	{
	  double dt, tmax;
	  
	  tmax = (max-min)/sqrt(2);
	  dt = 2*tmax/data->nrays;  
	  
	  for(k=0; k < raft_scan_get_nrays(data); k++)    
	    {
	      gsl_vector_set(data->t, k, -tmax + k*dt);
	    }
	  
	  data->raysid = RAFT_REG;
	}
    }

  /*-------*/
  /* views */

  subsection = cfg_getsec(section,"views");
  
  nviews = cfg_size(subsection,"array");
  filev = cfg_size(subsection,"file");
  
  if(filev!=0 || nviews!=0)
    {
      if(!cfg_size(subsection,"id"))
	return RAFT_CFG;
      else
	{
	  defined = 0;
	  
	  if(!defined && !strcmp("RAFT_REG",cfg_getstr(subsection,"id")))
	    {
	      defined = 1;
	      mesh = RAFT_REG;
	    }
	  
	  if(!defined && !strcmp("RAFT_NREG",cfg_getstr(subsection,"id")))
	    {
	      defined = 1;
	      mesh = RAFT_NREG;
	    }
	  
	  if(!defined && !strcmp("RAFT_CHEB",cfg_getstr(subsection,"id")))
	    {
	      defined = 1;
	      mesh = RAFT_CHEB;
	    }
	  
	  if(!defined)
	    return RAFT_CFG;
	}

      data->viewsid = mesh;
    }
  
  if(nviews!=0)
    {
      if(nviews > data->nviews)
	{
	  gsl_vector_free(data->theta);
	  data->theta = gsl_vector_alloc(nviews);
	}
      
      data->nviews = nviews;
      
      for(k = 0; k < nviews; k++)
	{
	  theta = cfg_getnfloat(subsection, "array", k);
	  gsl_vector_set(data->theta, k, theta);
	}
    }
  else
    {
      if(filev!=0)
	{
	  FILE *fp;
	  
	  fp = fopen(cfg_getstr(subsection,"file"),"r");
	  if(fp==NULL)
	    return RAFT_CFG;
	  else
	    {
	      for(k=0; k < raft_scan_get_nviews(data); k++)    
		{
		  ret = fscanf(fp,"%lf",&theta);
		  
		  cost = cos(theta);
		  sint = sin(theta);
		  tant = sint/cost;

		  gsl_vector_set(data->theta, k, theta);
		  gsl_vector_set(data->costheta, k, cost);
		  gsl_vector_set(data->sintheta, k, sint);
		  gsl_vector_set(data->tantheta, k, tant);
		}
	    }
	  fclose(fp);
	}
      else
	{
	  double dtheta;

	  if(data->scan == RAFT_CT)
	    dtheta = MPI/data->nviews;
	  else
	    dtheta = (2*MPI)/data->nviews;
	  
	  for(k=0; k < data->nviews; k++)
	    {
	      theta = k*dtheta;
	      cost  = cos(theta);
	      sint  = sin(theta);
	      tant  = sint/cost;

	      gsl_vector_set(data->theta, k, k*dtheta);
	      gsl_vector_set(data->costheta, k, cost);
	      gsl_vector_set(data->sintheta, k, sint);
	      gsl_vector_set(data->tantheta, k, tant);
	    }

	  data->viewsid = RAFT_REG;
	}
    }

  /*----------*/
  /* x,y mesh */
  
  for(k=0; k < size-1; k++)
    {
      z =  min + k*(max - min)/(size-1);
      gsl_vector_set(data->x, k, z);
      gsl_vector_set(data->y, k, z);      
    }
  
  gsl_vector_set(data->x, size-1, max);
  gsl_vector_set(data->y, size-1, max);  

  return RAFT_SUCCESS;
}
