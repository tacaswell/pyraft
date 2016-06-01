#include "raft_likelihood.h"
#include "raft_projection.h"
#include "raft_weight.h"
#include "projection.h"
#include "likelihood.h"
#include <gsl/gsl_math.h>

/*######################################################
  Title: Likelihood functions
  
  Header - <raft/raft_likelihood.h>
  Type - None
  $Id: likelihood.c,v 1.1 2008-09-29 16:12:09 miqueles Exp $ - Last Update
  ####################################################*/

/*+====================================================+
  
  FUNCTION  eval_loglikelihood_transmission

  Evaluate the transmission loglikelihood function. 
  
  Input  

  data - scan data . See <raft_scan_t>.
  workspace - workspace for iterative methods . See <raft_iterative_t>.
  f - given estimation
  ph - photon counts
   
  Return 

  Function value.
  
  +====================================================+
*/

double eval_loglikelihood_transmission(raft_scan_t *data,
				       raft_iterative_t *workspace,
				       gsl_vector *f,
				       gsl_vector *ph)
{
  double l,L;

  raft_projection_monodata(data, f, 
			   workspace->em.ephotons, 
			   workspace->em.radon);
  
  gsl_blas_ddot(workspace->em.ephotons, 
		workspace->em.ones, 
		&l);
  
  gsl_vector_mul(workspace->em.radon, 
		 ph);
  
  gsl_blas_ddot(workspace->em.radon, 
		workspace->em.ones, &L);
  
  L += l;
  
  return -L;
}

/*+====================================================+
  
  FUNCTION eval_loglikelihood_emission
  
  Evaluate emission-loglikelihood function. 
  
  Input
  
  proj - projection workspace. See <raft_proj_t>.
  radon - radon transform
  logradon - logarithm of the radon transform
  p - given data
  
  Return
  
  Likelihood value

  +====================================================+
*/

double 
eval_loglikelihood_emission(raft_proj_t *proj,
			    gsl_vector *radon,
			    gsl_vector *logradon,
			    gsl_vector *p)
{
  double L[2];
  
  gsl_blas_ddot(p, logradon, &L[0]);
  
  gsl_blas_ddot(radon, proj->likelihood.ones, &L[1]);
  
  return (L[0]-L[1]);
}

/*######################################################
  Section: Computing Likelihood functions 
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_loglikelihood_pet

  Evaluate the loglikelihood function for PET, 
  with known attenuation maps. 
  
  Input: 

  data - scan data . See <raft_scan_t>.
  workspace - projection workspace. See <raft_proj_t>.
  p - PET data
  act - activity map
  att - transmission attenuation map
   
  Return:

  Function value.
  
  +====================================================+
*/

double
raft_loglikelihood_pet(raft_scan_t *data,
		       raft_proj_t *workspace,
		       gsl_vector *p,
		       gsl_vector *act,
		       gsl_vector *att)
{
  double L;
  
  raft_projection_radon_pet(data, 
			    workspace,
			    act,
			    att,
			    workspace->likelihood.radon);

  L = eval_loglikelihood_emission(workspace,
				  workspace->likelihood.radon,
				  workspace->likelihood.lograd,
				  p);

  return L;
}


/*+====================================================+
  
  FUNCTION: raft_loglikelihood_spect
  
  Evaluate the loglikelihood function for SPECT, 
  with known attenuation maps. 
  
  Input: 

  data - scan data . See <raft_scan_t>.
  workspace - projection workspace. See <raft_proj_t>.
  p - SPECT data
  act - activity map
  att - attenuation map
    
  Return:

  Function value.
  
  +====================================================+
*/

double 
raft_loglikelihood_spect(raft_scan_t *data,
			 raft_proj_t *workspace,
			 gsl_vector *p,
			 gsl_vector *act,
			 gsl_vector *att)
{
  double L;
  
  raft_projection_radon_spect(data, 
			      workspace,
			      act,
			      att,
			      workspace->likelihood.radon);
  
  L = eval_loglikelihood_emission(workspace,
				  workspace->likelihood.radon,
				  workspace->likelihood.lograd,
				  p);

  return L;
}


/*+====================================================+
  
  FUNCTION: raft_loglikelihood_xfct
  
  Evaluate the loglikelihood function for XFCT, 
  with known attenuation maps. 
  
  Input: 

  data - scan data . See <raft_scan_t>.
  workspace - projection workspace. See <raft_proj_t>
  p - XFCT data
  act -  activity map
  attT - transmission attenuation map
  attF - XFCT attenuation map
  
  Return:

  Function value.
  
  +====================================================+
*/

double 
raft_loglikelihood_xfct(raft_scan_t *data,
			raft_proj_t *workspace,
			gsl_vector *p,
			gsl_vector *act,
			gsl_vector *attT,
			gsl_vector *attF)
{
  double L;
  
  raft_projection_radon_xfct(data, 
			     workspace,
			     act,
			     attF,
			     attT,
			     workspace->likelihood.radon);

  L = eval_loglikelihood_emission(workspace,
				  workspace->likelihood.radon,
				  workspace->likelihood.lograd,
				  p);

  return L;
}


