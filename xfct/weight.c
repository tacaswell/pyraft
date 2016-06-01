#include "raft_scan.h"
#include "raft_projection.h" 
#include "raft_cbf.h"
#include "raft_weight.h"
#include "projection.h"

/*######################################################
  Title: Weights
  
  Header - <raft/raft_weight.h>
  Type - None
  $Id: weight.c,v 1.28 2011-03-02 19:23:19 miqueles Exp $ - Last Update
  ####################################################*/

/*######################################################
  Section: Computing weight functions
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_weight_pet
  
  Compute the weight function for the PET radon transform.
  
  Input: 
  
  data - scan data. See <raft_scan_t>.
  att - attenuation map
  
  Output:
  
  proj - projection workspace. Weigth function is 
         computed within this workspace.  
  
  +====================================================+
*/

void raft_weight_pet(raft_scan_t *data,
	      	     raft_proj_t *proj,
		     gsl_vector *att)
{
  
  raft_projection_radon(data, proj, att, proj->pet.Radon);
  
  gsl_blas_dcopy(proj->pet.ExpRadon, proj->likelihood.exprad);
  
  gsl_vector_set_all(proj->pet.weight, 1.0);
      
  gsl_vector_div(proj->pet.weight, proj->pet.ExpRadon);
     
  proj->pet.defined = 1;      
}

/*+====================================================+
  
  FUNCTION: raft_weight_spect
  
  Compute the weight function for the SPECT radon transform.
  
  Input: 
  
  data - scan data. See <raft_scan_t>.
  att - attenuation map
  
  Output:
  
  proj - projection workspace. Weigth function is 
         computed within this workspace.
	      
  +====================================================+
*/

void raft_weight_spect(raft_scan_t *data,
		       raft_proj_t *proj,
		       gsl_vector *att)
{
  int j, nviews, npixels;
  
  nviews  = raft_scan_get_nviews(data);
  npixels = raft_scan_get_npixels(data);

  for(j = 0; j < nviews; j++)
    {      
      raft_projection_radon_dbt_view(data, proj, j, att, proj->dbt.Dbt[j]);
      
      gsl_vector_set_all(proj->spect.weight[j], 1.0);
      
      gsl_vector_div(proj->spect.weight[j], proj->dbt.ExpDbt[j]);
      
      gsl_blas_dcopy(proj->dbt.Dbt[j], proj->spect.Dbt[j]);
    }
  
  proj->AvAtt = gsl_blas_dasum(att)/npixels;  /* gsl_vector_max(att); */

  proj->spect.defined = 1;  
}

/*+====================================================+
  
  FUNCTION: raft_weight_xfct
  
  Compute the weight function for the XFCT
  Radon transform.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  attT - transmission attenuation map
  attF - XFCT attenuation map
  
  Output:
  
  proj - projection workspace. Weigth function is 
         computed within the projection workspace
	      
  +====================================================+
*/

void raft_weight_xfct(raft_scan_t *data,
		      raft_proj_t *proj,
		      gsl_vector *attT,
		      gsl_vector *attF)
{
  div_t D;
  int nrays, nviews, i, j, N, npixels;
  int jmin, jmax, nfrays, I, J, aperture;
  double delta, apert, dth, Log;

  nrays = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);
  npixels = raft_scan_get_npixels(data);

  N = (int)(nviews/2);
  
  aperture = data->xfctAperture;
  dth      = raft_scan_get_angle(data,1)-raft_scan_get_angle(data,0);
  apert    = MPI - 2 * aperture * dth;
  
  proj->xfct.apert = apert;  
  
  if(!proj->xfct.definedDbtT)
    {
      raft_projection_radon_dbt(data, 
				proj, 
				attT, 
				proj->xfct.DbtT);
      
      raft_projection_workspace_get_exp_dbt(proj, proj->xfct.ExpDbtT);
      
      proj->xfct.definedDbtT = 1;
    }
  
  raft_projection_radon_dbt(data, 
			    proj, 
			    attF, 
			    proj->xfct.DbtF);
  
  raft_projection_workspace_get_exp_dbt(proj, proj->xfct.ExpDbtF);
  
  
  nfrays = N - 2*aperture;
  delta  = apert/nfrays;
    
  for(j = 0; j < nviews; j++)
    {      
      jmin = j + aperture;
      jmax = j + N - aperture;
      
      /*------------------*/
      /* XFCT attenuation */

      gsl_vector_set_all(proj->xfct.weight[j], 0.0);
      
      if(jmax==jmin)
	{
	  delta = 1.0;
	  i     = jmin;
	  
	  I = i - floor(i/nviews) * nviews; 
	  
	  gsl_vector_set_all(proj->xfct.ones, 1.0);
	  
	  gsl_vector_div(proj->xfct.ones, proj->xfct.ExpDbtF[I]);
	  
	  gsl_vector_add(proj->xfct.weight[j], proj->xfct.ones);
	}
      else
	{
	  for(i = jmin; i < jmax; i++)
	    {
	      I = i - floor(i/nviews) * nviews; 
	      
	      gsl_vector_set_all(proj->xfct.ones, 1.0);
	      
	      gsl_vector_div(proj->xfct.ones, proj->xfct.ExpDbtF[I]);
	      
	      gsl_vector_add(proj->xfct.weight[j], proj->xfct.ones);	      
	    }
	}
     
      gsl_blas_dscal(delta, proj->xfct.weight[j]);
      
      /* averaging in the fluorescence fan */
      gsl_blas_dscal(1/apert, proj->xfct.weight[j]);
      /**/
  
      /*--------------------------*/
      /* transmission attenuation */
      
      D = div(j + N, nviews);
      J = D.rem;               
      
      gsl_vector_div(proj->xfct.weight[j], proj->xfct.ExpDbtT[J]);           

      /*------------------------------------------*/
      /* log of weightF/apert: for XFCT inversion */
      
      for(i=0; i<npixels; i++)
	{
	  Log = log( gsl_vector_get(proj->xfct.weight[j],i)/apert );
	  
	  gsl_vector_set(proj->xfct.LogWeightF[j], i, Log);
	}
    }
  
  proj->AvAtt = gsl_blas_dasum(attT)/npixels;

  proj->xfct.defined = 1;
}

/*+====================================================+
  
  FUNCTION: raft_weight_xfct_partial
  
  Compute the partial weight function for the XFCT
  Radon transform.
  
  Only for the perpendicular fluorescence ray.  

  Input: 
  
  data - scan data . See <raft_scan_t>.
  attT - transmission attenuation map
  attF - XFCT attenuation map
  
  Output:
  
  proj - projection workspace. Weigth function is 
         computed within the projection workspace
	      
  +====================================================+
*/

void raft_weight_xfct_partial(raft_scan_t *data,
			      raft_proj_t *proj,
			      gsl_vector *attT,
			      gsl_vector *attF)
{
  div_t D;
  int nrays, nviews, npixels, i,j, N, beta, J, aperture;
  double apert, dth, Log;
  
  nrays = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);
  npixels = raft_scan_get_npixels(data);

  N        = (int)(nviews/2);
  aperture = data->xfctAperture;
  dth      = raft_scan_get_angle(data,1)-raft_scan_get_angle(data,0);
  apert    = MPI - 2 * aperture * dth;

  if(!proj->xfct.definedDbtT)
    {
      raft_projection_radon_dbt(data, 
				proj, 
				attT, 
				proj->xfct.DbtT);
  
      raft_projection_workspace_get_exp_dbt(proj, proj->xfct.ExpDbtT);
      
      proj->xfct.definedDbtT = 1;
    }    
  
  raft_projection_radon_dbt(data, 
			    proj, 
			    attF, 
			    proj->xfct.DbtF);
  
  raft_projection_workspace_get_exp_dbt(proj, proj->xfct.ExpDbtF);

  
  /* perpendicular fluorescence ray*/
  
  beta = (int)(N/2) - data->xfctAperture;
  
  for(j = 0; j < nviews; j++)
    { 
      gsl_vector_set_all(proj->xfct.weightP[j], 1.0);
      
      /*------------------*/
      /* XFCT attenuation */
      
      D = div(j + beta, nviews);
      J = D.rem;   
      
      gsl_vector_div(proj->xfct.weightP[j], proj->xfct.ExpDbtF[J]);
      
      /*--------------------------*/
      /* transmission attenuation */

      gsl_vector_div(proj->xfct.weightP[j], proj->xfct.ExpDbtT[j]);

      /*--------------------------------------------*/
      /* log of weightP: for XFCT partial inversion */
      
      for(i=0; i<npixels; i++)
	{
	  Log = log( gsl_vector_get(proj->xfct.weightP[j],i)/apert );
	  
	  gsl_vector_set(proj->xfct.LogWeightP[j], i, Log);
	}
    }
  
  proj->AvAtt = gsl_blas_dasum(attT)/npixels;

  proj->xfct.apert   = apert; 
  proj->xfct.defined = 1;
}

/*+====================================================+
  
  FUNCTION: raft_weight_spect_view
  
  Compute the weight function for the SPECT radon 
  transform, at a given view number.
  
  Input: 
  
  data - scan data. See <raft_scan_t>.
  att - attenuation map
  j - view number
  
  Output:
  
  proj - projection workspace. Weigth function is 
         computed within this workspace.
	      
  +====================================================+
*/

void raft_weight_spect_view(raft_scan_t *data,
			    raft_proj_t *proj,
			    gsl_vector *att,
			    int j)
{
  raft_projection_radon_dbt_view(data, proj, j, att, proj->dbt.Dbt[j]);
  
  gsl_vector_set_all(proj->spect.weight[j], 1.0);
  
  gsl_vector_div(proj->spect.weight[j], proj->dbt.ExpDbt[j]);
  
  gsl_blas_dcopy(proj->dbt.Dbt[j], proj->spect.Dbt[j]);

  /*----------------------------------*/
  /* still undefined: only for a view */
  proj->spect.defined = 0;   
}


/*+====================================================+
  
  FUNCTION: raft_weight_xfct_view
  
  Compute the weight function for the XFCT
  Radon transform, at a given view number.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  attT - transmission attenuation map
  attF - XFCT attenuation map
  j - view number
  
  Output:
  
  proj - projection workspace. Weigth function is 
         computed within the projection workspace
	      
  +====================================================+
*/

void raft_weight_xfct_view(raft_scan_t *data,
			   raft_proj_t *proj,
			   gsl_vector *attT,
			   gsl_vector *attF,
			   int j)
{
  div_t D;
  int nrays, nviews, i, N, npixels;
  int jmin, jmax, nfrays, I, J, aperture;
  double delta, apert, dth, Log;

  nrays = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);
  npixels = raft_scan_get_npixels(data);

  N = (int)(nviews/2);
  
  aperture = data->xfctAperture;
  dth      = raft_scan_get_angle(data,1)-raft_scan_get_angle(data,0);
  apert    = MPI - 2 * aperture * dth;
  
  proj->xfct.apert = apert;  
  
  nfrays = N - 2*aperture;
  delta  = apert/nfrays;    
  
  jmin = j + aperture;
  jmax = j + N - aperture;
  
  /*------------------*/
  /* XFCT attenuation */
  
  gsl_vector_set_all(proj->xfct.weight[j], 0.0);
  
  if(jmax==jmin)
    {
      delta = 1.0;
      i     = jmin;
      
      I = i - floor(i/nviews) * nviews; 
      
      gsl_vector_set_all(proj->xfct.ones, 1.0);
      
      raft_projection_radon_dbt_view(data, proj, I, attF, proj->dbt.Dbt[I]);
      
      gsl_vector_div(proj->xfct.ones, proj->dbt.ExpDbt[I]);      
      
      gsl_vector_add(proj->xfct.weight[j], proj->xfct.ones);
    }
  else
    {
      for(i = jmin; i < jmax; i++)
	{
	  I = i - floor(i/nviews) * nviews; 
	  
	  gsl_vector_set_all(proj->xfct.ones, 1.0);
	  
	  raft_projection_radon_dbt_view(data, proj, I, attF, proj->dbt.Dbt[I]);
	  
	  gsl_vector_div(proj->xfct.ones, proj->dbt.ExpDbt[I]);  
	  
	  gsl_vector_add(proj->xfct.weight[j], proj->xfct.ones);	      
	}
    }
  
  gsl_blas_dscal(delta, proj->xfct.weight[j]);
  
  /*--------------------------*/
  /* transmission attenuation */

  D = div(j + N, nviews);
  J = D.rem;               
  
  raft_projection_radon_dbt_view(data, proj, J, attT, proj->dbt.Dbt[J]);
  
  gsl_vector_div(proj->xfct.weight[j], proj->dbt.ExpDbt[J]);           
  
  /*------------------------------------------*/
  /* log of weightF/apert: for XFCT inversion */
  
  for(i=0; i<npixels; i++)
    {
      Log = log( gsl_vector_get(proj->xfct.weight[j],i)/apert );
      
      gsl_vector_set(proj->xfct.LogWeightF[j], i, Log);
    }

  /*----------------------------------*/
  /* still undefined: only for a view */
  proj->xfct.defined = 0;
}


/*######################################################
  Section: Getting weight
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_weight_pet_get
  
  Get PET weight function from workspace
  
  Input: 
  
  proj - projection workspace
  
  Return:

  GSL view of the exponential of the Radon transform
  (PET).
   
  +====================================================+
*/

gsl_vector *
raft_weight_pet_get(raft_proj_t *proj)
{
  return proj->pet.weight;
}

/*+====================================================+
  
  FUNCTION: raft_weight_spect_get
  
  Get SPECT weight function from workspace
  
  Input: 
  
  proj - projection workspace
  
  Output:

  weight - weight matrix

  +====================================================+
*/

void
raft_weight_spect_get(raft_proj_t *proj,
		      gsl_vector **weight)
{
  int j;

  for(j = 0; j < proj->nviews; j++)
    {
      gsl_blas_dcopy(proj->spect.weight[j], weight[j]);
    }  
}

/*+====================================================+
  
  FUNCTION: raft_weight_xfct_get
  
  Get XFCT weight function from workspace
  
  Input: 
  
  proj - projection workspace
     
  Output:

  weight - weight matrix

  +====================================================+
*/

void
raft_weight_xfct_get(raft_proj_t *proj,
		     gsl_vector **weight)
{
  int j;

  for(j = 0; j < proj->nviews; j++)
    {
      gsl_blas_dcopy(proj->xfct.weight[j], weight[j]);
    } 
}

/*+====================================================+
  
  FUNCTION: raft_weight_spect_get_view
  
  Get SPECT weight function from workspace for a given
  angle.
  
  Input: 
  
  proj - projection workspace
  view - view angle
     
  Return:

  Whenever defined, it returns the weight vector. 

  +====================================================+
*/

gsl_vector *
raft_weight_spect_get_view(raft_proj_t *proj,
			   int view)
{
  return proj->spect.weight[view];
}

/*+====================================================+
  
  FUNCTION: raft_weight_xfct_get_view
  
  Get XFCT weight function from workspace for
  a given angle.
  
  Input: 
  
  proj - projection workspace
  view - view angle
     
  Return:

  Whenever defined, it returns the weight vector.
  
  +====================================================+
*/

gsl_vector *
raft_weight_xfct_get_view(raft_proj_t *proj,
			  int view)
{
  return proj->xfct.weight[view];
}


/*######################################################
  Section: Computing cbf weight functions
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_weight_pet_cbf
  
  Compute cbf-PET weight function.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  att - attenuation map
  proj - projection workspace. See <raft_proj_t>.

  Output:
  
  cbf - cbf workspace. See <raft_cbf_t>.  
  
  +====================================================+
*/

void raft_weight_pet_cbf(raft_scan_t *data,
			 raft_proj_t *proj,
			 raft_cbf_t *cbf,
			 gsl_vector *att)
{
  /*----------------*/
  /* current weight */
  
  raft_projection_radon(data, proj, att, proj->pet.Radon);

  gsl_blas_dcopy(proj->pet.ExpRadon, proj->likelihood.exprad);
  
  gsl_vector_set_all(cbf->projection.weight.pet.current, 1.0);
      
  gsl_vector_div(cbf->projection.weight.pet.current, 
		 proj->pet.ExpRadon);
     
  /*----------------*/
  /* forward weight */
  
  raft_projection_radon(data, proj, cbf->attplus, proj->pet.Radon);
  
  gsl_blas_dcopy(proj->pet.ExpRadon, proj->likelihood.exprad);
  
  gsl_vector_set_all(cbf->projection.weight.pet.forward, 1.0);
      
  gsl_vector_div(cbf->projection.weight.pet.forward, 
		 proj->pet.ExpRadon);

  /*-----------------*/
  /* backward weight */

  raft_projection_radon(data, proj, cbf->attminus, proj->pet.Radon);
  
  gsl_blas_dcopy(proj->pet.ExpRadon, proj->likelihood.exprad);

  gsl_vector_set_all(cbf->projection.weight.pet.backward, 1.0);
      
  gsl_vector_div(cbf->projection.weight.pet.backward, 
		 proj->pet.ExpRadon);

  /*--------------------------*/
  /* define on proj workspace */

  gsl_blas_dcopy(cbf->projection.weight.pet.current, proj->pet.weight);

  proj->pet.defined = 1;

}


/*+====================================================+
  
  FUNCTION: raft_weight_spect_cbf 
  
  Compute cbf-SPECT weight function.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  att - attenuation map
  proj - projection workspace. See <raft_proj_t>.
  
  Output:
  
  cbf - cbf workspace. See <raft_cbf_t>.
	      
  +====================================================+
*/

void raft_weight_spect_cbf(raft_scan_t *data,
			   raft_proj_t *proj,
			   raft_cbf_t *cbf,
			   gsl_vector *att)
{
  int j, nviews;
  
  nviews = raft_scan_get_nviews(data);
  
  raft_cbf_radon_dbt(data, proj, cbf, att);
  
  for(j = 0; j < nviews; j++)
    {      
      /*----------------*/
      /* current weight */

      gsl_vector_set_all(cbf->projection.weight.spect.current[j], 1.0);
      gsl_vector_div(cbf->projection.weight.spect.current[j], 
		     cbf->dbt.ExpDbt[j]);
      
      /*----------------*/
      /* forward weight */

      gsl_vector_set_all(cbf->projection.weight.spect.forward[j], 1.0);
      gsl_vector_div(cbf->projection.weight.spect.forward[j], 
		     cbf->dbt.ExpDbtForward[j]);
      
      /*-----------------*/
      /* backward weight */

      gsl_vector_set_all(cbf->projection.weight.spect.backward[j], 1.0);
      gsl_vector_div(cbf->projection.weight.spect.backward[j], 
		     cbf->dbt.ExpDbtBackward[j]);

      /*--------------------------*/
      /* define on proj workspace */
      
      gsl_blas_dcopy(cbf->projection.weight.spect.current[j], proj->spect.weight[j]);
    }

  proj->spect.defined = 1;  
}

/*+====================================================+
  
  FUNCTION: raft_weight_xfct_cbf
  
  Compute cbf-XFCT weight function.
  
  If the identifier 'id' is RAFT_SPECT, the cbf data will
  be computed with respect to the transmission attenuation
  map. Otherwise, if 'id' is RAFT_XFCT, they will be
  computed with the XFCT attenuation map.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  attT - transmission attenuation map
  attF - XFCT attenuation map
  id - identifier RAFT_XFCT or RAFT_SPECT
  proj - projection workspace
  
  Output:
  
  cbf - cbf workspace. See <raft_cbf_t>.

  +====================================================+
*/

void raft_weight_xfct_cbf(raft_scan_t *data,
			  raft_proj_t *proj,
			  raft_cbf_t *cbf,
			  gsl_vector *attT,
			  gsl_vector *attF,
			  int id)
{
  div_t D;
  int nrays, nviews, i, j, N;
  int jmin, jmax, nfrays, I, J, aperture;
  double delta, apert, dth;
    
  nrays = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);

  N        = (int)(nviews/2);
  aperture = data->xfctAperture;  
  dth      = raft_scan_get_angle(data,1)-raft_scan_get_angle(data,0);
  apert    = MPI - 2 * aperture * dth;
  
  if(id == RAFT_XFCT)
    {
      raft_projection_radon_dbt(data, 
				proj, 
				attT, 
				proj->xfct.DbtT);
  
      raft_projection_workspace_get_exp_dbt(proj,
					    proj->xfct.ExpDbtT);
      
      raft_cbf_radon_dbt(data, proj, cbf, attF);
      
      
      for(j = 0; j < nviews; j++)
	{
	  jmin = j + aperture;
	  jmax = j + N - aperture;
	  
	  nfrays = jmax - jmin;
	  
	  gsl_vector_set_all(cbf->projection.weight.xfct.current[j], 0.0);
	  gsl_vector_set_all(cbf->projection.weight.xfct.forward[j], 0.0);
	  gsl_vector_set_all(cbf->projection.weight.xfct.backward[j],0.0);	  
	  
	  for(i = jmin; i < jmax; i++)
	    {
	      I = i - floor(i/nviews) * nviews; 
	      
	      gsl_vector_set_all(proj->xfct.ones, 1.0);
	      gsl_vector_div(proj->xfct.ones, cbf->dbt.ExpDbt[I]);
	      gsl_vector_add(cbf->projection.weight.xfct.current[j], proj->xfct.ones);
	      
	      gsl_vector_set_all(proj->xfct.ones, 1.0);
	      gsl_vector_div(proj->xfct.ones, cbf->dbt.ExpDbtForward[I]);
	      gsl_vector_add(cbf->projection.weight.xfct.forward[j], proj->xfct.ones);

	      gsl_vector_set_all(proj->xfct.ones, 1.0);
	      gsl_vector_div(proj->xfct.ones, cbf->dbt.ExpDbtBackward[I]);
	      gsl_vector_add(cbf->projection.weight.xfct.backward[j], proj->xfct.ones);
	    }
	  
	  delta = apert/nfrays;
	  
	  gsl_vector_scale(cbf->projection.weight.xfct.current[j], delta);
	  gsl_vector_scale(cbf->projection.weight.xfct.forward[j], delta);
	  gsl_vector_scale(cbf->projection.weight.xfct.backward[j],delta);

	  D = div(j + N, nviews);
	  J = D.rem;
	  
	  gsl_vector_div(cbf->projection.weight.xfct.current[j], proj->xfct.ExpDbtT[J]);
	  gsl_vector_div(cbf->projection.weight.xfct.forward[j], proj->xfct.ExpDbtT[J]);
	  gsl_vector_div(cbf->projection.weight.xfct.backward[j], proj->xfct.ExpDbtT[J]);
	}
    }
  else
    {
      raft_projection_radon_dbt(data, 
				proj, 
				attF, 
				proj->xfct.DbtF);
      
      raft_projection_workspace_get_exp_dbt(proj,
					    proj->xfct.ExpDbtF);

      raft_cbf_radon_dbt(data, proj, cbf, attT);

      
      for(j = 0; j < nviews; j++)
	{	  
	  jmin = j + aperture;
	  jmax = j + N - aperture;
	  
	  nfrays = jmax - jmin;
	  
	  gsl_vector_set_all(cbf->projection.weight.xfct.current[j], 0.0);
	  gsl_vector_set_all(cbf->projection.weight.xfct.forward[j], 0.0);
	  gsl_vector_set_all(cbf->projection.weight.xfct.backward[j],0.0);	  
	  
	  for(i = jmin; i < jmax; i++)
	    {
	      I = i - floor(i/nviews) * nviews; 
	      
	      gsl_vector_set_all(proj->xfct.ones, 1.0);
	      gsl_vector_div(proj->xfct.ones, proj->xfct.ExpDbtF[I]);
	      
	      gsl_vector_add(cbf->projection.weight.xfct.current[j], proj->xfct.ones);
	      gsl_vector_add(cbf->projection.weight.xfct.forward[j], proj->xfct.ones);
	      gsl_vector_add(cbf->projection.weight.xfct.backward[j],proj->xfct.ones);
	    }
	  
	  delta = apert/nfrays;
	  
	  gsl_vector_scale(cbf->projection.weight.xfct.current[j], delta);
	  gsl_vector_scale(cbf->projection.weight.xfct.forward[j], delta);
	  gsl_vector_scale(cbf->projection.weight.xfct.backward[j],delta);
	  
	  D = div(j + N, nviews);
	  J = D.rem;	  

	  gsl_vector_div(cbf->projection.weight.xfct.current[j], 
			 cbf->dbt.ExpDbt[J]);
	  
	  gsl_vector_div(cbf->projection.weight.xfct.forward[j], 
			 cbf->dbt.ExpDbtForward[J]);
	  
	  gsl_vector_div(cbf->projection.weight.xfct.backward[j],
			 cbf->dbt.ExpDbtBackward[J]);
	}

  
      /*--------------------------*/
      /* define on proj workspace */
      
      gsl_blas_dcopy(cbf->projection.weight.xfct.current[j],
		     proj->xfct.weight[j]);
    }
  proj->xfct.defined = 1;
}

/*######################################################
  Section: Computing average weight functions
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_weight_average
  
  Compute average weight function.
  
  The average function is defined by

  a(x) = (1/2\pi) \int_0^{2\pi} w(x,\theta) d\theta 
  
  where w(x,\theta) is the weigth function.The resulting
  a = a(x) is a square phantom defining the average
  attenuation at a pixel, over all directions.

  Input: 
  
  data - scan data . See <raft_scan_t>.
  weight - weight function (for SPECT of XFCT)
    
  Output:
  
  a - phantom for the average function . See <raft_cbf_t>.

  +====================================================+
*/

void raft_weight_average(raft_scan_t *data,
			 gsl_vector **weight,
			 gsl_vector *a)
{
  int i,j,k,size,nviews;
  double sum, w;
  
  size   = raft_scan_get_size(data);
  nviews = raft_scan_get_nviews(data);
  
  for(i=0;i<size; i++)
    {
      for(j=0; j<size; j++)
	{
	  sum = 0;
	  
	  for(k=0; k < nviews; k++)
	    {
	      sum += raft_phantom_get(weight[k],i,j);
	    }
	  
	  w = sum/nviews;
	    	  
	  raft_phantom_set(a, i, j, w);
	}
    }
}

/*######################################################
  Section: Computing bounds for Kunyansky's algorithm
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_weight_boundKw
  
  Compute bounds (from Miqueles & De Pierro) for the 
  contraction operator on Kunyansky's weighted algorithm.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  weight - weight function (for SPECT or XFCT)
    
  Return:
  
  Contraction bound

  +====================================================+
*/

double raft_weight_boundKw(raft_scan_t *data,
			   gsl_vector **weight)
{
  int j,J,N,nviews;
  double max, maxa, h, bnd;
  div_t D;
  CBLAS_INDEX_t imax;
  
  raft_weight_average(data, weight, data->memory.Ave);
  
  maxa   = gsl_vector_max(data->memory.Ave); 
  h      = -0.50/maxa;
  nviews = raft_scan_get_nviews(data);
  N      = (int)(nviews/2); 

  gsl_vector_set_all(data->memory.Ave, 1.0);

  for (j=0; j< nviews; j++)
    {
      D = div(j + N, nviews);
      J = D.rem; 
      
      gsl_blas_dcopy(weight[j], data->memory.AveMinusWeight);	  
      gsl_blas_daxpy(1.0, weight[J], data->memory.AveMinusWeight);
      gsl_blas_dscal(h, data->memory.AveMinusWeight);
      gsl_blas_daxpy(1.0, data->memory.Ave, data->memory.AveMinusWeight);
      
      imax = gsl_blas_idamax(data->memory.AveMinusWeight);
      max  = gsl_vector_get(data->memory.AveMinusWeight, imax);
      
      gsl_vector_set(data->memory.MaxAveMinusWeight, j, fabs(max));
    } 
  
  bnd = gsl_vector_max(data->memory.MaxAveMinusWeight);  
  
  return bnd;  
}

/*+====================================================+
  
  FUNCTION: raft_weight_boundKuw
  
  Compute bounds (from Miqueles & De Pierro) for the 
  contraction operator on Kunyansky's unweighted algorithm.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  weight - weight function (for SPECT or XFCT)
    
  Return:
  
  Contraction bound

  +====================================================+
*/

double raft_weight_boundKuw(raft_scan_t *data,
			    gsl_vector **weight)
{
  int j,J,N,nviews;
  double max, h, bnd;
  div_t D;
  CBLAS_INDEX_t imax;

  nviews = raft_scan_get_nviews(data);
  
  N = (int)(nviews/2);
  h = -0.5;
  
  gsl_vector_set_all(data->memory.Ave, 1.0);
  
  for (j=0; j< nviews; j++)
    {
      D = div(j + N, nviews);
      J = D.rem; 
      
      gsl_blas_dcopy(weight[j], data->memory.AveMinusWeight);	  
      gsl_blas_daxpy(1.0, weight[J], data->memory.AveMinusWeight);
      gsl_blas_dscal(h, data->memory.AveMinusWeight);
      gsl_blas_daxpy(1.0, data->memory.Ave, data->memory.AveMinusWeight);
      
      imax = gsl_blas_idamax(data->memory.AveMinusWeight);
      max  = gsl_vector_get(data->memory.AveMinusWeight, imax);
      
      gsl_vector_set(data->memory.MaxAveMinusWeight, j, fabs(max));
    }
  
  bnd = gsl_vector_max(data->memory.MaxAveMinusWeight);

  return bnd;
}
