#include "cbf.h"
#include "projection.h"
#include "raft_weight.h"
#include "raft_projection.h"
#include <gsl/gsl_math.h>

/*######################################################
  Title: CBF 
  
  Header - <raft/raft_cbf.h>
  Type - <raft_cbf_t>
  $Id: cbf.c,v 1.6 2009-11-29 18:30:54 miqueles Exp $ - Last Update
  ####################################################*/

/*+====================================================+
  
  FUNCTION  eval_rayintegral_dbt_cbf
  
  Compute the ray integral and cbf data for 
  divergent beam transform
  
  Input  

  scan - scan data . See <raft_scan_t>.
  att - attenuation phantom vector
  j - angle index
  cbf - cbf workspace
  xi - x-coordinate index (column)
  yi - y-coordinate index (index)
  m - pixel number
  
  Output 

  current  - ray integral at current attenuation map
  forward  - ray integral at forward attenuation map
  backward - ray integral at backward attenuation map

  +====================================================+
*/

void eval_rayintegral_dbt_cbf(raft_scan_t *data, 
			      gsl_vector *att,
			      int j,
			      int xi,
			      int yi,
			      raft_cbf_t *cbf,
			      double *current,
			      double *forward,
			      double *backward)
{
  int i, begin, end, size, row, column;
  double sum[3], cost, sint, tant, x, y, X, Y;
  double fsint, fcost, tol, step;
  gsl_vector_view atten, attenp, attenm;
  
  tol = 1/sqrt(2);
  step = raft_scan_get_x(data,1)-raft_scan_get_x(data,0);
  size = raft_scan_get_size(data);

  cost  = gsl_vector_get(data->costheta, j);
  sint  = gsl_vector_get(data->sintheta, j);
  tant  = gsl_vector_get(data->tantheta, j);
  fsint = fabs(sint);
  fcost = fabs(cost);
  
  x = raft_scan_get_x(data, xi);
  y = raft_scan_get_y(data, size-1-yi);
  
  if(fcost<tol)
    {
      if(SIGN(sint)<0)
	{
	  begin = xi;
	  end   = size;
	}
      else
	{
	  begin = 0;
	  end   = xi+1;
	}
      
      sum[0] = 0;
      sum[1] = 0;
      sum[2] = 0;

      for(i=begin; i < end; i++)
	{
	  X   = raft_scan_get_x(data,i);
	  Y   = y + (x-X)/tant;
	  row = interp_nearest(data->y, Y, RAFT_REG);
	  
	  /* column of the attenuation phantom matrix*/
	  atten  = gsl_vector_subvector_with_stride(att, i, size, size);
	  attenp = gsl_vector_subvector_with_stride(cbf->attplus, i, size, size);
	  attenm = gsl_vector_subvector_with_stride(cbf->attminus, i, size, size);
	  
	  if(SIGN(row)>0)	  
	    {
	      sum[0] += interp_linear_reverse(&atten.vector, data->y, row, Y);
	      sum[1] += interp_linear_reverse(&attenp.vector, data->y, row, Y);
	      sum[2] += interp_linear_reverse(&attenm.vector, data->y, row, Y);
	    }
	}

      *current  = (sum[0]*step)/fsint;
      *forward  = (sum[1]*step)/fsint;
      *backward = (sum[2]*step)/fsint;      
    } 
  else
    {
      if(SIGN(cost)>0) 
	{
	  begin = 0;
	  end   = yi+1;
	}
      else
	{
	  begin = yi;
	  end   = size;
	}
      
      sum[0] = 0;
      sum[1] = 0;
      sum[2] = 0;

      for(i=begin; i < end; i++)
	{
	  Y      = raft_scan_get_y(data,size-1-i);
	  X      = x + (y-Y)*tant;
	  column = interp_nearest(data->x, X, RAFT_REG);
	  
	  /* row of the activity phantom matrix*/
	  atten   = gsl_vector_subvector(att, i*size, size);
	  attenp  = gsl_vector_subvector(cbf->attplus, i*size, size);
	  attenm  = gsl_vector_subvector(cbf->attminus, i*size, size);
	  
	  if(SIGN(column)>0)	  
	    {
	      sum[0] += interp_linear(&atten.vector, data->x, column, X);
	      sum[1] += interp_linear(&attenp.vector, data->x, column, X);
	      sum[2] += interp_linear(&attenm.vector, data->x, column, X);
	    }
	}

      *current  = (sum[0]*step)/fcost;
      *forward  = (sum[1]*step)/fcost;
      *backward = (sum[2]*step)/fcost;
    }
}

/*######################################################
  Section: Allocating
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_cbf_workspace_alloc
  
  Workspace allocation for cbf procedures.
  
  Input: 
  
  data - scan dat. See <raft_scan_t>

  Output:

  cbf - cbf workspace. See <raft_cbf_t>
  
  +====================================================+
*/

void
raft_cbf_workspace_alloc(raft_scan_t *data,
			 raft_cbf_t *cbf)
{
  int j;
  int ndata, npixels, nviews;

  ndata = raft_scan_get_ndata(data);
  npixels =  raft_scan_get_npixels(data);
  nviews  = raft_scan_get_nviews(data);

  cbf->nviews = nviews;
  
  cbf->attplus   = gsl_vector_alloc(npixels);
  cbf->attminus  = gsl_vector_alloc(npixels);
  cbf->actplus   = gsl_vector_alloc(npixels);
  cbf->actminus  = gsl_vector_alloc(npixels);
  cbf->canonical = gsl_vector_alloc(npixels);
  cbf->kernel    = gsl_vector_alloc(npixels);
   
  cbf->projection.FixAct.forward  = gsl_vector_alloc(ndata);
  cbf->projection.FixAct.backward = gsl_vector_alloc(ndata);
  cbf->projection.FixAct.current  = gsl_vector_alloc(ndata);
  cbf->projection.FixAct.kernelForward  = gsl_vector_alloc(npixels);
  cbf->projection.FixAct.kernelBackward = gsl_vector_alloc(npixels);
  
  cbf->projection.FixAtt.forward  = gsl_vector_alloc(ndata);
  cbf->projection.FixAtt.backward = gsl_vector_alloc(ndata);
  cbf->projection.FixAtt.current  = gsl_vector_alloc(ndata);
  cbf->projection.FixAtt.kernelForward  = gsl_vector_alloc(npixels);
  cbf->projection.FixAtt.kernelBackward = gsl_vector_alloc(npixels);
  
  cbf->projection.weight.pet.forward  = gsl_vector_alloc(ndata);
  cbf->projection.weight.pet.backward = gsl_vector_alloc(ndata);
  cbf->projection.weight.pet.current  = gsl_vector_alloc(ndata);
    
  cbf->dbt.Dbt      = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  cbf->dbt.forward  = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  cbf->dbt.backward = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  cbf->dbt.ExpDbt   = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  cbf->dbt.ExpDbtForward  = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  cbf->dbt.ExpDbtBackward = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);

  cbf->projection.weight.spect.forward  = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  cbf->projection.weight.spect.backward = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  cbf->projection.weight.spect.current  = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  cbf->projection.weight.xfct.forward   = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  cbf->projection.weight.xfct.backward  = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);
  cbf->projection.weight.xfct.current   = (gsl_vector **)malloc(sizeof(gsl_vector *)*nviews);

  for(j = 0; j < nviews; j++)
    {
      cbf->dbt.Dbt[j]            = gsl_vector_alloc(npixels);
      cbf->dbt.forward[j]        = gsl_vector_alloc(npixels);
      cbf->dbt.backward[j]       = gsl_vector_alloc(npixels);
      cbf->dbt.ExpDbt[j]         = gsl_vector_alloc(npixels);
      cbf->dbt.ExpDbtForward[j]  = gsl_vector_alloc(npixels);
      cbf->dbt.ExpDbtBackward[j] = gsl_vector_alloc(npixels);

      cbf->projection.weight.spect.forward[j]  = gsl_vector_alloc(npixels);
      cbf->projection.weight.spect.backward[j] = gsl_vector_alloc(npixels);
      cbf->projection.weight.spect.current[j]  = gsl_vector_alloc(npixels); 
      cbf->projection.weight.xfct.forward[j]   = gsl_vector_alloc(npixels);
      cbf->projection.weight.xfct.backward[j]  = gsl_vector_alloc(npixels);
      cbf->projection.weight.xfct.current[j]   = gsl_vector_alloc(npixels);
    }
}

/*+====================================================+
  
  FUNCTION: raft_cbf_workspace_free
  
  Frees cbf workspace.
  
  Input: 
  
  cbf - cbf workspace. See <raft_cbf_t>.
    
  +====================================================+
*/

void
raft_cbf_workspace_free(raft_cbf_t *cbf)
{
  int j;

  gsl_vector_free(cbf->attplus);
  gsl_vector_free(cbf->attminus);
  gsl_vector_free(cbf->actplus);
  gsl_vector_free(cbf->actminus);
  gsl_vector_free(cbf->canonical);
  gsl_vector_free(cbf->kernel);
  
  gsl_vector_free(cbf->projection.FixAct.forward);
  gsl_vector_free(cbf->projection.FixAct.backward);
  gsl_vector_free(cbf->projection.FixAct.current);
  gsl_vector_free(cbf->projection.FixAct.kernelForward);
  gsl_vector_free(cbf->projection.FixAct.kernelBackward);
  
  gsl_vector_free(cbf->projection.FixAtt.forward);
  gsl_vector_free(cbf->projection.FixAtt.backward);
  gsl_vector_free(cbf->projection.FixAtt.current);
  gsl_vector_free(cbf->projection.FixAtt.kernelForward);
  gsl_vector_free(cbf->projection.FixAtt.kernelBackward);

  gsl_vector_free(cbf->projection.weight.pet.forward);
  gsl_vector_free(cbf->projection.weight.pet.backward);
  gsl_vector_free(cbf->projection.weight.pet.current);  
  
  for(j = 0; j < cbf->nviews; j++)
    {
      gsl_vector_free(cbf->dbt.Dbt[j]);
      gsl_vector_free(cbf->dbt.forward[j]);
      gsl_vector_free(cbf->dbt.backward[j]); 
      gsl_vector_free(cbf->dbt.ExpDbt[j]);
      gsl_vector_free(cbf->dbt.ExpDbtForward[j]);
      gsl_vector_free(cbf->dbt.ExpDbtBackward[j]);
      
      gsl_vector_free(cbf->projection.weight.spect.forward[j]);
      gsl_vector_free(cbf->projection.weight.spect.backward[j]);
      gsl_vector_free(cbf->projection.weight.spect.current[j]);  
  
      gsl_vector_free(cbf->projection.weight.xfct.current[j]);
      gsl_vector_free(cbf->projection.weight.xfct.forward[j]);
      gsl_vector_free(cbf->projection.weight.xfct.backward[j]);      
    }

  free(cbf->dbt.Dbt);
  free(cbf->dbt.forward);
  free(cbf->dbt.backward);  
  free(cbf->dbt.ExpDbt);
  free(cbf->dbt.ExpDbtForward);
  free(cbf->dbt.ExpDbtBackward);
  
  free(cbf->projection.weight.spect.forward);
  free(cbf->projection.weight.spect.backward);
  free(cbf->projection.weight.spect.current);  

  free(cbf->projection.weight.xfct.current);
  free(cbf->projection.weight.xfct.forward);
  free(cbf->projection.weight.xfct.backward);  

}

/*######################################################
  Section: Setting cbf workspace
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_cbf_set_pet
  
  Set PET cbf data within workspace.
  
  Input: 
  
  h - step size
  k - pixel number
  att - attenuation phantom
  act - activity phantom
  
  Output:

  cbf - cbf workspace. See <raft_cbf_t>.
  
  +====================================================+
*/

void
raft_cbf_set_pet(raft_scan_t *data,
		 raft_cbf_t *cbf,
		 int k,
		 double h,
		 gsl_vector *att,
		 gsl_vector *act)
{
  raft_cbf_set_spect(data, cbf, k, h, att, act);
}


/*+====================================================+
  
  FUNCTION: raft_cbf_set_spect
  
  Set SPECT cbf data within workspace.
  
  Input: 
  
  h - step size
  k - pixel number
  att - attenuation phantom
  act - activity phantom
  
  Output:

  cbf - cbf workspace. See <raft_cbf_t>.
  
  +====================================================+
*/

void
raft_cbf_set_spect(raft_scan_t *data,
		   raft_cbf_t *cbf,
		   int k,
		   double h,
		   gsl_vector *att,
		   gsl_vector *act)
{
  double current, plus, minus;
  
  cbf->step = h;
  cbf->pixel = k;
  
  gsl_vector_set_all(cbf->canonical, 0.0);
  gsl_vector_set(cbf->canonical, k, 1.0);

  /**/

  current = gsl_vector_get(att, k);
  plus    = current + h;
  minus   = current - h;
  
  gsl_vector_memcpy(cbf->attplus, att);
  gsl_vector_memcpy(cbf->attminus, att);
  
  gsl_vector_set(cbf->attplus,  k, plus);
  gsl_vector_set(cbf->attminus, k, minus);
  
  current = gsl_vector_get(act, k);
  plus    = current + h;
  minus   = current - h;
  
  gsl_vector_memcpy(cbf->actplus, act);
  gsl_vector_memcpy(cbf->actminus, act);
  
  gsl_vector_set(cbf->actplus,  k, plus);
  gsl_vector_set(cbf->actminus, k, minus);
}

/*+===================================================+
  
  FUNCTION: raft_cbf_set_diff_xfct
  
  Set XFCT cbf data within workspace.

  Input: 
  
  h - step size
  k - pixel number
  att - attenuation map (either trasmission or fluorescence)
  act - activity map
  
  Output:

  cbf - cbf workspace. See <raft_cbf_t>.
  
  +====================================================+
*/

void
raft_cbf_set_xfct(raft_scan_t *data,
		  raft_cbf_t *cbf,
		  int k,
		  double h,
		  gsl_vector *att,
		  gsl_vector *act)
{
  raft_cbf_set_spect(data,cbf,k,h,att,act);
}

/*+====================================================+
  
  FUNCTION: raft_cbf_set_dbt
  
  Set DBT cbf data within workspace.
  
  Input: 
  
  h - step size
  k - pixel number
  att - attenuation phantom
  
  Output:

  cbf - cbf workspace. See <raft_cbf_t>.
  
  +====================================================+
*/

void
raft_cbf_set_dbt(raft_scan_t *data,
		 raft_cbf_t *cbf,
		 int k,
		 double h,
		 gsl_vector *att)
{
  double current, plus, minus;
    
  cbf->step = h;
  cbf->pixel = k;
  
  gsl_vector_set_all(cbf->canonical, 0.0);
  gsl_vector_set(cbf->canonical, k, 1.0);
  
  current = gsl_vector_get(att, k);
  plus    = current + h;
  minus   = current - h;
  
  gsl_vector_memcpy(cbf->attplus, att);
  gsl_vector_memcpy(cbf->attminus, att);
  
  gsl_vector_set(cbf->attplus,  k, plus);
  gsl_vector_set(cbf->attminus, k, minus);
}

/*######################################################
  Section: Computing cbf Radon transforms
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_cbf_radon_pet
  
  Compute cbf-PET Radon transform.
  
  Input: 
  
  data - scan data. See <raft_scan_t>.
  act - activity phantom vector
  att - attenuation phantom vector
  proj - projection workspace
  
  Output:
  
  cbf - cbf workspace. See <raft_cbf_t>.
  
  +====================================================+
*/

void raft_cbf_radon_pet(raft_scan_t *data, 
			raft_proj_t *proj,
			raft_cbf_t *cbf,
			gsl_vector *act,
			gsl_vector *att)
{
  raft_weight_spect_cbf(data, proj, cbf, att);

  /*----------------*/
  /* Fixed activity */
  
  raft_projection_radon(data,
			proj,
			act, 
			cbf->projection.FixAct.current);
  
  gsl_vector_memcpy(cbf->projection.FixAct.backward,
		    cbf->projection.FixAct.current);
      
  gsl_vector_memcpy(cbf->projection.FixAct.forward,
		    cbf->projection.FixAct.current);

  gsl_vector_mul(cbf->projection.FixAct.forward, 
		 cbf->projection.weight.pet.forward);

  gsl_vector_mul(cbf->projection.FixAct.backward, 
		 cbf->projection.weight.pet.backward);
  
  gsl_vector_mul(cbf->projection.FixAct.current, 
		 cbf->projection.weight.pet.current);
  
  /*-------------------*/
  /* Fixed attenuation */
  
  gsl_vector_memcpy(cbf->projection.FixAtt.current,
		    cbf->projection.FixAct.current);
  
  raft_projection_radon(data, 
			proj,
			cbf->actplus, 
			cbf->projection.FixAtt.forward);
  
  raft_projection_radon(data,
			proj,
			cbf->actminus, 
			cbf->projection.FixAtt.backward);
  
  gsl_vector_mul(cbf->projection.FixAtt.forward, 
		 cbf->projection.weight.pet.current);
  
  gsl_vector_mul(cbf->projection.FixAtt.backward, 
		 cbf->projection.weight.pet.current);  
}


/*+====================================================+
  
  FUNCTION: raft_cbf_radon_spect
  
  Compute cbf-SPECT Radon transform.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  act - activity phantom vector
  att - attenuation phantom vector
  proj - projection workspace

  Output:
  
  cbf - cbf workspace. See <raft_cbf_t>.
  
  +====================================================+
*/

void 
raft_cbf_radon_spect(raft_scan_t *data, 
		     raft_proj_t *proj,
		     raft_cbf_t *cbf,
		     gsl_vector *act,
		     gsl_vector *att)
{
  int j, k, I, nrays, nviews;
  double r, FixActrF, FixActrB, FixAttrF, FixAttrB;
  
  nrays = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);

  raft_weight_spect_cbf(data, proj, cbf, att);
  
  for(j = 0; j < nviews; j++)
    { 
      /*-----------------*/
      /*  Fixed activity */
      
      gsl_vector_memcpy(cbf->kernel, act);
      gsl_vector_mul(cbf->kernel, 
		     cbf->projection.weight.spect.current[j]);
      
      gsl_vector_memcpy(cbf->projection.FixAct.kernelForward, act);
      gsl_vector_mul(cbf->projection.FixAct.kernelForward, 
		     cbf->projection.weight.spect.forward[j]);
      
      gsl_vector_memcpy(cbf->projection.FixAct.kernelBackward, act);
      gsl_vector_mul(cbf->projection.FixAct.kernelBackward, 
		     cbf->projection.weight.spect.backward[j]);
      
      /*-------------------*/
      /* Fixed attenuation */
      
      gsl_vector_memcpy(cbf->projection.FixAtt.kernelForward, cbf->actplus);
      gsl_vector_mul(cbf->projection.FixAtt.kernelForward, 
		     cbf->projection.weight.spect.current[j]);
      
      gsl_vector_memcpy(cbf->projection.FixAtt.kernelBackward, cbf->actminus);
      gsl_vector_mul(cbf->projection.FixAtt.kernelBackward, 
		     cbf->projection.weight.spect.current[j]);
      
      for(k = 0; k < nrays; k++)
	{
	  I = k + j*nrays;
	  
	  r = eval_rayintegral_(data,cbf->kernel,j,k,
			       data->memory.a, 
			       &data->memory.norma,
			       &data->memory.dotp);
	  
	  FixActrF = eval_rayintegral_(data,cbf->projection.FixAct.kernelForward,j,k,
				      data->memory.a, 
				      &data->memory.norma,
				      &data->memory.dotp);
	  
	  FixActrB = eval_rayintegral_(data,cbf->projection.FixAct.kernelBackward,j,k,
				      data->memory.a, 
				      &data->memory.norma,
				      &data->memory.dotp);
	  
	  FixAttrF = eval_rayintegral_(data,cbf->projection.FixAtt.kernelForward,j,k,
				      data->memory.a, 
				      &data->memory.norma,
				      &data->memory.dotp);
	  
	  FixAttrB = eval_rayintegral_(data,cbf->projection.FixAtt.kernelBackward,j,k,
				      data->memory.a, 
				      &data->memory.norma,
				      &data->memory.dotp);
	  
	  gsl_vector_set(cbf->projection.FixAct.current, I, r);
	  gsl_vector_set(cbf->projection.FixAct.forward, I, FixActrF);
	  gsl_vector_set(cbf->projection.FixAct.backward, I, FixActrB);
	  
	  gsl_vector_set(cbf->projection.FixAtt.current,I,r);
	  gsl_vector_set(cbf->projection.FixAtt.forward,I,FixAttrF);
	  gsl_vector_set(cbf->projection.FixAtt.backward,I,FixAttrB);
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_cbf_radon_xfct
  
  Compute cbf-XFCT Radon transform.

  If the identifier 'id' is RAFT_SPECT, the cbf data will
  be computed with respect to the transmission attenuation
  map. Otherwise, if 'id' is RAFT_XFCT, they will be
  computed with the XFCT attenuation map.
  
  Input: 
  
  data - scan data. See <raft_scan_t>.
  act - activity map
  attT - transmission attenuation map
  attF - XFCT attenuation map
  id - identifier for RAFT_XFCT or RAFT_SPECT
  proj - projection workspace

  Output:

  cbf - cbf workspace. See <raft_cbf_t>.
  
  +====================================================+
*/

void 
raft_cbf_radon_xfct(raft_scan_t *data, 
		    raft_proj_t *proj,
		    raft_cbf_t *cbf,
		    gsl_vector *act,
		    gsl_vector *attT,
		    gsl_vector *attF,
		    int id)
{
  int j, k, I, nrays, nviews;
  double r, FixActrF, FixActrB, FixAttrF, FixAttrB;
    
  nrays = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);

  raft_weight_xfct_cbf(data, proj, cbf, attT, attF, id);
  
  for(j = 0; j < nviews; j++)
    { 
      /*----------------*/
      /* Fixed activity */
      
      gsl_vector_memcpy(cbf->kernel, act);
      gsl_vector_mul(cbf->kernel, 
		     cbf->projection.weight.xfct.current[j]);
      
      gsl_vector_memcpy(cbf->projection.FixAct.kernelForward, act);
      gsl_vector_mul(cbf->projection.FixAct.kernelForward, 
		     cbf->projection.weight.xfct.forward[j]);
      
      gsl_vector_memcpy(cbf->projection.FixAct.kernelBackward, act);
      gsl_vector_mul(cbf->projection.FixAct.kernelBackward, 
		     cbf->projection.weight.xfct.backward[j]);
      
      /*-------------------*/
      /* Fixed attenuation */
      
      gsl_vector_memcpy(cbf->projection.FixAtt.kernelForward, cbf->actplus);
      gsl_vector_mul(cbf->projection.FixAtt.kernelForward, 
		     cbf->projection.weight.xfct.current[j]);
      
      gsl_vector_memcpy(cbf->projection.FixAtt.kernelBackward, cbf->actminus);
      gsl_vector_mul(cbf->projection.FixAtt.kernelBackward, 
		     cbf->projection.weight.xfct.current[j]);
      
      for(k = 0; k < nrays; k++)
	{
	  I = k + j*nrays;
	  
	  r = eval_rayintegral_(data,cbf->kernel,j,k,
			       data->memory.a, 
			       &data->memory.norma,
			       &data->memory.dotp);
	  
	  FixActrF = eval_rayintegral_(data,cbf->projection.FixAct.kernelForward,j,k,
				      data->memory.a, 
				      &data->memory.norma,
				      &data->memory.dotp);
	  
	  FixActrB = eval_rayintegral_(data,cbf->projection.FixAct.kernelBackward,j,k,
				      data->memory.a, 
				      &data->memory.norma,
				      &data->memory.dotp);
	  
	  FixAttrF = eval_rayintegral_(data,cbf->projection.FixAtt.kernelForward,j,k,
				      data->memory.a, 
				      &data->memory.norma,
				      &data->memory.dotp);
	  
	  FixAttrB = eval_rayintegral_(data,cbf->projection.FixAtt.kernelBackward,j,k,
				      data->memory.a, 
				      &data->memory.norma,
				      &data->memory.dotp);
	  
	  gsl_vector_set(cbf->projection.FixAct.current, I, r);
	  gsl_vector_set(cbf->projection.FixAct.forward, I, FixActrF);
	  gsl_vector_set(cbf->projection.FixAct.backward, I, FixActrB);
	  
	  gsl_vector_set(cbf->projection.FixAtt.current,I,r);
	  gsl_vector_set(cbf->projection.FixAtt.forward,I,FixAttrF);
	  gsl_vector_set(cbf->projection.FixAtt.backward,I,FixAttrB);
	}
    }
}

/*+====================================================+
  
  FUNCTION: raft_cbf_radon_dbt
  
  Compute cbf divergent beam transform.
  
  Input: 
  
  data - scan data . See <raft_scan_t>.
  att - attenuation phantom vector
  proj - projection workspace
    
  Output:

  cbf - cbf workspace. See <raft_cbf_t>.
  
  +====================================================+
*/

void raft_cbf_radon_dbt(raft_scan_t *data, 
			raft_proj_t *proj,
			raft_cbf_t *cbf,
			gsl_vector *phantom)
{
  int I, i, j, k, size, nviews, nrays;
  double r, rF, rB, er, erF, erB; 
  
  size = raft_scan_get_size(data);
  nrays  = raft_scan_get_nrays(data);
  nviews = raft_scan_get_nviews(data);
    
  for(j=0; j < nviews; j++)
    {
      for(i = 0; i < size; i++)
	{
	  for(k = 0; k < size; k++)
	    {
	      I = k + i*size; 
	      
	      eval_rayintegral_dbt_cbf(data, 
				       phantom, j, k, i,
				       cbf,
				       &r, 
				       &rF, 
				       &rB);
	      er  = exp(r); 
	      erF = exp(rF);
	      erB = exp(rB);	      
	      
	      gsl_vector_set(cbf->dbt.Dbt[j],  I, r);
	      gsl_vector_set(cbf->dbt.forward[j], I, rF);
	      gsl_vector_set(cbf->dbt.backward[j], I, rB);
	      gsl_vector_set(cbf->dbt.ExpDbt[j], I, er);
	      gsl_vector_set(cbf->dbt.ExpDbtForward[j], I, erF);
	      gsl_vector_set(cbf->dbt.ExpDbtBackward[j], I, erB);	      
	    }
	}
    }
}

/*######################################################
  Section: Getting cbf data
  ####################################################*/

/*+====================================================+
  
  FUNCTION: raft_cbf_get_forward_fact
  
  Get forward attenuated Radon transform for
  a fixed activity map.
  
  Input: 
  
  cbf - cbf workspace. See <raft_cbf_t>.
  
  Output:

  rad - attenuated Radon transform

  +====================================================+
*/

void 
raft_cbf_get_forward_fact(raft_cbf_t *cbf,
			  gsl_vector *rad)
{
  gsl_vector_memcpy(rad, cbf->projection.FixAct.forward);
}

/*+====================================================+
  
  FUNCTION: raft_cbf_get_backward_fact
  
  Get backward attenuated Radon transform for a fixed
  activity map.
  
  Input: 
  
  cbf - cbf workspace
     
  Output:

  rad - attenuated Radon transform

  +====================================================+
*/

void 
raft_cbf_get_backward_fact(raft_cbf_t *cbf,
			   gsl_vector *rad)
{
  gsl_vector_memcpy(rad, cbf->projection.FixAct.backward);
}

/*+====================================================+
  
  FUNCTION: raft_cbf_get_current_fact
  
  Get current attenuated Radon transform for a fixed
  activity map.
  
  Input: 
  
  cbf - cbf workspace. See <raft_cbf_t>.
 
  Output:

  rad - attenuated Radon transform

  +====================================================+
*/

void 
raft_cbf_get_current_fact(raft_cbf_t *cbf,
			  gsl_vector *rad)
{
  gsl_vector_memcpy(rad, cbf->projection.FixAct.current);
}

/*+====================================================+
  
  FUNCTION: raft_cbf_get_forward_fatt
  
  Get forward attenuated Radon transform for
  a fixed attenuation map.
  
  Input: 
  
  cbf - cbf workspace. See <raft_cbf_t>.
     
  Output:

  rad - attenuated Radon transform

  +====================================================+
*/

void 
raft_cbf_get_forward_fatt(raft_cbf_t *cbf,
			  gsl_vector *rad)
{
  gsl_vector_memcpy(rad, cbf->projection.FixAtt.forward);
}

/*+====================================================+
  
  FUNCTION: raft_cbf_get_backward_fatt
  
  Get backward attenuated Radon transform for a fixed
  attenuation map.
  
  Input: 
  
  cbf - cbf workspace. See <raft_cbf_t>.
  
  Output:

  rad - attenuated Radon transform

  +====================================================+
*/

void 
raft_cbf_get_backward_fatt(raft_cbf_t *cbf,
			   gsl_vector *rad)
{
  gsl_vector_memcpy(rad, cbf->projection.FixAtt.backward);
}

/*+====================================================+
  
  FUNCTION: raft_cbf_get_current_fatt
  
  Get current attenuated Radon transform for a fixed
  attenuation map.
  
  Input: 
  
  cbf - cbf workspace. See <raft_cbf_t>.
     
  Output:

  rad - attenuated Radon transform

  +====================================================+
*/

void 
raft_cbf_get_current_fatt(raft_cbf_t *cbf,
			  gsl_vector *rad)
{
  gsl_vector_memcpy(rad, cbf->projection.FixAtt.current);
}

/*+====================================================+
  
  FUNCTION: raft_cbf_get_current_dbt
  
  Get current divergent beam transform.
  
  Input: 
  
  cbf - cbf workspace. See <raft_cbf_t>.
      
  Output:

  dbt - divergent beam transform

  +====================================================+
*/

void 
raft_cbf_get_current_dbt(raft_cbf_t *cbf,
			 gsl_vector **dbt)
{
  int j;

  for(j=0; j < cbf->nviews; j++)
    {
      gsl_vector_memcpy(dbt[j], cbf->dbt.Dbt[j]);
    }
}

/*+====================================================+
  
  FUNCTION: raft_cbf_get_forward_dbt
  
  Get forward divergent beam transform.
  
  Input: 
  
  cbf - cbf workspace. See <raft_cbf_t>.
      
  Output:

  dbt - divergent beam transform

  +====================================================+
*/

void 
raft_cbf_get_forward_dbt(raft_cbf_t *cbf,
			 gsl_vector **dbt)
{
  int j;

  for(j=0; j < cbf->nviews; j++)
    {
      gsl_vector_memcpy(dbt[j], cbf->dbt.forward[j]);
    }
}

/*+====================================================+
  
  FUNCTION: raft_cbf_get_backward_dbt
  
  Get backward divergent beam transform.
  
  Input: 
  
  cbf - cbf workspace. See <raft_cbf_t>.
      
  Output:

  dbt - divergent beam transform

  +====================================================+
*/

void 
raft_cbf_get_backward_dbt(raft_cbf_t *cbf,
			  gsl_vector **dbt)
{
  int j;
  
  for(j=0; j < cbf->nviews; j++)
    {
      gsl_vector_memcpy(dbt[j], cbf->dbt.backward[j]);
    }
}

