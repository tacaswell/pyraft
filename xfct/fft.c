#include "fft.h"
#include "raft_param.h"

/*######################################################
  Private
  
  $Id: fft.c,v 1.12 2008-11-24 16:13:07 miqueles Exp $ - Last Update
  ####################################################*/

/*+====================================================+
  
  FUNCTION fft_workspace_alloc
  
  FFT workspace allocation
  
  Input 

  workspace - fft workspace
  N - number of samples
  step - step size
  
  _Remark_

  To perform two-dimensional FFT's, the library assume
  that all matrices are square with size N.
  
  +====================================================+
*/

void fft_workspace_alloc(fft_t *workspace, 
			 int N,
			 double step)
{
  workspace->N    = N;
  workspace->step = step;
   
  workspace->gslwavetable = gsl_fft_complex_wavetable_alloc(N);
  workspace->gslworkspace = gsl_fft_complex_workspace_alloc(N);
  
  workspace->real = gsl_vector_alloc(N*N);
  workspace->imag = gsl_vector_alloc(N*N);
}

/*+====================================================+
  
  FUNCTION fft_workspace_free
  
  Frees FFT workspace
  
  Input 

  workspace - fft workspace
  
  +====================================================+
*/

void fft_workspace_free(fft_t *workspace)
{
  gsl_vector_free(workspace->real);
  gsl_vector_free(workspace->imag);

  gsl_fft_complex_wavetable_free(workspace->gslwavetable);
  gsl_fft_complex_workspace_free(workspace->gslworkspace);
} 

/*+====================================================+
  
  FUNCTION fft2_complex
  
  Complex two-dimensional FFT using GSL procedures.
  
  Input 

  req - real part of matrix q
  imq - imaginary part of matrix q
  workspace - FFT workspace

  Ouput

  reQ - real part of transformed matrix Q
  imQ - imaginary part of transformed matrix Q

  Remark 1

  Matrices are represented as vectors.

  Remark 2

  GSL doesn't shift the negative frequencies of 
  Fourier transform.
  
  +====================================================+
*/

void fft2_complex(gsl_vector *reQ,
		  gsl_vector *imQ,
		  gsl_vector *req,
		  gsl_vector *imq,
		  fft_t *workspace)
{
  int j, N;
  gsl_vector_view rew, imw, rev, imv, Re, Im;
  
  N = workspace->N;

  /*--------------------------------*/
  /* FFT of each column of matrix q */

  for(j=0; j < N; j++)
    {
      rev = gsl_vector_subvector(req, j*N, N);   
      imv = gsl_vector_subvector(imq, j*N, N);   
      
      Re = gsl_vector_subvector(workspace->real, j*N, N);
      Im = gsl_vector_subvector(workspace->imag, j*N, N);
      
      fft_complex(&Re.vector, 
		  &Im.vector, 
		  &rev.vector, 
		  &imv.vector,
		  workspace);
    }

  /*---------------------------------------------------*/
  /* FFT of resultant matrix transpose(real + i*imag) */

  for(j=0; j < N; j++)
    {
      rew = gsl_vector_subvector_with_stride(workspace->real, j, N, N);   
      imw = gsl_vector_subvector_with_stride(workspace->imag, j, N, N);   
  
      Re = gsl_vector_subvector(reQ, j*N, N);
      Im = gsl_vector_subvector(imQ, j*N, N);
      
      fft_complex(&Re.vector, 
		  &Im.vector, 
		  &rew.vector, 
		  &imw.vector,
		  workspace);      
    }
}


/*+====================================================+
  
  FUNCTION ifft2_complex
  
  Complex two-dimensional IFFT using GSL procedures.
  
  Input 

  req - real part of matrix q
  imq - imaginary part of matrix q
  workspace - FFT workspace

  Ouput

  reQ - real part of transformed matrix Q
  imQ - imaginary part of transformed matrix Q

  Remark 

  Matrices are represented as vectors. GSL doesn't shift the 
  negative frequencies of Fourier transform.
  
  +====================================================+
*/

void ifft2_complex(gsl_vector *req,
		   gsl_vector *imq,
		   gsl_vector *reQ,
		   gsl_vector *imQ,
		   fft_t *workspace)
{
  int j, N;
  gsl_vector_view rew, imw, rev, imv, Re, Im;
  
  N = workspace->N;

  /*---------------------------------*/
  /* IFFT of each column of matrix Q */

  for(j=0; j < N; j++)
    {
      rev = gsl_vector_subvector(reQ, j*N, N);   
      imv = gsl_vector_subvector(imQ, j*N, N);   
  
      Re = gsl_vector_subvector(workspace->real, j*N, N);
      Im = gsl_vector_subvector(workspace->imag, j*N, N);
      
      ifft_complex(&Re.vector, 
		   &Im.vector, 
		   &rev.vector, 
		   &imv.vector,
		   workspace);      
    }

  /*------------------------------------------------*/
  /* IFFT of resultant matrix transpose(re + i*im) */

  for(j=0; j < N; j++)
    {
      rew = gsl_vector_subvector_with_stride(workspace->real, j, N, N);   
      imw = gsl_vector_subvector_with_stride(workspace->imag, j, N, N);   
  
      Re = gsl_vector_subvector(req, j*N, N);
      Im = gsl_vector_subvector(imq, j*N, N);
      
      ifft_complex(&Re.vector, 
		   &Im.vector, 
		   &rew.vector, 
		   &imw.vector,
		   workspace);      
    }
}


/*+====================================================+
  
  FUNCTION fft_complex
  
  Complex FFT using GSL procedures.
  
  Input 

  rex - real part of signal 
  imx - imaginary part of signal
  workspace - FFT workspace

  Ouput

  reX - real part of transformed signal
  imX - imaginary part of transformed signal

  Remark 

  GSL doesn't shift the negative frequencies of 
  Fourier transform.
  
  +====================================================+
*/

void fft_complex(gsl_vector *reX,
		 gsl_vector *imX,
		 gsl_vector *rex,
		 gsl_vector *imx,
		 fft_t *workspace)
{
  int i;
  double data[2*workspace->N];
  
  for(i=0; i<workspace->N; i++){
    REAL(data,i) = gsl_vector_get(rex,i);
    IMAG(data,i) = gsl_vector_get(imx,i);
  }
  
  /*-----------------------------*/
  /* complex transform (forward) */

  gsl_fft_complex_forward(data, 
			  1, 
			  workspace->N, 
			  workspace->gslwavetable, 
			  workspace->gslworkspace);

  for(i=0;i<workspace->N;i++)
    {
      gsl_vector_set(reX, i, REAL(data,i));
      gsl_vector_set(imX, i, IMAG(data,i));
    }
}

/*+====================================================+
  
  FUNCTION ifft_complex
  
  Complex IFFT using GSL procedures.
  
  Input 

  reX - real part of transformed signal 
  imX - imaginary part of transformed signal
  workspace - FFT workspace

  Ouput

  rex - real part of original signal
  imx - imaginary part of original signal

  Remark 

  GSL doesn't shift the negative frequencies of 
  Fourier transform.
  
  +====================================================+
*/

void ifft_complex(gsl_vector *rex,
		  gsl_vector *imx,
		  gsl_vector *reX,
		  gsl_vector *imX,
		  fft_t *workspace)
{
  int i;
  double data[2*workspace->N];
  
  for(i=0; i<workspace->N; i++){
    REAL(data, i) = gsl_vector_get(reX, i);
    IMAG(data, i) = gsl_vector_get(imX, i);
  }

  /*-----------------------------*/
  /* complex transform (inverse) */

  gsl_fft_complex_inverse(data, 
			  1, 
			  workspace->N, 
			  workspace->gslwavetable, 
			  workspace->gslworkspace);
  
  for(i=0; i<workspace->N; i++)
    {
      gsl_vector_set(rex, i, REAL(data,i));
      gsl_vector_set(imx, i, IMAG(data,i));
    }
}

/*+====================================================+
  
  FUNCTION fft_shift
  
  Swaps the left and right halves of a vector.
  
  Input 

  x - vector 
  N - dimension
  
  +====================================================+
*/

void fft_shift(gsl_vector *x, int N)
{
  int N2;
  gsl_vector_view half1, half2;

  N2 = floor(N/2);

  gsl_vector_reverse(x);
  
  half1 = gsl_vector_subvector(x, 0, N2); 
  half2 = gsl_vector_subvector(x, N2, N-N2); 

  gsl_vector_reverse(&half1.vector);
  gsl_vector_reverse(&half2.vector);  
}
