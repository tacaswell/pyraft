#include "raft_math.h"
#include "fft.h"

/*######################################################
  Private
  
  $Id: mathematics.c,v 1.18 2009-05-19 19:09:02 miqueles Exp $ - Last Update
  ####################################################*/

/*+====================================================+
  
  FUNCTION intRandomGen
  
  GSL integer random generator between m and n
  
  Input 

  n - integer
    
  Return

  Generated integer.

  +====================================================+
*/

int
intRandomGen(int n)
{  
  const gsl_rng_type * T;
  gsl_rng * r;  
  int gn;
  
  gsl_rng_env_setup();
  
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  gn = (int) gsl_rng_uniform_int (r, n);
  
  gsl_rng_free(r);
  
  return gn;  
}


/*+====================================================+
  
  FUNCTION maxvector
  
  Computation of max(0,x) for a GSL vector
  
  Input 

  x - GSL vector
    
  Output

  x = max(0,x)

  +====================================================+
*/

void
maxVector(gsl_vector *x)
{
  int i, size;
  double xx;
  
  size = x->size;

  for(i = 0; i<size; i++)
    {
      xx = gsl_vector_get(x,i);
      
      gsl_vector_set(x,i,MAX(0,xx));
    }
    
}

/*+====================================================+
  
  FUNCTION eval_chebyshev2
  
  Evaluate a Chebyshev polynomial of the second kind.
  
  Input 

  x - real value
  k - polynomial order
  
  Output

  u - evaluated polynomial

  +====================================================+
*/

void eval_chebyshev2(int k, 
		     double x, 
		     double *u)
{
  double U, a;

  a = acos(x);
  U = sin(a);
  a = a * (k+1);
  U = sin(a)/U;

  *u = U;
}


/*+====================================================+
  
  FUNCTION eval_array_chebyshev2
  
  Evaluate a Chebyshev polynomial of the second kind.
  
  Input 

  x - array of real values
  k - polynomial order

  Output

  U - vector with function values    
  
  +====================================================+
*/

void eval_array_chebyshev2(int k, 
			   gsl_vector *x,
			   gsl_vector *U)
{
  int i;
  double u;

  for(i = 0; i < x->size; i++)
    {
      eval_chebyshev2(k, gsl_vector_get(x,i), &u);
      
      gsl_vector_set(U, i, u);
    }
}


/*+====================================================+
  
  FUNCTION svd_solve
  
  Solve a linear system using the SVD decomposition.
  
  Input 

  A - coefficient matrix
  b - coefficient vector
  n - number of columns

  Output

  V - lower triangular matrix of the decmposition
  S - array with singular values
  x - solution
  
  _Remark_

  GSL procedures are used to compute the singular 
  values. As a consequence, the matrix A is replaced 
  by matrix U (upper triangular matrix of the
  decomposition).   
  
  +====================================================+
*/

void svd_solve(gsl_matrix *A, 
	       gsl_vector *b,
	       int n,
	       gsl_matrix *V,
	       gsl_vector *S,
	       gsl_vector *x)
{
  int i;
  double si, div, tol;
  gsl_vector *aux;
  
  tol = 1e-06;

  aux = gsl_vector_alloc(n);

  /* Remark (least squares solution)

  x = pseudoinv(A) * b = (V * S^+ * U^T) * b
     
  1) x = 1 * U^T * b + 0*x 
  2) x = S^+ * x
  3) aux = 1 * V * x + 0*aux
  4) x = aux
  
  */
  
  gsl_linalg_SV_decomp_jacobi(A, V, S);
  
  gsl_blas_dgemv(CblasTrans, 1.0, A, b, 0.0, x);
  
  for(i=0; i < n; i++)
    {
      si = gsl_vector_get(S,i);
      
      if(fabs(si) < tol)
	{
	  gsl_vector_set(x,i,0.0);
	}
      else
	{
	  div = gsl_vector_get(x,i)/si;
	  
	  gsl_vector_set(x, i, div);
	}
    }
  
  gsl_blas_dgemv(CblasNoTrans, 1.0, V, x, 0.0, aux);
  
  gsl_vector_memcpy(x, aux);	
  
  gsl_vector_free(aux);
}


/*+====================================================+
  
  FUNCTION hilbert_workspace_alloc
  
  Workspace allocation to compute the Hilbert transform
  
  Input 

  N - number of samples
  filter - filter type
  step - step size
  workspace - hilbert transform workspace   
  
  +====================================================+
*/

void hilbert_workspace_alloc(hilb_t *workspace,
			     int N,
			     int filter,
			     double step)
{
  int i;
  double w,wc,dw,fun;
  
  workspace->fourier = (fft_t *)malloc(sizeof(fft_t));

  fft_workspace_alloc(workspace->fourier, N, step);
  
  workspace->kernel = gsl_vector_alloc(N);
  workspace->zero   = gsl_vector_alloc(N);
  workspace->aux    = gsl_vector_alloc(N);
  
  workspace->fftREsignal  = gsl_vector_alloc(N);
  workspace->fftIMsignal  = gsl_vector_alloc(N);
  workspace->fftREkernel  = gsl_vector_alloc(N);
  workspace->fftIMkernel  = gsl_vector_alloc(N);
  workspace->fftREproduct = gsl_vector_alloc(N);
  workspace->fftIMproduct = gsl_vector_alloc(N);
  
  wc = 1/(2*step);
  dw = (2*wc)/N;
  
  for(i=0;i<N; i++)
  {
    w   = -wc + i * dw;
    fun = SIGN(w);
    
    gsl_vector_set(workspace->fftREkernel, i, 0);
    gsl_vector_set(workspace->fftIMkernel, i, fun);
    /*nao e sinal negativo de fun, pelo fftshift*/
  }  
}

/*+====================================================+
  
  FUNCTION hilbert_workspace_free
  
  Frees Hilbert transform workspace
  
  Input 

  workspace - hilbert workspace
  
  +====================================================+
*/

void hilbert_workspace_free(hilb_t *workspace)
{
  fft_workspace_free(workspace->fourier);
  
  gsl_vector_free(workspace->zero);
  gsl_vector_free(workspace->kernel);
  gsl_vector_free(workspace->aux);
  
  gsl_vector_free(workspace->fftREsignal);
  gsl_vector_free(workspace->fftIMsignal);
  gsl_vector_free(workspace->fftREkernel);
  gsl_vector_free(workspace->fftIMkernel);
  gsl_vector_free(workspace->fftREproduct);
  gsl_vector_free(workspace->fftIMproduct);

  free(workspace->fourier);
}

/*+====================================================+
  
  FUNCTION hilbert_transform
  
  Compute the Hilbert transform.
  
  Input 
  
  workspace - hilbert transform workspace
  g - signal
    
  Output

  h - hilbert transform  
 
  +====================================================+
*/

void hilbert_transform(gsl_vector *g,
		       gsl_vector *h,
		       hilb_t *workspace)
{
  int N;
  double dt;

  dt = workspace->fourier->step;
  N  = workspace->fourier->N;
    
  gsl_vector_set_all(workspace->zero, 0.0);  
  
  fft_complex(workspace->fftREsignal, 
	      workspace->fftIMsignal,
	      g, 
	      workspace->zero, 
	      workspace->fourier);      

  /*
    fft_complex(workspace->fftREkernel,
    workspace->fftIMkernel,
    workspace->kernel,
    workspace->zero,
    workspace->fourier);
  */

  /* real part product */

  gsl_blas_dcopy(workspace->fftIMsignal, workspace->fftREproduct);
  
  gsl_vector_mul(workspace->fftREproduct, workspace->fftIMkernel);

  gsl_vector_scale(workspace->fftREproduct, -1.0);

  /* imaginary part product */

  gsl_blas_dcopy(workspace->fftREsignal, workspace->fftIMproduct);
  
  gsl_vector_mul(workspace->fftIMproduct, workspace->fftIMkernel);

  /* inverse */

  ifft_complex(h,
	       workspace->zero,
	       workspace->fftREproduct, 
	       workspace->fftIMproduct, 
	       workspace->fourier);  
  
  
  /*
    fft_shift(h, N);  
    gsl_vector_scale(h, dt); 
  */
    
}
