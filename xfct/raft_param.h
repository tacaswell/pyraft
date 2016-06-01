#ifndef _RAFT_PARAM_H_
#define _RAFT_PARAM_H_

#ifndef MAX 
#define MAX(a,b) ((a) > (b) ? (a):(b))
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a):(b))
#endif

#ifndef SIGN
#define SIGN(x) ((x) >= 0 ? 1 : -1) 
#endif 

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

#ifndef AVRG
#define AVRG(x,y) (((x)+(y))/2)
#endif

#ifndef CVXCOMB
#define CVXCOMB(a,x,b,y) ((a)*(x)+(b)*(y))
#endif

#ifndef REAL
#define REAL(z,i) ((z)[2*(i)])
#endif

#ifndef IMAG
#define IMAG(z,i) ((z)[2*(i)+1])
#endif

#ifndef ZERO
#define ZERO 1e-10
#endif

#ifndef MPI
#define MPI 3.14159265358979323846
#endif

#ifndef XYMIN 
#define XYMIN -0.70710678
#endif

#ifndef XYMAX 
#define XYMAX 0.70710678
#endif

/*######################################################
  Title: RAFT parameters
  ####################################################*/

/*+====================================================+
  
  CONSTANT: enum

  RAFT parameters.

  RAFT_RAMLAK -	Ram-Lak filter for FBP algorithm
  RAFT_COSINE -	Cosine filter for FBP algorithm
  RAFT_SHEPP -	Sheep filer for FBP algorithm
  RAFT_EMCONV -	Convex EM transmission algorithm
  RAFT_EMGRAD -	Gradient type transmission EM algorithm
  RAFT_EMSTD -	Standard transmission EM algorithm
  RAFT_KUNSTD -	Standard acceleration for Kunyansky algorithm (i.e., none!)
  RAFT_KUNGSJB - GSJB acceleration for Kunyansky algorithm
  RAFT_MONO -	Define monochromatic projections
  RAFT_POLY -	Define polychromatic projections
  RAFT_REG -	Regular mesh of rays (equispaced)
  RAFT_NREG -	Non regular mesh of rays (nonequispaced)
  RAFT_CT -	Define CT scan geometry
  RAFT_SPECT -	Define SPECT scan geometry
  RAFT_XFCT -	Define FLUOR scan geometry
  RAFT_YKTBND - Compute bounds for Kunyansky's method
  RAFT_NKTBND - Do not compute bounds for Kunyansky's method
  
  +====================================================+
*/

enum{ 
  WHITE = 1  
};

enum{
  RAFT_RAMLAK = 0,
  RAFT_COSINE = 1,
  RAFT_SHEPP = 2
};

enum{
  RAFT_EMCONV = 0,
  RAFT_EMGRAD = 1,
  RAFT_EMSTD = 2,
  RAFT_KUNSTD = 3,
  RAFT_KUNGSJB = 4,  
  RAFT_YKTBND = 5,
  RAFT_NKTBND = 6
};

enum{
  RAFT_MONO = 1,
  RAFT_POLY = 0
};


enum{
  RAFT_REG  = 0,
  RAFT_NREG = 1,
  RAFT_CHEB = 2,
  RAFT_CT   = 3,
  RAFT_SPECT= 4,
  RAFT_XFCT = 5
};



#endif
