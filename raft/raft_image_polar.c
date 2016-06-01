#include "raft_image_polar.h"
#include <pthread.h>
#include <math.h>

#define MAX_NTHREADS 20000
#define PI 3.1415926535897932384626433832795
#define MAX( x, y ) ( ( ( x ) > ( y ) ) ? ( x ) : ( y ) )
#define SIGN( x ) ( ( ( x ) > 0.0 ) ? 1.0 : ( ( ( x ) < 0.0 ) ? -1.0 : 0.0 ) )

#define LOGEPS 0.01

typedef struct{

  double x0,y0,dx,dy;
  double theta, dth;
  double r, r0, rf, dr; 

}mesh_t;

typedef struct{

  raft_image cartesian;
  raft_image polar;
  int size, views, rays;
  int nthread;
  int colIndex[2];  
  mesh_t mesh;
  
}param_t;

void *polar_inverse_loop(void *t);

void *logpolar_inverse_loop(void *t);

void *polar_direct_loop(void *t);

void *logpolar_direct_loop(void *t);

void *sinogram_direct_loop(void *t);

void *sinogramlp_direct_loop(void *t);

/*------------------------------------------------------*/

void raft_image_c2p(raft_image cartesian,
		    raft_image polar,
		    int nthreads)
{
  pthread_t thread[MAX_NTHREADS];
  pthread_attr_t attr;
  int R, V, size, e, n, rc;    
  param_t param[MAX_NTHREADS];  
  mesh_t mesh;

  void *status;

  size = raft_matrix_nlines(cartesian.data);
  V    = raft_matrix_ncolumns(polar.data);
  R    = raft_matrix_nlines(polar.data);

  //Initialize mesh;

  mesh.x0 = -0.707106781;
  mesh.y0 = -0.707106781; 
  mesh.dx = (double) (-2*mesh.x0)/(size-1); 
  mesh.dy = mesh.dx;

  mesh.r0 = 0;
  mesh.rf = 1;
  mesh.dr = (double) (mesh.rf-mesh.r0)/(R-1);
  mesh.dth = (double) (2*PI)/(V-1);
      
  //Initialize and set thread detached attribute
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  e = (int) floor((V*R)/nthreads);
      
  for(n = 0; n < nthreads; n++) 
    {
      param[n].mesh  = mesh;
      param[n].size  = size;
      param[n].views = V;
      param[n].rays  = R;
      
      param[n].cartesian = cartesian;
      param[n].polar     = polar; 
      param[n].nthread   = n;      
      
      param[n].colIndex[0] = e * n;
      param[n].colIndex[1] = (n+1) * e - 1;     
      
      rc = pthread_create(&thread[n], 
			  &attr, 
			  polar_direct_loop, 
			  (void *)&param[n]);
    }

  // Free attribute and wait for the other threads
  pthread_attr_destroy(&attr);
  for(n = 0; n < nthreads; n++) 
    {
      rc = pthread_join(thread[n], &status);
    }  
}


/*-------------------------------*/

void raft_image_p2c(raft_image polar,
		    raft_image cartesian,
		    int nthreads)
{
  pthread_t thread[MAX_NTHREADS];
  pthread_attr_t attr;
  int R, V, size, e, n, rc;    
  param_t param[MAX_NTHREADS];  
  mesh_t mesh;

  void *status;

  size = raft_matrix_nlines(cartesian.data);
  V    = raft_matrix_ncolumns(polar.data);
  R    = raft_matrix_nlines(polar.data);

  //Initialize mesh;
  
  mesh.x0 = -1.0; //-0.707106781;
  mesh.y0 = -1.0; //-0.707106781; 
  mesh.dx = (double) (-2*mesh.x0)/(size-1); 
  mesh.dy = mesh.dx;

  mesh.r0 = 0;
  mesh.rf = 1;
  mesh.dr = (double) (mesh.rf-mesh.r0)/(R-1);
  mesh.dth = (double) (2*PI)/(V-1);
      
  //Initialize and set thread detached attribute
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  e = (int) floor((size*size)/nthreads);
      
  for(n = 0; n < nthreads; n++) 
    {
      param[n].mesh  = mesh;
      param[n].size  = size;
      param[n].views = V;
      param[n].rays  = R;
      
      param[n].cartesian = cartesian;
      param[n].polar     = polar; 
      param[n].nthread   = n;      
      
      param[n].colIndex[0] = e * n;
      param[n].colIndex[1] = (n+1) * e - 1;     
      
      rc = pthread_create(&thread[n], 
			  &attr, 
			  polar_inverse_loop, 
			  (void *)&param[n]);
    }

  // Free attribute and wait for the other threads
  pthread_attr_destroy(&attr);
  for(n = 0; n < nthreads; n++) 
    {
      rc = pthread_join(thread[n], &status);
    }    
}

/*-------------------------------*/

void raft_image_lp2c(raft_image logpolar,
		     raft_image cartesian,
		     int nthreads)
{
  pthread_t thread[MAX_NTHREADS];
  pthread_attr_t attr;
  int R, V, size, e, n, rc;    
  param_t param[MAX_NTHREADS];  
  mesh_t mesh;

  void *status;

  size = raft_matrix_nlines(cartesian.data);
  V    = raft_matrix_ncolumns(logpolar.data);
  R    = raft_matrix_nlines(logpolar.data);
  
  //Initialize mesh;
  
  mesh.x0 = -0.707106781;
  mesh.y0 = -0.707106781; 
  mesh.dx = (double) (-2*mesh.x0)/(size-1); 
  mesh.dy = mesh.dx;

  mesh.r0 = log(LOGEPS);
  mesh.rf = 0;
  mesh.dr = (double) (mesh.rf-mesh.r0)/(R-1);
  mesh.dth = (double) (2*PI)/(V-1);
      
  //Initialize and set thread detached attribute
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  e = (int) floor((size*size)/nthreads);
      
  for(n = 0; n < nthreads; n++) 
    {
      param[n].mesh  = mesh;
      param[n].size  = size;
      param[n].views = V;
      param[n].rays  = R;
      
      param[n].cartesian = cartesian;
      param[n].polar     = logpolar; 
      param[n].nthread   = n;      
      
      param[n].colIndex[0] = e * n;
      param[n].colIndex[1] = (n+1) * e - 1;     
      
      rc = pthread_create(&thread[n], 
			  &attr, 
			  logpolar_inverse_loop, 
			  (void *)&param[n]);
    }

  // Free attribute and wait for the other threads
  pthread_attr_destroy(&attr);
  for(n = 0; n < nthreads; n++) 
    {
      rc = pthread_join(thread[n], &status);
    }    
}

/*--------------------------------*/

void raft_image_c2lp(raft_image cartesian,
		     raft_image logpolar,
		     int nthreads)
{
  pthread_t thread[MAX_NTHREADS];
  pthread_attr_t attr;
  int R, V, size, e, n, rc;    
  param_t param[MAX_NTHREADS];  
  mesh_t mesh;
  double eps;

  void *status;

  size = raft_matrix_nlines(cartesian.data);
  V    = raft_matrix_ncolumns(logpolar.data);
  R    = raft_matrix_nlines(logpolar.data);

  //Initialize mesh;

  mesh.x0 = -0.707106781;
  mesh.y0 = -0.707106781; 
  mesh.dx = (double) (-2*mesh.x0)/(size-1); 
  mesh.dy = mesh.dx;
  
  mesh.r0 = log(LOGEPS);
  mesh.rf = 0;
  mesh.dr = (double) (mesh.rf-mesh.r0)/(R-1);
  mesh.dth = (double) (2*PI)/(V-1);
      
  //Initialize and set thread detached attribute
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  e = (int) floor((V*R)/nthreads);
      
  for(n = 0; n < nthreads; n++) 
    {
      param[n].mesh  = mesh;
      param[n].size  = size;
      param[n].views = V;
      param[n].rays  = R;
      
      param[n].cartesian = cartesian;
      param[n].polar     = logpolar; 
      param[n].nthread   = n;      
      
      param[n].colIndex[0] = e * n;
      param[n].colIndex[1] = (n+1) * e - 1;     
      
      rc = pthread_create(&thread[n], 
			  &attr, 
			  logpolar_direct_loop, 
			  (void *)&param[n]);
    }

  // Free attribute and wait for the other threads
  pthread_attr_destroy(&attr);
  for(n = 0; n < nthreads; n++) 
    {
      rc = pthread_join(thread[n], &status);
    }    
}

/*------------------------------------------------------*/

void raft_image_s2p(raft_image sinogram,
		    raft_image polar,
		    int nthreads)
{
  pthread_t thread[MAX_NTHREADS];
  pthread_attr_t attr;
  int R, V, size, e, n, rc;    
  param_t param[MAX_NTHREADS];  
  mesh_t mesh;

  void *status;

  V    = raft_matrix_ncolumns(sinogram.data)*2;
  R    = raft_matrix_nlines(sinogram.data)/2;

  //Initialize and set thread detached attribute
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  e = (int) floor((V*R)/nthreads);
      
  for(n = 0; n < nthreads; n++) 
    {
      param[n].mesh  = mesh;
      param[n].size  = 0;
      param[n].views = V;
      param[n].rays  = R;
      
      param[n].cartesian = sinogram;
      param[n].polar     = polar; 
      param[n].nthread   = n;      
      
      param[n].colIndex[0] = e * n;
      param[n].colIndex[1] = (n+1) * e - 1;     
      
      rc = pthread_create(&thread[n], 
			  &attr, 
			  sinogram_direct_loop, 
			  (void *)&param[n]);
    }

  // Free attribute and wait for the other threads
  pthread_attr_destroy(&attr);
  for(n = 0; n < nthreads; n++) 
    {
      rc = pthread_join(thread[n], &status);
    }    
}

/*------------------------------------------------------*/

void raft_image_s2lp(raft_image sinogram,
		     raft_image logpolar,
		     int nthreads)
{
  pthread_t thread[MAX_NTHREADS];
  pthread_attr_t attr;
  int R, V, size, e, n, rc;    
  param_t param[MAX_NTHREADS];  
  mesh_t mesh;

  void *status;

  V    = raft_matrix_ncolumns(sinogram.data)*2;
  R    = raft_matrix_nlines(sinogram.data)/2;

   //Initialize mesh;

  mesh.dx = 2.0/(2*R-1);
  
  mesh.r0 = log(LOGEPS);
  mesh.rf = 0;
  mesh.dr = (double) (mesh.rf-mesh.r0)/(R-1);
  mesh.dth = (double) (2*PI)/(V-1);

  //Initialize and set thread detached attribute
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  e = (int) floor((V*R)/nthreads);
      
  for(n = 0; n < nthreads; n++) 
    {
      param[n].mesh  = mesh;
      param[n].size  = 0;
      param[n].views = V;
      param[n].rays  = R;
      
      param[n].cartesian = sinogram;
      param[n].polar     = logpolar; 
      param[n].nthread   = n;      
      
      param[n].colIndex[0] = e * n;
      param[n].colIndex[1] = (n+1) * e - 1;     
      
      rc = pthread_create(&thread[n], 
			  &attr, 
			  sinogramlp_direct_loop, 
			  (void *)&param[n]);
    }

  // Free attribute and wait for the other threads
  pthread_attr_destroy(&attr);
  for(n = 0; n < nthreads; n++) 
    {
      rc = pthread_join(thread[n], &status);
    }    
}


/*--------------------------------------------------------*/

void *polar_inverse_loop(void *t)
{
  int w, n, R, V;
  param_t *param;
  
  param = (param_t *)t;
  
  n = param->size;
  R = param->rays;
  V = param->views;
  
  for(w = param->colIndex[0]; w <= param->colIndex[1]; w++)
    {
      int i, j, q, m, k;
      double aS, S, N, Z, d[4]={1,-1,1,-1};
      double x, y, A[4]={0,PI,PI,2*PI};
      double r, th, rrr, tth, C[2], a, b;
            
      i = w/n;
      j = w%n;
      
      x = param->mesh.x0 + i*param->mesh.dx;  
      y = param->mesh.y0 + j*param->mesh.dy;  

      // find (x,y) quadrant
      
      S = SIGN(x*y);
      aS = fabs(S);
      
      N = SIGN(-x*S) + (3-S)/2;
      
      Z = (1-SIGN(x-0.001))*(1-SIGN(y))/2;
      
      q = (int) (aS*N + (1-aS)*Z);

      r  = sqrt(x*x + y*y);
	
      th = A[q] + d[q] * atan(fabs(y/x)); 		   

      m = (int) floor((r-param->mesh.r0)/param->mesh.dr);
      k = (int) floor((th)/param->mesh.dth);	
      
      rrr = param->mesh.r0 + m * param->mesh.dr;
      tth = k*param->mesh.dth;
      
      if(k < V && k>-1 && m < R && m>-1)
	{
	  // Nearest
	  // raft_image_sample(param->cartesian, j, i) = raft_image_sample(param->polar, k, m);
	  
	  // Bilinear
	  
	  C[0] = raft_image_sample(param->polar, m+1, k+1);
	  C[1] = raft_image_sample(param->polar, m, k+1);

	  a = (C[0]-C[1])*(th-tth)/param->mesh.dth + C[1];

	  C[0] = raft_image_sample(param->polar, m, k+1);
	  C[1] = raft_image_sample(param->polar, m, k);

	  b = (C[0]-C[1])*(th-tth)/param->mesh.dth + C[1];	
	  
	  raft_image_sample(param->cartesian, j, i) = (a-b)*(r-rrr)/param->mesh.dr + b;
	}
      else
	raft_image_sample(param->cartesian, j, i) = 0; 
    }

  pthread_exit(NULL);

}

/*---------*/

void *logpolar_inverse_loop(void *t)
{
  int w, n, R, V;
  param_t *param;
  
  param = (param_t *)t;
  
  n = param->size;
  R = param->rays;
  V = param->views;
  
  for(w = param->colIndex[0]; w <= param->colIndex[1]; w++)
    {
      int i, j, q, m, k;
      double aS, S, N, Z, d[4]={1,-1,1,-1};
      double x, y, A[4]={0,PI,PI,2*PI};
      double r, th, rrr, tth, C[2], a, b;
            
      i = w/n;
      j = w%n;
      
      x = param->mesh.x0 + i*param->mesh.dx;  
      y = param->mesh.y0 + j*param->mesh.dy;  

      // find (x,y) quadrant
      
      S = SIGN(x*y);
      aS = fabs(S);
      
      N = SIGN(-x*S) + (3-S)/2;
      
      Z = (1-SIGN(x-0.001))*(1-SIGN(y))/2;
      
      q = (int) (aS*N + (1-aS)*Z);

      r  = log(sqrt(x*x + y*y));
	
      th = A[q] + d[q] * atan(fabs(y/x)); 		   

      m = (int) floor((r-param->mesh.r0)/param->mesh.dr);
      k = (int) floor((th)/param->mesh.dth);	
      
      rrr = param->mesh.r0 + m * param->mesh.dr;
      tth = k*param->mesh.dth;
      
      if(k < V && k>-1 && m < R && m>-1)
	{
	  // Nearest
	  // raft_image_sample(param->cartesian, j, i) = raft_image_sample(param->polar, k, m);
	  
	  // Bilinear
	  
	  C[0] = raft_image_sample(param->polar, m+1, k+1);
	  C[1] = raft_image_sample(param->polar, m, k+1);

	  a = (C[0]-C[1])*(th-tth)/param->mesh.dth + C[1];

	  C[0] = raft_image_sample(param->polar, m, k+1);
	  C[1] = raft_image_sample(param->polar, m, k);

	  b = (C[0]-C[1])*(th-tth)/param->mesh.dth + C[1];	
	  
	  raft_image_sample(param->cartesian, j, i) = (a-b)*(r-rrr)/param->mesh.dr + b;
	}
      else
	raft_image_sample(param->cartesian, j, i) = 0; 
    }

  pthread_exit(NULL);

}

/*---------*/

void *polar_direct_loop(void *t)
{
  int w, n, R;
  param_t *param;
  
  param = (param_t *)t;
  
  n = param->size;
  R = param->rays;
  
  for(w = param->colIndex[0]; w <= param->colIndex[1]; w++)
    {
      int i,j,x,y;
      double xx,yy,xxx,yyy;
      double theta, r;
      double A,B,C[2];
      
      j = w/R;
      i = w%R;

      theta = j*param->mesh.dth;
	  
      r= i*param->mesh.dr;

      xx = r*cos(theta);
      yy = r*sin(theta);
      
      x = (int) floor((xx-param->mesh.x0)/param->mesh.dx);
      y = (int) floor((yy-param->mesh.y0)/param->mesh.dy);	 	
      
      xxx = param->mesh.x0 + x*param->mesh.dx;
      yyy = param->mesh.y0 + y*param->mesh.dy;
      
      if (x<n && x>-1 && y<n && y>-1)
	{
	  // Nearest
	  // raft_image_sample(polar, j, i) = raft_image_sample(cartesian, y, x);
	  
	  // Bilinear
	  
	  C[0] = raft_image_sample(param->cartesian, y+1, x+1);
	  C[1] = raft_image_sample(param->cartesian, y, x+1);
	  
	  A = (C[0]-C[1])*(yy-yyy)/param->mesh.dy + C[1];
	  
	  C[0] = raft_image_sample(param->cartesian, y+1, x);
	  C[1] = raft_image_sample(param->cartesian, y, x);
	  
	  B = (C[0]-C[1])*(yy-yyy)/param->mesh.dy + C[1];	
	  
	  raft_image_sample(param->polar, i, j) = (A-B)*(xx-xxx)/param->mesh.dx + B;
	}
      else
	raft_image_sample(param->polar, i, j) = 0;
      
    }
  
  pthread_exit(NULL);
}

/*---------*/

void *logpolar_direct_loop(void *t)
{
  int w, n, R, V;
  param_t *param;
  
  param = (param_t *)t;
  
  n = param->size;
  R = param->rays;
  V = param->views;
  
  for(w = param->colIndex[0]; w <= param->colIndex[1]; w++)
    {
      int i,j,x,y;
      double xx,yy,xxx,yyy;
      double theta, r;
      double A,B,C[2];
      
      j = w/R;
      i = w%R;

      theta = j*param->mesh.dth;
	  
      r= param->mesh.r0 + i*param->mesh.dr;

      xx = exp(r)*cos(theta);
      yy = exp(r)*sin(theta);
      
      x = (int) floor((xx-param->mesh.x0)/param->mesh.dx);
      y = (int) floor((yy-param->mesh.y0)/param->mesh.dy);	 	
      
      xxx = param->mesh.x0 + x*param->mesh.dx;
      yyy = param->mesh.y0 + y*param->mesh.dy;
      
      if (x<n && x>-1 && y<n && y>-1)
	{
	  // Nearest
	  // raft_image_sample(polar, j, i) = raft_image_sample(cartesian, y, x);
	  
	  // Bilinear
	  
	  C[0] = raft_image_sample(param->cartesian, y+1, x+1);
	  C[1] = raft_image_sample(param->cartesian, y, x+1);
	  
	  A = (C[0]-C[1])*(yy-yyy)/param->mesh.dy + C[1];
	  
	  C[0] = raft_image_sample(param->cartesian, y+1, x);
	  C[1] = raft_image_sample(param->cartesian, y, x);
	  
	  B = (C[0]-C[1])*(yy-yyy)/param->mesh.dy + C[1];	
	  
	  raft_image_sample(param->polar, i, j) = (A-B)*(xx-xxx)/param->mesh.dx + B;
	}
      else
	raft_image_sample(param->polar, i, j) = 0;
      
    }
  
  pthread_exit(NULL);
}

/*---------*/

void *sinogram_direct_loop(void *t)
{
  int w, R, V;
  param_t *param;
  
  param = (param_t *)t;
  
  R = param->rays;
  V = param->views;
  
  for(w = param->colIndex[0]; w <= param->colIndex[1]; w++)
    {
      int i,j;
            
      j = w/R;
      i = w%R;
      
      if( j < V/2 )
	{
	  raft_image_sample(param->polar, i, j) = raft_image_sample(param->cartesian, R + i, j);
	}
      else
	{
	  raft_image_sample(param->polar, i, j) = raft_image_sample(param->cartesian, R + 1 - i, j - V);
	}
    }

  pthread_exit(NULL);
}

/*---------*/

void *sinogramlp_direct_loop(void *t)
{
  int w, R, V, m;
  param_t *param;
  double z;

  param = (param_t *)t;
  
  R = param->rays;
  V = param->views;
  
  for(w = param->colIndex[0]; w <= param->colIndex[1]; w++)
    {
      int i,j;
            
      j = w/R;
      i = w%R;
      
      z = exp(param->mesh.r0 + i*param->mesh.dr);

      m = (int) floor((z-0)/param->mesh.dx);   //porque dx? malandragem da um tempo ...

      if( j < V/2 )
	{
	  raft_image_sample(param->polar, i, j) = raft_image_sample(param->cartesian, R + m, j);
	}
      else
	{
	  raft_image_sample(param->polar, i, j) = raft_image_sample(param->cartesian, R - m, j - V);
	}
    }
  
  pthread_exit(NULL);
}

