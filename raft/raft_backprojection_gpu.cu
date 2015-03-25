#include "raft_backprojection_gpu.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159265358979
#define INIC -1.0

#define TPBXb 32
#define TPBYb 32

texture<float, cudaTextureType2D, cudaReadModeElementType> texSino;

#define CUDA_CHECK_RETURN(value) {											\
	cudaError_t _m_cudaStat = value;										\
	if (_m_cudaStat != cudaSuccess) {										\
		fprintf(stderr, "Error %s at line %d in file %s\n",					\
				cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);		\
		exit(1);															\
	} }

__global__ void backp_kernel(float *image, int wdI, int nrays, int nangles, float delta, float dt, float dth)
{
  int i, j, T;
  float t, cs1, cs2, cs3, cs4, k;
  float	x, y;
  float cosk, sink;

  i = 2*(blockDim.x * blockIdx.x + threadIdx.x);
  j = 2*(blockDim.y * blockIdx.y + threadIdx.y);

  if ( ((i+1)<wdI) && ((j+1) < wdI) ){

	  cs1 = 0;
	  cs2 = 0;
	  cs3 = 0;
	  cs4 = 0;

	  for(k=0; k < (nangles-1); k++)
	  {
		  sincosf(k * dth, &sink, &cosk);
		  ///////////////////////////
		  x = (float)INIC + i * delta;
		  y = (float)INIC + j * delta;
	      t = x*cosk + y*sink;
	      T = (float)((t + 1)/dt);
	      //if(T >= 0 && T <= (nrays-1))
		  	  cs1 = cs1 + tex2D(texSino, k + 0.5f, T + 0.5f);
		  //////////////////////////
		  x = (float)INIC + (i+1) * delta;
		  y = (float)INIC + j * delta;
	      t = x*cosk + y*sink;
	      T = (float)((t + 1)/dt);
	      //if(T >= 0 && T <= (nrays-1))
		  	  cs2 = cs2 + tex2D(texSino, k + 0.5f, T + 0.5f);
		  //////////////////////////
		  x = (float)INIC + i * delta;
		  y = (float)INIC + (j+1) * delta;
	      t = x*cosk + y*sink;
	      T = (float)((t + 1)/dt);
	      //if(T >= 0 && T <= (nrays-1))
		  	  cs3 = cs3 + tex2D(texSino, k + 0.5f, T + 0.5f);
		  //////////////////////////
		  x = (float)INIC + (i+1) * delta;
		  y = (float)INIC + (j+1) * delta;
	      t = x*cosk + y*sink;
	      T = (float)((t + 1)/dt);
	      //if(T >= 0 && T <= (nrays-1))
		  	  cs4 = cs4 + tex2D(texSino, k + 0.5f, T + 0.5f);
	  }

	  image[(j)*wdI + (wdI-1-i)] 			= (cs1*dth);
	  image[(j)*wdI + (wdI-1-i-1)] 		= (cs2*dth);
	  image[(j+1)*wdI + (wdI-1-i)]		= (cs3*dth);
	  image[(j+1)*wdI + (wdI-1-i-1)] 		= (cs4*dth);
  }
}




extern "C" {
void raft_backprojection_slantstack_gpu(float *image, float *sino, int sizeImage, int nrays, int nangles){
	float* d_output;
	int sizeSino = nrays * nangles * sizeof(float);
	float dt  = 2.0/(nrays-1);
	float dth = PI/(nangles-1);
	float delta = (float) 2*fabsf(INIC)/(sizeImage-1);

	/////// KERNEL EXECUTION TIME TEST
			cudaEvent_t start, stop;
			cudaEventCreate(&start);
			cudaEventCreate(&stop);
			float milliseconds = 0;
			cudaEventRecord(start);
	//////////////////////////////////

    // Allocate CUDA array in device memory (sinogram matrix)
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0,cudaChannelFormatKindFloat);
    cudaArray* cuArray;
    cudaMallocArray(&cuArray, &channelDesc, nangles, nrays);

    // Copy to device memory the sinogram matrix
    cudaMemcpyToArray(cuArray, 0, 0, sino, sizeSino , cudaMemcpyHostToDevice);

    // Set texture parameters
    texSino.addressMode[0] = cudaAddressModeBorder;
    texSino.addressMode[1] = cudaAddressModeBorder;
    texSino.filterMode     = cudaFilterModeLinear;
    /*texSino.normalized     = true; */

    // Bind the array to the texture reference
    cudaBindTextureToArray(texSino, cuArray, channelDesc);

    // Allocate GPU memory for the output image
    cudaMalloc(&d_output, sizeof(float) * sizeImage *sizeImage);

    //GRID and BLOCKS SIZE
    dim3 threadsPerBlock(TPBXb,TPBYb);
    dim3 grid((sizeImage/threadsPerBlock.x)/2 + 1, (sizeImage/threadsPerBlock.y)/2 + 1);

    //KERNEL EXECUTION
    backp_kernel<<<grid, threadsPerBlock>>>(d_output, sizeImage, nrays, nangles, delta, dt, dth);
    cudaDeviceSynchronize();

	/////// PRINT KERNEL EXECUTION TIME
			cudaEventRecord(stop);
			cudaEventSynchronize(stop);
			cudaEventElapsedTime(&milliseconds, start, stop);
			//fprintf(stderr, "backp: %f ms\n", milliseconds);
	///////////////////////////////////

	//Copy the output image from device memory to host memory
	cudaMemcpy (image , d_output , sizeImage*sizeImage*sizeof(float) , cudaMemcpyDeviceToHost);

    cudaUnbindTexture(texSino);
    cudaFreeArray(cuArray);
    cudaFree(d_output);
    cudaDeviceReset();
}
}
