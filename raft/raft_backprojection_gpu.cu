#include "raft_backprojection_gpu.h"
#include "raft_backprojection_gpu_function.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159265358979
#define INIC -1.0

/*
Autor: Joao Carlos Cerqueira	email: jc.cerqueira13@gmail.com
*/

extern "C" {
void raft_backprojection_slantstack_gpu(float *image, float *sino, int sizeImage, int nrays, int nangles){

	float *d_output, *d_input;

	// Allocate GPU memory for the output image
	cudaMalloc(&d_output, sizeof(float) * sizeImage *sizeImage);

	// Allocate GPU memory for input image and copy
	cudaMalloc(&d_input, sizeof(float) * nrays * nangles);
	cudaMemcpy(d_input, sino, sizeof(float) * nrays * nangles, cudaMemcpyHostToDevice);	

	//KERNEL EXECUTION
	//raft_backprojection_gpu_kernel<<<grid, threadsPerBlock>>>(d_output, sizeImage, nrays, nangles, delta, dt, dth);
	raft_backprojection_gpu_function(d_output, d_input, sizeImage, nrays, nangles);

	//Copy the output image from device memory to host memory
	cudaMemcpy (image , d_output , sizeImage*sizeImage*sizeof(float) , cudaMemcpyDeviceToHost);

	cudaFree(d_output);
	cudaFree(d_input);
	cudaDeviceReset();
}
}
