#include "raft_radon_gpu.h"
#include "raft_radon_gpu_function.h"
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define INIC -1.0

/*
Autor: Joao Carlos Cerqueira	email: jc.cerqueira13@gmail.com
*/

extern "C" {
void raft_radon_slantstack_gpu(float* h_output, float* h_input, int sizeImage, int nrays, int nangles)
{
	float *d_output, *d_input;

	// Allocate GPU buffers for the output sinogram
	cudaMalloc(&d_output, sizeof(float) * nrays * nangles);
	
	// Allocate GPU memory for input image and copy
	cudaMalloc(&d_input, sizeof(float) * sizeImage * sizeImage);
	cudaMemcpy(d_input, h_input, sizeof(float) * sizeImage * sizeImage, cudaMemcpyHostToDevice);	

	raft_radon_gpu_function(d_output, d_input, sizeImage, nrays, nangles, 0.0);

	// Copy output vector from GPU buffer to host memory.
	cudaMemcpy(h_output, d_output, sizeof(float) * nrays * nangles, cudaMemcpyDeviceToHost);

	cudaFree(d_input);
	cudaFree(d_output);
	cudaDeviceReset();
}
}
