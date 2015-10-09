#include "raft_cuda_aux.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TB_CUDA_AUX 256

/*
Autor: Joao Carlos Cerqueira	email: jc.cerqueira13@gmail.com

Funcoes auxiliares para os programas em GPU utilizando CUDA
*/


extern "C" {
__global__ void gpu_mtx_elementwise_div(float* output, float* input1, float* input2, int dim){
	
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	
	if((k<dim))
		output[k] = input1[k]/input2[k];

	return;

}
}

extern "C" {
__global__ void gpu_mtx_elementwise_mult(float* output, float* input1, float* input2, int dim){
	
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(k<dim)
		output[k] = input1[k]*input2[k];

	return;

}
}

extern "C" {
__global__ void gpu_mtx_elementwise_sum(float* output, float* input1, float* input2, int dim){
	
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(k<dim)
		output[k] = input1[k]+input2[k];

	return;

}
}

extern "C" {
__global__ void gpu_mtx_elementwise_minus(float* output, float* input1, float* input2, int dim){
	
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(k<dim)
		output[k] = input1[k] - input2[k];

	return;

}
}

extern "C" {
void mtx_elementwise_div(float* output, float* input1, float* input2, int dim){
	int TPB = (int)TB_CUDA_AUX;
	int grid;
	grid = (int) ceil(	(float)dim/(float)TPB	);
	gpu_mtx_elementwise_div<<< grid, TPB >>>(output, input1, input2, dim);
}
}

extern "C" {
void mtx_elementwise_mult(float* output, float* input1, float* input2, int dim){
	int TPB = (int)TB_CUDA_AUX;
	int grid;
	grid = (int) ceil(	(float)dim/(float)TPB	);
	gpu_mtx_elementwise_mult<<< grid, TPB >>>(output, input1, input2, dim);
}
}

extern "C" {
void mtx_elementwise_sum(float* output, float* input1, float* input2, int dim){
	int TPB = (int)TB_CUDA_AUX;
	int grid;
	grid = (int) ceil(	(float)dim/(float)TPB	);
	gpu_mtx_elementwise_sum<<< grid, TPB >>>(output, input1, input2, dim);
}
}

