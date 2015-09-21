#include "raft_backprojection_gpu_function.h"
#include "raft_radon_gpu_function.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define THREADS_BLOCK_EM 256

extern "C" {
__global__ void tr_gpu_mtx_elementwise_div(float* output, float* input1, float* input2, int dim){
	
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	
	if((k<dim))
		output[k] = input1[k]/input2[k];

	return;

}
}

extern "C" {
__global__ void tr_gpu_mtx_elementwise_mult(float* output, float* input1, float* input2, int dim){
	
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(k<dim)
		output[k] = input1[k]*input2[k];

	return;

}
}

extern "C" {
__global__ void tr_gpu_mtx_elementwise_sum(float* output, float* input1, float* input2, int dim){
	
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(k<dim)
		output[k] = input1[k]+input2[k];

	return;

}
}

extern "C" {
void tr_mtx_elementwise_div(float* output, float* input1, float* input2, int dim){
	int TPB = (int)THREADS_BLOCK_EM;
	int grid;
	grid = (int) ceil(	(float)dim/(float)TPB	);
	tr_gpu_mtx_elementwise_div<<< grid, TPB >>>(output, input1, input2, dim);
}
}

extern "C" {
void tr_mtx_elementwise_mult(float* output, float* input1, float* input2, int dim){
	int TPB = (int)THREADS_BLOCK_EM;
	int grid;
	grid = (int) ceil(	(float)dim/(float)TPB	);
	tr_gpu_mtx_elementwise_mult<<< grid, TPB >>>(output, input1, input2, dim);
}
}

extern "C" {
void tr_mtx_elementwise_sum(float* output, float* input1, float* input2, int dim){
	int TPB = (int)THREADS_BLOCK_EM;
	int grid;
	grid = (int) ceil(	(float)dim/(float)TPB	);
	tr_gpu_mtx_elementwise_sum<<< grid, TPB >>>(output, input1, input2, dim);
}
}



extern "C"{
void raft_tr_em_gpu(float *output, float *rawsino, float *flat, float *dark, int sizeImage, int nrays, int nangles, int niter)
{
	int i;	

	int SinoMem = sizeof(float)*nrays*nangles;
	int ImageMem = sizeof(float)*sizeImage*sizeImage;	

	float *d_image1, *d_image2, *d_image3, *d_flat, *d_dark, *d_sino1, *d_sino2, *d_sino3, *d_Sino;
	float *img_inicial;


// CUDA LINEAR MEMORY
	cudaMalloc(&d_sino1, SinoMem);
	cudaMalloc(&d_sino2, SinoMem);
	cudaMalloc(&d_sino3, SinoMem);	
	cudaMalloc(&d_Sino, SinoMem);	
	cudaMalloc(&d_flat, SinoMem);
	cudaMalloc(&d_dark, SinoMem);
	cudaMalloc(&d_image1, ImageMem);
	cudaMalloc(&d_image2, ImageMem);
	cudaMalloc(&d_image3, ImageMem);


// INITIAL IMAGE
	img_inicial = (float *)malloc(ImageMem);
	for (i = 0; i < (sizeImage*sizeImage); i++)		
		img_inicial[i] = 1.0;

	cudaMemcpy(d_image1, img_inicial, ImageMem, cudaMemcpyHostToDevice);
	free(img_inicial);

// COPIES
	cudaMemcpy(d_flat, flat, SinoMem, cudaMemcpyHostToDevice);
	cudaMemcpy(d_dark, dark, SinoMem, cudaMemcpyHostToDevice);
	cudaMemcpy(d_Sino, rawsino, SinoMem, cudaMemcpyHostToDevice);

	raft_backprojection_gpu_function(d_image3, d_Sino, sizeImage, nrays, nangles);

	for(i=0; i<niter; i++){
		raft_exp_radon_gpu_function(d_sino1, d_image1, sizeImage, nrays, nangles);
		tr_mtx_elementwise_mult(d_sino1, d_sino1, d_flat, nrays*nangles);
		//tr_mtx_elementwise_mult(d_sino2, d_sino1, d_Sino, nrays*nangles);
		//tr_mtx_elementwise_sum(d_sino3, d_sino1, d_dark, nrays*nangles);
		//tr_mtx_elementwise_div(d_sino3, d_sino2, d_sino3, nrays*nangles);
		
		raft_backprojection_gpu_function(d_image2, d_sino1, sizeImage, nrays, nangles);
		tr_mtx_elementwise_div(d_image2, d_image2, d_image3, sizeImage*sizeImage);
		tr_mtx_elementwise_mult(d_image1, d_image1, d_image2, sizeImage*sizeImage);
	}
	


	cudaMemcpy(output , d_image1 , ImageMem , cudaMemcpyDeviceToHost);


}// END VOID
}
