#include "raft_cuda_aux.h"
#include "raft_backprojection_gpu_function.h"
#include "raft_radon_gpu_function.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define THREADS_BLOCK_EM 256

/*Autor: Joao Carlos Cerqueira	email: jc.cerqueira13@gmail.com

A funcao 'void raft_em_gpu' implementa o metodo de reconstrucao iterativo
EM, onde 'niter' eh o numero de iteracoes. Este metodo de reconstrucao é
limitado a tomografias baseadas em emissao. Para inversoes de tomografias
baseadas em transmissao, utilize a funcao 'void raft_tr_em_gpu'
*/


extern "C" {
void raft_em_gpu(float *output, float *sino, int sizeImage, int nrays, int nangles, int niter){
	int i, nit;	
	
	int SinoMem = sizeof(float)*nrays*nangles;
	int ImageMem = sizeof(float)*sizeImage*sizeImage;	

	float *d_image1, *d_image2, *d_image3, *d_sino1, *d_sino2, *d_Sino;
	float *h_sino_ones, *d_bp_ones;
	float *img_inicial;
	
	
	//	SET CUDA LINEAR MEMORY
	cudaMalloc(&d_sino1, SinoMem);
	cudaMalloc(&d_sino2, SinoMem);	
	cudaMalloc(&d_Sino, SinoMem);	
	cudaMalloc(&d_image1, ImageMem);
	cudaMalloc(&d_image2, ImageMem);
	cudaMalloc(&d_image3, ImageMem);
	cudaMalloc(&d_bp_ones, ImageMem);

	
	//	Imagem inicial (provisório)
	img_inicial	= (float *)malloc(ImageMem);
	for (i = 0; i < (sizeImage*sizeImage); i++)		
		img_inicial[i] = 1.0;

	cudaMemcpy(d_image1, img_inicial, ImageMem, cudaMemcpyHostToDevice);

	// Sinogram of ones
	h_sino_ones = (float *)malloc(SinoMem);
	for(i = 0; i < (nrays*nangles); i++)
		h_sino_ones[i] = 1;

	cudaMemcpy(d_sino1, h_sino_ones, SinoMem, cudaMemcpyHostToDevice);
	


	//	BP OF ONES and INPUT SINOGRAM COPY
	raft_backprojection_gpu_function(d_bp_ones, d_sino1, sizeImage, nrays, nangles);

	cudaMemcpy(d_Sino, sino, SinoMem, cudaMemcpyHostToDevice);	

	//	THE FOR LOOP
	for(nit = 0; nit < niter; nit++){
		raft_radon_gpu_function(d_sino1, d_image1, sizeImage, nrays, nangles, 0.0);
		mtx_elementwise_div(d_sino2, d_Sino, d_sino1, nrays*nangles);
		raft_backprojection_gpu_function(d_image2, d_sino2, sizeImage, nrays, nangles);
		mtx_elementwise_div(d_image3, d_image2, d_bp_ones, sizeImage*sizeImage);
		mtx_elementwise_mult(d_image1, d_image3, d_image1, sizeImage*sizeImage);
	}

	cudaMemcpy(output , d_image1 , ImageMem , cudaMemcpyDeviceToHost);

	
	cudaFree(d_image1);
	cudaFree(d_image2);
	cudaFree(d_image3);
	cudaFree(d_sino1);
	cudaFree(d_sino2);
	cudaFree(d_Sino);
	cudaFree(d_bp_ones);

	
	free(img_inicial);
	free(h_sino_ones);

	return;
}
}
