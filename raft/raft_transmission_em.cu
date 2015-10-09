#include "raft_cuda_aux.h"
#include "raft_backprojection_gpu_function.h"
#include "raft_radon_gpu_function.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
Autor: Joao Carlos Cerqueira	email: jc.cerqueira13@gmail.com

Implementacao de um algoritmo iterativo de reconstrucao restrito
a problemas de tomografia baseados em transmissao.
*/

extern "C"{
void raft_tr_em_gpu(float *output, float *rawsino, float *flat, int sizeImage, int nrays, int nangles, int niter)
{
	int i;	

	int SinoMem = sizeof(float)*nrays*nangles;
	int ImageMem = sizeof(float)*sizeImage*sizeImage;	

	float *d_image1, *d_image2, *d_image3, *d_flat, *d_sino1, *d_sino2, *d_sino3, *d_Sino;
	float *img_inicial;

// CUDA LINEAR MEMORY
	cudaMalloc(&d_sino1, SinoMem);
	cudaMalloc(&d_sino2, SinoMem);
	cudaMalloc(&d_sino3, SinoMem);	
	cudaMalloc(&d_Sino, SinoMem);	
	cudaMalloc(&d_flat, SinoMem);
	cudaMalloc(&d_image1, ImageMem);
	cudaMalloc(&d_image2, ImageMem);
	cudaMalloc(&d_image3, ImageMem);

// IMAGEM INICIAL
/* Futuramente, pode ser interessante ter a imagem inicial como uma variável */
	img_inicial = (float *)malloc(ImageMem);
	for (i = 0; i < (sizeImage*sizeImage); i++)		
		img_inicial[i] = 1.0;

	cudaMemcpy(d_image1, img_inicial, ImageMem, cudaMemcpyHostToDevice);


// COPIAS DA CPU PARA A GPU
	cudaMemcpy(d_flat, flat, SinoMem, cudaMemcpyHostToDevice);
	cudaMemcpy(d_Sino, rawsino, SinoMem, cudaMemcpyHostToDevice);


	raft_backprojection_gpu_function(d_image3, d_Sino, sizeImage, nrays, nangles);

	for(i=0; i<niter; i++){
		raft_radon_gpu_function(d_sino1, d_image1, sizeImage, nrays, nangles, -1.0);
		mtx_elementwise_mult(d_sino1, d_sino1, d_flat, nrays*nangles);
		raft_backprojection_gpu_function(d_image2, d_sino1, sizeImage, nrays, nangles);
		mtx_elementwise_div(d_image2, d_image2, d_image3, sizeImage*sizeImage);
		mtx_elementwise_mult(d_image1, d_image1, d_image2, sizeImage*sizeImage);

		mtx_elementwise_mult(d_image1, d_image1, d_image2, sizeImage*sizeImage);
	
	}

// COPIA DO RESULTADO DA GPU PARA A CPU
	cudaMemcpy(output , d_image1 , ImageMem , cudaMemcpyDeviceToHost);

	cudaFree(d_sino1);
	cudaFree(d_sino2);
	cudaFree(d_sino3);
	cudaFree(d_Sino);
	cudaFree(d_flat);
	cudaFree(d_image1);
	cudaFree(d_image2);
	cudaFree(d_image3);
	free(img_inicial);

	return;
}// END VOID
}
