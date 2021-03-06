#include "raft_radon_gpu_function.h"
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.141592653589793238462643383279502884
#define INIC -1.0

#define TPBXr 16
#define TPBYr 16
#define TPBE 256

/*
Autor: Joao Carlos Cerqueira	email: jc.cerqueira13@gmail.com

A funcao 'void raft_radon_gpu_function' implementa a transformada de
radon classica em GPU. A saida desta funcao sera a exponenciacao da
transformada de radon da imagem de entrada multiplicada pelo valor 
da variavel 'expo': saida = exp[expo*radon]
Para obter a radon classica sem exponenciacao, deve-se selecionar o
valor: expo = 0.0
*/



texture<float, cudaTextureType2D, cudaReadModeElementType> texImage;


extern "C" {
__global__ void radon_gpu_exponential(float* output, float* input1, float input2, int dim){
	
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	
	if((k<dim))
		output[k] = expf( input1[k] * input2 );

	return;

}
}

extern "C" {
__global__ void raft_radon_gpu_kernel(float* output, float dt, float dtheta, int sizeImage, int nrays, int nangles, float delta, float idelta)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
     int j = 2*(blockIdx.y * blockDim.y + threadIdx.y);

	float cond;
	float ini, x, y, linesum1, linesum2, ctheta, stheta, ttheta, theta, t, X, Y;

    if((i < nangles) && (j < nrays))
    {
    	theta = PI - i*dtheta;

    	ini = (float)INIC;

	cond = sqrtf(2)/2;

    	ctheta =cosf(theta);
    	stheta =sinf(theta);
    	ttheta =tanf(theta);

    	if(stheta < cond){
    		linesum1 = 0;
    		linesum2 = 0;
    	for(Y = 0; Y < sizeImage; Y++){
				y = ini + Y*delta;
				t = -1.0 + j*dt;
				x = (t/ctheta - y*ttheta);
				X = (float)((x - ini)*idelta);
				//if(X >= 0 && X <= (sizeImage-1))
					linesum1 += tex2D(texImage, X + 0.5f, Y + 0.5f);
				////////////////////////////
				t = -1.0 + (j+1)*dt;
				x = (t/ctheta - y*ttheta);
				X = (float)((x - ini)*idelta);
				//if(X >= 0 && X < (sizeImage-1))
					linesum2 += tex2D(texImage, X + 0.5f, Y + 0.5f);
    		}
    		output[(j)*nangles + (i)] = delta*linesum1/fabsf(ctheta);
    		output[(j+1)*nangles + (i)] = delta*linesum2/fabsf(ctheta);
    	}

    	else{
    		linesum1 = 0;
    		linesum2 = 0;
    	for(X = 0; X < sizeImage; X++){
    			x = ini + X*delta;

    			t = -1.0 + j*dt;
    			y = (t/stheta - x/ttheta);
    			Y = (float)((y - ini)*idelta);
    			//if(Y >= 0 && Y <= (sizeImage-1))
    				linesum1 += tex2D(texImage, X + 0.5f, Y + 0.5f);
    			/////////////////////
    			t = -1.0 + (j+1)*dt;
    			y = (t/stheta - x/ttheta);
    			Y = (float)((y - ini)*idelta);
    			//if(Y >= 0 && Y <= (sizeImage-1))
    				linesum2 += tex2D(texImage, X + 0.5f, Y + 0.5f);
    		}
    		output[(j)*nangles + (i)] = delta*linesum1/fabsf(stheta);
    		output[(j+1)*nangles + (i)] = delta*linesum2/fabsf(stheta);
    	}

    }
}
}


extern "C" {
void raft_radon_gpu_function(float* d_output, float* d_input, int sizeImage, int nrays, int nangles, float expo){
	
	float dt = 2.0/(nrays-1);
	float dtheta = PI/(nangles);
	float delta = (float) 2*fabsf(INIC)/(sizeImage-1);
	float idelta = 1/delta;

    // Allocate CUDA array in device memory (phantom matrix)
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0,cudaChannelFormatKindFloat);
    cudaArray* cuArray;
    cudaMallocArray(&cuArray, &channelDesc, sizeImage, sizeImage);

    // Transfer device input data to texture
    cudaMemcpyToArray(cuArray, 0, 0, d_input, sizeImage*sizeImage*sizeof(float) , cudaMemcpyDeviceToDevice);

    // Set texture parameters
    texImage.addressMode[0] = cudaAddressModeBorder;
    texImage.addressMode[1] = cudaAddressModeBorder;
    texImage.filterMode     = cudaFilterModeLinear;
    /*texImage.normalized     = true;*/

    // Bind the array to the texture reference
    cudaBindTextureToArray(texImage, cuArray, channelDesc);

	dim3 threadsPerBlock(TPBXr,TPBYr);
    	dim3 grid((nangles/threadsPerBlock.x) + 1, (nrays/threadsPerBlock.y)/2 + 1);
	raft_radon_gpu_kernel<<<grid, threadsPerBlock>>>(d_output, dt, dtheta, sizeImage, nrays, nangles, delta, idelta);
	
	// OUTPUT = e^[expo * std_output]
	if( expo != 0.0 ){
		int TPB = (int)TPBE;
		int grid2;
		int dim;
		dim = nrays*nangles;
		grid2 = (int) ceil(	(float)dim/(float)TPB	);
		radon_gpu_exponential<<<grid2, TPB>>>(d_output, d_output, expo, dim);
	}

	cudaDeviceSynchronize();
	
	cudaUnbindTexture(texImage);
	cudaFreeArray(cuArray);
	return;

}
}

