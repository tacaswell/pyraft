#include "raft_radon_gpu.h"
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159265358979
#define INIC -1.0

#define TPBXr 16
#define TPBYr 16

texture<float, cudaTextureType2D, cudaReadModeElementType> texImage;

__global__ void radon_kernel(float* output, float dt, float dtheta, int sizeImage, int nrays, int nangles, float delta, float idelta)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
     int j = 2*(blockIdx.y * blockDim.y + threadIdx.y);

	float cond;
	float ini, x, y, cumsum1, cumsum2, ctheta, stheta, ttheta, theta, t, X, Y;

    if((i < nangles) && (j < nrays))
    {
    	theta = i*dtheta;

    	ini = (float)INIC;

	cond = sqrtf(2)/2;

    	ctheta =cosf(theta);
    	stheta =sinf(theta);
    	ttheta =tanf(theta);

    	if(stheta < cond){
    		cumsum1 = 0;
    		cumsum2 = 0;
    	for(Y = 0; Y < sizeImage; Y++){
				y = ini + Y*delta;
				t = -1.0 + j*dt;
				x = (t/ctheta - y*ttheta);
				X = (float)((x - ini)*idelta);
				//if(X >= 0 && X <= (sizeImage-1))
					cumsum1 += tex2D(texImage, X + 0.5f, Y + 0.5f);
				////////////////////////////
				t = -1.0 + (j+1)*dt;
				x = (t/ctheta - y*ttheta);
				X = (float)((x - ini)*idelta);
				//if(X >= 0 && X < (sizeImage-1))
					cumsum2 += tex2D(texImage, X + 0.5f, Y + 0.5f);
    		}
    		output[(j)*nangles + (nangles-1-i)] = delta*cumsum1/fabsf(ctheta);
    		output[(j+1)*nangles + (nangles-1-i)] = delta*cumsum2/fabsf(ctheta);
    	}

    	else{
    		cumsum1 = 0;
    		cumsum2 = 0;
    	for(X = 0; X < sizeImage; X++){
    			x = ini + X*delta;

    			t = -1.0 + j*dt;
    			y = (t/stheta - x/ttheta);
    			Y = (float)((y - ini)*idelta);
    			//if(Y >= 0 && Y <= (sizeImage-1))
    				cumsum1 += tex2D(texImage, X + 0.5f, Y + 0.5f);
    			/////////////////////
    			t = -1.0 + (j+1)*dt;
    			y = (t/stheta - x/ttheta);
    			Y = (float)((y - ini)*idelta);
    			//if(Y >= 0 && Y <= (sizeImage-1))
    				cumsum2 += tex2D(texImage, X + 0.5f, Y + 0.5f);
    		}
    		output[(j)*nangles + (nangles-1-i)] = delta*cumsum1/fabsf(stheta);
    		output[(j+1)*nangles + (nangles-1-i)] = delta*cumsum2/fabsf(stheta);
    	}

    }
}


extern "C" {
void raft_radon_slantstack_gpu(float* h_output, float* h_input, int sizeImage, int nrays, int nangles)
{
    float* d_output;
    float dt = 2.0/(nrays-1);
    float dtheta = PI/(nangles-1);
    int size = sizeImage*sizeImage*sizeof(float);
    float delta = (float) 2*fabsf(INIC)/(sizeImage-1);
    float idelta = 1/delta;




    /////// KERNEL EXECUTION TIME TEST
			cudaEvent_t start, stop;
			cudaEventCreate(&start);
			cudaEventCreate(&stop);
			float milliseconds2 = 0;

			cudaEventRecord(start);
	//////////////////////////////////

    // Allocate CUDA array in device memory (phantom matrix)
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0,cudaChannelFormatKindFloat);
    cudaArray* cuArray;
    cudaMallocArray(&cuArray, &channelDesc, sizeImage, sizeImage);

    // Copy to device memory the phantom matrix
    cudaMemcpyToArray(cuArray, 0, 0, h_input, size , cudaMemcpyHostToDevice);

    // Set texture parameters
    texImage.addressMode[0] = cudaAddressModeBorder;
    texImage.addressMode[1] = cudaAddressModeBorder;
    texImage.filterMode     = cudaFilterModeLinear;
    /*texImage.normalized     = true;*/

    // Bind the array to the texture reference
    cudaBindTextureToArray(texImage, cuArray, channelDesc);

    // Allocate GPU buffers for the output image
    cudaMalloc(&d_output, sizeof(float) * nrays * nangles);

    // GRID and BLOCKS SIZE
    dim3 threadsPerBlock(TPBXr,TPBYr);
    dim3 grid((nangles/threadsPerBlock.x) + 1, (nrays/threadsPerBlock.y)/2 + 1);

    //KERNEL EXECUTION
    radon_kernel<<<grid, threadsPerBlock>>>(d_output, dt, dtheta, sizeImage, nrays, nangles, delta, idelta);
    cudaDeviceSynchronize();

    /////// PRINT KERNEL EXECUTION TIME
			cudaEventRecord(stop);
			cudaEventSynchronize(stop);
			cudaEventElapsedTime(&milliseconds2, start, stop);
			//fprintf(stderr, "radon: %f ms\n", milliseconds2);
    ///////////////////////////////////

    // Copy output vector from GPU buffer to host memory.
    cudaMemcpy(h_output, d_output, sizeof(float) * nrays * nangles, cudaMemcpyDeviceToHost);

    // unbind texture from buffer
    cudaUnbindTexture(texImage);
    cudaFreeArray(cuArray);
    cudaFree(d_output);
    cudaDeviceReset();
}
}
