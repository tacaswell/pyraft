#include "raft_backprojection_gpu.h"
#include "raft_radon_gpu.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.141592
#define INIC -1.0

#define TPBXb 32
#define TPBYb 32

#define TPBXr 16
#define TPBYr 16


extern "C" {
void raft_em_gpu(float *output, float *sino, int sizeImage, int nrays, int nangles, int niter){
	int i, k;
	float avg;

	float *sino_ones, *img_ones, *it_image1, *it_image2, *it_sino1, *it_sino2;
	div_t D;
	int angle, ray;
	float t;

	it_image1	= (float *)malloc(sizeImage*sizeImage*sizeof(float));
	it_image2	= (float *)malloc(sizeImage*sizeImage*sizeof(float));
	it_sino1	= (float *)malloc(nrays*nangles*sizeof(float));
	it_sino2	= (float *)malloc(nrays*nangles*sizeof(float));

	///// BACKPROJECTION OF ONES/////
	img_ones	= (float *)malloc(sizeImage*sizeImage*sizeof(float));
	sino_ones	= (float *)malloc(nrays*nangles*sizeof(float));

	for(i = 0; i < (nrays*nangles); i++) 	sino_ones[i] = 10;

	raft_backprojection_slantstack_gpu(img_ones, sino_ones, sizeImage, nrays, nangles);

	free(sino_ones);
	////////////////////////////////
	avg = 0;
	for(i=0;i<(nrays*nangles); i++){
		avg +=sino[i];
	}
	avg = avg/(nrays*nangles);

	//// INITIAL IMAGE//////////////
	for (i = 0; i < (sizeImage*sizeImage); i++)		it_image1[i] = avg+1;
	///////////////////////////////

for(k = 0; k < niter ; k++){


	raft_radon_slantstack_gpu(it_sino1, it_image1, sizeImage, nrays, nangles);

	for(i = 0; i < (nrays*nangles); i++){
		D = div(i, nrays);
		ray = D.quot;
		angle = D.rem;
		t = -1.0 + ray*(2.0/(nrays-1));

		////if (fabs(t) > sqrt(2)/2){
		//	it_sino2[i] = 0.0;
		//}
		//else{
			if(it_sino1[i]!=0) it_sino2[i] = sino[i]/it_sino1[i];
		//}
	}


	raft_backprojection_slantstack_gpu(it_image2, it_sino2, sizeImage, nrays, nangles);

	for (i = 0; i < (sizeImage*sizeImage); i++)
		it_image1[i] = (it_image2[i]*it_image1[i]) / img_ones[i];

}
	for (i = 0; i < (sizeImage*sizeImage); i++)
		output[i] = it_image1[sizeImage*sizeImage-i-1];

	return;
}
}
