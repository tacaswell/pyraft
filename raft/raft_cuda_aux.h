#ifndef RAFT_CUDA_AUX_H
#define RAFT_CUDA_AUX_H

#ifdef __cplusplus
extern "C" {
#endif

	void mtx_elementwise_div(float* output, float* input1, float* input2, int dim);
  
	void mtx_elementwise_mult(float* output, float* input1, float* input2, int dim);
	
	void mtx_elementwise_sum(float* output, float* input1, float* input2, int dim);
  
  	void mtx_elementwise_minus(float* output, float* input1, float* input2, int dim);
  

#ifdef __cplusplus
} 
#endif

#endif 
