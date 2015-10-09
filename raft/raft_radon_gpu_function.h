#ifndef RAFT_RADON_GPU_FUNCTION_H
#define RAFT_RADON_GPU_FUNCTION_H

#ifdef __cplusplus
extern "C" {
#endif

void raft_radon_gpu_function(float* d_output, float* d_input, int sizeImage, int nrays, int nangles, float expo);
void raft_exp_radon_gpu_function(float* d_output, float* d_input, int sizeImage, int nrays, int nangles);
  

#ifdef __cplusplus
} // extern "C" {
#endif

#endif // #ifndef RAFT_RADON_GPU_FUNCTION_H
