#ifndef RAFT_BACKPROJECTION_GPU_FUNCTION_H
#define RAFT_BACKPROJECTION_GPU_FUNCTION_H

#ifdef __cplusplus
extern "C" {
#endif

void raft_backprojection_gpu_function(float *d_output, float *h_input, int sizeImage, int nrays, int nangles);  

#ifdef __cplusplus
} // extern "C" {
#endif

#endif // #ifndef RAFT_BACKPROJECTION_GPU_FUNCTION_H
