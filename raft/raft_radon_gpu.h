#ifndef RAFT_RADON_GPU_H
#define RAFT_RADON_GPU_H

#ifdef __cplusplus
extern "C" {
#endif

void raft_radon_slantstack_gpu(float* h_output, float* h_input, int sizeImage, int nrays, int nangles);

#ifdef __cplusplus
}
#endif

#endif // #ifndef RAFT_RADON_GPU_H
