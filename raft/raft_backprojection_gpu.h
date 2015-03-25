#ifndef RAFT_IMAGE_BACKPROJECTION_GPU_H
#define RAFT_IMAGE_BACKPROJECTION_GPU_H

#ifdef __cplusplus
extern "C" {
#endif

void raft_backprojection_slantstack_gpu(float *image, float *sino, int sizeImage, int nrays, int nangles);

#ifdef __cplusplus
}
#endif

#endif // #ifndef RAFT_IMAGE_BACKPROJECTION_GPU_H



