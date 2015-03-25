#ifndef RAFT_SHEAR_LINEAR_H
#define RAFT_SHEAR_LINEAR_H

#include "raft_image.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Apply horizontal shear.
 *
 * This function applies a horizontal shear to an image. If you wish
 * to apply vertical shear use the transpose image:
 * raft_shear_linear( raft_image_transpose( image ), alpha, nthreads )
 *
 * \param image Image to be sheared;
 * \param alpha Shearing factor;
 * \param nthreads Number of worker threads.
 */
void raft_shear_linear( raft_image image, double alpha, int nthreads );

#ifdef __cplusplus
} // extern "C" {
#endif

#endif // #ifndef RAFT_SHEAR_LINEAR_H
