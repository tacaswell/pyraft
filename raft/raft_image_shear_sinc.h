#ifndef RAFT_IMAGE_SHEAR_SINC_H
#define RAFT_IMAGE_SHEAR_SINC_H

#include "raft_image.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Apply horizontal shear.
 *
 * This function applies a horizontal shear to an image. If you wish
 * to apply vertical shear use the transpose image:
 * raft_sep3_sinc_rotate( raft_image_transpose( image ), alpha, nthreads )
 *
 * \param image Image to be sheared;
 * \param alpha Shearing factor;
 * \param nthreads Number of worker threads.
 */
void raft_image_shear_sinc( raft_image image, double alpha, int nthreads );

#ifdef __cplusplus
} // extern "C" {
#endif

#endif // #ifndef RAFT_IMAGE_SHEAR_SINC_H
