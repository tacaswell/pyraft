#ifndef RAFT_IMAGE_ROTATE_SEP3SINC_H
#define RAFT_IMAGE_ROTATE_SEP3SINC_H

#include "raft_image.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Rotate an image using sep3 with sinc interpolation.
 *
 * \TODO Make shear relative to the origin, not the center of the image!
 *
 * \param image Image to be rotated;
 * \param theta Angle of rotation;
 * \param nthreads Number of worker threads.
 */
void raft_image_rotate_sep3sinc( raft_image image, double theta, int nthreads );

#ifdef __cplusplus
} // extern "C" {
#endif

#endif // #ifndef RAFT_IMAGE_ROTATE_SEP3SINC_H
