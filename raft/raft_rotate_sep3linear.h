#ifndef RAFT_ROTATE_SEP3LINEAR_H
#define RAFT_ROTATE_SEP3LINEAR_H

#include "raft_image.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Rotate an image using sep3 with linear interpolation.
 *
 * \param image Image to be rotated;
 * \param theta Angle of rotation;
 * \param nthreads Number of worker threads.
 */
void raft_rotate_sep3linear( raft_image image, double theta, int nthreads );

#ifdef __cplusplus
} // extern "C" {
#endif

#endif // #ifndef RAFT_ROTATE_SEP3LINEAR_H
