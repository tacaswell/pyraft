#ifndef RAFT_IMAGE_PHANTOM_FROMDESC_H
#define RAFT_IMAGE_PHANTOM_FROMDESC_H

#include "raft_image.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Regular sample of the phantom.
 *
 * This function computes the value of a phantom at the points
 * sampled by the raft_image image.
 *
 * \param image Image to contain the phantom. This parameter
 * determines, through the viewport and size of storage, the
 * points where the image will be sampled,
 * \param desc Description of the phantom. Each line describes an
 * ellipse. It must contain 6 collumns with the following parameters:
 * value of elliptic region; horizontal axis; vertical axis;
 * x-coordinate of center; y-coordinate of center; rotation angle
 * 
 */

void raft_image_phantom_fromdesc( raft_image image, raft_matrix desc );

#ifdef __cplusplus
} // extern "C" {
#endif

#endif // #ifndef RAFT_IMAGE_PHANTOM_FROMDESC_H
