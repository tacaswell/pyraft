#ifndef RAFT_RADON_FROMDESC_H
#define RAFT_RADON_FROMDESC_H

#include "raft_image.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Regular sample of the sinogram.
 *
 * This function computes the Radon Transform at the points
 * sampled by the raft_image sino. The result is the exact Radon
 * Transform of the phantom described by desc at such points,
 * up to rounding errors.
 *
 * \param sino Image to contain the sinogram. This parameter
 * determines, through the viewport and size of storage, the
 * points where the RT will be sampled,
 * \param desc Description of the phantom.
 */
void raft_radon_fromdesc( raft_image sino,
                          raft_matrix desc
                        );

#ifdef __cplusplus
} // extern "C" {
#endif

#endif // #ifndef RAFT_RADON_FROMDESC_H
