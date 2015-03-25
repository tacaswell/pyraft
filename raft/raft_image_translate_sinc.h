#ifndef RAFT_IMAGE_TRANSLATE_SINC_H
#define RAFT_IMAGE_TRANSLATE_SINC_H

#include "raft_image.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Periodic horizontal translation.
 *
 * Periodically translate an image in the horizontal direction using sinc
 * interpolation. If you wish to translate vertically, apply to the transpose.
 *
 * \param image Image to be translated;
 * \param delta Amount of translation;
 * \param nthreads Number of threads.
 */
void raft_image_translate_sinc( raft_image image, double delta, int nthreads );

#ifdef __cplusplus
} // extern "C" {
#endif

#endif // #ifndef RAFT_IMAGE_TRANSLATE_SINC_H
