#ifndef RAFT_IMAGE_BACKPROJECTION_BRESENHAM_H
#define RAFT_IMAGE_BACKPROJECTION_BRESENHAM_H

#include "raft_image.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Backprojection operator.
 *
 * This function backprojects a sinogram at an image. It uses a ray
 * ray-tracing algorithm inspired in Bresenham's pioneering work as
 * implemented in raft_ddaiterator.h.
 *
 * \param image Image which will receive the backprojection;
 * \param sino Sinogram to be backprojected;
 * \param nthreads Number of working threads.
 */
void raft_backprojection_bresenham( raft_image sino,
                                    raft_image image,
                                    int nthreads
                                  );

#ifdef __cplusplus
}
#endif

#endif // #ifndef RAFT_IMAGE_BACKPROJECTION_BRESENHAM_H
