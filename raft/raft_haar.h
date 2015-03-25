#ifndef RAFT_HAAR_H
#define RAFT_HAAR_H

#include "raft_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Haar transform.
 *
 * This function computes the Haar transform of a matrix.
 *
 * \param data Matrix to be transformed;
 * \param nlevels Number of levels for transformation;
 * \param nthreads Number of working threads.
 */
void raft_haar( raft_matrix data,
                int nlevels,
                int nthreads
              );

/**
 * \brief Inverse Haar transform.
 *
 * This function computes the inverse Haar transform of a matrix.
 *
 * \param data Matrix to be transformed;
 * \param nlevels Number of levels for transformation;
 * \param nthreads Number of working threads.
 */
void raft_ihaar( raft_matrix data,
                 int nlevels,
                 int nthreads
               );

#ifdef __cplusplus
}
#endif

#endif // #ifndef RAFT_IMAGE_BACKPROJECTION_BRESENHAM_H
