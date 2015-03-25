#ifndef RAFT_TRANSLATE_LINEAR_H
#define RAFT_TRANSLATE_LINEAR_H

#include "raft_vector.h"
#include "raft_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
   * \brief Translates a signal.
   *
   * This function will translate a sampled signal using linear interpolation.
   *
   * \param delta Amount of translation;
   * \param x Vector with samples of the signal to be translated.
   */
void raft_translate_linear( double delta, raft_vector x );

/**
   * \brief Translates many signals.
   *
   * This function will translate many sampled signal using linear interpolation.
   *
   * \param deltas Vector with horizontal displacements;
   * \param image Matrix with samples of the signal to be translated;
   * \param nthreads Number of computational threads to be used.
   */
void raft_translate_many_linear( raft_vector deltas,
                                 raft_matrix mat,
                                 int         nthreads
                               );

#ifdef __cplusplus
} // extern "C" {
#endif

#endif //#ifndef RAFT_TRANSLATE_LINEAR_H
