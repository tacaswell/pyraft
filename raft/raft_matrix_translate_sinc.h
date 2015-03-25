#ifndef RAFT_MATRIX_TRANSLATE_SINC_H
#define RAFT_MATRIX_TRANSLATE_SINC_H

#include "raft_vector.h"
#include "raft_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
   * \brief Translates a signal.
   *
   * This function will translate a sampled signal using sinc interpolation
   * by the fft method.
   *
   * Attention: this function creates fftw3 plans. While its own planning is
   * synchronized, using it in parallel with other functions which may call
   * fftw planning-related functions is not safe.
   *
   * \param delta Amount of translation;
   * \param p_x Pointer to vector with samples of the signal to be translated.
   */
void raft_matrix_translate_sinc( double delta, raft_vector x );

/**
   * \brief Translates many signals.
   *
   * This function will translate many sampled signal using sinc interpolation
   * by the fft method.
   *
   * Attention: this function creates fftw3 plans. While its own planning is
   * synchronized, using it in parallel with other functions which may call
   * fftw planning-related functions is not safe.
   *
   * \param deltas Vector with horizontal displacements;
   * \param image Matrix with samples of the signal to be translated;
   * \param nthreads Number of computational threads to be used.
   */
void raft_matrix_translate_many_sinc( raft_vector deltas,
                                      raft_matrix mat,
                                      int         nthreads
                                    );

#ifdef __cplusplus
} // extern "C" {
#endif

#endif //#ifndef RAFT_MATRIX_TRANSLATE_SINC_H
