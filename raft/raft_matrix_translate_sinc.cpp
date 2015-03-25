#include "raft_matrix_translate_sinc.h"

// Implementation is in C++, because of std::thread and associated utilities.

#include "raft_mutex_fftw.h" // For synchronization routines.
#include <fftw3.h>   // FFT routines;
#include <math.h>    // for cos() and sin();
#include <thread>    // for std::thread;
#include <vector>    // for std::vector (will hold the threads).

// Create plans for translation:
inline void raft_matrix_translate_sinc_create_plans( raft_vector    x,
                                                     fftw_complex * p_out,
                                                     fftw_plan    * p_plan,
                                                     fftw_plan    * p_inverse_plan
                                                   )
{
   //Locks resources:
   raft_mutex_fftw_lock();

      // FFT plan. We need to use fftw's advanced interface because of variable stride.
      *p_plan = fftw_plan_many_dft_r2c( 1, // Number of dimensions - only one;
                                        &(x.size), // Size of each dimension - only length;
                                        1, // How many? 1;
                                        x.p_data, // Address of input vector;
                                        NULL, // Size of the "host" input array - not used;
                                        x.stride, // Stride of the input vector;
                                        0, // Distance between successive input arrays;
                                        p_out, // Output vector;
                                        NULL, // Size of the "host" output array - not used;
                                        1, // Stride of the output vector;
                                        0, // Distance between successive output arrays;
                                        FFTW_ESTIMATE // Flags.
                                      );
      // IFFT plan.
      *p_inverse_plan = fftw_plan_many_dft_c2r( 1,
                                                &(x.size),
                                                1,
                                                p_out,
                                                NULL,
                                                1,
                                                0,
                                                x.p_data,
                                                NULL,
                                                x.stride,
                                                0,
                                                FFTW_ESTIMATE
                                              );

      // Unlock resources:
      raft_mutex_fftw_unlock();
}

// Perform translation given the resources (plans and space):
inline void raft_matrix_translate_sinc_with_resources( double         delta,
                                                       raft_vector    x,
                                                       fftw_complex * p_out,
                                                       fftw_plan      plan,
                                                       fftw_plan      inverse_plan
                                                     )
{
   // Excecutes FFT:
   fftw_execute_dft_r2c( plan, x.p_data, p_out );

   // Translation in the Fourier domain:
   double delta_angle = -6.283185307179586476925286766559 * delta / x.size;
   unsigned osize =  x.size / 2 + 1;
   unsigned i = 0;
   for ( ; i < osize; ++i )
   {
      double angle = delta_angle * i;
      double s = sin( angle );
      double c = cos( angle );
      double tmp_real = ( ( p_out[ i ][ 0 ] * c ) - ( p_out[ i ][ 1 ] * s ) ) / x.size;
      p_out[ i ][ 1 ] = ( ( p_out[ i ][ 0 ] * s ) + ( p_out[ i ][ 1 ] * c ) ) / x.size;
      p_out[ i ][ 0 ] = tmp_real;
   }

   // Excecutes IFFT:
   fftw_execute_dft_c2r( inverse_plan, p_out, x.p_data );
}

void raft_matrix_translate_sinc( double delta, raft_vector x )
{
   // Space for transformation:
   // TODO: Test for allocation failure!
   fftw_complex * p_out = fftw_alloc_complex( x.size / 2 + 1 );

   // Get things planned:
   fftw_plan plan, inverse_plan;
   raft_matrix_translate_sinc_create_plans( x, p_out, &plan, &inverse_plan );

   // Perform translation:
   raft_matrix_translate_sinc_with_resources( delta, x, p_out, plan, inverse_plan );

   // Release resources:
   fftw_destroy_plan( inverse_plan );
   fftw_destroy_plan( plan );
   fftw_free( p_out );
}

// To be run by the working threads. Thanks to std::thread
// it can receive multiple arguments.
void raft_matrix_translate_many_sinc_worker( raft_vector    deltas,
                                             raft_matrix    mat,
                                             fftw_complex * p_out,
                                             fftw_plan      plan,
                                             fftw_plan      inverse_plan,
                                             int            starting_line,
                                             int            n_lines
                                           )
{
   for ( int i( 0 ); i < n_lines; ++i )
      raft_matrix_translate_sinc_with_resources( raft_vector_element( deltas, i + starting_line ),
                                                 raft_matrix_get_line( mat, i + starting_line ),
                                                 p_out,
                                                 plan,
                                                 inverse_plan
                                               );
}

void raft_matrix_translate_many_sinc( raft_vector deltas,
                                      raft_matrix mat,
                                      int nthreads
                                    )
{
   // Number of translations to be performed:
   int size = ( mat.lines < deltas.size ) ? mat.lines : deltas.size;
   // Make sure we do not have too many or too little threads:
   nthreads = ( nthreads <= size ) ? nthreads : size;
   nthreads = ( nthreads > 0 ) ? nthreads : 1;

   // Space for transformation:
   // TODO: Test for allocation failure!
   int osize( mat.columns / 2 + 1 );
   fftw_complex * p_out = fftw_alloc_complex( osize * nthreads );

   // Get things planned:
   fftw_plan plan, inverse_plan;
   raft_matrix_translate_sinc_create_plans( raft_matrix_get_line( mat, 0 ),
                                            p_out,
                                            &plan,
                                            &inverse_plan
                                          );

   // Base number of lines per thread:
   int base_nlines( size / nthreads );
   // Remainder, i.e., number of threads with an extra line:
   int remainder_nlines( size % nthreads );
   // Current starting_line for worker thread:
   int cur_starting_line( 0 );

   // Create working threads:
   std::vector< std::thread > threads;
   threads.reserve( nthreads );
   int cur_thread = 0;
   for ( ; cur_thread < nthreads; ++cur_thread )
   {
      int cur_nlines( base_nlines + ( cur_thread < remainder_nlines ) );

      threads.push_back( std::thread( raft_matrix_translate_many_sinc_worker,
                                      deltas,
                                      mat,
                                      p_out + cur_thread * osize,
                                      plan,
                                      inverse_plan,
                                      cur_starting_line,
                                      cur_nlines
                                    )
                       );
      cur_starting_line += cur_nlines;
   }

   // Wait for threads to finish the job:
   for ( auto& thread : threads )
      if ( thread.joinable() )
         thread.join();

   // Release resources:
   fftw_destroy_plan( inverse_plan );
   fftw_destroy_plan( plan );
   fftw_free( p_out );
}

//g++ -std=c++0x -malign-double -O3 fft_translate.cpp -lfftw3
