#include "raft_filter.h"

// Implementation is in C++, because of std::thread and associated utilities.

#include "raft_mutex_fftw.h"  // For synchronization routines.
#include <fftw3.h>            // FFT routines;
#include <thread>             // for std::thread;
#include <vector>             // for std::vector (will hold the threads).
#include <cmath>              // for std::fabs

extern "C" {

// BLAS dcopy:
void dcopy_( int const *, double const *, int const *, double *, int const * );

}

// Create plans for filtering:
inline void raft_filter1d_create_plans( raft_vector    x,
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

// Perform filtering given the resources (plans and space):
inline void raft_filter1d_with_resources( raft_vector    filter,
                                          raft_vector    x,
                                          fftw_complex * p_out,
                                          raft_vector    x_zp,
                                          fftw_plan      plan,
                                          fftw_plan      inverse_plan
                                        )
{
   // Copy data to zero-padded vector:
   dcopy_( &(x.size), x.p_data, &(x.stride), x_zp.p_data, &(x_zp.stride) );
   for ( int i( x.size ); i < x_zp.size; ++i  )
      raft_vector_element( x_zp, i ) = 0.0;

   // Excecutes FFT:
   fftw_execute_dft_r2c( plan, x_zp.p_data, p_out );

   // Filtering in the Fourier domain:
   unsigned osize =  x.size / 2 + 1;
   unsigned i = 0;
   for ( ; i < osize; ++i )
   {
      double factor = raft_vector_element( filter, i ) / x_zp.size;
      p_out[ i ][ 1 ] *= factor;
      p_out[ i ][ 0 ] *= factor;
   }

   // Excecutes IFFT:
   fftw_execute_dft_c2r( inverse_plan, p_out, x_zp.p_data );

   // Copy back to original vector:
    dcopy_( &(x.size), x_zp.p_data, &(x_zp.stride), x.p_data, &(x.stride) );
}

// To be run by the working threads. Thanks to std::thread
// it can receive multiple arguments.
void raft_filter1d_worker( raft_vector    filter,
                           raft_matrix    mat,
                           fftw_complex * p_out,
                                          raft_vector    x_zp,
                           fftw_plan      plan,
                           fftw_plan      inverse_plan,
                           int            starting_column,
                           int            end_column
                         )
{
   int cur_column = starting_column;
   for ( ; cur_column < end_column; ++cur_column )
      raft_filter1d_with_resources( filter,
                                    raft_matrix_get_column( mat, cur_column ),
                                    p_out,
                                    x_zp,
                                    plan,
                                    inverse_plan
                                  );
}

void raft_filter1d( raft_vector filter,
                    raft_matrix data,
                    int nthreads
                  )
{
   // Make sure we do not have too many or too little threads:
   nthreads = ( nthreads <= data.columns ) ? nthreads : data.columns;
   nthreads = ( nthreads > 0 ) ? nthreads : 1;

   // Space for transformation:
   raft_matrix zpdata = raft_matrix_create( 2 * data.lines - 1, nthreads );
   int osize( zpdata.lines / 2 + 1 );
   fftw_complex * p_out = fftw_alloc_complex( osize * nthreads );
   // Allocation failure, nothing else to do:
   if ( !p_out )
      return;

   // Get things planned:
   fftw_plan plan, inverse_plan;
   raft_filter1d_create_plans( raft_matrix_get_column( zpdata, 0 ),
                               p_out,
                               &plan,
                               &inverse_plan
                             );

   // Base number of columns per thread:
   int base_ncolumns( data.columns / nthreads );
   // Remainder, i.e., number of threads with an extra column:
   int remainder_ncolumns( data.columns % nthreads );
   // Current starting_line for worker thread:
   int cur_starting_column( 0 );

   // Create working threads:
   std::vector< std::thread > threads;
   threads.reserve( nthreads );
   int cur_thread = 0;
   for ( ; cur_thread < nthreads; ++cur_thread )
   {
      int cur_ncolumns( base_ncolumns + ( cur_thread < remainder_ncolumns ) );

      threads.push_back( std::thread( raft_filter1d_worker,
                                      filter,
                                      data,
                                      p_out + cur_thread * osize,
                                      raft_matrix_get_column( zpdata, cur_thread ),
                                      plan,
                                      inverse_plan,
                                      cur_starting_column,
                                      cur_starting_column + cur_ncolumns
                                    )
                       );
      cur_starting_column += cur_ncolumns;
   }

   // Wait for threads to finish the job:
   for ( auto& thread : threads )
      if ( thread.joinable() )
         thread.join();

   // Release resources:
   fftw_destroy_plan( inverse_plan );
   fftw_destroy_plan( plan );
   fftw_free( p_out );
   raft_matrix_destroy( &zpdata );
}

void raft_filter1d_ramp( double cutoff,
                         raft_image sino,
                         int nthreads
                       )
{
   // Vector to hold the filter:
   raft_vector filter = raft_vector_create( sino.data.lines / 2 + 1 );
   // Allocation failed, quit without doing nothing:
   if ( !( filter.p_data ) )
      return;

   // Make sure cutoff values is not too large:
   cutoff = ( cutoff <= 1.0 ) ? cutoff : 1.0;

   // Creation of the filter:
   double delta = 3.14159265358979323846264338327 / ( filter.size - 1 );
   delta /= std::fabs( raft_image_vsamplingdistance( sino ) );
   int i = 0;
   int length = filter.size * cutoff;
   for ( ; i < length; ++i )
      raft_vector_element( filter, i ) = delta * i;
   for ( ; i < filter.size; ++i )
      raft_vector_element( filter, i ) = 0.0;

   // Perform filtering:
   raft_filter1d( filter, sino.data, nthreads );

   // Release resources:
   raft_vector_destroy( &filter );
}
