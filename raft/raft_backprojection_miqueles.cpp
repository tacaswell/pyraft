#include "raft_backprojection_miqueles.h"
#include "raft_image_functions.h"

#include <fftw3.h>
#include <algorithm>
#include <iostream>
#include <thread>
#include <string.h>

////////////////////// polar.cp:
void to_pos(raft_matrix &x)
{
	int N1 = x.lines;
	int N2 = x.columns;
	for(int i=0; i<N1; i++) {
		for(int j=0; j<N2; j++) {
			raft_matrix_element(x, i, j) = fabs(raft_matrix_element(x, i, j));	
		}
	}
}

////////////////////// Getting the HG ////////////////////////
extern "C" {

// BLAS dcopy:
void dcopy_( int const *, double const *, int const *, double *, int const * );

}

inline void raft_getg_with_resources(raft_vector sigma, 
				raft_vector x, 
				fftw_complex *p_out, 
				raft_vector x_zp, 
				raft_vector r_r,
				raft_vector r_i,
				fftw_plan plan) 
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
	for ( ; i < osize; ++i ) {
// 		double factor = 1/(raft_vector_element( sigma, i ) * x_zp.size);
		double factor = 1/(raft_vector_element( sigma, i ));
// 		double factor = 1/( x_zp.size);
		p_out[ i ][ 1 ] *= factor;
		p_out[ i ][ 0 ] *= factor;
	}

	// Copy back to r_r and r_i:
	for(int i=0; i<osize; i++)
	{
		raft_vector_element(r_r, i) = p_out[i][0];
		raft_vector_element(r_i, i) = p_out[i][1];
	}
}

void raft_getg_worker( raft_vector   		sigma,
				raft_matrix    	mat,
				raft_matrix	mat_out_r,
				raft_matrix	mat_out_i,
				fftw_complex * 	p_out,
				raft_vector    	x_zp,
				fftw_plan      	plan,
				int            	starting_column,
				int            	end_column)
{
	int cur_column = starting_column;
	for ( ; cur_column < end_column; ++cur_column )
		raft_getg_with_resources( sigma,
					    raft_matrix_get_column( mat, cur_column ),
					    p_out,
					    x_zp,
					    raft_matrix_get_column( mat_out_r, cur_column ),
					    raft_matrix_get_column( mat_out_i, cur_column ),
					    plan
					  );
}

void raft_getg( raft_vector sigma,
			raft_matrix data,
			raft_matrix &rmap_r, 
			raft_matrix &rmap_i,
			int nthreads)
{
	rmap_r = raft_matrix_create( data.lines, data.columns );
	rmap_i = raft_matrix_create( data.lines, data.columns );
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
	fftw_plan plan;
	raft_vector x = raft_matrix_get_column( zpdata, 0 );
	plan = fftw_plan_many_dft_r2c( 1, // Number of dimensions - only one;
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
		std::thread thread =  std::thread( raft_getg_worker,
					      sigma,
					      data,
					      rmap_r, 
					      rmap_i,
					      p_out + cur_thread * osize,
					      raft_matrix_get_column( zpdata, cur_thread ),
					      plan,
					      cur_starting_column,
					      cur_starting_column + cur_ncolumns
					    );

		threads.push_back(move(thread));
		cur_starting_column += cur_ncolumns;
	}

	// Wait for threads to finish the job:
	for ( auto& thread : threads )
		if ( thread.joinable() )
			thread.join();

	// Release resources:
	fftw_destroy_plan( plan );
	fftw_free( p_out );
	raft_matrix_destroy( &zpdata );
}

void get_hg_ramp(raft_image sino, raft_matrix &rmap_r, raft_matrix &rmap_i, int nthreads)
{
	int Ns = sino.data.lines;
	double ds = 1.0/Ns;
	raft_vector sigma = raft_vector_create(sino.data.lines / 2 + 1);
// 	double sigma_0 = sino.tl_y;
// 	double sigma_1 = sino.br_y/(2*pi);
	double sigma_0 = 0;
	double sigma_1 = 1.0/(2*ds); 

	double d_sigma = sigma_1/Ns; 

	raft_vector_element(sigma, 0) = Ns;
	for(int i=1; i<sigma.size; i++) {
		double s = sigma_0 + i*d_sigma;
		raft_vector_element(sigma, i) = s; 
	}

	raft_getg(sigma, sino.data, rmap_r, rmap_i, nthreads);
	raft_vector_destroy(&sigma);
}

//////////////////////// Backprojection function //////////////////////////////
void raft_backprojection_miqueles(raft_image sino, raft_image res_, int nthreads)
{
	int Nx = res_.data.columns;
	int Ny = res_.data.lines;
	// Transform from sinogram coordinates to semi-polar
	int Nt = sino.data.columns;
	int Ns = sino.data.lines;
	//std::cout<<"ORIG: Nt="<<Nt<<"; Ns="<<Ns<<std::endl;
	raft_image res = sino2sp(sino);
	// The smallest errors occurs when we are using square matrixes
	Nt = res.data.columns;
	Ns = res.data.lines;
	//std::cout<<"SP: Nt="<<Nt<<"; Ns="<<Ns<<std::endl;
	Ns = std::max(Nt, res.data.lines);
	res = zero_padding_on_s(res, 2*Ns);
	Nt = res.data.columns;
	Ns = res.data.lines;
	//std::cout<<"SP ZP: Nt="<<Nt<<"; Ns="<<Ns<<std::endl;
	// Initializing the arrays for Fourier image
	raft_image res_r = res;
	raft_image res_i = res;
	// Calculating the values of Fourier image of sinogram and dividing on sigma
	get_hg_ramp(res, res_r.data, res_i.data, nthreads);
	// Transforming the frequency domain to Cartezian coordinates
	raft_image tmp_r, tmp_i;
	sp2c_miqueles(res_r, res_i, tmp_r, tmp_i, Nx, Ny, nthreads);
	// doing the 2D inverse FFT to obtain the result  
	ifftw_2d(tmp_r.data, tmp_i.data, nthreads);
	res = tmp_r;
	// Shifting
	fft_shift_2d(res.data);
	// Some postprocessing...
	// Mirroring to positive only values
	to_pos(res.data);
	// Interpolation to the given size Nx, Ny
	res = bl_interpolate_mt(res, Nx, Ny, nthreads);
	// Just set mean nearest values in center
	int r = Ns/200;
	for(int j=-r; j<r; j++) {
		for(int i=-r; i<r; i++) {
			raft_matrix_element(res.data, Nx/2 + j, Ny/2+i) = 
				raft_matrix_element(res.data, Nx/2 - r - 1, Ny/2 - r - 1);
		}
	}

	memcpy((double*)res_.data.p_data, (double*)res.data.p_data, sizeof(double)*Nx*Ny); 
}


