#include "raft_backprojection_bst.h"
#include "raft_image_functions2.h"
#include "raft_image.h"

#include <fftw3.h>
#include <algorithm>
#include <iostream>
#include <thread>
#include <string.h>

#include <boost/math/special_functions/bessel.hpp>


double const pi = 3.1415926535897932384626433832795;

raft_bst raft_bst_plan_create( raft_image sino, double padding_coeff )
{
	int nrays, nviews;
	raft_bst result;

	nrays  = sino.data.lines;
	nviews = sino.data.columns; 
	/* sinogram @ polar coordinates */

	result.polarsino = raft_image_create( nrays/2, 2*nviews );

	if ( raft_image_is_empty( result.polarsino ) )
	{
		result.polarsino.tl_x = result.polarsino.tl_y = 
			result.polarsino.br_x = result.polarsino.br_y = 0.0;
		return result;
	}

	result.polarsino.tl_x = -pi;
	result.polarsino.br_x = pi;
	result.polarsino.tl_y = 1;
	result.polarsino.br_y = 0;

	// Zero-padded sinogram
	int nrays_zp = snapTransformSize_bst(padding_coeff*nrays/2 - 1);

	result.zero_padded_sino = raft_image_create( nrays_zp, 2*nviews );

	std::cout<<"NRAYS_ZP="<<nrays_zp<<std::endl;
	if ( raft_image_is_empty( result.zero_padded_sino ) )
	{
		result.zero_padded_sino.tl_x = result.zero_padded_sino.tl_y = 
			result.zero_padded_sino.br_x = result.zero_padded_sino.br_y = 0.0;
		return result;
	}

	double ds = 1.0/result.polarsino.data.lines;
	double s_max = nrays_zp*ds;

	result.zero_padded_sino.br_x = result.polarsino.br_x;
	result.zero_padded_sino.tl_x = result.polarsino.tl_x;
	result.zero_padded_sino.br_y = result.polarsino.br_y;
	result.zero_padded_sino.tl_y = s_max;

	// Fourier image of zero padded	sinogram (polar representation)
	result.fftp_re = raft_image_create( nrays_zp, 2*nviews );
	result.fftp_im = raft_image_create( nrays_zp, 2*nviews );
	
	if ( raft_image_is_empty( result.fftp_re ) )
	{
		result.fftp_re.tl_x = result.fftp_re.tl_y = 
			result.fftp_re.br_x = result.fftp_re.br_y = 0.0;
		return result;
	}
	if ( raft_image_is_empty( result.fftp_im ) )
	{
		result.fftp_im.tl_x = result.fftp_im.tl_y = 
			result.fftp_im.br_x = result.fftp_im.br_y = 0.0;
		return result;
	}

	result.fftp_re.tl_x = result.zero_padded_sino.tl_x;
	result.fftp_re.br_x = result.zero_padded_sino.br_x;
	result.fftp_re.tl_y = result.zero_padded_sino.tl_y;
	result.fftp_re.br_y = result.zero_padded_sino.br_y;

	result.fftp_im.tl_x = result.zero_padded_sino.tl_x;
	result.fftp_im.br_x = result.zero_padded_sino.br_x;
	result.fftp_im.tl_y = result.zero_padded_sino.tl_y;
	result.fftp_im.br_y = result.zero_padded_sino.br_y;

	// Fourier image, interpolateed to the Cartesian coordinates
		
	result.fftc_re = raft_image_create( nrays_zp, nrays_zp );  // FIXME think about sizes here
	result.fftc_im = raft_image_create( nrays_zp, nrays_zp );
	
	if ( raft_image_is_empty( result.fftc_re ) )
	{
		result.fftc_re.tl_x = result.fftc_re.tl_y = 
			result.fftc_re.br_x = result.fftc_re.br_y = 0.0;
		return result;
	}
	if ( raft_image_is_empty( result.fftc_im ) )
	{
		result.fftc_im.tl_x = result.fftc_im.tl_y = 
			result.fftc_im.br_x = result.fftc_im.br_y = 0.0;
		return result;
	}

	double s1 = std::max(result.fftp_re.br_y, result.fftp_re.tl_y);
	
	double x0 = -s1;
	double x1 =  s1;
	double y0 = -s1;
	double y1 =  s1;
	result.fftc_re.tl_x = x0;
	result.fftc_re.br_x = x1;
	result.fftc_re.tl_y = y1;
	result.fftc_re.br_y = y0;

	result.fftc_im.tl_x = x0;
	result.fftc_im.br_x = x1;
	result.fftc_im.tl_y = y1;
	result.fftc_im.br_y = y0;
	// Postprocessing and interpolation
	result.int1 = raft_image_create(result.fftc_re.data.lines/2, 
			result.fftc_re.data.columns/2);

	double coeff = (double)result.polarsino.data.lines/(double)result.zero_padded_sino.data.lines;
	std::cout<<"coeff="<<coeff<<std::endl;
	int Cx = coeff*result.fftc_re.data.lines;
	int Cy = coeff*result.fftc_re.data.columns;

	result.cut  = raft_image_create(Cx, Cy);
	
	result.int1.tl_x = result.fftc_re.tl_x;
	result.int1.br_x = result.fftc_re.br_x;
	result.int1.tl_y = result.fftc_re.tl_y;
	result.int1.br_y = result.fftc_re.br_y;

	result.cut.tl_x = coeff*result.fftc_re.tl_x;
	result.cut.br_x = coeff*result.fftc_re.br_x;
	result.cut.tl_y = coeff*result.fftc_re.tl_y;
	result.cut.br_y = coeff*result.fftc_re.br_y;

	return result;
}

void raft_bst_plan_destroy( raft_bst * plan )
{
	/* destroy sinogram @ polar coordinates */

	raft_image_destroy( &( plan->polarsino ) );
	raft_image_destroy( &( plan->zero_padded_sino ) );
	raft_image_destroy( &( plan->fftp_re ) );
	raft_image_destroy( &( plan->fftp_im ) );
	raft_image_destroy( &( plan->fftc_re ) );
	raft_image_destroy( &( plan->fftc_im ) );
	raft_image_destroy( &( plan->int1 ) );
	raft_image_destroy( &( plan->cut ) );

	/* */
}

void raft_bst_plan_set_corners( raft_image sino, raft_bst * plan )
{
	plan->polarsino.tl_x = -pi;
	plan->polarsino.br_x = pi;
	plan->polarsino.tl_y = 1.0;
	plan->polarsino.br_y = 0;

	int nrays = sino.data.lines;
	double ds = 2.0/nrays;
	int nrays_zp = snapTransformSize_bst(plan->padding_coeff*nrays/2 - 1);
	double s_max = nrays_zp*ds;
	plan->zero_padded_sino.br_x = plan->polarsino.br_x;
	plan->zero_padded_sino.tl_x = plan->polarsino.tl_x;
	plan->zero_padded_sino.br_y = plan->polarsino.br_y;
	plan->zero_padded_sino.tl_y = s_max;

	plan->fftp_re.tl_x = plan->zero_padded_sino.tl_x;
	plan->fftp_re.br_x = plan->zero_padded_sino.br_x;
	plan->fftp_re.tl_y = plan->zero_padded_sino.tl_y;
	plan->fftp_re.br_y = plan->zero_padded_sino.br_y;

	plan->fftp_im.tl_x = plan->zero_padded_sino.tl_x;
	plan->fftp_im.br_x = plan->zero_padded_sino.br_x;
	plan->fftp_im.tl_y = plan->zero_padded_sino.tl_y;
	plan->fftp_im.br_y = plan->zero_padded_sino.br_y;

	double s1 = std::max(plan->fftp_re.br_y, plan->fftp_re.tl_y);
	
	double x0 = -s1;
	double x1 =  s1;
	double y0 = -s1;
	double y1 =  s1;
	plan->fftc_re.tl_x = x0;
	plan->fftc_re.br_x = x1;
	plan->fftc_re.tl_y = y1;
	plan->fftc_re.br_y = y0;

	plan->fftc_im.tl_x = x0;
	plan->fftc_im.br_x = x1;
	plan->fftc_im.tl_y = y1;
	plan->fftc_im.br_y = y0;

	plan->int1.tl_x = plan->fftc_re.tl_x;
	plan->int1.br_x = plan->fftc_re.br_x;
	plan->int1.tl_y = plan->fftc_re.tl_y;
	plan->int1.br_y = plan->fftc_re.br_y;

	double coeff = (double)plan->polarsino.data.lines/(double)plan->zero_padded_sino.data.lines;
	plan->cut.tl_x = coeff*plan->fftc_re.tl_x;
	plan->cut.br_x = coeff*plan->fftc_re.br_x;
	plan->cut.tl_y = coeff*plan->fftc_re.tl_y;
	plan->cut.br_y = coeff*plan->fftc_re.br_y;
}


void raft_getg_with_resources(raft_vector sigma, 
				raft_vector w,
				raft_vector x, 
				fftw_complex *p_out, 
				raft_vector x_zp, 
				raft_vector r_r,
				raft_vector r_i,
				fftw_plan plan 
				) 
{
	// Copy data to zero-padded vector:
	for ( int i( 0 ); i < x_zp.size; ++i  ) {
		raft_vector_element( x_zp, i ) = (i<x.size) ? raft_vector_element(x, i)*raft_vector_element(w, i) : 0.0;
	}

	// Excecutes FFT:
	fftw_execute_dft_r2c( plan, x_zp.p_data, p_out );

	// Filtering in the Fourier domain:
	unsigned osize =  x_zp.size / 2 + 1;

	for( int i=0; i<r_r.size; i++ ) {
		double factor = 1.0/(x_zp.size*raft_vector_element( sigma, i ));
// 		double factor = 1.0/(x_zp.size*raft_vector_element( sigma, i ));
		if(i<osize) {
			raft_vector_element(r_r, i) = p_out[i][0]*factor;
			raft_vector_element(r_i, i) = p_out[i][1]*factor;
		} 
// 		else {
// 			raft_vector_element(r_r, i) = p_out[osize - 1 - (i-osize)][0]*factor;
// 			raft_vector_element(r_i, i) = p_out[osize - 1 - (i-osize)][1]*factor;
// 		}
	}
}

void raft_getg_worker( raft_vector   		sigma,
				raft_vector 	w,
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
						w,
					    raft_matrix_get_column( mat, cur_column ),
					    p_out,
					    x_zp,
					    raft_matrix_get_column( mat_out_r, cur_column ),
					    raft_matrix_get_column( mat_out_i, cur_column ),
					    plan 
					  );
}

void raft_getg( raft_vector sigma,
		raft_vector w,
			raft_matrix data,
			raft_matrix &rmap_r, 
			raft_matrix &rmap_i,
			int nthreads)
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
					      w,
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

void get_hg_ramp(raft_image source, raft_image rmap_r, raft_image rmap_i, int nthreads, double d)
{
	int Ns = source.data.lines;
	raft_vector sigma = raft_vector_create( Ns );
	double d_sigma =  2.0/d;//2*pi/Ns;
	
	raft_vector_element(sigma, 0) = Ns;
	for(int i=1; i<sigma.size; i++) {
		double s = i*d_sigma;
		raft_vector_element(sigma, i) = s; 
	}

	//window
	raft_vector w = raft_vector_create(Ns);
	double beta = 1.8;
	for(int i =0; i<Ns; i++ ) {
		double arg_bessel = beta*sqrt(1 - ((2*i - Ns + 1)/(Ns - 1))*((2*i - Ns + 1)/(Ns - 1)));
		double b1 = boost::math::cyl_bessel_i(0, arg_bessel);
		double b2 = boost::math::cyl_bessel_i(0, beta);
		raft_vector_element(w, i) = fabs(b1)/fabs(b2);
// 		raft_vector_element(w, i) = 1.0;
	}

	raft_getg(sigma, w, source.data, rmap_r.data, rmap_i.data, nthreads);

	raft_vector_destroy(&sigma);
	raft_vector_destroy(&w);
}














double meandc_fft(raft_matrix data)
{
	double res = 0;
	int size = data.columns;
	for(int i=0; i<size; i++) 
		res += raft_matrix_element(data, 0, i);
	res /= size;
	return res;
}

//////////////////////// Backprojection function //////////////////////////////
// #include "g_output.h"
void raft_backprojection_bst(raft_image sino,
			     raft_image res, 
			     raft_bst plan,
			     int nthreads)
{
	
	raft_bst_plan_set_corners( sino, &plan );
	sino2sp_bst(sino, plan.polarsino);
	// Subtracting mean DC component
// 	int i,j;
// 	double meandc = mean_dc(plan.polarsino.data);
// 	printf("MEAN DC BEFORE THE SUBTRACTION: %f\n", meandc);
// 	double smin = raft_matrix_minimum_value(plan.polarsino.data, i,j);
// 	double smax = raft_matrix_maximum_value(plan.polarsino.data, i,j);
// 	double scale = (meandc - smin)/(smax - smin);
// 	subtract_mean_dc(plan.polarsino, meandc);
// 	printf("MEAN DC AFTER THE SUBTRACTION: %f\n", mean_dc(plan.polarsino.data));

	zero_padding_bst(plan.polarsino, plan.zero_padded_sino);
	get_hg_ramp(plan.zero_padded_sino, plan.fftp_re, plan.fftp_im, nthreads, sino.data.lines);
	
// 	printf("fftp_re max: %f, fftp_re min: %f \n", raft_matrix_maximum_value(plan.fftp_re.data, i,j), 
// 			raft_matrix_minimum_value(plan.fftp_re.data, i,j));
	sp2c_miqueles_bst(plan.fftp_re, plan.fftp_im, plan.fftc_re, plan.fftc_im, nthreads);
// 	printf("fftc_re max: %f, fftc_re min: %f \n", raft_matrix_maximum_value(plan.fftc_re.data, i,j), 
// 			raft_matrix_minimum_value(plan.fftc_re.data, i,j));
// 	printf("MEAN DC FFT polar: %f\n", meandc_fft(plan.fftp_re.data));
	ifftw_2d_bst(plan.fftc_re, plan.fftc_im, nthreads);
	fft_shift_2d_bst(plan.fftc_re);
	bl_interpolate_mt_bst(plan.fftc_re, plan.int1, nthreads);

	// Adding mean DC back
// 	double rmin = raft_matrix_minimum_value(plan.int1.data, i,j);
// 	double rmax = raft_matrix_maximum_value(plan.int1.data, i,j);
// 	meandc = scale*(rmax-rmin) + rmin;
// 	printf("MEAN DC BEFORE ADDING: %f\n", mean_dc(plan.int1.data));
// 	add_mean_dc(plan.int1, meandc);
// 	printf("MEAN DC AFTER ADDING: %f\n", mean_dc(plan.fftp_re.data));
// 	printf("int1 max: %f, int1 min: %f \n", raft_matrix_maximum_value(plan.int1.data, i,j), 
	for(int i=0; i<plan.int1.data.lines*plan.int1.data.columns; i++) 
		plan.int1.data.p_data[i] /= plan.polarsino.data.lines;
		
	symmetric_copy_bst(plan.int1, plan.cut);
	bl_interpolate_mt_bst(plan.cut, res, nthreads);

// 	printf("int1 max: %f in point (%d,%d); int1 min: %f in point (%d,%d) \n", raft_matrix_maximum_value(res.data, i,j), i, j, 
// 			raft_matrix_minimum_value(res.data, i,j));
// 
}



