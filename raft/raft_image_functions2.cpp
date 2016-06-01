#include "raft_image_functions2.h"
#include "raft_matrix.h"
#include "raft_image.h"
#include <iostream>
#include <thread>
#include <fftw3.h>
#include <algorithm>
#include <string.h>

double const pi = 3.1415926535897932384626433832795;

int iDivUp_bst(int a, int b)
{
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

//Align a to nearest higher multiple of b
int iAlignUp_bst(int a, int b)
{
    return (a % b != 0) ? (a - a % b + b) : a;
}

int snapTransformSize_bst(int dataSize)
{
    int hiBit;
    unsigned int lowPOT, hiPOT;

    dataSize = iAlignUp_bst(dataSize, 16);

    for (hiBit = 31; hiBit >= 0; hiBit--)
        if (dataSize & (1U << hiBit))
        {
            break;
        }

    lowPOT = 1U << hiBit;

    if (lowPOT == (unsigned int)dataSize)
    {
        return dataSize;
    }

    hiPOT = 1U << (hiBit + 1);

    if (hiPOT <= 1024)
    {
        return hiPOT;
    }
    else
    {
        return iAlignUp_bst(dataSize, 512);
    }
}

void sino2sp_bst(raft_image source, raft_image res)
{
	int t, s, Ns, Nt, Ns2;

	Ns = source.data.lines;
	Nt = source.data.columns;

	memset(res.data.p_data, 0, sizeof(double)*res.data.lines*res.data.columns);
// 	fprintf(stdout,"(1) Ns:%d - Nt:=%d\n",res.data.lines,res.data.columns);
	Ns2 = (int) Ns/2;

	for(t = 0; t < Nt; ++t) {
		for(s = 0; s < Ns2; ++s) {
			raft_matrix_element(res.data, s, t) =    raft_matrix_element(source.data, Ns2 + s, t);
			raft_matrix_element(res.data, s, t+Nt) = raft_matrix_element(source.data, Ns2 - s, t);
		}
	}

	//res.data.lines = (int) Ns/2;
	//res.data.columns = 2*Nt;

// 	fprintf(stdout,"(2) Ns:%d - Nt:=%d => %lf \n",res.data.lines,res.data.columns, raft_matrix_element(res.data, 0,0));
}

void zero_padding_bst(raft_image source, raft_image res)
{
	int Ns = res.data.lines;
	int ns = source.data.lines;
	int nt = source.data.columns;
	if(Ns < ns) {
		return;
	}
	
// 	std::cout<<"Zerro padding on s : ("<<ns<<","<<nt<<" ) -> (";
// 	std::cout<<Ns<<","<<nt<<")."<<std::endl;
	memset(res.data.p_data, 0, sizeof(double)*Ns*nt);
	for(int j=0; j<nt; j++) {
		for(int i=0; i<ns; i++) {
			raft_matrix_element(res.data, i, j) = raft_matrix_element(source.data, i, j);
		}	
	}
}

inline void sp2c_miqueles_worker_bst( 	raft_image source_r, 
					raft_image source_i,
					raft_image res_r, 
					raft_image res_i, 
					double t0, 
					double s0, 
					double x0, 
					double y0,
					double dx, 
					double dy, 
					double ds, 
					double dt, 
					int col1, 
					int col2
 )
{
	int Nt = source_r.data.columns;
	int Ns = source_r.data.lines;
	int Nx = res_r.data.lines;
	int idxt_min, idxt_maj, idxs_min, idxs_maj;
	double s, t, t_min, t_maj, s_min, s_maj, 
	       f11, f12, f21, f22, fr,
	       f11i, f12i, f21i, f22i, fi;
	double y = y0;
	double x = x0;

	double factor = 1.0/(ds*dt);
	for(int j=col1; j<col2; j++) {
		for(int i=0; i<Nx; i++) {
			s = sqrt(x*x + y*y);
			t = acos(x/s);
			t = (y>0) ? t :  - t;

			if(t==pi) {
				t = (x>=0) ? t : -t;
			}

			idxt_min = floor((t-t0)/dt);
			idxt_maj = idxt_min + 1;
			t_min = t0 + dt*idxt_min;
			t_maj = t_min + dt;

			idxs_min = floor((s - s0)/ds);
			idxs_maj = idxs_min + 1;
			s_min = s0 + ds*idxs_min;
			s_maj = s_min + ds;
			
			if(idxt_min == -1) idxt_min = Nt-1;
			if(idxt_maj == Nt) idxt_maj = 0;
			if(idxs_maj > Ns-1 || idxt_maj > Nt-1 ||
					idxs_min < 0 || idxt_min < 0) {
				x += dx;
				continue;
			}

			f11 = raft_matrix_element(source_r.data, idxs_min, idxt_min);
			f12 = raft_matrix_element(source_r.data, idxs_maj, idxt_min);
			f21 = raft_matrix_element(source_r.data, idxs_min, idxt_maj);
			f22 = raft_matrix_element(source_r.data, idxs_maj, idxt_maj);
		
			f11i = raft_matrix_element(source_i.data, idxs_min, idxt_min);
			f12i = raft_matrix_element(source_i.data, idxs_maj, idxt_min);
			f21i = raft_matrix_element(source_i.data, idxs_min, idxt_maj);
			f22i = raft_matrix_element(source_i.data, idxs_maj, idxt_maj);

			fr =  f11*(t_maj - t)*(s_maj - s) + 
					f21*(t - t_min)*(s_maj - s) + 
					f12*(t_maj - t)*(s - s_min) + 
					f22*(t - t_min)*(s - s_min) ;
			fi =  f11i*(t_maj - t)*(s_maj - s) + 
					f21i*(t - t_min)*(s_maj - s) + 
					f12i*(t_maj - t)*(s - s_min) + 
					f22i*(t - t_min)*(s - s_min) ;
			raft_matrix_element(res_r.data, i, j) = fr*factor;
			raft_matrix_element(res_i.data, i, j) = fi*factor;
// 			double b1 = f11 + (f12 - f11)*(s - s_min)/ds;
// 			double b2 = f21 + (f22 - f21)*(s - s_min)/ds;
// 			double fr = b1  + (b2 - b1)*(t - t_min)/dt;
// 			
// 			b1 = f11i + (f12i - f11i)*(s - s_min)/ds;
// 			b2 = f21i + (f22i - f21i)*(s - s_min)/ds;
// 			double fi = b1 + (b2 - b1)*(t - t_min)/dt;
// 			raft_matrix_element(res_r.data, i, j) = fr;//*factor;
// 			raft_matrix_element(res_i.data, i, j) = fi;//*factor;
			x += dx;
		}
		x = x0;
		y += dy;
	}
}

void sp2c_miqueles_bst(raft_image source_r, raft_image source_i,
		raft_image res_r, raft_image res_i, int nthreads)
{
	int Nx = res_r.data.columns;
	int Ny = res_r.data.lines;
	int Nt = source_r.data.columns;
	double t0 = std::min(source_r.br_x, source_r.tl_x);
	double t1 = std::max(source_r.br_x, source_r.tl_x);
	double dt = (t1-t0)/Nt;
	
	int Ns = source_r.data.lines;
	double s0 = std::min(source_r.br_y, source_r.tl_y);
	double s1 = std::max(source_r.br_y, source_r.tl_y);
	double ds = (s1-s0)/Ns;
	
	double x0 = -s1;
	double x1 =  s1;
	double y0 = -s1;
	double y1 =  s1;

	double dx = s1*2/Nx;
	double dy = s1*2/Ny;

	memset(res_r.data.p_data, 0, sizeof(double)*Nx*Ny);
	memset(res_i.data.p_data, 0, sizeof(double)*Nx*Ny);

// 	std::cout<<"___________ semi-polar 2 cartesian interpolation ________"<<std::endl;
// 	std::cout<<"old size: Ns="<<Ns<<"; Nt="<<Nt<<std::endl;
// 	std::cout<<"new size: Nx="<<Nx<<"; Ny="<<Ny<<std::endl;
// 	std::cout<<"s0="<<s0<<"; s1="<<s1<<"; ds="<<ds<<std::endl;
// 	std::cout<<"t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
// 	std::cout<<"x0="<<x0<<"; x1="<<x1<<"; dx="<<dx<<std::endl;
// 	std::cout<<"y0="<<y0<<"; y1="<<y1<<"; dy="<<dy<<std::endl;
// 
// 	res_r.data = raft_matrix_create(Ny, Nx);
// 	res_i.data = raft_matrix_create(Ny, Nx);

	// Make sure we do not have too many or too little threads:
	nthreads = ( nthreads <= Ny ) ? nthreads : Ny;
	nthreads = ( nthreads > 0 ) ? nthreads : 1;
	
	// Base number of columns per thread:
	int base_ncolumns( Ny/ nthreads );
	// Remainder, i.e., number of threads with an extra column:
	int remainder_ncolumns( Ny % nthreads );
	// Current starting_line for worker thread:
	int cur_starting_column( 0 );

	// Create working threads:
	std::vector< std::thread > threads;
	threads.reserve( nthreads );
	int cur_thread = 0;
// 	std::cout<<"ntreads="<< nthreads<<"; base_ncolumns="<<base_ncolumns<<"; remainder_ncolumns="<<remainder_ncolumns<<std::endl;
	for ( ; cur_thread < nthreads; ++cur_thread )
	{
		double y0_cur = y0 + cur_starting_column*dy; 
		int cur_ncolumns( base_ncolumns + ( cur_thread < remainder_ncolumns ) );
// 		std::cout<<"add thr for cur_starting_column="<<cur_starting_column<<"; cur_ncolumns="<<cur_ncolumns<<std::endl;
		std::thread thread =  std::thread( sp2c_miqueles_worker_bst,
					     source_r, 
					     source_i, 
					     res_r,
					     res_i,
					     t0, s0, x0, y0_cur, 
					     dx, dy, ds, dt, 
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

// 	res_r.tl_x = x0;
// 	res_r.tl_y = y1;
// 	res_r.br_x = x1;
// 	res_r.br_y = y0;
// 	res_i.tl_x = x0;
// 	res_i.tl_y = y1;
// 	res_i.br_x = x1;
// 	res_i.br_y = y0;
}

// Semi polar -> log-polar interpolation
inline void sp2lp_worker(raft_image source, 
			raft_image res, 
			double r0, 
			double s0, 
			double s1,
			double dr,
			double ds, 
			int col1, 
			int col2
			)
{
	int Nt = source.data.columns;
	int Ns = source.data.lines;
	int Nr = res.data.lines;

	double r = r0;
	double s, s_min, s_maj, f1, f2, f_res;
	int i_min, i_maj;
	for(int j = col1; j<col2; j++) {
		for(int i=0; i<Nr; ++i) {
			s = exp(r);
			if(s>s1) {
// 				std::cout<<"@!!!!!"<<std::endl;
				r+=dr;
				continue;
			}
			i_min = floor((s - s0)/ds);
			i_maj = i_min + 1; 
			s_min = s0 + i_min*ds;
			s_maj = s_min + ds; 
			
			if(i_maj > Ns-1 || i_min < 0) {
				r += dr;
				continue;
			}
			
			f1 = raft_matrix_element(source.data, i_min, j);
			f2 = raft_matrix_element(source.data, i_maj, j);
		
			f_res = f2*(s - s_min) + f1*(s_maj - s);
			f_res = f_res/ds;
			
			raft_matrix_element(res.data, i, j) = f_res;
			r+=dr;
		}
		r = r0;
	}
}

raft_image sp2lp(raft_image source, raft_image res, double r0, int nthreads)
{
	int Nt = source.data.columns;
	int Ns = source.data.lines;
	int Nr = res.data.lines;

	double s0 = std::min(source.tl_y, source.br_y);
	double s1 = std::max(source.tl_y, source.br_y);
	double ds = (s1 - s0)/Ns; 
	double r1 = log(s1);
	double dr = (r1 - r0)/Nr;

// 	std::cout<<"__________Interpolating from semi-polar to log-polar system________"<<std::endl;
// 	std::cout<<"old sizes: Nt="<<Nt<<"; Ns="<<Ns<<std::endl;
// 	std::cout<<"new sizes: Nt="<<Nt<<"; Nr="<<Nr<<std::endl; 
// 	std::cout<<"s1="<<s1<<"; ds="<<ds<<std::endl;
// 	std::cout<<"r0="<<r0<<"; r1="<<r1;
// 	std::cout<<"; dr="<<dr<<std::endl;
// 
	// Make sure we do not have too many or too little threads:
	nthreads = ( nthreads <= Nt ) ? nthreads : Nt;
	nthreads = ( nthreads > 0 ) ? nthreads : 1;
	
	// Base number of columns per thread:
	int base_ncolumns( Nt/ nthreads );
	// Remainder, i.e., number of threads with an extra column:
	int remainder_ncolumns( Nt % nthreads );
	// Current starting_line for worker thread:
	int cur_starting_column( 0 );

	// Create working threads:
	std::vector< std::thread > threads;
	threads.reserve( nthreads );
	int cur_thread = 0;
	for ( ; cur_thread < nthreads; ++cur_thread )
	{
		int cur_ncolumns( base_ncolumns + ( cur_thread < remainder_ncolumns ) );
		std::thread thread =  std::thread( sp2lp_worker,
					     source, 
					     res,
					     r0, s0, s1, 
					     dr, ds, 
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
}


// Log-Polar -> Cartezian interpolation
inline void lp2c_worker(raft_image source, 
			raft_image res, 
			double t0, 
			double r0, 
			double x0, 
			double y0,
			double dx, 
			double dy, 
			double dr, 
			double dt, 
			int col1, 
			int col2
			)
{
	int Nt = source.data.columns;
	int Nr = source.data.lines;
	int Nx = res.data.lines;
	int idxt_min, idxt_maj, idxr_min, idxr_maj;
	double s, r, t, t_min, t_maj, r_min, r_maj, 
	       f11, f12, f21, f22, factor, fr;


	double x = x0;
	double y = y0;

	factor = 1.0/dt*dr/2;

	for(int j=col1; j<col2; j++) {
		for(int i=0; i<Nx; i++) {
			s = sqrt(x*x + y*y);
			r = log(s);
			t = acos(x/s);
			t = (y>0) ? t :  - t;
			
			idxt_min = floor((t-t0)/dt);
			idxt_maj = idxt_min + 1;
			t_min = t0 + dt*idxt_min;
			t_maj = t_min + dt;

			idxr_min = floor((r - r0)/dr);
			idxr_maj = idxr_min + 1;
			r_min = r0 + dr*idxr_min;
			r_maj = r_min + dr;


			if(idxr_maj > Nr-1 || idxt_maj > Nt ||
					idxr_min < 0 || idxt_min < 0) {
				x += dx;
				continue;
			}

			f11 = raft_matrix_element(source.data, idxr_min, idxt_min);
			f12 = raft_matrix_element(source.data, idxr_maj, idxt_min);
			f21 = raft_matrix_element(source.data, idxr_min, idxt_maj);
			f22 = raft_matrix_element(source.data, idxr_maj, idxt_maj);
		
			fr =  f11*(t_maj - t)*(r_maj - r) + 
					f21*(t - t_min)*(r_maj - r) + 
					f12*(t_maj - t)*(r - r_min) + 
					f22*(t - t_min)*(r - r_min) ;
			fr = fr*factor;
			raft_matrix_element(res.data, i, j) = fr;
			x += dx;
		}
		x = x0;
		y += dy;
	}
}

raft_image lp2c(raft_image source, raft_image res, int nthreads)
{
	int Nt = source.data.columns;
	double t0 = std::min(source.br_x, source.tl_x);
	double t1 = std::max(source.br_x, source.tl_x);
	double dt = (t1-t0)/Nt;
	
	int Nr = source.data.lines;
	double r0 = std::min(source.br_y, source.tl_y);
	double r1 = std::max(source.br_y, source.tl_y);
	double dr = (r1-r0)/Nr;


	int Nx = res.data.lines;
	int Ny = res.data.columns;

	double s1 = exp(r1); 
	double x0 =  -s1;
	double x1 =   s1;
	double y0 =  -s1;
	double y1 =   s1;

	double dx = s1*2/Nx;
	double dy = s1*2/Ny;

// 	std::cout<<"___________ log-polar 2 cartesian interpolation ________"<<std::endl;
// 	std::cout<<"old size: Nr="<<Nr<<"; Nt="<<Nt<<std::endl;
// 	std::cout<<"new size: Nx="<<Nx<<"; Ny="<<Ny<<std::endl;
// 	std::cout<<"r0="<<r0<<"; r1="<<r1<<"; dr="<<dr<<std::endl;
// 	std::cout<<"t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
// 	std::cout<<"x0="<<x0<<"; x1="<<x1<<"; dx="<<dx<<std::endl;
// 	std::cout<<"y0="<<y0<<"; y1="<<y1<<"; dy="<<dy<<std::endl;

	// Make sure we do not have too many or too little threads:
	nthreads = ( nthreads <= Ny ) ? nthreads : Ny;
	nthreads = ( nthreads > 0 ) ? nthreads : 1;
	
	// Base number of columns per thread:
	int base_ncolumns( Ny/ nthreads );
	// Remainder, i.e., number of threads with an extra column:
	int remainder_ncolumns( Ny % nthreads );
	// Current starting_line for worker thread:
	int cur_starting_column( 0 );

	// Create working threads:
	std::vector< std::thread > threads;
	threads.reserve( nthreads );
	int cur_thread = 0;
	for ( ; cur_thread < nthreads; ++cur_thread )
	{
		double y0_cur = y0 + cur_starting_column*dy; 
		int cur_ncolumns( base_ncolumns + ( cur_thread < remainder_ncolumns ) );
		std::thread thread =  std::thread( lp2c_worker,
					     source, 
					     res,
					     t0, r0, x0, y0_cur, 
					     dx, dy, dr, dt, 
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

	int r = 2;
	for(int i=0; i<Nx; i++) {
		for(int j=0; j<Ny; j++) {
			if((i-Nx/2)*(i-Nx/2) + (j-Ny/2)*(j-Ny/2)< r*r ) {
				raft_matrix_element(res.data, i,j) =  
					raft_matrix_element(res.data, Nx/2-r, Ny/2-r);
			}
		}
	}
}

void convolution_2d(raft_image source, raft_image kernel, raft_image res, int nthreads) 
{
	int Nx = source.data.lines;
	int Ny = source.data.columns;
	if(Nx != kernel.data.lines || Ny != kernel.data.columns ) {
		std::cout<<"Size error: source = <"<<Nx<<","<<Ny<<">";
		std::cout<<"; kernel = <"<<kernel.data.lines<<","<<kernel.data.columns<<">"<<std::endl;
		return;
		
	}
// 	std::cout<<"conv_2d starts: Mx="<<Nx<<"; Ny="<<Ny<<std::endl;
	int data_size = Nx*Ny;
	int spectra_size = Nx*(Ny/2+1);

	double scale = 1.0/data_size;
	fftw_complex *spectrum_data = 	(fftw_complex*) fftw_malloc(spectra_size*sizeof(fftw_complex));
	fftw_complex *spectrum_kernel = (fftw_complex*) fftw_malloc(spectra_size*sizeof(fftw_complex));

	fftw_init_threads();
	fftw_plan_with_nthreads(nthreads);

	fftw_plan p = fftw_plan_dft_r2c_2d(Nx, Ny, source.data.p_data, spectrum_data, FFTW_ESTIMATE);
	fftw_plan pinv = fftw_plan_dft_c2r_2d(Nx, Ny, spectrum_data, res.data.p_data, FFTW_ESTIMATE);

	fftw_execute_dft_r2c(p, source.data.p_data, spectrum_data);
	fftw_execute_dft_r2c(p, kernel.data.p_data, spectrum_kernel);

// 	std::cout<<"conv2d: straight FFT(x) succeed"<<std::endl;
	double xr, xi, kr, ki;
	for(int i=0; i<spectra_size; i++) {
		xr = spectrum_data[i][0];
		xi = spectrum_data[i][1];

		kr = spectrum_kernel[i][0];
		ki = spectrum_kernel[i][1];
			
		spectrum_data[i][0] = 	(xr*kr - xi*ki)*scale;
		spectrum_data[i][1] = 	(xr*ki + xi*kr)*scale;
	}
// 	std::cout<<"conv2d: Multiplication succeed"<<std::endl;
	fftw_execute_dft_c2r(pinv, spectrum_data, res.data.p_data);
// 	std::cout<<"conv2d: inverse fourier procedure succeed"<<std::endl;
	fftw_destroy_plan(p);
	fftw_destroy_plan(pinv);
	fftw_free(spectrum_data);
	fftw_free(spectrum_kernel);
	fftw_cleanup_threads();

// 	std::cout<<"conv2d: SUCCEED"<<std::endl;
}




















void ifftw_2d_bst(raft_image xr, raft_image xi, int threads)
{
// 	std::cout<<"iFFTW 2D started"<<std::endl;
	int N = xr.data.lines;
	int M = xr.data.columns;
	fftw_complex *in;
	in =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*M);
	int idx; 
	for(int j=0; j<M; j++) {
		for(int i=0; i<N; i++) {
			idx = j*N + i;
			in[idx][0] = raft_matrix_element(xr.data, i, j);	
			in[idx][1] = raft_matrix_element(xi.data, i, j);	

		}
	}
	fftw_init_threads();
	fftw_plan_with_nthreads(threads);
	fftw_plan p = fftw_plan_dft_2d(N, M, in, in, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p); 
	fftw_destroy_plan(p);
	fftw_cleanup_threads();
// 	double factor = 1.0/(M*N);
	double factor = 1.0;
	for(int j=0; j<M; j++) {
		for(int i=0; i<N; i++) {
			idx = j*N + i;
			raft_matrix_element(xr.data, i, j) = in[idx][0]*factor;
			raft_matrix_element(xi.data, i, j) = in[idx][1]*factor;
// 			raft_matrix_element(xr, i, j) = out[idx][0]/(2*pi);
// 			raft_matrix_element(xi, i, j) = out[idx][1]/(2*pi);
		}
	}
// 	std::cout<<"iFFT_2D succeed"<<std::endl;
}

void fft_shift_2d_bst(raft_image x)
{
	int Nx = x.data.lines;
	int Ny = x.data.columns;
	double tmp;
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny/2;j++){
			tmp = raft_matrix_element(x.data, i,j);
			raft_matrix_element(x.data, i,j) = raft_matrix_element(x.data, i,Ny/2+j);
			raft_matrix_element(x.data, i,Ny/2+j) = tmp;
		}
	}

	for(int i=0;i<Nx/2;i++){
		for(int j=0;j<Ny;j++){
			tmp = raft_matrix_element(x.data, i,j);
			raft_matrix_element(x.data, i,j) = raft_matrix_element(x.data, i+Nx/2,j);
			raft_matrix_element(x.data, i+Nx/2,j) = tmp;
		}
	}
}

void ifftw_2d_C2R_bst(raft_image source_r, raft_image source_i, raft_image res, int nthreads)
{
	int old_Nx = source_r.data.columns;
	int old_Ny = source_r.data.lines;

	double scale = 1.0/(old_Nx*old_Ny);
	int Nx = old_Nx;
	int Ny = 2*old_Ny-1;

	fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex*)*old_Nx*old_Ny);

	for(int j=0; j<old_Nx; j++) {
		for(int i=0; i<old_Ny; i++) {
			int idx = j*old_Nx + i;
			in[idx][0] = (1.0/2*pi)*raft_matrix_element(source_r.data, i, j);	
			in[idx][1] = (1.0/2*pi)*raft_matrix_element(source_i.data, i, j);	

		}
	}


	fftw_init_threads();
	fftw_plan_with_nthreads(nthreads);

	fftw_plan pinv = fftw_plan_dft_c2r_2d(Nx, Ny, in, res.data.p_data, FFTW_ESTIMATE);
	fftw_execute_dft_c2r(pinv, in, res.data.p_data);
	fftw_destroy_plan(pinv);
	fftw_free(in);
	
}

// Bilinear interpolation to any sizes Nx, Ny
inline void bl_interpolate_worker_bst(raft_image source, 
			raft_image res, 
			double x0, 
			double y0,
			double y0_cur,
			double old_dx,
			double old_dy,
			double dx,
			double dy, 
			int col1, 
			int col2
			)
{
	int old_Nx = source.data.lines;
	int old_Ny = source.data.columns;
	int Nx = res.data.lines;

	double x = x0;
	double y = y0_cur; 

	int i_min, i_maj, j_min, j_maj;
	double x_min, x_maj, y_min, y_maj, f_res, f11, f12, f21, f22, tmp;
	tmp = old_dx*old_dy;
	for(int j=col1; j<=col2; ++j) {
		for(int i=0; i<Nx; ++i) {
			i_min = floor((x - x0)/old_dx);
			j_min = floor((y - y0)/old_dy);
			i_maj = i_min + 1;
			j_maj = j_min + 1;

			x_min = x0 + i_min*old_dx;
			y_min = y0 + j_min*old_dy;

			x_maj = x_min+old_dx;
			y_maj = y_min+old_dy;

			if((j_maj > old_Ny-1) || (j_min < 0) || 
					(i_maj > old_Nx - 1) || (i_min<0)) {
				x+=dx;
				continue;
			}

			int i1=i_min, i2=i_maj, j1=j_min, j2=j_maj;
			j1 = (j_min < 0 ) ? 0 : j_min;
			j2 = (j_maj > old_Ny-1 ) ? old_Ny-1 : j_maj;
			i1 = (i_min < 0 ) ? 0 : i_min;
			i2 = (i_maj > old_Nx-1 ) ? old_Nx-1 : i_maj;

			if(j_maj > old_Ny - 1) j_maj=old_Ny-1;
			if(j_min < 0 ) j_min=0;
			if(i_maj > old_Nx - 1) i_maj = old_Nx - 1;
			if(i_min < 0 ) i_min = 0;

			f11 = raft_matrix_element(source.data, i1, j1);
			f12 = raft_matrix_element(source.data, i1, j2);
			f21 = raft_matrix_element(source.data, i2, j1);
			f22 = raft_matrix_element(source.data, i2, j2);

// 			if(x-x_min == 0) printf("x-x_min = 0; i=%d, j=%d, i_min=%d, i_max=%d, x_min=%f, x=%f, x_maj=%f \n",
// 				       i, j, i_min, i_maj, x_min, x, x_maj);
// 			if(y-y_min == 0) printf("y-y_min = 0; i=%d, j=%d, i_min=%d, i_may=%d, y_min=%f, y=%f, y_maj=%f \n",
// 				       i, j, i_min, i_maj, y_min, y, y_maj);
// 			if(x_maj-x == 0) printf("x_maj-x = 0; i=%d, j=%d, i_min=%d, i_max=%d, x_min=%f, x=%f, x_maj=%f \n",
// 				       i, j, i_min, i_maj, x_min, x, x_maj);
// 			if(y_maj-y == 0) printf("y_maj-y = 0; i=%d, j=%d, i_min=%d, i_may=%d, y_min=%f, y=%f, y_maj=%f \n",
// 					i, j, i_min, i_maj, y_min, y, y_maj);
			f_res =  f11*(x_maj - x)*(y_maj - y) + 
					f21*(x - x_min)*(y_maj - y) + 
					f12*(x_maj - x)*(y - y_min) + 
					f22*(x - x_min)*(y - y_min) ;
			f_res = f_res/tmp;
			
			raft_matrix_element(res.data, i, j) = f_res;
			x+=dx;
		}
		x = x0;
		y+=dy;
	}
}

void bl_interpolate_mt_bst(raft_image source, raft_image res, int nthreads)
{
	int old_Nx = source.data.lines;
	int old_Ny = source.data.columns;
	int Nx = res.data.lines;
	int Ny = res.data.columns;

	memset(res.data.p_data, 0, sizeof(double)*Nx*Ny);
	double x0 =std::min(source.tl_x, source.br_x); 
	double y0 =std::min(source.tl_y, source.br_y); 
	double x1 =std::max(source.tl_x, source.br_x); 
	double y1 =std::max(source.tl_y, source.br_y); 
	double old_dx = (x1-x0)/old_Nx;
	double old_dy = (y1-y0)/old_Ny;
	double dx = (x1-x0)/Nx;
	double dy = (y1-y0)/Ny;
// 	raft_image res;
// 	res.data = raft_matrix_create(Nx, Ny);

// 	std::cout<<"Multithread bilinear interpolation: resizing ("<<old_Nx<<","<<old_Ny<<") -> (";
// 	std::cout<<Nx<<","<<Ny<<")."<<std::endl;
// 
	// Make sure we do not have too many or too little threads:
	nthreads = ( nthreads <= Ny ) ? nthreads : Ny;
	nthreads = ( nthreads > 0 ) ? nthreads : 1;
	
	// Base number of columns per thread:
	int base_ncolumns( Ny/ nthreads );
	// Remainder, i.e., number of threads with an extra column:
	int remainder_ncolumns( Ny % nthreads );
	// Current starting_line for worker thread:
	int cur_starting_column( 0 );

	// Create working threads:
	std::vector< std::thread > threads;
	threads.reserve( nthreads );
	int cur_thread = 0;
	for ( ; cur_thread < nthreads; ++cur_thread )
	{
		double y0_cur = y0 + cur_starting_column*dy; 
		int cur_ncolumns( base_ncolumns + ( cur_thread < remainder_ncolumns ) );
		std::thread thread =  std::thread( bl_interpolate_worker_bst,
					     source, 
					     res,
					     x0, y0, y0_cur, 
					     old_dx, old_dy, dx, dy, 
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


// 	std::cout<<"________bl_interpolation succeed___________"<<std::endl;
}

void symmetric_copy_bst(raft_image source, raft_image res)
{
	int old_nx = source.data.lines;
	int old_ny = source.data.columns;
	
	
	int nx = res.data.lines;
	int ny = res.data.columns;
// 	std::cout<<"Symmetric copy: ("<<old_nx<<","<<old_ny<<") -> ("<<nx<<","<<ny<<")"<<std::endl;


	memset(res.data.p_data, 0, sizeof(double)*nx*ny);
	int idx_x_start = (old_nx - nx)/2;
	int idx_y_start = (old_ny - ny)/2;


// 	std::cout<<"Symmetric copy: idx_x="<<idx_x_start<<" ,idx_y="<<idx_y_start<<std::endl;
	int idx = 0;
	if(idx_x_start >= 0 && idx_y_start >= 0) {
		for(int j=0; j<ny; j++) {
			for(int i=0; i<nx; i++) {
				if(idx_x_start + i < 0 || idx_y_start < 0) {
					raft_matrix_element(res.data, i, j) = 0; 	
					continue;
				}
				raft_matrix_element(res.data, i, j) = 
					raft_matrix_element(source.data, i+ idx_x_start, 
							j+idx_y_start);
			}
		
		}
	} else if(idx_x_start < 0 && idx_y_start >= 0) {
		for(int j=0; j<ny; j++) {
			for(int i=0; i<old_nx; i++) {
				raft_matrix_element(res.data, i - idx_x_start, j) = 
					raft_matrix_element(source.data, i, 
							j+idx_y_start);
			}
		
		}
		
	} else if(idx_x_start >= 0 && idx_y_start < 0) {
		for(int j=0; j<old_ny; j++) {
			for(int i=0; i<nx; i++) {
				raft_matrix_element(res.data, i, j-idx_y_start) = 
					raft_matrix_element(source.data, i + idx_x_start, 
							j);
			}
		
		}
		
	} else if(idx_x_start < 0 && idx_y_start < 0) {
		printf("case3");
		for(int j=0; j<old_ny; j++) {
			for(int i=0; i<old_nx; i++) {
				raft_matrix_element(res.data, i-idx_x_start, j-idx_y_start) = 
					raft_matrix_element(source.data, i, j);
			}
		
		}
		
	}
}





















double mean_dc(raft_matrix data)
{
	double res = 0;
// 	int Ns = data.lines;
// 	int Nt = data.columns;
// 	for(int i=0; i<Ns; i++) {
// 		res += raft_matrix_element(data, i, 0);
// 	}
// 	return res/Ns;
	int size = data.columns*data.lines;
	for(int j=0; j<size; j++) {
		res += data.p_data[j]; 
	}
	return res/(size);

}

void subtract_mean_dc(raft_image img, double meandc)
{
	int size = img.data.columns*img.data.lines;
	for(int i=0; i<size; i++) {
		img.data.p_data[i] -= meandc;
	}
}

void add_mean_dc(raft_image img, double meandc)
{
	int size = img.data.columns*img.data.lines;
	for(int i=0; i<size; i++) {
		img.data.p_data[i] += meandc;
	}
	
}

double raft_matrix_maximum_value(raft_matrix data, int &i_max, int &j_max)
{
	int Ns = data.lines;
	int Nt = data.columns;
	double res = fabs(raft_matrix_element(data, 0,0)); 
	double el;
	for(int j=0; j<Nt; j++) {
		for(int i=0; i<Ns; i++) {
			el = fabs(raft_matrix_element(data, i, j));
			if(el > res) {
				res = el;
				i_max = i;
				j_max = j;
			}
		}
	}
	return res;
}	

double raft_matrix_minimum_value(raft_matrix data, int &i_min, int &j_min)
{
	int Ns = data.lines;
	int Nt = data.columns;
	double res = raft_matrix_element(data, 0,0); 
	double el;
	for(int j=0; j<Nt; j++) {
		for(int i=0; i<Ns; i++) {
			el = raft_matrix_element(data, i, j);
			if(el < res) {
				res = el;
				i_min = i;
				j_min = j;
			}
		}
	}
	return res;
}	




