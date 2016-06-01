#include "raft_functions_mt.h"
#include "raft_image_functions.h"
#include <fftw3.h>
#include <algorithm>
#include <iostream>
#include <thread>
// #include "g_output.h"


// Cartesian -> Semi-polar interpolation
inline void c2sp_mt_worker( 		raft_image source, 
					raft_image res, 
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
	int Ns = res.data.lines;
	int Nx = source.data.lines;
	int Ny = source.data.columns;
	
	double x = x0; 
	double y = y0;
	double t = t0;
	double s = s0;

	int i_min, i_maj, j_min, j_maj;
	double x_min, x_maj, y_min, y_maj, f11, f12, f21, f22, f, dv;

	for(int j=col1; j<col2; ++j) {
		double _cos = cos(t);
		double _sin = sin(t);
		for(int i=0; i<Ns; ++i) {
			x = s*_cos; 
			y = s*_sin;
			y = (fabs(y)>=dy) ? y : 0;

			j_min = floor((x-x0)/dx);
			j_maj = j_min + 1;
			i_min = (fabs(y)>=dy) ? floor((y-y0)/dy) : Ny/2-1;
			i_maj = i_min + 1;
			
			x_min = x0 + j_min*dx;
			x_maj = x_min + dx;
			y_min = y0 + i_maj*dy;
			y_maj = y_min + dy;
		
			if((j_maj > Nx-1) || (j_min < 0) || (i_maj > Ny-1) || (i_min < 0)) {
				s += ds;
				continue;
			}
		
			f11 = raft_matrix_element(source.data, i_min, j_min);
			f12 = raft_matrix_element(source.data, i_min, j_maj);
			f21 = raft_matrix_element(source.data, i_maj, j_min);
			f22 = raft_matrix_element(source.data, i_maj, j_maj);
			
			dv = dx*dy;
// 			f =  f11*(x_maj - x)*(y_maj - y) + 
// 					f21*(x - x_min)*(y_maj - y) + 
// 					f12*(x_maj - x)*(y - y_min) + 
// 					f22*(x - x_min)*(y - y_min) ;

			f = (f11+f12+f21+f22)/4; //FIXME!!! Negative values with bilinear approximation
			raft_matrix_element(res.data, i, j) = f;
			s+=ds;
		}
		s = s0;
		t += dt;
	}	

}

void c2sp_mt(raft_image source, raft_image &res, int Nt, int Ns, int nthreads)
{
	int Nx = source.data.columns;
	double x0 = std::min(source.br_x, source.tl_x);
	double x1 = std::max(source.br_x, source.tl_x);
	double dx = (x1-x0)/Nx;
	
	int Ny = source.data.lines;
	double y0 = std::min(source.br_y, source.tl_y);
	double y1 = std::max(source.br_y, source.tl_y);
	double dy   = (y1-y0)/Ny;

	double s0 = 0;
	double s1 = sqrt(x1*x1 + y1*y1); //FIXME
	double ds = s1/Ns;

	double t0 = -pi;
	double t1 = pi;
	double dt = 2*pi/Nt;

// 	std::cout<<"___________ cartesian 2 semi-polar interpolation ________"<<std::endl;
// 	std::cout<<"old size: Nx="<<Nx<<"; Ny="<<Ny<<std::endl;
// 	std::cout<<"new size: Nt="<<Nt<<"; Ns="<<Ns<<std::endl;
// 	std::cout<<"x0="<<x0<<"; x1="<<x1<<"; dx="<<dx<<std::endl;
// 	std::cout<<"y0="<<y0<<"; y1="<<y1<<"; dy="<<dy<<std::endl;
// 	std::cout<<"s0="<<s0<<"; s1="<<s1<<"; ds="<<ds<<std::endl;
// 	std::cout<<"t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
// 
	// Create the outpout image
	res.data = raft_matrix_create(Ns, Nt);
	res.tl_x = t0;
	res.br_x = t1;
	res.br_y = s0;
	res.tl_y = s1;
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
		double t0_cur = t0 + cur_starting_column*dt;
		int cur_ncolumns( base_ncolumns + ( cur_thread < remainder_ncolumns ) );
		std::cout<<"adding thread for t0="<<t0_cur<<"; t1="<<t0_cur + cur_ncolumns*dt<<std::endl;
		std::thread thread =  std::thread( c2sp_mt_worker,
					     source, res,
					     t0_cur, s0, x0, y0, 
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


}

// Cartesian -> log-polar interpolation
inline void c2lp_mt_worker( 		raft_image source, 
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
	int Nr = res.data.lines;
	int Nx = source.data.lines;
	int Ny = source.data.columns;
	
	double x = x0; 
	double y = y0;
	double t = t0;
	double r = r0;
	double s;

	int i_min, i_maj, j_min, j_maj;
	double x_min, x_maj, y_min, y_maj, f11, f12, f21, f22, f, dv;

	for(int j=col1; j<col2; ++j) {
		double _cos = cos(t);
		double _sin = sin(t);
		for(int i=0; i<Nr; ++i) {
			s = exp(r);
			x = s*_cos; 
			y = s*_sin;
			y = (fabs(y)>=dy) ? y : 0;

			j_min = floor((x-x0)/dx);
			j_maj = j_min + 1;
			i_min = (fabs(y)>=dy) ? floor((y-y0)/dy) : Ny/2-1;
			i_maj = i_min + 1;
			
			x_min = x0 + j_min*dx;
			x_maj = x_min + dx;
			y_min = y0 + i_maj*dy;
			y_maj = y_min + dy;
	
			if((j_maj > Nx-1) || (j_min < 0) || (i_maj > Ny-1) || (i_min < 0)) {
				r += dr;
				continue;
			}
		
			f11 = raft_matrix_element(source.data, i_min, j_min);
			f12 = raft_matrix_element(source.data, i_min, j_maj);
			f21 = raft_matrix_element(source.data, i_maj, j_min);
			f22 = raft_matrix_element(source.data, i_maj, j_maj);
			
			dv = dx*dy;

			f = (f11+f12+f21+f22)/4; //FIXME!!! Negative values with bilinear approximation
			raft_matrix_element(res.data, i, j) = f;
			r+=dr;
		}
		r = r0;
		t += dt;
	}	

}

void c2lp_mt(raft_image source, raft_image &res, int Nt, int Nr, double r0, int nthreads)
{
	int Nx = source.data.columns;
	double x0 = std::min(source.br_x, source.tl_x);
	double x1 = std::max(source.br_x, source.tl_x);
	double dx = (x1-x0)/Nx;
	
	int Ny = source.data.lines;
	double y0 = std::min(source.br_y, source.tl_y);
	double y1 = std::max(source.br_y, source.tl_y);
	double dy = (y1-y0)/Ny;
	double r1 = 0; //FIXME
// 	double r1 = log(std::max(source.tl_y, source.br_y));
// 	double r1 = log(sqrt(x1*x1 + y1*y1)); //FIXME
	double dr = (r1-r0)/Nr;

	double t0 = -pi;
	double t1 = pi;
	double dt = 2*pi/Nt;

// 	std::cout<<"___________ cartesian 2 log-polar interpolation ________"<<std::endl;
// 	std::cout<<"old size: Nx="<<Nx<<"; Ny="<<Ny<<std::endl;
// 	std::cout<<"new size: Nt="<<Nt<<"; Nr="<<Nr<<std::endl;
// 	std::cout<<"x0="<<x0<<"; x1="<<x1<<"; dx="<<dx<<std::endl;
// 	std::cout<<"y0="<<y0<<"; y1="<<y1<<"; dy="<<dy<<std::endl;
// 	std::cout<<"r0="<<r0<<"; r1="<<r1<<"; dr="<<dr<<std::endl;
// 	std::cout<<"t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
// 
	// Create the outpout image
	res.data = raft_matrix_create(Nr, Nt);
	res.tl_x = t0;
	res.br_x = t1;
	res.br_y = r0;
	res.tl_y = r1;
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
		double t0_cur = t0 + cur_starting_column*dt;
		int cur_ncolumns( base_ncolumns + ( cur_thread < remainder_ncolumns ) );
		std::cout<<"adding thread for t0="<<t0_cur<<"; t1="<<t0_cur + cur_ncolumns*dt<<std::endl;
		std::thread thread =  std::thread( c2lp_mt_worker,
					     source, res,
					     t0_cur, r0, x0, y0, 
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


}

// Semi-polar -> Cartesian interpolation
// Simple interpolation
inline void sp2c_mt_worker( 		raft_image source, 
					raft_image res, 
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
	int Nt = source.data.columns;
	int Ns = source.data.lines;
	int Nx = res.data.lines;
	int idxt_min, idxt_maj, idxs_min, idxs_maj;
	double s, t, t_min, t_maj, s_min, s_maj, 
	       f11, f12, f21, f22, f;
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

			f11 = raft_matrix_element(source.data, idxs_min, idxt_min);
			f12 = raft_matrix_element(source.data, idxs_maj, idxt_min);
			f21 = raft_matrix_element(source.data, idxs_min, idxt_maj);
			f22 = raft_matrix_element(source.data, idxs_maj, idxt_maj);
		
			f =  f11*(t_maj - t)*(s_maj - s) + 
					f21*(t - t_min)*(s_maj - s) + 
					f12*(t_maj - t)*(s - s_min) + 
					f22*(t - t_min)*(s - s_min) ;

			f *= factor;
			raft_matrix_element(res.data, j, i) = f;
			x += dx;
		}
		x = x0;
		y += dy;
	}
}

void sp2c_mt(	raft_image source, 
		raft_image &res, 
		int Nx, 
		int Ny, 
		int nthreads)
{
	// Defining sizes of input array and steps
	int Nt = source.data.columns;
	double t0 = std::min(source.br_x, source.tl_x);
	double t1 = std::max(source.br_x, source.tl_x);
	double dt = (t1-t0)/Nt;
	
	int Ns = source.data.lines;
	double s0 = std::min(source.br_y, source.tl_y);
	double s1 = std::max(source.br_y, source.tl_y);
	double ds = (s1-s0)/Ns;

	// Defining borders of output and steps
	double x0 = -s1;
	double x1 =  s1;
	double y0 = -s1;
	double y1 =  s1;
	double dx = s1*2/Nx;
	double dy = s1*2/Ny;

// 	std::cout<<"___________ semi-polar 2 cartesian interpolation ________"<<std::endl;
// 	std::cout<<"old size: Ns="<<Ns<<"; Nt="<<Nt<<std::endl;
// 	std::cout<<"new size: Nx="<<Nx<<"; Ny="<<Ny<<std::endl;
// 	std::cout<<"s0="<<s0<<"; s1="<<s1<<"; ds="<<ds<<std::endl;
// 	std::cout<<"t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
// 	std::cout<<"x0="<<x0<<"; x1="<<x1<<"; dx="<<dx<<std::endl;
// 	std::cout<<"y0="<<y0<<"; y1="<<y1<<"; dy="<<dy<<std::endl;

	// Creating output image
	res.data = raft_matrix_create(Nx, Ny);
	res.tl_x = x0;
	res.tl_y = y1;
	res.br_x = x1;
	res.br_y = y0;

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
	std::cout<<"ntreads="<< nthreads<<"; base_ncolumns="<<base_ncolumns<<"; remainder_ncolumns="<<remainder_ncolumns<<std::endl;
	for ( ; cur_thread < nthreads; ++cur_thread )
	{
		double y0_cur = y0 + cur_starting_column*dy; 
		int cur_ncolumns( base_ncolumns + ( cur_thread < remainder_ncolumns ) );
		std::cout<<"add thr for cur_starting_column="<<cur_starting_column<<"; cur_ncolumns="<<cur_ncolumns<<std::endl;
		std::thread thread =  std::thread( sp2c_mt_worker,
					     source, 
					     res,
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

}

// The same, but for using two arrays for real and imaginary part
// A  little bit faster for the algorithm of BST
inline void sp2c_miqueles_worker( 	raft_image source_r, 
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
	double factor = 1.0;//(ds*dt);
	
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

// 			fr = (f11+f12+f21+f22)/4;
// 			fi = (f11i+f12i+f21i+f22i)/4;
			fr =  f11*(t_maj - t)*(s_maj - s) + 
					f21*(t - t_min)*(s_maj - s) + 
					f12*(t_maj - t)*(s - s_min) + 
					f22*(t - t_min)*(s - s_min) ;
			fi =  f11i*(t_maj - t)*(s_maj - s) + 
					f21i*(t - t_min)*(s_maj - s) + 
					f12i*(t_maj - t)*(s - s_min) + 
					f22i*(t - t_min)*(s - s_min) ;
			raft_matrix_element(res_r.data, i, j) = fr;
			raft_matrix_element(res_i.data, i, j) = fi;
			x += dx;
		}
		x = x0;
		y += dy;
	}
}

void sp2c_miqueles(raft_image source_r, raft_image source_i,
		raft_image &res_r, raft_image &res_i, int Nx, int Ny, int nthreads)
{
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

// 	std::cout<<"___________ semi-polar 2 cartesian interpolation ________"<<std::endl;
// 	std::cout<<"old size: Ns="<<Ns<<"; Nt="<<Nt<<std::endl;
// 	std::cout<<"new size: Nx="<<Nx<<"; Ny="<<Ny<<std::endl;
// 	std::cout<<"s0="<<s0<<"; s1="<<s1<<"; ds="<<ds<<std::endl;
// 	std::cout<<"t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
// 	std::cout<<"x0="<<x0<<"; x1="<<x1<<"; dx="<<dx<<std::endl;
// 	std::cout<<"y0="<<y0<<"; y1="<<y1<<"; dy="<<dy<<std::endl;
// 
	res_r.data = raft_matrix_create(Ny, Nx);
	res_i.data = raft_matrix_create(Ny, Nx);

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
		std::thread thread =  std::thread( sp2c_miqueles_worker,
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

	res_r.tl_x = x0;
	res_r.tl_y = y1;
	res_r.br_x = x1;
	res_r.br_y = y0;
	res_i.tl_x = x0;
	res_i.tl_y = y1;
	res_i.br_x = x1;
	res_i.br_y = y0;
}
//////////////////////////////////////////
void sp2c_miqueles2(raft_image source_r, raft_image source_i,
		raft_image &res_r, raft_image &res_i, int nthreads)
{
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

	// Calculating new size...
	int Nx = 2*Ns;
	int Ny = 2*Ns;
	double dx = s1*2/Nx;
	double dy = s1*2/Ny;

// 	std::cout<<"___________ semi-polar 2 cartesian interpolation ________"<<std::endl;
// 	std::cout<<"old size: Ns="<<Ns<<"; Nt="<<Nt<<std::endl;
// 	std::cout<<"new size: Nx="<<Nx<<"; Ny="<<Ny<<std::endl;
// 	std::cout<<"s0="<<s0<<"; s1="<<s1<<"; ds="<<ds<<std::endl;
// 	std::cout<<"t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
// 	std::cout<<"x0="<<x0<<"; x1="<<x1<<"; dx="<<dx<<std::endl;
// 	std::cout<<"y0="<<y0<<"; y1="<<y1<<"; dy="<<dy<<std::endl;
//



	res_r.data = raft_matrix_create(Ny, Nx);
	res_i.data = raft_matrix_create(Ny, Nx);

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
		std::thread thread =  std::thread( sp2c_miqueles_worker,
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

	res_r.tl_x = x0;
	res_r.tl_y = y1;
	res_r.br_x = x1;
	res_r.br_y = y0;
	res_i.tl_x = x0;
	res_i.tl_y = y1;
	res_i.br_x = x1;
	res_i.br_y = y0;
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

raft_image lp2c_mt(raft_image source, int Nx, int Ny, int nthreads)
{
	int Nt = source.data.columns;
	double t0 = std::min(source.br_x, source.tl_x);
	double t1 = std::max(source.br_x, source.tl_x);
	double dt = (t1-t0)/Nt;
	
	int Nr = source.data.lines;
	double r0 = std::min(source.br_y, source.tl_y);
	double r1 = std::max(source.br_y, source.tl_y);
	double dr = (r1-r0)/Nr;

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


	raft_image res;
	res.data = raft_matrix_create(Ny, Nx);
	
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
//

	res.tl_x = x0;
	res.tl_y = y1;
	res.br_x = x1;
	res.br_y = y0;

	return res;
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

raft_image sp2lp_mt(raft_image source, double Nr, double r0, int nthreads)
{
	int Nt = source.data.columns;
	int Ns = source.data.lines;
	raft_image res = source;
	res.data = raft_matrix_create(Nr, Nt);
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
	res.tl_y = r0;
	res.br_y = r1;
// 	std::cout<<"interpolation p2lp succeed"<<std::endl;
	return res; 
}

// Log-polar -> semi polar interpolation
inline void lp2sp_worker(raft_image source, 
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
	int Nr = source.data.lines;
	int Ns = res.data.lines;

	double s = s0;
	double r, r_min, r_maj, f1, f2, f_res;
	int i_min, i_maj;
	for(int j = col1; j<col2; j++) {
		for(int i=0; i<Nr; ++i) {
			r = log(s); 
			i_min = floor((r - r0)/dr);
			i_maj = i_min + 1; 
			r_min = r0 + i_min*dr;
			r_maj = r_min + dr; 
			
			if(i_maj > Nr-1 || i_min < 0) {
				s += ds;
				continue;
			}
			
			f1 = raft_matrix_element(source.data, i_min, j);
			f2 = raft_matrix_element(source.data, i_maj, j);
		
			f_res = f2*(r - r_min) + f1*(r_maj - r);
			
			raft_matrix_element(res.data, i, j) = f_res;
			s+=ds;
		}
		s = s0;
	}
}

raft_image lp2sp_mt(raft_image source, double Ns, int nthreads)
{
	int Nt = source.data.columns;
	int Nr = source.data.lines;
	raft_image res = source;
	res.data = raft_matrix_create(Ns, Nt);
	// Old mesh
	double r0 = std::min(source.tl_y, source.br_y);
	double r1 = std::max(source.tl_y, source.br_y);
	double dr = (r1 - r0)/Nr;
	// New mesh
	double s1 = exp(r1);
	double s0 = exp(r0);
	double ds = (s1 - s0)/Ns;

//	std::cout<<"__________Interpolating from log-polar to semi-polar system________"<<std::endl;
//	std::cout<<"old sizes: Nt="<<Nt<<"; Ns="<<Ns<<std::endl;
//	std::cout<<"new sizes: Nt="<<Nt<<"; Nr="<<Nr<<std::endl; 
//	std::cout<<"s1="<<s1<<"; ds="<<ds<<std::endl;
//	std::cout<<"r0="<<r0<<"; r1="<<r1;
//	std::cout<<"; dr="<<dr<<std::endl;

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
		std::thread thread =  std::thread( lp2sp_worker,
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
	res.tl_y = s0;
	res.br_y = s1;
// 	std::cout<<"interpolation p2lp succeed"<<std::endl;
	return res; 
}

// Bilinear interpolation to any sizes Nx, Ny
inline void bl_interpolate_worker(raft_image source, 
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

			f11 = raft_matrix_element(source.data, i_min, j_min);
			f12 = raft_matrix_element(source.data, i_min, j_maj);
			f21 = raft_matrix_element(source.data, i_maj, j_min);
			f22 = raft_matrix_element(source.data, i_maj, j_maj);
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

raft_image bl_interpolate_mt(raft_image source, int Nx, int Ny, int nthreads)
{
	int old_Nx = source.data.lines;
	int old_Ny = source.data.columns;

	double x0 =std::min(source.tl_x, source.br_x); 
	double y0 =std::min(source.tl_y, source.br_y); 
	double x1 =std::max(source.tl_x, source.br_x); 
	double y1 =std::max(source.tl_y, source.br_y); 
	double old_dx = (x1-x0)/old_Nx;
	double old_dy = (y1-y0)/old_Ny;
	double dx = (x1-x0)/Nx;
	double dy = (y1-y0)/Ny;
	raft_image res;
	res.data = raft_matrix_create(Nx, Ny);

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
		std::thread thread =  std::thread( bl_interpolate_worker,
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


	res.tl_x = source.tl_x;
	res.tl_y = source.tl_y;
	res.br_x = source.br_x;
	res.br_y = source.br_y;
// 	std::cout<<"________bl_interpolation succeed___________"<<std::endl;
	return res;
}

// 1D vector linear interpolation to new size N
raft_vector interpolate_1d_mt(raft_vector f, int N, int nthreads)
{
	int old_N = f.size;
	raft_vector res = raft_vector_create(N);


	double old_dx = 1.0/old_N;
	double dx = 1.0/N;


	int i_min, i_maj;
	double x, x_min, x_maj, f_min, f_maj;
	for(int i=0; i<N; i++) {
		// getting old elements numbers
		x = i*dx;
		i_min = floor(x/old_dx);
		i_maj = i_min + 1;
		x_min = i_min*old_dx;
		x_maj = x_min + old_dx;
		f_min = raft_vector_element(f, i_min);
		f_maj = raft_vector_element(f, i_maj);
		
		raft_vector_element(res, i) = f_min + (f_maj - f_min)*
			(x - x_min)/(x_maj - x_min);
	}
	return res;
}

