#include "raft_image_functions.h"
#include <fftw3.h>
#include <algorithm>
#include <iostream>
#include <thread>

extern "C" {

// BLAS dcopy:
void dcopy_( int const *, double const *, int const *, double *, int const * );

}

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

double mean(raft_matrix data)
{
	int Ns = data.lines;
	int Nt = data.columns;

	double res;
	for(int j = 0; j<Nt; j++) {
		for(int i=0; i<Ns; i++) {
			res += raft_matrix_element(data, i, j);
		}
	}
	return res/(Nt*Ns);
}

double mean_dc(raft_matrix data)
{
	int Ns = data.lines;
	int Nt = data.columns;

	double res;
	for(int j = 0; j<Nt; j++) {
		res += raft_matrix_element(data, 0, j);
	}
	return res/(Nt);
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

double desc_l2(raft_image img1, raft_image img2) 
{
	int Nx = img1.data.lines;
	int Ny = img1.data.columns;
	if(Nx != img2.data.lines || Ny != img2.data.columns) {
		//std::cout<<"shopa with sizes"<<std::endl;
	}
	
	double res = 0;	
	for(int i = 0; i<Nx; i++) {
		for(int j=0; j<Ny; j++) {
			double a = raft_matrix_element(img1.data, i, j) - raft_matrix_element(img2.data, i, j);
			res += a*a/(Nx*Ny);
		}
	}
	return res;


}

void raft_image_normalize(raft_image &source)
{
	int Ns = source.data.lines;
	int Nt = source.data.columns;

	int im, jm;
	double max = fabs(raft_matrix_maximum_value(source.data, im, jm));
	//std::cout<<"maximum in the matrix: "<<max<<"with indexes "<<im<<", "<<jm<<std::endl;
	double el;
	for(int j=0; j<Nt; j++) {
		for(int i=0; i<Ns; i++) {
			raft_matrix_element(source.data, i, j) = 
				fabs(raft_matrix_element(source.data, i, j))/max;
		}
	}

// 	source.tl_x = -1.0;
// 	source.br_x = 1.0;
// 	source.tl_y = 1.0;
// 	source.br_x = -1.0;
}

raft_image sino2sp(raft_image source)
{
	int Ns = source.data.lines;
	int Nt = source.data.columns;

	raft_image res;
	res.data = raft_matrix_create(Ns/2, 2*Nt);
	res.tl_x = -pi;
	res.br_x = pi;
	res.tl_y = 1;
	res.br_y = 0;

	for(int t=0;t<Nt;++t) {
		for(int s=0;s<Ns/2;++s) {
			raft_matrix_element(res.data, s, t) =    raft_matrix_element(source.data, Ns/2 + s, t);
			raft_matrix_element(res.data, s, t+Nt) = raft_matrix_element(source.data, Ns/2 - s, t);
		}
	}
// 	std::cout<<"sino2sp succeed"<<std::endl;
	return res;
}

raft_image get_sector(raft_image source, double t0, double t1)
{
	int Ns = source.data.lines;
	int Nt = source.data.columns;
	
	double old_t0 = std::min(source.br_x, source.tl_x);
	double old_t1 = std::max(source.br_x, source.tl_x);
	double dt = (old_t1 - old_t0)/Nt;

	//std::cout<<"_____________ GET SECTOR _____________"<<std::endl;
	//std::cout<<"Old angles: t0="<<old_t0<<"; t1="<<source.tl_x<<std::endl;
	//std::cout<<"New angles: t0="<<t0<<"; t1="<<t1<<std::endl;
	//std::cout<<"dt="<<dt<<std::endl;


	raft_image res;
	if(t0 > -pi && t1 < pi) {
		int idx_start = floor((t0 - old_t0)/dt);
		int idx_end = floor((t1 - old_t0)/dt);
		int N = idx_end - idx_start;
		//std::cout<<"!!!!!!!!!!!!!!!!   "<<idx_start<<", "<<idx_end<<" !!!!!!"<<//std::cout;
		res.data = raft_matrix_create(Ns, N);
		for(int j=0;j<N;j++) {
			for(int i=0;i<Ns; ++i) {
				if(j+idx_start < Nt) {
					raft_matrix_element(res.data, i,  j) = raft_matrix_element(source.data, i, j+idx_start);
				}
			}
		}
	} else if(t1 > pi && t0 < pi && t0 > -pi) {
		t1 -= 2*pi;
		int idx1 = floor((t0 - old_t0)/dt);
		int idx2 = floor((t1 - old_t0)/dt);
		//std::cout<<"2 case: idx1="<<idx1<<"; idx2="<<idx2<<std::endl;
		int N1 = Nt - idx1 - 1;
		int N2 = idx2;
		int N = N1 + N2; 
		res.data = raft_matrix_create(Ns, N);
		for(int i=0; i<Ns; i++) {
			for(int j=idx1; j<Nt-1; j++) {
				raft_matrix_element(res.data, i, j-idx1) = raft_matrix_element(source.data, i, j);	
			}
			for(int j=0; j<idx2-1; j++) {
				raft_matrix_element(res.data, i, j+N1) = raft_matrix_element(source.data, i, j);
			}
		}
		t1 += 2*pi;	
	} else if(t0 < -pi && t1> -pi && t1 < pi) {
		t0 += 2*pi;
		int idx1 = floor((t0 - old_t0)/dt);
		int idx2 = floor((t1 - old_t0)/dt);
		//std::cout<<"3 case: idx1="<<idx1<<"; idx2="<<idx2<<std::endl;
		int N1 = Nt - idx1 - 1;
		int N2 = idx2;
		int N = N1 + N2; 
		res.data = raft_matrix_create(Ns, N);
		for(int i=0; i<Ns; i++) {
			for(int j=idx1; j<Nt-1; j++) {
				raft_matrix_element(res.data, i, j-idx1) = raft_matrix_element(source.data, i, j);	
			}
			for(int j=0; j<N2; j++) {
				raft_matrix_element(res.data, i, j+N1) = raft_matrix_element(source.data, i, j);
			}
		}
		t0 -= 2*pi;
	} else if(t0 < -pi && t1 > pi) {
		t0 += 2*pi;
		t1 -= 2*pi;
		int idx1 = floor((t0 - old_t0)/dt);
		int idx2 = floor((t1 - old_t0)/dt);
		//std::cout<<"4 case: idx1="<<idx1<<"; idx2="<<idx2<<std::endl;
		int N1 = Nt - 1 - idx1;
		int N2 = idx2;
		int N = Nt + N1 + N2;
		res.data = raft_matrix_create(Ns, N);
		for(int i=0; i<Ns; i++) {
			for(int j=idx1; j<Nt-1; j++) {
				raft_matrix_element(res.data, i, j-idx1) = raft_matrix_element(source.data, i, j);
			}
			for(int j=0; j<Nt-1; j++) {
				raft_matrix_element(res.data, i, j+N1) = raft_matrix_element(source.data, i, j);
			}
			for(int j=0; j<N2; j++) {
				raft_matrix_element(res.data, i, N - N2 - 1 + j) = raft_matrix_element(source.data, i, j);
			}
		}
		t0 -= 2*pi;
		t1 += 2*pi;
	}
	res.tl_x = t1;
	res.br_x = t0;
	res.tl_y = source.tl_y;
	res.br_y = source.br_y;
	//std::cout<<"get_sector suceed"<<std::endl;
	return res;
}

raft_image filter_sector(raft_image source, double t0, double t1)
{
	int Ns = source.data.lines;
	int Nt = source.data.columns;
	
	double old_t0 = std::min(source.br_x, source.tl_x);
	double old_t1 = std::max(source.br_x, source.tl_x);
	double dt = (old_t1 - old_t0)/Nt;

	//std::cout<<"_____________ GET SECTOR _____________"<<std::endl;
	//std::cout<<"Old angles: t0="<<old_t0<<"; t1="<<source.tl_x<<std::endl;
	//std::cout<<"New angles: t0="<<t0<<"; t1="<<t1<<std::endl;
	//std::cout<<"dt="<<dt<<std::endl;

	raft_image res;
	res.data = raft_matrix_create(Ns, Nt);

	int idx_start = floor((t0 - old_t0)/dt);
	int idx_end = floor((t1 - old_t0)/dt);
	int N = idx_end - idx_start;
	for(int j=idx_start;j<idx_end;j++) {
		for(int i=0;i<Ns; ++i) {
			raft_matrix_element(res.data, i,  j) = raft_matrix_element(source.data, i, j);
		}
	}

	res.tl_x = source.tl_x;
	res.br_x = source.br_x;
	res.tl_y = source.tl_y;
	res.br_y = source.br_y;
	//std::cout<<"get_sector suceed"<<std::endl;
	return res;
}

raft_image rotate(raft_image source, double t)
{
	int Ns = source.data.lines;
	int Nt = source.data.columns;
	
	raft_image res;
	
	double t0 = source.tl_x;
	double t1 = source.br_x;


	double dt = (t1 - t0)/Nt;
	int idxt =floor(t/dt);

	raft_matrix tmp = raft_matrix_create(Ns, Nt);
	for(int j=0;j<Nt;j++) {
		for(int i=0;i<Ns; ++i) {
			if(idxt>0) {
				if(j+idxt<Nt) {
					raft_matrix_element(tmp, i,  j + idxt) = raft_matrix_element(source.data, i, j);
				} else {
					raft_matrix_element(tmp, i,  j + idxt - Nt) = raft_matrix_element(source.data, i, j);
				}
			} else {
				if(j+idxt>0) {
					raft_matrix_element(tmp, i,  j + idxt) = raft_matrix_element(source.data, i, j);
				} else {
					raft_matrix_element(tmp, i,  j + idxt + Nt -1) = raft_matrix_element(source.data, i, j);
				}
			}
		}
	}
	res.data = tmp;
	res.tl_x = source.tl_x;
	res.br_x = source.br_x;
	res.tl_y = source.tl_y;
	res.br_y = source.br_y;
	return res;
}

raft_image cut(raft_image source, double tl_x, double br_x, 
		double br_y, double tl_y)
{

	// x0 -> br_x
	// y0 -> tl_y
	// x1 -> tl_x
	// y1 -> br_y
	int Nx = source.data.columns;
	int Ny = source.data.lines;

	double x0 = source.tl_x;
	double x1 = source.br_x;

	double y0 = source.br_y;
	double y1 = source.tl_y;

	double dx = (x1-x0)/Nx;	
	double dy = (y1-y0)/Ny;	

	//std::cout<<"________cutting______________"<<std::endl;
	//std::cout<<"old tl_x="<<x0<<"; old br_x="<<x1<<"; new tl_x="<<tl_x<<"; new br_x="<<br_x<<std::endl;
	//std::cout<<"old br_y="<<y0<<"; old tl_y="<<y1<<"; new br_y="<<br_y<<"; new tl_y="<<tl_y<<std::endl;
	if(tl_x < x0) tl_x = x0;
	if(br_x > x1) br_x = x1;
	
	if(br_y < y0) br_y = y0;
	if(tl_y > y1) br_y = y1;

	int idx_x0 = floor((tl_x - x0)/dx);
	int idx_x1 = floor((br_x - x0)/dx);

	int idx_y0 = floor((br_y - y0)/dy);
	int idx_y1 = floor((tl_y - y0)/dy);	

	int nx = idx_x1 - idx_x0;
	int ny = idx_y1 - idx_y0;

	//std::cout<<"idx_x0="<<idx_x0<<"; idx_x1="<<idx_x1<<" idx_y0="<<idx_y0<<"; idx_y1="<<idx_y1<<std::endl;
	//std::cout<<"old size: <"<<Nx<<" x "<<Nx<<">; dx="<<dx<<"; dy="<<dy<<std::endl;
	//std::cout<<"new size: <"<<nx<<" x "<<ny<<">"<<std::endl; 

	raft_image res;
	res.data = raft_matrix_create(ny, nx);
	for(int j=0; j<nx; ++j) {
		for(int i=0; i<ny; ++i) {
			raft_matrix_element(res.data, i, j) = raft_matrix_element(source.data, i+idx_y0, j+idx_x0);	
		}
	}
	res.tl_x = tl_x; 
	res.tl_y = tl_y;
	res.br_x = br_x;
	res.br_y = br_y;
	return res;
}

raft_image raft_image_flip_x(raft_image source)
{
	raft_image res;
	int Nt = source.data.columns;
	int Ns = source.data.lines;
	res.data = raft_matrix_create(Nt, Ns);
	for(int j=0; j<Nt; j++) {
		for(int i=0; i<Ns; i++) {
			double tmp = raft_matrix_element(source.data, i, j);
			raft_matrix_element(res.data, i, j) = 
				raft_matrix_element(source.data, i, Nt - j - 1);
			raft_matrix_element(res.data, i, Nt -j - 1) = tmp;
		}
	}
	res.tl_x = source.tl_x;
	res.tl_y = source.tl_y;
	res.br_x = source.br_x;
	res.br_y = source.br_y;
	return res;
}

raft_image raft_image_flip_y(raft_image &source)
{
	raft_image res;
	int Nt = source.data.columns;
	int Ns = source.data.lines;
	res.data = raft_matrix_create(Nt, Ns);
	for(int j=0; j<Nt; j++) {
		for(int i=0; i<Ns; i++) {
			double tmp = raft_matrix_element(source.data, i, j);
			raft_matrix_element(res.data, i, j) = 
				raft_matrix_element(source.data, Ns - i - 1, j);
			raft_matrix_element(res.data, Ns - i - 1, j) = tmp;
		}
	}
	res.tl_x = source.tl_x;
	res.tl_y = source.tl_y;
	res.br_x = source.br_x;
	res.br_y = source.br_y;
	return res;
}

raft_image zero_padding(raft_image source, int Nl, int Nc)
{
	int min_l = (source.data.lines < Nl) ? source.data.lines : Nl;
	int min_c = (source.data.columns < Nl) ? source.data.columns : Nl;

// 	double dl = source.data.line_stride;
// 	double dc = source.data.column_stride;

	double dc = (source.tl_x - source.br_x)/source.data.columns;
	double dl = (source.br_y - source.tl_y)/source.data.lines;


// 	std::cout<<"Zerro padding started:"<<std::endl;
// 	std::cout<<"dc="<<dc<<"; dl="<<dl<<std::endl;;
// 	std::cout<<"OLD size: lines="<<source.data.lines<<"; columns="<<source.data.columns<<std::endl;
// 	std::cout<<"NEW size: lines="<<Nl<<"; columns="<<Nc<<std::endl;
	raft_image res = source;
	res.data = raft_matrix_create(Nl, Nc);
	int start_j = (Nc - source.data.columns)/2;
	int start_i = (Nl - source.data.lines)/2;
	for(int j=0; j<min_c; j++) {
		for(int i=0; i<min_l; i++) {
			raft_matrix_element(res.data, i + start_i, j + start_j) = raft_matrix_element(source.data, i, j);	
		}
	}	

	double new_max_x = source.br_x + Nc*dc;
	double new_max_y = source.tl_y + Nl*dl;

// 	std::cout<<"OLD range ("<<source.br_x<<", "<<source.tl_x<<") X ("<<source.tl_y<<", "<<source.br_y<<")."<<std::endl;
// 	std::cout<<"NEW range ("<<source.br_x<<", "<<new_max_x<<") X ("<<source.tl_y<<", "<<new_max_y<<")."<<std::endl;
	res.tl_x = new_max_x;
	res.br_y = new_max_y;
// 	std::cout<<"zerro padding succeed"<<std::endl;

	return res;
}

raft_image zero_padding_on_s(raft_image source, int Ns)
{
	int nt = source.data.columns;
	int ns = source.data.lines;
	if(Ns <= ns) return source;
	
	double ds = fabs(source.br_y - source.tl_y)/ns;

	double s_max = source.br_y + Ns*ds;

	raft_image res;
	res.br_x = source.br_x;
	res.tl_x = source.tl_x;
	res.br_y = source.br_y;
	res.tl_y = s_max;

	//std::cout<<"Zerro padding started:"<<std::endl;
	//std::cout<<"ds="<<ds<<std::endl;;
	//std::cout<<"OLD size: lines="<<source.data.lines<<"; columns="<<source.data.columns<<std::endl;
	//std::cout<<"NEW size: lines="<<Ns<<"; columns="<<nt<<std::endl;

	res.data = raft_matrix_create(Ns, nt);
	for(int j=0; j<ns; j++) {
		for(int i=0; i<nt; i++) {
			raft_matrix_element(res.data, j, i) = raft_matrix_element(source.data, j, i);
		}	
	}
	return res;

}

raft_image zero_padding_sino(raft_image source, int Ns)
{
	int nt = source.data.columns;
	int ns = source.data.lines;
	if(Ns <= ns) return source;
	
	double ds = fabs(source.br_y - source.tl_y)/ns;

	double s_max = Ns*ds/2;

	raft_image res;
	res.br_x = source.br_x;
	res.tl_x = source.tl_x;
	res.tl_y = s_max;
	res.br_y = - s_max;

	//std::cout<<"Zerro padding started:"<<std::endl;
	//std::cout<<"; ds="<<ds<<std::endl;;
	//std::cout<<"OLD size: lines="<<source.data.lines<<"; columns="<<source.data.columns<<std::endl;
	//std::cout<<"NEW size: lines="<<Ns<<"; columns="<<nt<<std::endl;

	res.data = raft_matrix_create(Ns, nt);
	int idx_start = (Ns - ns)/2;
	for(int j=0; j<ns; j++) {
		for(int i=0; i<nt; i++) {
			raft_matrix_element(res.data, j + idx_start, i) = 
				raft_matrix_element(source.data, j, i);
		}	
	}
	return res;

}

raft_image zero_padding_on_s2(raft_image source, int Ns)
{
	int nt = source.data.columns;
	int ns = source.data.lines;
	if(Ns <= ns) return source;
	
	double ds = fabs(source.br_y - source.tl_y)/ns;
	double s_max = source.tl_y + Ns*ds;

	raft_image res;
	res.br_x = source.br_x;
	res.tl_x = source.tl_x;
	res.tl_y = source.tl_y;
	res.br_y = s_max;

// 	std::cout<<"Zerro padding started:"<<std::endl;
// 	std::cout<<"; ds="<<ds<<std::endl;;
// 	std::cout<<"OLD size: lines="<<source.data.lines<<"; columns="<<source.data.columns<<std::endl;
// 	std::cout<<"NEW size: lines="<<Ns<<"; columns="<<nt<<std::endl;

	res.data = raft_matrix_create(Ns, nt);
	for(int j=0; j<Ns; j++) {
		for(int i=0; i<nt; i++) {
			if(j<ns)
				raft_matrix_element(res.data, j, i) = raft_matrix_element(source.data, j, i);
			if(j+ns < Ns) {
				raft_matrix_element(res.data, j+ns, i) = 
					raft_matrix_element(source.data, ns-j-1, i);
			}
		}	
	}
	return res;

}

raft_image bl_interpolate(raft_image source, int Nx, int Ny)
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

	double x = x0;
	double y = y0; 

	int i_min, i_maj, j_min, j_maj;
	double x_min, x_maj, y_min, y_maj, f_res, f11, f12, f21, f22, tmp;
	tmp = old_dx*old_dy;
// 	std::cout<<"_________bl_interpolate starts____________"<<std::endl;
// 	std::cout<<"old_size: "<<old_Nx<<" X "<<old_Ny;
// 	std::cout<<"; new size: "<<Nx<<" X "<<Ny<<std::endl;
	raft_image res;
	res.data = raft_matrix_create(Nx, Ny);
	for(int j=0; j<Ny; ++j) {
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
	res.tl_x = source.tl_x;
	res.tl_y = source.tl_y;
	res.br_x = source.br_x;
	res.br_y = source.br_y;
// 	std::cout<<"________bl_interpolation succeed___________"<<std::endl;
	return res;
}

raft_image sp2c(raft_image source, int Nx, int Ny)
{
	int Nt = source.data.columns;
	double t0 = std::min(source.br_x, source.tl_x);
	double t1 = std::max(source.br_x, source.tl_x);
	double dt = (t1-t0)/Nt;
	
	int Ns = source.data.lines;
	double s0 = std::min(source.br_y, source.tl_y);
	double s1 = std::max(source.br_y, source.tl_y);
	double ds   = (s1-s0)/Ns;
	
	double x0 =  -s1;
	double x1 =   s1;
	double y0 =  -s1;
	double y1 =   s1;

	double dx = s1*2/Nx;
	double dy = s1*2/Ny;

// 	std::cout<<"___________ semi-polar 2 cartesian interpolation ________"<<std::endl;
// 	std::cout<<"old size: Ns="<<Ns<<"; Nt="<<Nt<<std::endl;
// 	std::cout<<"new size: Nx="<<Nx<<"; Ny="<<Ny<<std::endl;
// 	std::cout<<"s0="<<s0<<"; s1="<<s1<<"; ds="<<ds<<std::endl;
// 	std::cout<<"t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
// 	std::cout<<"x0="<<x0<<"; x1="<<x1<<"; dx="<<dx<<std::endl;
// 	std::cout<<"y0="<<y0<<"; y1="<<y1<<"; dy="<<dy<<std::endl;

	double x = x0; 
	double y = y0;
	raft_image res;
	res.data = raft_matrix_create(Ny, Nx);


	int idxt_min, idxt_maj, idxs_min, idxs_maj;
	double s, t, t_min, t_maj, s_min, s_maj, 
	       f11, f12, f21, f22, fr;
	double factor = 1.0/dt*ds;
	for(int j=0; j<Ny; j++) {
		for(int i=0; i<Nx; i++) {
			s = sqrt(x*x + y*y);
			t = 180*acos(x/s)/pi;
			t = (y>0) ? t :  - t;

			idxt_min = floor((t-t0)/dt);
			idxt_maj = idxt_min + 1;
			t_min = t0 + dt*idxt_min;
			t_maj = t_min + dt;

			idxs_min = floor((s - s0)/ds);
			idxs_maj = idxs_min + 1;
			s_min = s0 + ds*idxs_min;
			s_maj = s_min + ds;

			if(idxs_maj > Ns-1 || idxt_maj > Nt-1 ||
					idxs_min < 0 || idxt_min < 0) {
				x += dx;
				continue;
			}

			f11 = raft_matrix_element(source.data, idxs_min, idxt_min);
			f12 = raft_matrix_element(source.data, idxs_maj, idxt_min);
			f21 = raft_matrix_element(source.data, idxs_min, idxt_maj);
			f22 = raft_matrix_element(source.data, idxs_maj, idxt_maj);
		
			fr =  f11*(t_maj - t)*(s_maj - s) + 
					f21*(t - t_min)*(s_maj - s) + 
					f12*(t_maj - t)*(s - s_min) + 
					f22*(t - t_min)*(s - s_min) ;
			fr *= factor;
			raft_matrix_element(res.data, i, j) = fr;
			x += dx;
		}
		x = x0;
		y += dy;
	}
	res.br_x = x0;
	res.tl_y = y0;
	res.tl_x = x1;
	res.br_y = y1;
	return res;
}

raft_image c2sp(raft_image source, int Nt, int Ns)
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
	double ds = (s1 - s0)/Ns;

	double t0 = -180;
	double t1 = 180;
	double dt = 360.0/Nt;

// 	//std::cout<<"___________ cartesian 2 semi-polar interpolation ________"<<std::endl;
// 	//std::cout<<"old size: Nx="<<Nx<<"; Ny="<<Ny<<std::endl;
// 	//std::cout<<"new size: Nt="<<Nt<<"; Ns="<<Ns<<std::endl;
// 	std::cout<<"s0="<<s0<<"; s1="<<s1<<"; ds="<<ds<<std::endl;
// 	std::cout<<"t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
// 	std::cout<<"x0="<<x0<<"; x1="<<x1<<"; dx="<<dx<<std::endl;
// 	std::cout<<"y0="<<y0<<"; y1="<<y1<<"; dy="<<dy<<std::endl;

	double x = x0; 
	double y = y0;

	double t = t0;
	double s = s0;

	raft_image res;
	res.data = raft_matrix_create(Ns, Nt);
	res.br_x = t0;
	res.tl_x = t1; 
	res.tl_y = s0;
	res.br_y = s1;

	int i_min, i_maj, j_min, j_maj;
	double x_min, x_maj, y_min, y_maj, f11, f12, f21, f22, f, dv;

	for(int j=0; j<Nt; ++j) {
		double _cos = cos(pi*t/180);
		double _sin = sin(pi*t/180);
		for(int i=0; i<Ns; ++i) {
			x = s*_cos; 
			y = s*_sin;

			j_min = floor((x-x0)/dx);
			j_maj = j_min + 1;
			i_min = floor((y-y0)/dy);
			i_maj = i_min + 1;
			
			x_min = x0 + j_min*dx;
			x_maj = x_min + dx;
			y_min = y0 + i_maj*dy;
			y_maj = y_min + dy;
			
			if((j_maj > Nx-1) || (j_min < 0) || (i_maj > Ny-1) || (i_min < 0)) {
				s+=ds;
				continue;
			}
		
			f11 = raft_matrix_element(source.data, i_min, j_min);
			f12 = raft_matrix_element(source.data, i_min, j_maj);
			f21 = raft_matrix_element(source.data, i_maj, j_min);
			f22 = raft_matrix_element(source.data, i_maj, j_maj);
			
			dv = dx*dy;
			f =  f11*(x_maj - x)*(y_maj - y) + 
					f21*(x - x_min)*(y_maj - y) + 
					f12*(x_maj - x)*(y - y_min) + 
					f22*(x - x_min)*(y - y_min) ;
			f = f/dv;
			raft_matrix_element(res.data, i, j) = f;
			s+=ds;
		}
		s = s0;
		t+=dt;
	}	

	return res;
}

raft_image sp2lp(raft_image source, double Nr, double r0)
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

	double r = r0;
	double s, s_min, s_maj, f1, f2, f_res;
	int i_min, i_maj;
	for(int j = 0; j<Nt; j++) {
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
	res.tl_y = r0;
	res.br_y = r1;
// 	std::cout<<"interpolation p2lp succeed"<<std::endl;
	return res; 
}

void lp2sp(raft_image source, raft_image &res, int Ns)
{
	int Nt = source.data.columns;
	int Nr = source.data.lines;
	double r0 = std::min(source.tl_y, source.br_y);
	double r1 = std::max(source.tl_y, source.br_y);
	double dr = (r1 - r0)/Nr;

	double s0 = exp(r0);
	double s1 = exp(r1); 
	double ds = (s1 - s0)/Nr;

// 	std::cout<<"__________Interpolating from semi-polar to log-polar system________"<<std::endl;
// 	std::cout<<"old sizes: Nt="<<Nt<<"; Ns="<<Ns<<std::endl;
// 	std::cout<<"new sizes: Nt="<<Nt<<"; Nr="<<Nr<<std::endl; 
// 	std::cout<<"s1="<<s1<<"; ds="<<ds<<std::endl;
// 	std::cout<<"r0="<<r0<<"; r1="<<r1;
// 	std::cout<<"; dr="<<dr<<std::endl;

	double s = s0;
	double r, r_min, r_maj, f1, f2, f_res;
	int i_min, i_maj;
	res = source;
	res.data = raft_matrix_create(Nr, Nt);
	for(int j = 0; j<Nt; j++) {
		for(int i=0; i<Ns; ++i) {
			r = log(s);
			if(r>r1) {
// 				std::cout<<"@!!!!!"<<std::endl;
				s+=ds;
				continue;
			}
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
			f_res = f_res/dr;
			
			raft_matrix_element(res.data, i, j) = f_res;
			s+=ds;
		}
		s = s0;
	}
	res.tl_y = s0;
	res.br_y = s1;
// 	std::cout<<"interpolation p2lp succeed"<<std::endl;
}

raft_image lp2c(raft_image source, int Nx, int Ny)
{
	int Nt = source.data.columns;
	double t0 = std::min(source.br_x, source.tl_x);
	double t1 = std::max(source.br_x, source.tl_x);
	double dt = (t1-t0)/Nt;
	
	int Nr = source.data.lines;
	double r0 = std::min(source.br_y, source.tl_y);
	double r1 = std::max(source.br_y, source.tl_y);
	double dr   = (r1-r0)/Nr;

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

	double x = x0; 
	double y = y0;

	raft_image res;
	res.data = raft_matrix_create(Ny, Nx);
	int idxt_min, idxt_maj, idxr_min, idxr_maj;
	double s, r, t, t_min, t_maj, r_min, r_maj, 
	       f11, f12, f21, f22, factor, fr;
	for(int j=0; j<Ny; j++) {
		for(int i=0; i<Nx; i++) {
			s = sqrt(x*x + y*y);
			r = log(s);
			t = 180*acos(x/s)/pi;
			t = (y>0) ? t :  - t;
			
			idxt_min = floor((t-t0)/dt);
			idxt_maj = idxt_min + 1;
			t_min = t0 + dt*idxt_min;
			t_maj = t_min + dt;

			idxr_min = floor((r - r0)/dr);
			idxr_maj = idxr_min + 1;
			r_min = r0 + dr*idxr_min;
			r_maj = r_min + dr;

// 			if(idxt_maj = Nt && fabs(t) > 179.9) idxt_maj=0;
			if(idxr_maj > Nr-1 || idxt_maj > Nt ||
					idxr_min < 0 || idxt_min < 0) {
				x += dx;
				continue;
			}

			f11 = raft_matrix_element(source.data, idxr_min, idxt_min);
			f12 = raft_matrix_element(source.data, idxr_maj, idxt_min);
			f21 = raft_matrix_element(source.data, idxr_min, idxt_maj);
			f22 = raft_matrix_element(source.data, idxr_maj, idxt_maj);
		
			factor = 1.0/dt*dr/2;
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
	res.tl_x = x0;
	res.tl_y = y0;
	res.br_x = x1;
	res.br_y = y1;

// 	std::cout<<"lp2c suceed"<<std::endl;
	return res;
}

raft_image c2lp(raft_image source, int Nt, int Nr, double r0)
{
	int Nx = source.data.columns;
	double x0 = std::min(source.tl_x, source.br_x);
	double x1 = std::max(source.tl_x, source.br_x);
	double dx = (x1-x0)/Nx;
	
	int Ny = source.data.lines;
	double y0 = std::min(source.br_y, source.tl_y);
	double y1 = std::max(source.br_y, source.tl_y);
	double dy   = (y1-y0)/Ny;

// 	double r1 = log(sqrt(x1*x1 + y1*y1)); 
	double r1 = 0; 
	double dr = (r1 - r0)/Nr;

	double t0 = -180;
	double t1 = 180;
	double dt = 360.0/Nt;

	//std::cout<<"___________ cartesian 2 semi-polar interpolation ________"<<std::endl;
	//std::cout<<"old size: Nx="<<Nx<<"; Ny="<<Ny<<std::endl;
	//std::cout<<"new size: Nt="<<Nt<<"; Nr="<<Nr<<std::endl;
	//std::cout<<"r0="<<r0<<"; r1="<<r1<<"; dr="<<dr<<std::endl;
	//std::cout<<"t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
	//std::cout<<"x0="<<x0<<"; x1="<<x1<<"; dx="<<dx<<std::endl;
	//std::cout<<"y0="<<y0<<"; y1="<<y1<<"; dy="<<dy<<std::endl;

	double x = x0; 
	double y = y0;

	double t = t0;
	double r = r0;

	raft_image res;
	res.data = raft_matrix_create(Nr, Nt);
	res.br_x = t0;
	res.tl_x = t1; 
	res.tl_y = r0;
	res.br_y = r1;

	int i_min, i_maj, j_min, j_maj;
	double x_min, x_maj, y_min, y_maj, f11, f12, f21, f22, f, dv;

	for(int j=0; j<Nt; ++j) {
		double _cos = cos(pi*(t-90)/180);
		double _sin = sin(pi*(t-90)/180);
		for(int i=0; i<Nr; ++i) {
			x = exp(r)*_cos; 
			y = exp(r)*_sin;

			j_min = floor((x-x0)/dx);
			j_maj = j_min + 1;
			i_min = floor((y-y0)/dy);
			i_maj = i_min + 1;
			
			x_min = x0 + j_min*dx;
			x_maj = x_min + dx;
			y_min = y0 + i_maj*dy;
			y_maj = y_min + dy;
			
			if((j_maj > Nx-1) || (j_min < 0) || (i_maj > Ny-1) || (i_min < 0)) {
				r+=dr;
				continue;
			}
		
			f11 = raft_matrix_element(source.data, i_min, j_min);
			f12 = raft_matrix_element(source.data, i_min, j_maj);
			f21 = raft_matrix_element(source.data, i_maj, j_min);
			f22 = raft_matrix_element(source.data, i_maj, j_maj);
			
			dv = dx*dy;
			f =  f11*(x_maj - x)*(y_maj - y) + 
					f21*(x - x_min)*(y_maj - y) + 
					f12*(x_maj - x)*(y - y_min) + 
					f22*(x - x_min)*(y - y_min) ;
			f = f/dv;
			raft_matrix_element(res.data, i, j) = f;
			r+=dr;
		}
		r = r0;
		t+=dt;
	}
	return res;
}

////////////////////// FFT's /////////////////////////// 
void fft_shift_2d(raft_matrix &x)
{
	int Nx = x.lines;
	int Ny = x.columns;
	double tmp;
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny/2;j++){
			tmp = raft_matrix_element(x, i,j);
			raft_matrix_element(x, i,j) = raft_matrix_element(x, i,Ny/2+j);
			raft_matrix_element(x, i,Ny/2+j) = tmp;
		}
	}

	for(int i=0;i<Nx/2;i++){
		for(int j=0;j<Ny;j++){
			tmp = raft_matrix_element(x, i,j);
			raft_matrix_element(x, i,j) = raft_matrix_element(x, i+Nx/2,j);
			raft_matrix_element(x, i+Nx/2,j) = tmp;
		}
	}
}

void convolution_2d(raft_matrix x, raft_matrix k, raft_matrix &res, int threads)
{
	int Nx = x.lines;
	int Ny = x.columns;

	double scale = 1.0/(Nx*Ny);
// 	std::cout<<"conv_2d starts: Mx="<<Nx<<"; Ny="<<Ny<<std::endl;
	fftw_complex *x_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	fftw_complex *k_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	fftw_complex *x_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	fftw_complex *k_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	for(int i=0; i<Nx; i++) {
		for(int j=0; j<Ny; j++) {
			int idx = j*Nx+i;
			x_in[idx][0] = 
				raft_matrix_element(x, i, j);
			x_in[idx][1] = 0.0; 
			k_in[idx][0] = 
				raft_matrix_element(k, i, j);	
			k_in[idx][1] = 0.0; 
		}
	}
	fftw_init_threads();
	fftw_plan_with_nthreads(threads);

	fftw_plan p = fftw_plan_dft_2d(Nx, Ny, x_in, x_out, FFTW_FORWARD, 
			FFTW_ESTIMATE);
	fftw_execute_dft(p, x_in, x_out); 
// 	std::cout<<"conv2d: straight FFT(x) succeed"<<std::endl;
	fftw_plan pk = fftw_plan_dft_2d(Nx, Ny, k_in, k_out, FFTW_FORWARD, 
			FFTW_ESTIMATE);
	fftw_execute(pk); 
// 	std::cout<<"conv2d: straight FFT(k) succeed"<<std::endl;
	double xr, xi, kr, ki;
	for(int j=0; j<Ny; j++) {
		for(int i=0; i<Nx; i++) { 
			int ij = j*Nx + i;
			xr = x_out[ij][0];
			xi = x_out[ij][1];

			kr = k_out[ij][0];
			ki = k_out[ij][1];
			
			x_in[ij][0] = (xr*kr - xi*ki)*scale;
			x_in[ij][1] = (xr*ki + xi*kr)*scale;
		}
	}

	fftw_free(k_in);
	fftw_free(k_out);

// 	std::cout<<"conv2d: start inverse fourier procedure"<<std::endl;
	fftw_plan p_r = fftw_plan_dft_2d(Nx, Ny, x_in ,x_out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p_r); 

	res = raft_matrix_create(Nx, Ny);
	for(int i=0; i<Nx; i++) {
		for(int j=0; j<Ny; j++) {
			raft_matrix_element(res, i, j) = x_out[j*Nx+i][0];
		}
	}

// 	fft_shift_2d(res);
	fftw_destroy_plan(p);
	fftw_destroy_plan(p_r);
	fftw_free(x_in);
	fftw_free(x_out);
}

void convolution_2d_2(raft_matrix x, raft_matrix k, raft_matrix &res, int threads)
{
	int Nx = x.lines;
	int Ny = x.columns;

	double scale = 1.0/(Nx*Ny);
// 	std::cout<<"conv_2d starts: Mx="<<Nx<<"; Ny="<<Ny<<std::endl;
	fftw_complex *x_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	fftw_complex *k_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	fftw_complex *x_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	fftw_complex *k_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	for(int i=0; i<Nx; i++) {
		for(int j=0; j<Ny; j++) {
			int idx = j*Nx+i;
			x_in[idx][0] = 
				raft_matrix_element(x, i, j);
			x_in[idx][1] = 0.0; 
			k_in[idx][0] = 
				raft_matrix_element(k, i, j);	
			k_in[idx][1] = 0.0; 
		}
	}
	fftw_init_threads();
	fftw_plan_with_nthreads(threads);

	fftw_plan p = fftw_plan_dft_2d(Nx, Ny, x_in, x_out, FFTW_FORWARD, 
			FFTW_ESTIMATE);
	fftw_execute_dft(p, x_in, x_out); 
	fftw_execute_dft(p, k_in, k_out); 
	double xr, xi, kr, ki;
	// convolution
	for(int j=0; j<Ny; j++) {
		for(int i=0; i<Nx; i++) { 
			int ij = j*Nx + i;
			xr = x_out[ij][0]*(i*i+j*j)/(Nx*Ny);
			xi = x_out[ij][1]*(i*i+j*j)/(Nx*Ny);

			kr = k_out[ij][0];
			ki = k_out[ij][1];
			
			x_in[ij][0] = (xr*kr - xi*ki)*scale;
			x_in[ij][1] = (xr*ki + xi*kr)*scale;
			x_out[ij][0] = 0;
			x_out[ij][1] = 0;
		}
	}
	fftw_plan p_r = fftw_plan_dft_2d(Nx, Ny, x_in ,x_out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p_r); 
	// it is a result in x_out





// // 	//deconvolution (TEST)
// 
// 
// 	std::cout<<"FFT of the result done"<<std::endl;
// 	for(int i=0; i<Nx; i++) {
// 		for(int j=0; j<Ny; j++) {
// 			int idx = j*Nx+i;
// 			x_in[idx][0] = x_out[idx][0];
// 				
// 			x_in[idx][1] = 0.0; 
// 			k_in[idx][0] = 
// 				raft_matrix_element(k, i, j);	
// 			k_in[idx][1] = 0.0; 
// 		}
// 	}
// 	fftw_execute_dft(p, x_in, x_out);
// 	fftw_execute_dft(p, k_in, k_out); 
// 	std::cout<<"FFT of the kernel done"<<std::endl;
// 	for(int j=0; j<Ny; j++) {
// 		for(int i=0; i<Nx; i++) { 
// 			int ij = j*Nx + i;
// 			xr = x_out[ij][0];
// 			xi = x_out[ij][1];
// 
// 			kr = k_out[ij][0];
// 			ki = k_out[ij][1];
// 			
// 			scale = 1.0/(Nx*Ny*2*pi*((kr*kr+ki*ki)));
// 			x_in[ij][0] = (xr*kr + xi*ki)*scale;
// 			x_in[ij][1] = (xi*kr - xr*ki)*scale;
// 			x_out[ij][0] = 0;
// 			x_out[ij][1] = 0;
// 		}
// 	}
// 
// 
// 	fftw_execute(p_r); 
// 

	res = raft_matrix_create(Nx, Ny);
	for(int i=0; i<Nx; i++) {
		for(int j=0; j<Ny; j++) {
			raft_matrix_element(res, i, j) = x_out[j*Nx+i][0];
		}
	}

	fft_shift_2d(res);
	fftw_destroy_plan(p);
	fftw_destroy_plan(p_r);
	fftw_free(x_in);
	fftw_free(x_out);
	fftw_free(k_in);
	fftw_free(k_out);
}

void deconvolution_2d(raft_matrix x, raft_matrix k, raft_matrix &res, double a, int threads)
{
	int Nx = x.lines;
	int Ny = x.columns;

	//std::cout<<"deconv_2d starts: Mx="<<Nx<<"; Ny="<<Ny<<std::endl;
	fftw_complex *x_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	fftw_complex *k_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	fftw_complex *x_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	fftw_complex *k_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);


	fftw_init_threads();
	fftw_plan_with_nthreads(threads);
	fftw_plan p = fftw_plan_dft_2d(Nx, Ny, x_in, x_out, FFTW_FORWARD, 
			FFTW_ESTIMATE);


	for(int i=0; i<Nx; i++) {
		for(int j=0; j<Ny; j++) {
			int idx = j*Nx+i;
			x_in[idx][0] = raft_matrix_element(x, i, j);	
			k_in[idx][0] = raft_matrix_element(k, i, j);	
			k_in[idx][1] = 0.0; 
			x_in[idx][1] = 0.0; 
		}
	}
	fftw_execute_dft(p, x_in, x_out);
	fftw_execute_dft(p, k_in, k_out); 
	//std::cout<<"FFT of the kernel done"<<std::endl;
	double xr, xi, kr, ki, scale;
	double cx, cy, d;
	double asize = 1.0; 
	d = asize/Nx;
	for(int j=0; j<Ny; j++) {
		for(int i=0; i<Nx; i++) { 
			int ij = j*Nx + i;
			xr = x_out[ij][0];
			xi = x_out[ij][1];

			kr = k_out[ij][0];
			ki = k_out[ij][1];
		
// 			cx = -Nx/2+i;
// 			cy = -Ny/2+j;
			cx = i;
			cy = j;
// 			scale = 1.0/(Nx*Ny*2*pi*((kr*kr+ki*ki) + a*(cx*cx+cy*cy)*(cx*cx+cy*cy)));
			scale = d*d/((kr*kr+ki*ki) + 
						a*pi*pi*d*d*d*d*(cx*cx+cy*cy)*(cx*cx+cy*cy));
			x_in[ij][0] = (xr*kr + xi*ki)*scale;
			x_in[ij][1] = (xi*kr - xr*ki)*scale;
			x_out[ij][0] = 0;
			x_out[ij][1] = 0;
		}
	}

	fftw_plan p_r = fftw_plan_dft_2d(Nx, Ny, x_in ,x_out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p_r); 

	res = raft_matrix_create(Nx, Ny);
	for(int i=0; i<Nx; i++) {
		for(int j=0; j<Ny; j++) {
			raft_matrix_element(res, i, j) = x_out[j*Nx+i][0];
		}
	}

	fft_shift_2d(res);
	fftw_destroy_plan(p);
// 	fftw_destroy_plan(p_r);
	fftw_free(x_in);
	fftw_free(x_out);
	fftw_free(k_in);
	fftw_free(k_out);
}

void ifftw_2d(raft_matrix &xr, raft_matrix &xi, int threads)
{
// 	std::cout<<"iFFT_2D started"<<std::endl;
	int Nx = xr.columns;
	int Ny = xr.lines;
// 	std::cout<<"Nx="<<Nx<<"; N="<<Ny<<std::endl;

	fftw_complex *in, *out;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	int idx; 
	for(int i=0; i<Nx; i++) {
		for(int j=0; j<Ny; j++) {
			idx = j*Nx + i;
			in[idx][0] = raft_matrix_element(xr, i, j);	
			in[idx][1] = raft_matrix_element(xi, i, j);	

		}
	}
	fftw_init_threads();
	fftw_plan_with_nthreads(threads);
	fftw_plan p = fftw_plan_dft_2d(Nx, Ny, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p); /* repeat as needed */
	fftw_destroy_plan(p);
	fftw_free(in);
	for(int i=0; i<Nx; i++) {
		for(int j=0; j<Ny; j++) {
			idx = j*Nx + i;
			raft_matrix_element(xr, i, j) = out[idx][0];
			raft_matrix_element(xi, i, j) = out[idx][1];
		}
	}
	fftw_free(out);
// 	std::cout<<"iFFT_2D succeed"<<std::endl;
}

void fftw_2d(raft_matrix &xr, raft_matrix &xi, int threads)
{
// 	std::cout<<"iFFT_2D started"<<std::endl;
	int Nx = xr.columns;
	int Ny = xr.lines;
// 	std::cout<<"Nx="<<Nx<<"; N="<<Ny<<std::endl;

	fftw_complex *in, *out;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	int idx; 
	for(int i=0; i<Nx; i++) {
		for(int j=0; j<Ny; j++) {
			idx = j*Nx + i;
			in[idx][0] = raft_matrix_element(xr, i, j);	
			in[idx][1] = raft_matrix_element(xi, i, j);	

		}
	}
	fftw_init_threads();
	fftw_plan_with_nthreads(threads);
	fftw_plan p = fftw_plan_dft_2d(Nx, Ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p); /* repeat as needed */
	fftw_destroy_plan(p);
	fftw_free(in);
	for(int i=0; i<Nx; i++) {
		for(int j=0; j<Ny; j++) {
			idx = j*Nx + i;
			raft_matrix_element(xr, i, j) = out[idx][0];
			raft_matrix_element(xi, i, j) = out[idx][1];
		}
	}
	fftw_free(out);
// 	std::cout<<"iFFT_2D succeed"<<std::endl;
}

raft_image raft_kernel_lp_create(raft_image source, double acc)
{
	int Nt = source.data.columns;
	int Ns = source.data.lines;
	raft_image res = source;
	res.data = raft_matrix_create(Ns, Nt);

	double t0 = std::min(source.tl_x, source.br_x);
	double t1 = std::max(source.tl_x, source.br_x);
	double dt = (t1-t0)/Nt;

	double r0 = std::min(source.tl_y, source.br_y);
	double r1 = std::max(source.tl_y, source.br_y);
	double dr = (r1-r0)/Nt;

	double t=t0;
	r0 = -std::max(source.tl_y, source.br_y);
	double r  = r0;

// 	std::cout<<"Creating kernel: t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
// 	std::cout<<"r0="<<r0<<"; r1="<<r1<<"; dr="<<dr<<std::endl;


	for(int i=0; i<Ns; i++) {
		for(int j=0; j<Nt; j++) {
			double arg = fabs(exp(r)*cos(t) - 1);
			if(arg< acc) raft_matrix_element(res.data, j, i) = 1/acc;
			r+=dr;
		}
		r=r0;
		t+=dt;
	}
	return res;
}

raft_image raft_straight_kernel_lp_create(raft_image source, double acc)
{
	int Nt = source.data.columns;
	int Ns = source.data.lines;
	raft_image res = source;
	res.data = raft_matrix_create(Ns, Nt);

	double t0 = std::min(source.tl_x, source.br_x);
	double t1 = std::max(source.tl_x, source.br_x);
	double dt = (t1-t0)/Nt;

	double r0 = std::min(source.tl_y, source.br_y);
	double r1 = std::max(source.tl_y, source.br_y);
	double dr = (r1-r0)/Nt;

	double t  = t0;
	double r  = r0;

	//std::cout<<"Creating kernel: t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
	//std::cout<<"r0="<<r0<<"; r1="<<r1<<"; dr="<<dr<<std::endl;
	for(int i=0; i<Ns; i++) {
		for(int j=0; j<Nt; j++) {
			double arg = (cos(pi*t/180) - exp(r))*(cos(pi*t/180) - exp(r));
			if(arg< acc*acc) raft_matrix_element(res.data, j, i) = 1/acc;
			r+=dr;
		}
		r=r0;
		t+=dt;
	}
	return res;
}


/////////////// MULTITHREAD IMPLEMNENTATIONS //////////////////////////

// Semi polar -> Cartesian (For miqueles algorithm it takes two arrays at the time 
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
// 	double factor = 1.0/(ds*dt);
	double factor = 1;
	for(int j=col1; j<col2; j++) {
		for(int i=0; i<Nx; i++) {
			s = sqrt(x*x + y*y);
			t = acos(x/s);
			t = (y>0) ? t :  - t;

			idxt_min = floor((t-t0)/dt);
			idxt_maj = idxt_min + 1;
			t_min = t0 + dt*idxt_min;
			t_maj = t_min + dt;

			idxs_min = floor((s - s0)/ds);
			idxs_maj = idxs_min + 1;
			s_min = s0 + ds*idxs_min;
			s_maj = s_min + ds;

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
			fr *= factor;
			fi *= factor;
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
	//std::cout<<"ntreads="<< nthreads<<"; base_ncolumns="<<base_ncolumns<<"; remainder_ncolumns="<<remainder_ncolumns<<std::endl;
	for ( ; cur_thread < nthreads; ++cur_thread )
	{
		double y0_cur = y0 + cur_starting_column*dy; 
		int cur_ncolumns( base_ncolumns + ( cur_thread < remainder_ncolumns ) );
		//std::cout<<"add thr for cur_starting_column="<<cur_starting_column<<"; cur_ncolumns="<<cur_ncolumns<<std::endl;
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

// Log-polar -> Cartesian
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

	//std::cout<<"___________ log-polar 2 cartesian interpolation ________"<<std::endl;
	//std::cout<<"old size: Nr="<<Nr<<"; Nt="<<Nt<<std::endl;
	//std::cout<<"new size: Nx="<<Nx<<"; Ny="<<Ny<<std::endl;
	//std::cout<<"r0="<<r0<<"; r1="<<r1<<"; dr="<<dr<<std::endl;
	//std::cout<<"t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
	//std::cout<<"x0="<<x0<<"; x1="<<x1<<"; dx="<<dx<<std::endl;
	//std::cout<<"y0="<<y0<<"; y1="<<y1<<"; dy="<<dy<<std::endl;


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

	res.tl_x = x0;
	res.tl_y = y0;
	res.br_x = x1;
	res.br_y = y1;

// 	std::cout<<"lp2c suceed"<<std::endl;
	return res;
}

// Semi polar -> log-polar

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
// 	std::cout<<"_________bl_interpolate starts____________"<<std::endl;
// 	std::cout<<"old_size: "<<old_Nx<<" X "<<old_Ny;
// 	std::cout<<"; new size: "<<Nx<<" X "<<Ny<<std::endl;
	for(int j=col1; j<col2; ++j) {
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









raft_image mo_sino(raft_image source, double dx, double dy, int Ns)
{
	double k = sqrt(dx*dx+dy*dy);
	double old_s0 = source.br_y;
	double old_s1 = source.tl_y;
	double s0 = 0; //FIXME
	double s1 = k + old_s1;
	double ds = (s1 - s0)/Ns;

	int Nt = source.data.columns;
	double t0 = source.tl_x;
	double t1 = source.br_x;
	double dt = (t1 - t0)/Nt;

	//std::cout<<"__________ Moving the origin on the sinogram ________"<<std::endl;
	//std::cout<<"sizes: Nt="<<Nt<<"; Ns="<<Ns<<std::endl;
	//std::cout<<"k="<<k<<"; old_s1="<<old_s1<<"; s1="<<s1<<"; ds="<<ds<<std::endl;
	//std::cout<<"t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
	double s = s0;
	double t = t0;
	int i_min, i_maj;
	double s_old, s_min, s_maj, f1, f2, f;

	raft_image res;
	res.data = raft_matrix_create(Ns, Nt);
	for(int j=0;j<Nt;j++) {
		for(int i=0;i<Ns;i++) {
			s_old = s - dx*cos(t) - dy*sin(t);
			if(s_old<old_s0|| s_old>old_s1) {
				s+=ds;
				continue;
			}

			i_min = floor((s_old - s0)/ds);
			i_maj = i_min + 1;

			if(i_min < 0 || i_maj > source.data.lines - 1) {
				s+=ds;
				continue;
			}

			s_min = s0 + i_min*ds;
			s_maj = s_min + ds;

			f1 = raft_matrix_element(source.data, i_min, j);
			f2 = raft_matrix_element(source.data, i_maj, j);

			f = f2*(s_old - s_min) + f1*(s_maj - s_old);
			f = f/ds;
			
			raft_matrix_element(res.data, i, j) = f;
			s+=ds;
		}
		s = s0;
		t+=dt;
	}
	res.tl_x = t0;
	res.br_x = t1;
	res.tl_y = s1;
	res.br_y = s0;
	//std::cout<<"_________mo_sino succeed__________"<<std::endl;
	return res;

}
