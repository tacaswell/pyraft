#include "raft_image_functions.h"
#include <iostream>
#include <thread>
// #include "g_output.h"

double const pi = 3.1415926535897932384626433832795;
extern "C" {

// BLAS dcopy:
void dcopy_( int const *, double const *, int const *, double *, int const * );

}

double norm(double *data, int size)
{
	double res = 0;
	for(int i=0; i<size; i++) {
		res += data[i]/size;
	}
	return res;
}

int iDivUp(int a, int b)
{
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

//Align a to nearest higher multiple of b
int iAlignUp(int a, int b)
{
    return (a % b != 0) ? (a - a % b + b) : a;
}

int snapTransformSize(int dataSize)
{
    int hiBit;
    unsigned int lowPOT, hiPOT;

    dataSize = iAlignUp(dataSize, 16);

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
        return iAlignUp(dataSize, 512);
    }
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

	std::cout<<"Transforming sinogram to polar : ("<<Ns<<","<<Nt<<" ) -> (";
	std::cout<<Ns/2<<","<<2*Nt<<")."<<std::endl;
	for(int t=0;t<Nt;++t) {
		for(int s=0;s<Ns/2;++s) {
			raft_matrix_element(res.data, s, t) =    raft_matrix_element(source.data, Ns/2 + s, t);
			raft_matrix_element(res.data, s, t+Nt) = raft_matrix_element(source.data, Ns/2 - s, t);
		}
	}
// 	std::cout<<"sino2sp succeed"<<std::endl;
	return res;
}

raft_image sp2sino(raft_image source)
{
	int Ns = source.data.lines;
	int Nt = source.data.columns;

	raft_image res;
	res.data = raft_matrix_create(2*Ns, Nt/2);
	res.tl_x = -pi;
	res.br_x = pi;
	res.tl_y = 1;
	res.br_y = -1;

	std::cout<<"Transforming polar to sinogram : ("<<Ns<<","<<Nt<<" ) -> (";
	std::cout<<2*Ns<<","<<Nt/2<<")."<<std::endl;
	for(int t=0;t<Nt/2;++t) {
		for(int s=0;s<Ns;++s) {
			raft_matrix_element(res.data, Ns + s, t) = raft_matrix_element(source.data, s, t);
			raft_matrix_element(res.data, Ns - s, t) = raft_matrix_element(source.data, s, t + Nt/2);
		}
	}
// 	std::cout<<"sino2sp succeed"<<std::endl;
	return res;
}

raft_image raft_image_rescale(raft_image source, double scale)
{
	raft_image res = source;
	int Nt = res.data.columns;
	int Ns = res.data.lines;
	for(int j=0; j<Nt; j++) {
		for(int i=0; i<Ns; i++) {
			raft_matrix_element(res.data, i, j) = 
				scale*raft_matrix_element(res.data, i, j);
		}
	}
	return res;
}	

raft_image pad_sector(raft_image sector)
{
	int Ns = sector.data.lines;
	int N = sector.data.columns;
	
	double t0 = sector.tl_x;
	double t1 = sector.br_x;
// 	if(t0 < 0) t0 += 2*pi;
// 	if(t1 < 0) t1 += 2*pi;
// 
// 	if(t0 > 2*pi) t0 -= 2*pi;
// 	if(t1 > 2*pi) t1 -= 2*pi;

	double dt = std::abs(t1 - t0)/N;
	int Nt = 2*pi/dt; 
	raft_image res = raft_image_create(Ns, Nt);

	res.tl_x = -pi;
	res.br_x = pi;
	res.br_y = sector.br_y;
	res.tl_y = sector.tl_y;

	int idx_start = (t0 + pi)/dt;
	int pos;

	std::cout<<"Padding the sector: ("<<Ns<<","<<N<<") -> (";
	std::cout<<Ns<<","<<Nt<<")."<<std::endl;
	for(int i=0; i<Ns; i++) {
		for(int j=0; j<N; j++) {
			pos = j+idx_start;
		
			if(pos< Nt-1 && pos >=0) {
				raft_matrix_element(res.data, i, pos) =
					raft_matrix_element(sector.data, i, j);
			} else if(pos < 0) {
				raft_matrix_element(res.data, i, Nt-pos) =
					raft_matrix_element(sector.data, i, j);
			} else if(pos > Nt-1) {
				raft_matrix_element(res.data, i, pos - Nt) =
					raft_matrix_element(sector.data, i, j);
			}
		}
	}
	return res;

}	

raft_image pad_sector(raft_image sector, int Nt, double t0, double t1)
{
	int Ns = sector.data.lines;
	int N = sector.data.columns;
	raft_image res = raft_image_create(Ns, Nt);
	res.tl_x = -pi;
	res.br_x = pi;

	res.br_y = sector.br_y;
	res.tl_y = sector.tl_y;
	double dt = (t1 - t0)/N;

	int idx_start = (t0 + pi)/dt;
	int pos;

	std::cout<<"Padding the sector: t0="<<t0<<"t1="<<t1<<std::endl;;
	for(int i=0; i<Ns; i++) {
		for(int j=0; j<N; j++) {
			pos = j+idx_start;
		
			if(pos< Nt-1 && pos >=0) {
				raft_matrix_element(res.data, i, pos) =
					raft_matrix_element(sector.data, i, j);
			} else if(pos < 0) {
				raft_matrix_element(res.data, i, Nt-pos) =
					raft_matrix_element(sector.data, i, j);
			} else if(pos > Nt-1) {
				raft_matrix_element(res.data, i, pos - Nt) =
					raft_matrix_element(sector.data, i, j);
			}
		}
	}
	return res;

}	

raft_image get_sector(raft_image source, double t0, double t1)
{
	int Ns = source.data.lines;
	int Nt = source.data.columns;
	
	double old_t0 = std::min(source.br_x, source.tl_x);
	double old_t1 = std::max(source.br_x, source.tl_x);
	double dt = (old_t1 - old_t0)/Nt;

	std::cout<<"_____________ GET SECTOR _____________"<<std::endl;
	std::cout<<"Old angles: t0="<<old_t0<<"; t1="<<source.tl_x<<std::endl;
	std::cout<<"New angles: t0="<<t0<<"; t1="<<t1<<std::endl;
	std::cout<<"dt="<<dt<<std::endl;


	raft_image res;
	if(t0 > -pi && t1 < pi) {
		int idx_start = floor((t0 - old_t0)/dt);
		int idx_end = floor((t1 - old_t0)/dt);
		int N = idx_end - idx_start;
		std::cout<<"!!!!!!!!!!!!!!!!   "<<idx_start<<", "<<idx_end<<" !!!!!!"<<std::cout;
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
		std::cout<<"2 case: idx1="<<idx1<<"; idx2="<<idx2<<std::endl;
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
		std::cout<<"3 case: idx1="<<idx1<<"; idx2="<<idx2<<std::endl;
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
		std::cout<<"4 case: idx1="<<idx1<<"; idx2="<<idx2<<std::endl;
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
	std::cout<<"get_sector suceed"<<std::endl;
	return res;
}

raft_image filter_sector(raft_image source, double t0, double t1)
{
	int Ns = source.data.lines;
	int Nt = source.data.columns;
	
	double old_t0 = std::min(source.br_x, source.tl_x);
	double old_t1 = std::max(source.br_x, source.tl_x);
	double dt = (old_t1 - old_t0)/Nt;

	std::cout<<"_____________ GET SECTOR _____________"<<std::endl;
	std::cout<<"Old angles: t0="<<old_t0<<"; t1="<<source.tl_x<<std::endl;
	std::cout<<"New angles: t0="<<t0<<"; t1="<<t1<<std::endl;
	std::cout<<"dt="<<dt<<std::endl;

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
	std::cout<<"get_sector suceed"<<std::endl;
	return res;
}

raft_image rotate(raft_image source, double t)
{
	int Ns = source.data.lines;
	int Nt = source.data.columns;

	std::cout<<"Rotating image on "<<t<<"; Ns="<<Ns<<"; Nt="<<Nt<<std::endl;
	raft_image res;
	
	double t0 = source.tl_x;
	double t1 = source.br_x;


	double dt = fabs(t1 - t0)/Nt;
	int idxt =floor(t/dt);

	raft_matrix tmp = raft_matrix_create(Ns, Nt);
	for(int j=0;j<Nt;j++) {
		for(int i=0;i<Ns; ++i) {
// 			std::cout<<"("<<i<<","<<j<<")"<<std::endl;
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

// 			if(j+idxt < Ns && j+idxt>=0) {
// 				raft_matrix_element(tmp, i, j + idxt) = raft_matrix_element(source.data, i, j);
// 			} else if(j+idxt>Nt-1) {
// 				raft_matrix_element(tmp, i, j+idxt - Nt) = raft_matrix_element(source.data, i, j);
// 			} else if(j+idxt<0) {
// 				raft_matrix_element(tmp, i, j+idxt + Nt - 1) = raft_matrix_element(source.data, i, j);
// 			}
// 		 
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

	std::cout<<"________cutting______________"<<std::endl;
	std::cout<<"old tl_x="<<x0<<"; old br_x="<<x1<<"; new tl_x="<<tl_x<<"; new br_x="<<br_x<<std::endl;
	std::cout<<"old br_y="<<y0<<"; old tl_y="<<y1<<"; new br_y="<<br_y<<"; new tl_y="<<tl_y<<std::endl;
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

	std::cout<<"idx_x0="<<idx_x0<<"; idx_x1="<<idx_x1<<" idx_y0="<<idx_y0<<"; idx_y1="<<idx_y1<<std::endl;
	std::cout<<"old size: <"<<Nx<<" x "<<Nx<<">; dx="<<dx<<"; dy="<<dy<<std::endl;
	std::cout<<"new size: <"<<nx<<" x "<<ny<<">"<<std::endl; 

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

raft_image zero_padding(raft_image source, int Ns, int Nt)
{
	if(Ns < source.data.lines || Nt < source.data.columns) return source;
	double dt = fabs(source.tl_x - source.br_x)/source.data.columns;
	double ds = fabs(source.br_y - source.tl_y)/source.data.lines;


	std::cout<<"Zerro padding started:"<<std::endl;
	std::cout<<"dt="<<dt<<"; ds="<<ds<<std::endl;;
	std::cout<<"OLD size: lines="<<source.data.lines<<"; columns="<<source.data.columns<<std::endl;
	std::cout<<"NEW size: lines="<<Ns<<"; columns="<<Nt<<std::endl;
	int start_j = (Nt - source.data.columns);
	int start_i = (Ns - source.data.lines);
	
	double new_min_x = source.tl_x - start_j*dt/2;
	double new_min_y = source.br_y - start_i*ds/2;
	double new_max_x = new_min_x + Nt*dt ;
	double new_max_y = new_min_y + Ns*ds;

	raft_image res;
	res.tl_x = new_min_x;
	res.br_x = new_max_x;
	res.br_y = new_min_y;
	res.tl_y = new_max_y;
	res.data = raft_matrix_create(Ns, Nt);
	for(int j=0; j<Nt - start_j; j++) {
		for(int i=0; i<Ns - start_i; i++) {
			raft_matrix_element(res.data, i + start_i/2, j + start_j/2) = raft_matrix_element(source.data, i, j);	
		}
	}	


// 	std::cout<<"OLD range ("<<source.br_x<<", "<<source.tl_x<<") X ("<<source.tl_y<<", "<<source.br_y<<")."<<std::ends;
// 	std::cout<<"NEW range ("<<source.br_x<<", "<<new_max_x<<") X ("<<source.tl_y<<", "<<new_max_y<<")."<<std::ends;
	res.br_x = new_max_x;
	res.tl_y = new_max_y;
	std::cout<<"zerro padding succeed"<<std::ends;

	return res;
}

raft_image zero_padding_on_s(raft_image source, int Ns)
{
	int ns = source.data.lines;
	int nt = source.data.columns;
	if(Ns <= ns) return source;
	
	double ds = fabs(source.br_y - source.tl_y)/ns;

	double s_max = source.br_y + Ns*ds;

	raft_image res;
	res.br_x = source.br_x;
	res.tl_x = source.tl_x;
	res.br_y = source.br_y;
	res.tl_y = s_max;


	res.data = raft_matrix_create(Ns, nt);
	std::cout<<"Zerro padding on s : ("<<ns<<","<<nt<<" ) -> (";
	std::cout<<Ns<<","<<nt<<")."<<std::endl;
	for(int j=0; j<nt; j++) {
		for(int i=0; i<ns; i++) {
			raft_matrix_element(res.data, i, j) = raft_matrix_element(source.data, i, j);
		}	
	}
	return res;

}

#include <string.h>
void zero_padding_on_s(raft_image source, raft_image res_, int Ns)
{
	int ns = source.data.lines;
	int nt = source.data.columns;
	if(Ns <= ns) return;
	
	double ds = fabs(source.br_y - source.tl_y)/ns;

	double s_max = source.br_y + Ns*ds;

	raft_image res;
	res.br_x = source.br_x;
	res.tl_x = source.tl_x;
	res.br_y = source.br_y;
	res.tl_y = s_max;


	res.data = raft_matrix_create(Ns, nt);
	std::cout<<"Zerro padding on s : ("<<ns<<","<<nt<<" ) -> (";
	std::cout<<Ns<<","<<nt<<")."<<std::endl;
	for(int j=0; j<nt; j++) {
		for(int i=0; i<ns; i++) {
			raft_matrix_element(res.data, i, j) = raft_matrix_element(source.data, i, j);
		}	
	}
	memcpy((double*)res_.data.p_data, (double*)res.data.p_data, sizeof(double)*Ns*nt);
	raft_image_destroy(&res);
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

	std::cout<<"Zerro padding started:"<<std::endl;
	std::cout<<"; ds="<<ds<<std::endl;;
	std::cout<<"OLD size: lines="<<source.data.lines<<"; columns="<<source.data.columns<<std::endl;
	std::cout<<"NEW size: lines="<<Ns<<"; columns="<<nt<<std::endl;

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

	std::cout<<"___________ semi-polar 2 cartesian interpolation ________"<<std::endl;
	std::cout<<"old size: Ns="<<Ns<<"; Nt="<<Nt<<std::endl;
	std::cout<<"new size: Nx="<<Nx<<"; Ny="<<Ny<<std::endl;
	std::cout<<"s0="<<s0<<"; s1="<<s1<<"; ds="<<ds<<std::endl;
	std::cout<<"t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
	std::cout<<"x0="<<x0<<"; x1="<<x1<<"; dx="<<dx<<std::endl;
	std::cout<<"y0="<<y0<<"; y1="<<y1<<"; dy="<<dy<<std::endl;

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

	double t0 = -pi;
	double t1 = pi;
	double dt = 2*pi/Nt;

// 	std::cout<<"___________ cartesian 2 semi-polar interpolation ________"<<std::endl;
// 	std::cout<<"old size: Nx="<<Nx<<"; Ny="<<Ny<<std::endl;
// 	std::cout<<"new size: Nt="<<Nt<<"; Ns="<<Ns<<std::endl;
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
		double _cos = cos(t);
		double _sin = sin(t);
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

	std::cout<<"___________ cartesian 2 semi-polar interpolation ________"<<std::endl;
	std::cout<<"old size: Nx="<<Nx<<"; Ny="<<Ny<<std::endl;
	std::cout<<"new size: Nt="<<Nt<<"; Nr="<<Nr<<std::endl;
	std::cout<<"r0="<<r0<<"; r1="<<r1<<"; dr="<<dr<<std::endl;
	std::cout<<"t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
	std::cout<<"x0="<<x0<<"; x1="<<x1<<"; dx="<<dx<<std::endl;
	std::cout<<"y0="<<y0<<"; y1="<<y1<<"; dy="<<dy<<std::endl;

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

void fft_shift_2d(double *x, int Nx, int Ny)
{
	double tmp;
	int idx;
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny/2;j++){
			idx = j*Nx+i;
			tmp = x[idx];
			x[idx] = x[(Ny/2+j)*Nx + i];
			x[(Ny/2+j)*Nx+i] = tmp;
		}
	}

	for(int i=0;i<Nx/2;i++){
		for(int j=0;j<Ny;j++){
			idx = j*Nx+i;
			tmp = x[idx];
			x[idx] = x[j*Nx + i+Nx/2];
			x[j*Nx+i+Nx/2] = tmp;
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

double *convolution_2d_C2R(double *x, double *k, int Nx, int Ny, int threads )
{
	std::cout<<"conv_2d starts: Mx="<<Nx<<"; Ny="<<Ny<<std::endl;
	int data_size = Nx*Ny;
	int spectra_size = Nx*(Ny/2+1);

	double scale = 1.0/spectra_size;
	fftw_complex *spectrum_data = 	(fftw_complex*) fftw_malloc(spectra_size*sizeof(fftw_complex));
	fftw_complex *spectrum_kernel = (fftw_complex*) fftw_malloc(spectra_size*sizeof(fftw_complex));

	fftw_init_threads();
	fftw_plan_with_nthreads(threads);

	fftw_plan p = fftw_plan_dft_r2c_2d(Nx, Ny, x, spectrum_data, 
			FFTW_ESTIMATE);
	fftw_execute_dft_r2c(p, x, spectrum_data);
	fftw_execute_dft_r2c(p, k, spectrum_kernel);

	std::cout<<"conv2d: straight FFT(x) succeed"<<std::endl;
	double xr, xi, kr, ki;
	for(int i=0; i<spectra_size; i++) {
		xr = spectrum_data[i][0];
		xi = spectrum_data[i][1];

		kr = spectrum_kernel[i][0];
		ki = spectrum_kernel[i][1];
			
		spectrum_data[i][0] = 	(xr*kr - xi*ki)*scale;
		spectrum_data[i][1] = 	(xr*ki + xi*kr)*scale;
	}
	std::cout<<"conv2d: Multiplication succeed"<<std::endl;
	fftw_plan pinv = fftw_plan_dft_c2r_2d(Nx, Ny, spectrum_data, x, FFTW_ESTIMATE);
	fftw_execute_dft_c2r(pinv, spectrum_data, x);
	std::cout<<"conv2d: inverse fourier procedure succeed"<<std::endl;
	fftw_destroy_plan(p);
	fftw_destroy_plan(pinv);
	fftw_free(spectrum_data);
	fftw_free(spectrum_kernel);

	std::cout<<"conv2d: SUCCEED"<<std::endl;
	return x;
}

fftw_complex *fftw_2d_R2C(double *x, int Nx, int Ny, int threads)
{
	std::cout<<"ffw_2d_r2c starts: Mx="<<Nx<<"; Ny="<<Ny<<std::endl;
	double scale = 1.0/(Nx*Ny);
	int data_size = Nx*Ny;
	int spectra_size = Nx*(Ny/2+1);

	fftw_complex *spectrum_data = 	(fftw_complex*) fftw_malloc(spectra_size*sizeof(fftw_complex));
	fftw_init_threads();
	fftw_plan_with_nthreads(threads);

	fftw_plan p = fftw_plan_dft_r2c_2d(Nx, Ny, x, spectrum_data, 
			FFTW_ESTIMATE);
	fftw_execute_dft_r2c(p, x, spectrum_data);
	return spectrum_data;

}

double *ifftw_2d_C2R(fftw_complex *sp, int Nx, int Ny, int nthreads)
{
	double scale = 1.0/(Nx*Ny);
	int data_size = Nx*Ny;
	int spectra_size = Nx*(Ny/2+1);

	double *res = new double[data_size];

	fftw_init_threads();
	fftw_plan_with_nthreads(nthreads);

	fftw_plan pinv = fftw_plan_dft_c2r_2d(Nx, Ny, sp, res, FFTW_ESTIMATE);
	fftw_execute_dft_c2r(pinv, sp, res);
	fftw_destroy_plan(pinv);
	return res;
	
}

double *semi_convolution_2d_C2R(double *x, 
		fftw_complex *spectrum_kernel, int Nx, int Ny, int threads)
{
	std::cout<<"conv_2d starts: Mx="<<Nx<<"; Ny="<<Ny<<std::endl;
	double scale = 1.0/(Nx*Ny);
	int data_size = Nx*Ny;
	int spectra_size = Nx*(Ny/2+1);

	fftw_complex *spectrum_data = 	(fftw_complex*) fftw_malloc(spectra_size*sizeof(fftw_complex));

	fftw_init_threads();
	fftw_plan_with_nthreads(threads);

	fftw_plan p = fftw_plan_dft_r2c_2d(Nx, Ny, x, spectrum_data, 
			FFTW_ESTIMATE);
	fftw_execute_dft_r2c(p, x, spectrum_data);

	std::cout<<"conv2d: straight FFT(x) succeed"<<std::endl;
	double xr, xi, kr, ki;
	for(int i=0; i<spectra_size; i++) {
		xr = spectrum_data[i][0];
		xi = spectrum_data[i][1];

		kr = spectrum_kernel[i][0];
		ki = spectrum_kernel[i][1];
			
		spectrum_data[i][0] = 	(xr*kr - xi*ki)*scale;
		spectrum_data[i][1] = 	(xr*ki + xi*kr)*scale;
	}
	std::cout<<"conv2d: Multiplication succeed"<<std::endl;
	fftw_plan pinv = fftw_plan_dft_c2r_2d(Nx, Ny, spectrum_data, x, FFTW_ESTIMATE);
	fftw_execute_dft_c2r(pinv, spectrum_data, x);
	std::cout<<"conv2d: inverse fourier procedure succeed"<<std::endl;
	fftw_destroy_plan(p);
	fftw_destroy_plan(pinv);
	fftw_free(spectrum_data);

	std::cout<<"conv2d: SUCCEED"<<std::endl;
	return x;
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

// 	fft_shift_2d(res);
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

	std::cout<<"deconv_2d starts: Mx="<<Nx<<"; Ny="<<Ny<<std::endl;
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
	std::cout<<"FFT of the kernel done"<<std::endl;
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
	std::cout<<"iFFTW 2D started"<<std::endl;
	int N = xr.lines;
	int M = xr.columns;
	fftw_complex *in, *out;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*M);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*M);
	int idx; 
	for(int j=0; j<M; j++) {
		for(int i=0; i<N; i++) {
			idx = j*N + i;
			in[idx][0] = raft_matrix_element(xr, i, j);	
			in[idx][1] = raft_matrix_element(xi, i, j);	

		}
	}
	fftw_init_threads();
	fftw_plan_with_nthreads(threads);
	fftw_plan p = fftw_plan_dft_2d(N, M, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p); /* repeat as needed */
	fftw_destroy_plan(p);
	fftw_free(in);
// 	double factor = 1.0/(N*M);
	double factor = 1.0;
	for(int j=0; j<M; j++) {
		for(int i=0; i<N; i++) {
			idx = j*N + i;
			raft_matrix_element(xr, i, j) = out[idx][0]*factor;
			raft_matrix_element(xi, i, j) = out[idx][1]*factor;
// 			raft_matrix_element(xr, i, j) = out[idx][0]/(2*pi);
// 			raft_matrix_element(xi, i, j) = out[idx][1]/(2*pi);
		}
	}
        //fftw_cleanup_threads();
	fftw_free(out);
	std::cout<<"iFFT_2D succeed"<<std::endl;
}

raft_matrix ifftw_2d_c2r(raft_matrix &xr, raft_matrix &xi, int threads)
{
	std::cout<<"iFFT_2D_c2r started"<<std::endl;
	int Nx = xr.columns;
	int Ny = xr.lines;
	std::cout<<"Nx="<<Nx<<"; N="<<Ny<<std::endl;

	fftw_complex *in;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	
	
	double *out = (double*)fftw_malloc(sizeof(double)*Nx*(2*Ny+1));
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
	fftw_plan p = fftw_plan_dft_c2r_2d(Nx, Ny, in, out, FFTW_ESTIMATE);
	fftw_execute(p); /* repeat as needed */
	fftw_destroy_plan(p);
	fftw_free(in);
	raft_matrix res_r = raft_matrix_create(Nx, 2*Ny+1);
	for(int i=0; i<Nx; i++) {
		for(int j=0; j<2*Ny+1; j++) {
			idx = j*Nx + i;
			raft_matrix_element(res_r, i, j) = out[idx];
		}
	}
	fftw_free(out);
// 	std::cout<<"iFFT_2D succeed"<<std::endl;
}

void fftw_2d(raft_matrix &xr, raft_matrix &xi, int threads)
{
	std::cout<<"FFT_2D started"<<std::endl;
	int N = xr.lines;
	int M = xr.columns;
// 	std::cout<<"N="<<N<<"; N="<<M<<std::endl;

	fftw_complex *in, *out;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*M);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*M);
	int idx; 
	for(int j=0; j<M; j++) {
		for(int i=0; i<N; i++) {
			idx = j*N + i;
			in[idx][0] = raft_matrix_element(xr, i, j);	
			in[idx][1] = raft_matrix_element(xi, i, j);	

		}
	}
	fftw_init_threads();
	fftw_plan_with_nthreads(threads);
	fftw_plan p = fftw_plan_dft_2d(N, M, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p); /* repeat as needed */
	fftw_destroy_plan(p);
	fftw_free(in);
	double factor = 1.0/(N*M);
	for(int j=0; j<M; j++) {
		for(int i=0; i<N; i++) {
			idx = j*N + i;
			raft_matrix_element(xr, i, j) = out[idx][0]*factor;
			raft_matrix_element(xi, i, j) = out[idx][1]*factor;
		}
	}
	fftw_free(out);
	std::cout<<"FFT_2D succeed"<<std::endl;
}





 	

// Sinogram transformation with moving the source (Radon shift property, polar coords)
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

	std::cout<<"__________ Moving the origin on the sinogram ________"<<std::endl;
	std::cout<<"sizes: Nt="<<Nt<<"; Ns="<<Ns<<std::endl;
	std::cout<<"k="<<k<<"; old_s1="<<old_s1<<"; s1="<<s1<<"; ds="<<ds<<std::endl;
	std::cout<<"t0="<<t0<<"; t1="<<t1<<"; dt="<<dt<<std::endl;
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
	std::cout<<"_________mo_sino succeed__________"<<std::endl;
	return res;

}

// Simple moving in Cartesian coordinates
raft_image mo_dec(raft_image source, double dist_x, double dist_y)
{
	// Size 
	int Nx = source.data.lines;  	
	int Ny = source.data.columns;
	// Corners
	double x0 = source.tl_x;
	double x1 = source.br_x;
	double y0 = source.br_y;
	double y1 = source.tl_y;
	// Mesh step
	double dx = (x1 - x0)/Nx;
	double dy = (y1 - y0)/Ny;
	// Number of pixels for moving 
	int di = dist_x/dx;
	int dj = dist_y/dy;
	// Result image
	raft_image res;
	res.data = raft_matrix_create(Nx, Ny);

	res.tl_x = x0;
	res.br_x = x1;
	res.tl_y = y1;
	res.br_y = y0;

	std::cout<<"_____ mo_dec: x0="<<x0<<"; x1="<<x1<<"; y0="<<y0<<"; y1="<<y1<<std::endl;
	std::cout<<"dx="<<dx<<"; dy="<<dy<<"; di="<<di<<"; dj="<<dj<<std::endl;

	int old_i, old_j;
	double val;
	for(int i=0; i<Nx; i++) {
		for(int j=0; j<Ny; j++) {
			old_i = i-di;
			old_j = j-dj;
			if(old_i > Nx-1 || old_j >Ny-1 || old_i<0 || old_j < 0 ) {
				val = 0;
			} else {
				val = raft_matrix_element(source.data, old_i, old_j);
			}
			raft_matrix_element(res.data, i, j) = val;
		}
	}

	return res;
}












