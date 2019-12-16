
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include <helper_cuda.h>

#define CCE checkCudaErrors

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include "pimage_local.h"

#include "../pimage.h"

using namespace std;

const int gl_win_size = 800;

const double phong_exponent = 1.5;

const double zeye_global = 25;

__device__ int akt_deg_global;


__device__ int gl_max_expo[3];


#include "ppoly.cu"

 __device__ int cudaPoly3::degree()
{
	int i;
	int deg = 0;
	int deg_akt;
	for(i= 0; i < len; ++i) {
		deg_akt = xexpo[i] + yexpo[i] + zexpo[i];
		if (deg_akt > deg) {
			deg = deg_akt;
		}
	}
	return deg;
}

__device__ cudaPoly3 & cudaPoly3::diff(int iv)
{
	int i;
	for(i = 0; i < len; ++i) {
		switch (iv) {
		case 0: coefs[i] *= xexpo[i];
				if (xexpo[i] > 0)--xexpo[i];
				break;
		case 1: coefs[i] *= yexpo[i];
				if (yexpo[i] > 0)--yexpo[i];
				break;
		case 2: coefs[i] *= zexpo[i];
				if (zexpo[i] > 0) --zexpo[i];
				break;
		default: break;
		}
	}

	return *this;

}

__device__ double cudaPoly3::eval(double* sl)
{
	double val = 0;
	int i;

	const int max_expo = 15;

	double xx, yy, zz;
	double xl[max_expo], yl[max_expo], zl[max_expo];

	xx = 1;
	yy = 1;
	zz = 1;

	for(i = 0; i <= gl_max_expo[0]; ++i) {
		xl[i] = xx;
		xx *= sl[0];
	}

	for(i = 0; i <= gl_max_expo[1]; ++i) {
		yl[i] = yy;
		yy *= sl[1];
	}

	for(i = 0; i <= gl_max_expo[2]; ++i) {
		zl[i] = zz;
		zz *= sl[2];
	}

	for(i = 0; i < len; ++i) {
		if (coefs[i] != 0) {
			val += coefs[i] * xl[xexpo[i]] * yl[yexpo[i]] * zl[zexpo[i]];

			//val += coefs[i] * pow(sl[0], xexpo[i]) * pow(sl[1], yexpo[i]) * pow(sl[2], zexpo[i]);
		}
	}
	return val;
}

__device__ void print_cudaPoly3(cudaPoly3 & pol)
{
	for (int i = 0; i < pol.len; ++i) {
		printf("%f * x^%d y^%d z^%d + ", pol.coefs[i], pol.xexpo[i], pol.yexpo[i], pol.zexpo[i]);
	}
	printf("\n deg f = %d \n", pol.degree());
}






inline __device__ int & mref(int* mat, int i, int j)
{
	return mat[2 * i * gl_win_size + j];
}


#define EPS 1e-12
#define M_INF -1e38

__constant__ __device__ cudaPoly3 f3;

static __device__ cudaPoly3 f3x, f3y, f3z;

static __device__ double m_euler[3 * 3];

const int max_deg = 20;

static __device__ int akt_deg;

static __device__ double akt_xbase[max_deg+1] = { 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 };

static __device__ double lagrange_basis[max_deg+1][max_deg+1];


__device__ double cuda_skalp_3(double* x, double* y)
{
	return x[0]*y[0] + x[1] * y[1] + x[2] * y[2];
}


__device__ void cuda_skalmult_3(double c, double *x)
{
	x[0] *= c;
	x[1] *= c;
	x[2] *= c;
}

__device__ void cuda_norml_3(double* x)
{
	double len = sqrt(cuda_skalp_3(x,x));
	if (len > 0) {
		cuda_skalmult_3(1/len, x);
	}
}

__device__ void cuda_vecadd_3(double * z, double *x, double *y)
{
	z[0] = x[0] + y[0];
	z[1] = x[1] + y[1];
	z[2] = x[2] + y[2];

}

__device__ void cuda_vecsub_3(double * z, double *x, double *y)
{
	z[0] = x[0] - y[0];
	z[1] = x[1] - y[1];
	z[2] = x[2] - y[2];

}



__device__ double cuda_zscale_factor(double z)
{
	return ((z - zeye_global)/(0 - zeye_global));
}


__device__ void cuda_lb_poly_mult(int deg, double * coefs, double a)
{
	// multiply coefs poly coefs[0] * t^deg + ... with (t + a)
	coefs[deg+1] = 0;
	for(int i = deg + 1; i >=  1; --i) {
		coefs[i] += a * coefs[i-1];
	}
}

__device__ void cuda_lb_poly_scal_mult(int deg, double * coefs, double a)
{
	for(int i = 0; i <= deg; ++i) {
		coefs[i] *= a;
	}
}


__device__ void cuda_lb_gen_lagrange_basis(int deg, double* xbase, double lagr_basis[][max_deg + 1])
{
	// xbase stützstellen, deg+1 Stück

	for(int i = 0; i <= deg; ++i) {
		double aux_poly[max_deg + 1];

		aux_poly[0] = 1;
		double pre_coef = 1;
		int curr_deg = 0;

		for(int j = 0; j <= deg; ++j) {
			if (i == j) {
				continue;
			}
			pre_coef *= (xbase[i] - xbase[j]);
			cuda_lb_poly_mult(curr_deg, aux_poly, -xbase[j]);
			++curr_deg;
		}
		cuda_lb_poly_scal_mult(deg, aux_poly, 1/pre_coef);

		memcpy(lagr_basis[i], aux_poly, sizeof(aux_poly));

	}
}

__device__ void cuda_m_mult_vec_mat(double* vec, double *mat, double * vec_res) {
	vec_res[0] = vec[0] * mat[0] + vec[1] * mat[3] + vec[2] * mat[6];
	vec_res[1] = vec[0] * mat[1] + vec[1] * mat[4] + vec[2] * mat[7];
	vec_res[2] = vec[0] * mat[2] + vec[1] * mat[5] + vec[2] * mat[8];
}

__device__ void cuda_m_mult_mat_vec(double* vec, double *mat, double * vec_res) {
	vec_res[0] = vec[0] * mat[0] + vec[1] * mat[1] + vec[2] * mat[2];
	vec_res[1] = vec[0] * mat[3] + vec[1] * mat[4] + vec[2] * mat[5];
	vec_res[2] = vec[0] * mat[6] + vec[1] * mat[7] + vec[2] * mat[8];
}


__device__ void cuda_m_rot_z(double a, double * mat)
{
	const int ndim = 3;
	for(int i = 0; i < ndim; ++i) {
		for(int j = 0; j < ndim; ++j) {
			mat[3*i+j] = 0;
		}
	}
	mat[0] = cos(a);
	mat[1] = sin(a);
	mat[3] = -sin(a);
	mat[4] = cos(a);
	mat[8] = 1;
}

__device__ void cuda_m_rot_x(double a, double * mat)
{
	const int ndim = 3;
	for(int i = 0; i < ndim; ++i) {
		for(int j = 0; j < ndim; ++j) {
			mat[3*i+j] = 0;
		}
	}
	mat[4] = cos(a);
	mat[5] = sin(a);
	mat[7] = -sin(a);
	mat[8] = cos(a);
	mat[0] = 1;

}

__device__ void cuda_m_mult(double * m1, double * m2, double * m3)
{
	const int ndim = 3;
	for(int i = 0; i < ndim; ++i) {
		for(int j = 0; j < ndim; ++j) {
			double sum = 0;
			for(int k = 0; k < ndim; ++k) {
				sum += m1[i*3 + k] * m2[k*3 + j];
			}
			m3[i*3 + j] = sum;
		}
	}
}

__device__ void cuda_m_print(double* m)
{
	for(int i = 0; i < 3; ++i) {
		for(int j = 0; j < 3; ++j) {
			printf("%f , ", m[i*3+j]);
		}
		printf("\n");
	}
}

__device__ double cuda_eval_f3_mat(cudaPoly3 & f3, double *sl)
{
	double sl1[3];
	cuda_m_mult_mat_vec(sl, m_euler, sl1);

	double f = f3.eval(sl1);

	return f;

}


__device__ int cuda_eval_poly_poly(cudaPoly3 & f3,
		double x, double y, double z, double & f, double & fx, double & fy, double & fz)
{

	// hier muss
	// 1) sl mit m_euler substituiert werden vor dem Einsetzen
	// 2) (fx fy fz) * m_euler multipliziert werden

	double sl[3];
	double sl1[3];
	sl[0] = x;
	sl[1] = y;
	sl[2] = z;

	cuda_m_mult_mat_vec(sl, m_euler, sl1);

	f = f3.eval(sl1);

	sl[0] = f3x.eval(sl1);
	sl[1] = f3y.eval(sl1);
	sl[2] = f3z.eval(sl1);

	cuda_m_mult_vec_mat(sl, m_euler, sl1);

	fx = sl1[0];
	fy = sl1[1];
	fz = sl1[2];

	return 0;
}

#if 0

__device__ int cuda_eval_coefs_poly(cudaPoly3 & f3, double x, double y, double & a0, double & a1, double & a2, double & a3,
		double & a4)
{
	// hier die Lagrangeinterpolation einfügen
	// z wird nacheinander z0, z1, z2,....zn gesetzt mit n = deg f3
	// m_euler wird zur Substitution der x y z benutzt

	double sl[3];
	sl[0] = x;
	sl[1] = y;
	double coefs[max_deg + 1];
	for(int i = 0; i <= max_deg; ++i) {
		coefs[i] = 0;
	}
	for(int i = 0; i <= akt_deg; ++i) {
		double z = akt_xbase[i];
		sl[2] = z;
		double fi = cuda_eval_f3_mat(f3, sl);
		for(int j = 0; j <= akt_deg; ++j) {
			coefs[j] += fi * lagrange_basis[i][j];
		}
	}

	a0 = coefs[4];
	a1 = coefs[3];
	a2 = coefs[2];
	a3 = coefs[1];
	a4 = coefs[0];
	return 0;
}

#endif

__device__ int cuda_eval_coefs_poly(cudaPoly3 & f3, double x, double y, double* coefs_lis)
{
	// hier die Lagrangeinterpolation einfügen
	// z wird nacheinander z0, z1, z2,....zn gesetzt mit n = deg f3
	// m_euler wird zur Substitution der x y z benutzt

	double sl[3];
	sl[0] = x;
	sl[1] = y;
	double coefs[max_deg + 1];
	for(int i = 0; i <= max_deg; ++i) {
		coefs[i] = 0;
	}
	for(int i = 0; i <= akt_deg_global; ++i) {
		double z = akt_xbase[i];
		sl[2] = z;
		double fi = cuda_eval_f3_mat(f3, sl);
		for(int j = 0; j <= akt_deg_global; ++j) {
			coefs[j] += fi * lagrange_basis[i][j];
		}
	}

	memcpy(coefs_lis, coefs, sizeof(coefs));

	return 0;
}


__device__ int cuda_eval_coefs_poly_centr_persp(cudaPoly3 & f3, double x, double y, double* coefs_lis)
{
	// hier die Lagrangeinterpolation einfügen
	// z wird nacheinander z0, z1, z2,....zn gesetzt mit n = deg f3
	// m_euler wird zur Substitution der x y z benutzt

	double sl[3];
	double coefs[max_deg + 1];
	for(int i = 0; i <= max_deg; ++i) {
		coefs[i] = 0;
	}
	for(int i = 0; i <= akt_deg_global; ++i) {
		double z = akt_xbase[i];

		double zscale_factor = cuda_zscale_factor(z);

		sl[0] = x * zscale_factor;
		sl[1] = y * zscale_factor;
		sl[2] = z;
		double fi = cuda_eval_f3_mat(f3, sl);
		for(int j = 0; j <= akt_deg_global; ++j) {
			coefs[j] += fi * lagrange_basis[i][j];
		}
	}

	memcpy(coefs_lis, coefs, sizeof(coefs));

	return 0;
}





__global__ void cuda_rotate_mat(double phi, double theta, double psi)
{
	const int ndim = 3;
	const int nsize = ndim * ndim;

	double m_z_phi[nsize];
	double m_x_theta[nsize];
	double m_z_psi[nsize];

	double m_aux[nsize];

	cuda_m_rot_z(phi, m_z_phi);
	cuda_m_rot_x(theta, m_x_theta);
	cuda_m_rot_z(psi, m_z_psi);

	cuda_m_mult(m_x_theta, m_z_phi, m_aux);
	cuda_m_mult(m_z_psi, m_aux, m_euler);


}

__global__ void cuda_init_f3_diff()
{
	f3x = f3;
	f3x.diff(0);
	f3y = f3;
	f3y.diff(1);
	f3z = f3;
	f3z.diff(2);

/*
	printf("f3x = ");
	print_cudaPoly3(f3x);
	printf("f3y = ");
	print_cudaPoly3(f3y);
	printf("f3z = ");
	print_cudaPoly3(f3z);
*/

	akt_deg = f3.degree();

	cuda_lb_gen_lagrange_basis(akt_deg, akt_xbase, lagrange_basis);

}




const double clip_radius = 8; //20;

__device__ inline bool cuda_in_clip_radius(double x, double y, double z)
{
	return x * x + y * y + z * z <= clip_radius * clip_radius;
}

__device__ void cuda_normalize_poly_coefs(double* poly_coefs, int akt_deg, int & deg_new)
{

	// in poly_coefs poly_coefs[0] is coefficient in term of highest degree

	int deg;
	int i_tst;

	for(i_tst = 0; i_tst < akt_deg; ++i_tst) {
		if (fabs(poly_coefs[i_tst]) > 1e-8)
			break;
	}

	deg = akt_deg - i_tst;

	deg_new = deg;

	if (deg < akt_deg) {
		for(int i = 0; i < max_deg - (akt_deg - deg); ++i ) {
			poly_coefs[i] = poly_coefs[i + akt_deg - deg];
		}
	}

}


__device__ int cuda_get_z_intersect_poly(double x, double y, double *z, double *n_xyz, bool & disc_zero)
{

#if 0
	double a0 = 0;
	double a1 = 0;
	double a2 = 0;
	double a3 = 0;
	double a4 = 0;
	double poly_coefs[11];

	double z_erg;
	double z_erg_new = M_INF;
	double z_erg_list[11];

	cuda_eval_coefs_poly(f3, x, y, a0, a1, a2, a3, a4);

	poly_coefs[0] = a4;
	poly_coefs[1] = a3;
	poly_coefs[2] = a2;
	poly_coefs[3] = a1;
	poly_coefs[4] = a0;

	int deg = (a4 != 0) ? 4 : (a3 != 0) ? 3 : (a2 != 0) ? 2 : 1;
	int num_z_erg;
#endif

	double poly_coefs[max_deg + 1];

	double z_erg;
	double z_erg_new = M_INF;
	double z_erg_list[max_deg + 1];

	//cuda_eval_coefs_poly(f3, x, y, poly_coefs);
	cuda_eval_coefs_poly_centr_persp(f3, x, y, poly_coefs);


	// coefficient of leading monomial is in poly_coefs[0]
	// akt_deg is intended degree

	int deg;
	int num_z_erg;

	cuda_normalize_poly_coefs(poly_coefs, akt_deg_global, deg);


#if 0
	if (deg < 4) {
		for(int i = 0; i < 6; ++i ) {
			poly_coefs[i] = poly_coefs[i + 4 - deg];
		}
	}

	double disc_poly = x_comp_disc(deg, poly_coefs);

#endif

	if (1) {

		x_root_final_list(deg, poly_coefs, 20, z_erg_list, num_z_erg);

		int j;

		j = num_z_erg - 1;


		while (j >= 0) {

			double zsf = cuda_zscale_factor(z_erg_list[j]);

			double x1 = x * zsf;
			double y1 = y * zsf;

			if (cuda_in_clip_radius(x1,y1,z_erg_list[j])) {
				//cout << "num_z_erg = " << num_z_erg << " j = " << j << " z_erg_new = " << z_erg_list[j] << endl;
				z_erg_new = z_erg_list[j];
				break;
			}
			--j;
		}

	} else {
		z_erg_new = M_INF;
	}

	z_erg = z_erg_new;

	double zsf = cuda_zscale_factor(z_erg);

	double x1 = x * zsf;
	double y1 = y * zsf;



	if (! cuda_in_clip_radius(x1, y1, z_erg)) {
		z_erg = M_INF;
	}

	*z = z_erg;

	if (z_erg > M_INF) {
		double f, fx, fy, fz;
		cuda_eval_poly_poly(f3, x1, y1, z_erg, f, fx, fy, fz);
		n_xyz[0] = fx;
		n_xyz[1] = fy;
		n_xyz[2] = fz;
		cuda_norml_3(n_xyz);
		//*n_z = fz/sqrt(fx*fx+fy*fy+fz*fz);
	} else {
		n_xyz[0] = 0;
		n_xyz[1] = 0;
		n_xyz[2] = 0;
	}

	return 0;
}


const int win_size = gl_win_size;

#define RAND_COL ((int)(255.0*((float)rand())/RAND_MAX))

/*
double z_buf[win_size][win_size];
double n_buf[win_size][win_size];
*/
//#define SCALE 10.0
#define SCALE 50.0



__global__ void compute_colmat(double a, double b, int xmax, int ymax,
		double euler_phi, double euler_theta, double euler_psi,
		int *colmat_r_d, int *colmat_g_d, int *colmat_b_d)
{

	int xx = threadIdx.x + blockIdx.x * blockDim.x;
	int yy = threadIdx.y + blockIdx.y * blockDim.y;

	double n;
	double x1, y1;
	double z;
	double n_xyz[3];
	double local_scale;

	int win_offset;
	bool disc_zero;

	double color_red, color_green, color_blue;
	double phong_kernel, spec_coef;
	int col_red, col_green, col_blue;

	int col_z;

	win_offset = xmax/2;


	local_scale = gl_win_size/xmax;

	//y1 = (yy - win_offset)/SCALE * local_scale;
    //x1 = (xx - win_offset)/SCALE * local_scale;

	y1 = (yy - win_offset)/(0.8 * SCALE) * local_scale;
    x1 = (xx - win_offset)/(0.8 * SCALE) * local_scale;



    //printf("x1 = %f, y1 = %f ", x1, y1);

    cuda_get_z_intersect_poly(x1, y1, &z, n_xyz, disc_zero);


/*
    z_buf[x][y] = z;

    if (z > M_INF) {
    	n_buf[x][y] = n_z;
    } else {
    	n_buf[x][y] = 0;
    }
*/



    // n is cos(i) for light incident along z axis from infinity
	n = (z > M_INF) ? n_xyz[2] : 0;

	double Lin[3] = {0,0,1};

	double Lin_p[3];

	memcpy(Lin_p, n_xyz, sizeof(Lin_p));

	cuda_skalmult_3(2 * cuda_skalp_3(Lin, n_xyz), Lin_p);

	double Lout[3];

	cuda_vecsub_3(Lout, Lin_p, Lin);

	double zsf = cuda_zscale_factor(z);

	double x11 = x1 * zsf;
	double y11 = y1 * zsf;

	double Peye[3] = {0, 0, zeye_global};

	double Psurf[3] = {x11, y11, z};

	double Pvec[3];

	cuda_vecsub_3(Pvec, Peye, Psurf);

	cuda_norml_3(Pvec);

	double coss = cuda_skalp_3(Pvec, Lout);

	assert(fabs(coss) <= 1.0);

	if (z > M_INF) {
		color_red = 0.0;
		color_green = 0.0;
		color_blue = 0.0;

		if (n < 0) {
			color_red = -n/2 + 0.1;
			color_green = 0;
			color_blue = 0.0;
		} else if (n >= 0) {
			color_green = n/2 + 0.1;
			color_red = 0;
			color_blue = 0.0;
		}

		//phong_kernel = 2 * n * n - 1;
		//spec_coef = 0.3 * pow(phong_kernel, phong_exponent);

		spec_coef = coss > 0 ? 0.3 * pow(coss, phong_exponent) : 0.0;

		color_red += spec_coef;
		color_green += spec_coef ;
		color_blue += spec_coef;

		col_red = (int) (250 * color_red);
		col_green = (int) (250 * color_green);
		col_blue = (int) (250 * color_blue);

#if 0

#define MAX_COL_Z (1 << 20)
#define NUM_STRIPES 32
#define STRIPE_PART 8

		col_z = (int)((z/15 + 1.0) * MAX_COL_Z);
		col_z %= MAX_COL_Z/NUM_STRIPES;
		if (0 <= col_z  && col_z <= MAX_COL_Z/(NUM_STRIPES * STRIPE_PART)) {
			col_blue = ::max(col_red, col_green);
			col_red = 0;
			col_green = 0;
		};

#endif

		mref(colmat_r_d, xx, yy) = col_red;
		mref(colmat_g_d, xx, yy) = col_green;
		mref(colmat_b_d, xx, yy) = col_blue;
	} else {
		// background color
		mref(colmat_r_d, xx, yy) = 250; //64;
		mref(colmat_g_d, xx, yy) = 250; //32;
		mref(colmat_b_d, xx, yy) = 250; //64;

	}

	__syncthreads();

}

#define THREAD_NUMXY 16

void gpu_compute_colmat(double a, double b, int xmax, int ymax, const cudaPoly3 & f3_h,
		double euler_phi, double euler_theta, double euler_psi,
		int *colmat_r, int *colmat_g, int *colmat_b) {

	int *colmat_r_d;
	int *colmat_g_d;
	int *colmat_b_d;
	
	int akt_deg_global_host = -1;
	
	int gl_max_expo_host[3] = {-1, -1, -1};

	for(int i = 0; i < 40; ++i) {
		int degi = f3_h.xexpo[i] + f3_h.yexpo[i] + f3_h.zexpo[i];
		if (degi > akt_deg_global_host) {
			akt_deg_global_host = degi;
		}
		if (f3_h.xexpo[i] > gl_max_expo_host[0]) {
			gl_max_expo_host[0] = f3_h.xexpo[i];
		}
		if (f3_h.yexpo[i] > gl_max_expo_host[1]) {
			gl_max_expo_host[1] = f3_h.yexpo[i];
		}
		if (f3_h.zexpo[i] > gl_max_expo_host[2]) {
			gl_max_expo_host[2] = f3_h.zexpo[i];
		}
	}

	const int N = 4 * gl_win_size * gl_win_size;

	cout << "N = " << N << endl;

	cout << "xmax = " << xmax << " ymax = " << ymax << endl;

	//print_cudaPoly3(f3_h);

	printf("sizeof(cudaPoly3) = %ld\n ", sizeof(cudaPoly3));

	x_prepare_binom<<<1,1>>>();

	cudaDeviceSynchronize();

	//x_print_binom<<<1,1>>>();

	CCE(cudaDeviceSetLimit(cudaLimitStackSize , 128 * 1024));

	CCE( cudaMemcpyToSymbol(akt_deg_global, &akt_deg_global_host, sizeof(int),
					0, cudaMemcpyHostToDevice ) );

	CCE( cudaMemcpyToSymbol(gl_max_expo, &gl_max_expo_host[0], 3 * sizeof(int) ) );



	CCE( cudaMemcpyToSymbol(f3, &f3_h, sizeof(cudaPoly3), 0, cudaMemcpyHostToDevice) );

	cudaDeviceSynchronize();


    cuda_rotate_mat<<<1,1>>>(euler_phi, euler_theta, euler_psi);

	cuda_init_f3_diff<<<1,1>>>();

	cudaDeviceSynchronize();


	CCE( cudaMalloc((void**) &colmat_r_d, N * sizeof(int) ) );
	CCE( cudaMalloc((void**) &colmat_g_d, N * sizeof(int) ) );
	CCE( cudaMalloc((void**) &colmat_b_d, N * sizeof(int) ) );

	dim3 grids(gl_win_size/THREAD_NUMXY, gl_win_size/THREAD_NUMXY);
	dim3 threads(THREAD_NUMXY, THREAD_NUMXY);

	compute_colmat<<<grids, threads>>>(a, b, xmax, ymax, euler_phi, euler_theta, euler_psi,
			colmat_r_d, colmat_g_d, colmat_b_d);

	cout << "gpu computation done." << endl;


	CCE(cudaMemcpy(colmat_r, colmat_r_d, N * sizeof(int), cudaMemcpyDeviceToHost));
	CCE(cudaMemcpy(colmat_g, colmat_g_d, N * sizeof(int), cudaMemcpyDeviceToHost));
	CCE(cudaMemcpy(colmat_b, colmat_b_d, N * sizeof(int), cudaMemcpyDeviceToHost));

	CCE(cudaFree(colmat_r_d));
	CCE(cudaFree(colmat_g_d));
	CCE(cudaFree(colmat_b_d));


}

