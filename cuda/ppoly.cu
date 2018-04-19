


#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include <helper_cuda.h>

#define CCE checkCudaErrors




#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>



#include "pimage_local.h"


//#define DEBUG(x) cerr << (x) << endl
#define DEBUG(x) {}
//#define KEYWAIT KEYWAIT
#define KEYWAIT {}
//#define PRINT_COEFS(x,y) print_coefs((x),(y))
#define PRINT_COEFS(x,y) {};

#define POLY_MAX_DEG 10

// poly an x^n + ... + a0
// list: an an-1 .... a0 a0 on highest address

using namespace std;

__device__ int x_root_interv_list(int deg, const double* poly_coefs, double b_bound, double sign, double *alist, double *blist, int & root_list_len);

__device__ int x_core_alg_01(int deg, const double* poly_coefs, double* palist, double *pblist,
							double* & palist_next, double* & pblist_next);


__device__ double binom[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1];

__device__ double coefs_save[POLY_MAX_DEG + 1];

__device__ void x_print_coefs(int deg, double *poly_coefs)
{
	for(int i = 0; i <= deg; ++i) {
		printf("%f x^%d + ", poly_coefs[i], (deg- i));
	}
	printf("\n");
}


__global__ void x_prepare_binom()
{
	for(int i = 0; i <= POLY_MAX_DEG; ++i) {
		binom[0][i] = 0;
	}
	binom[0][0] = 1;
	for(int j = 1; j <= POLY_MAX_DEG; ++j) {
		binom[j][0] = 1;
		for(int i = 1; i < j + 1; ++i) {
			binom[j][i] = binom[j-1][i-1] + binom[j-1][i];
		}
	}
}

__global__ void x_print_binom()
{
	for(int i = 0; i <= POLY_MAX_DEG; ++i) {
		for(int j = 0; j < i + 1; ++j) {
			printf(" %f ", binom[i][j]);
		}
		printf("\n");
	}
}

__device__ int x_make_trafo_1(int deg, double* coefs_in)
{
	DEBUG("start_trafo1");
	PRINT_COEFS(deg, coefs_in);

	int i = 0;
	int j = deg;
	while (i < j) {
		double m = coefs_in[i];
		coefs_in[i] = coefs_in[j];
		coefs_in[j] = m;
		++i;
		--j;
	}
	i = deg - 1;
	while (i >= 0) {
		j = i + 1;
		int cnt = 1;
		while (j <= deg) {
			coefs_in[j] += binom[deg-i][cnt] * coefs_in[i];
			++j;
			++cnt;
		}
		--i;
	}

	PRINT_COEFS(deg, coefs_in);
	DEBUG("end trafo1. press key.");
	KEYWAIT;
	return 0;
}

__device__ int x_make_trafo_2(int deg, double* coefs_in)
{
	double mult = 1;
	int i = 0;

	while (i <= deg) {
		coefs_in[i] *= mult;
		mult *= 2;
		++i;
	}
	return 0;
}

__device__ int x_make_trafo_3(int deg, double* coefs_in)
{
	int i = deg - 1;
	int j;

	DEBUG("start trafo 3");
	PRINT_COEFS(deg, coefs_in);
	while (i >= 0) {
		j = i + 1;
		int cnt = 1;
		while (j <= deg) {
			coefs_in[j] += binom[deg-i][cnt] * coefs_in[i];
			++j;
			++cnt;
		}
		--i;
	}
	PRINT_COEFS(deg, coefs_in);
	DEBUG("end trafo 3. press key");
	KEYWAIT;

	return 0;

}

__device__ int x_make_trafo_divide(int deg, double * coefs_in)
{
	double coefs_aux[POLY_MAX_DEG + 1];
	const int copy_size = (POLY_MAX_DEG+1) * sizeof(double);
	memset(coefs_aux, 0, copy_size);

	for(int i = 0; i < deg; ++i) {
		double mult = 0.5;
		for(int j = i; j <= deg-1; ++j) {
			coefs_aux[j] += mult * coefs_in[i];
			mult *= 0.5;
		}
	}

	memcpy(coefs_in, coefs_aux, copy_size);

	return 0;

}

__device__ double x_eval(int deg, const double *coefs, double x)
{
	double sum = 0;
	for(int i = 0; i <= deg; ++i) {
		sum = sum * x + coefs[i];
	}
	return sum;
}

// compute p(x/b)
__device__ int x_poly_scale(int deg, double* coefs_in, double b)
{
	double mult = 1;
	for(int i = deg; i >= 0; --i) {
		coefs_in[i] *= mult;
		mult *= b;
	}
	return 0;
}

__device__ int x_sign_variations(int deg, double* coefs)
{
	double start = coefs[0];
	int i = 1;
	int n_sign_vars = 0;

	DEBUG("enter sign_variations");
	PRINT_COEFS(deg, coefs);

	while (i <= deg) {
		double next = coefs[i];
		if (next != 0) {
			if (start * next < 0) {
				++n_sign_vars;
				start = next;
			}
		}
		i++;
	}

	// cerr << "sign vars = " << n_sign_vars << endl;
	DEBUG("leave sign_variations. press key.");
	KEYWAIT;
	return n_sign_vars;
}

__device__ double x_comp_disc(int deg, double* a)
{
	double disc = -999;
	if (deg > 3) {
		return disc;
	}
	if (deg == 2) {
		disc = -4 * a[2] * a[0] + a[1]*a[1];
	} else if (deg == 3) {
		disc = -27 *pow(a[0],2) * pow(a[3],2)  + 18 *  a[0] * a[3] * a[2] * a[1] + pow(a[1],2)  * pow(a[2],2)
				- 4 * pow(a[1],3) *  a[3] - 4 * pow(a[2],3) * a[0];
	}
	return disc;

}

__device__ int x_bisect(int deg, const double* poly_coefs, double x0, double x1, double & x_root)
{
	double y0, y1;
	double ymid;
	double xmid = 0.5*(x1 + x0);
	const double eps = 1e-6;

	//cout << "bisect:" << x0 << " " << x1 << endl;

	y0 = x_eval(deg, poly_coefs, x0);
	y1 = x_eval(deg, poly_coefs, x1);

	ymid = y0;

	if (!(y0 * y1 <= 0)) {
		//cout << "x0 = " << x0 << " x1 = " << x1 << " y0 = " << y0 << " y1 = " << y1 << endl;
		//print_coefs(deg, coefs_save);

		//cout << "disc = " << comp_disc(deg, coefs_save) << endl;

		x_root = 0;

		printf(".");

		return 0;

		//cout << "Press key." << endl;
		//getchar();
	}

	while (fabs(x1 - x0) > eps) {
		//cout << "x0 = " << x0 << " x1 = " << x1 << endl;
		xmid = 0.5*(x0 + x1);
		ymid = x_eval(deg, poly_coefs, xmid);
		if (fabs(ymid) < eps) {

			x_root = xmid;
			goto ende;
		}
		if (y0 * ymid >= 0) {
			x0 = xmid;
			y0 = ymid;
		} else if (y1 * ymid >= 0) {
			x1 = xmid;
			y1 = ymid;
		}
	}
	x_root = xmid;
ende:

	//cout << "x_root = " << x_root << endl;
	return 0;
}

__device__ int x_find_root_list(int deg, const double* poly_1, double *alist, double* blist, int n_len, double* x0_list,
		double * & x0_list_next)
{
	int j = 0;
	double x0 = -DBL_MAX;
	while (j < n_len) {
		x_bisect(deg, poly_1, alist[j], blist[j], x0);
		*x0_list = x0;
		++x0_list;
		++j;
	}
	x0_list_next = x0_list;
	return 0;
}

#define SWAP_DOUBLE(a,b) do { double m = (a); (a) = (b); (b) = m; } while(0)

// sorts into ascending order
__device__ int x_sort_list(int len, double* list)
{
	for(int i = 0; i < len - 1; ++i) {
		double min = list[i];
		int jmin = -1;
		for(int j = i+1; j < len; ++j ) {
			if (list[j] < min) {
				jmin = j;
				min = list[j];
			}
		}
		if (jmin > 0) {
			SWAP_DOUBLE(list[i], list[jmin]);
		}
	}
	return 0;
}

__device__ int x_root_final_list(int deg, const double* poly_coefs, double b_bound, double * x0_list, int & x0_len)
{
	double alist[POLY_MAX_DEG + 1];
	double blist[POLY_MAX_DEG + 1];
	double poly_1[POLY_MAX_DEG + 1];

	double *x0_list_next;
	double *x0_list_0 = x0_list;

	memcpy(coefs_save, poly_coefs, sizeof(double) * (POLY_MAX_DEG + 1));

	x0_len = 0;
	x0_list_next = x0_list;

	int deg0 = deg;

	while (poly_coefs[deg0] == 0) {
		deg0 = deg0 - 1;
		*x0_list_next = 0;
		++x0_list_next;
		++x0_len;
	}

	x0_list = x0_list_next;

	memcpy(poly_1, poly_coefs, (deg + 1)* sizeof(double));

	int ablist_len = -1;

	x_root_interv_list(deg, poly_1, b_bound, +1, alist, blist, ablist_len);

	//cout << "ablist_len = " << ablist_len << endl;

	if (ablist_len > 0) {
		int ret = x_find_root_list(deg, poly_coefs, alist, blist, ablist_len, x0_list, x0_list_next);

		x0_len += x0_list_next - x0_list;

		x0_list = x0_list_next;
	};

	//cout << "negative case" << endl;


	x_poly_scale(deg, poly_1, -1);
	x_root_interv_list(deg, poly_1, 200, -1, alist, blist, ablist_len);

	if (ablist_len > 0) {
		//cout << "negative ablist_len = " << ablist_len << endl;
		int ret = x_find_root_list(deg, poly_coefs, alist, blist, ablist_len, x0_list, x0_list_next);

		x0_len += x0_list_next - x0_list;

	}

	// sorts into ascending order
	x_sort_list(x0_len, x0_list_0);

	return 0;
}

__device__ int x_root_interv_list(int deg, const double* poly_coefs, double b_bound, double sign,
		double *alist, double *blist, int & root_list_len)
{
	double * palist_next;
	double * pblist_next;
	double poly_1[POLY_MAX_DEG + 1];

	memcpy(poly_1, poly_coefs, (deg + 1) * sizeof(double));

	x_poly_scale(deg, poly_1, b_bound);

	DEBUG( "poly scaled");

	KEYWAIT;

	x_core_alg_01(deg, poly_1, alist, blist, palist_next, pblist_next);

	DEBUG( "core alg done " );

	root_list_len = palist_next - alist;
	for(int i = 0; i < root_list_len; ++i) {
		alist[i] *= b_bound * sign;
		blist[i] *= b_bound * sign;
		if (sign < 0) {
			double m = alist[i];
			alist[i] = blist[i];
			blist[i] = m;
		}
	}

	return 0;
}



#define NO_ROOTS_FOUND -1
#define ROOT_01_FOUND 1


__device__ int x_core_alg_01(int deg, const double* poly_coefs, double* palist, double *pblist,
		double* & palist_next, double* & pblist_next)
{
	double poly_1[POLY_MAX_DEG + 1];
	double poly_2[POLY_MAX_DEG + 1];

	const int copy_size = (deg+1) * sizeof(double);
	memcpy(poly_1, poly_coefs, copy_size);

	x_make_trafo_1(deg, poly_1);

	DEBUG("trafo 1 done");

	int vars = x_sign_variations(deg, poly_1);
    if (vars == 0) {
    	palist_next = palist;
    	pblist_next = pblist;
    	DEBUG("no roots found");
    	return NO_ROOTS_FOUND;
    }
    if (vars == 1) {
    	*palist = 0;
    	palist_next = palist + 1;
    	*pblist = 1;
    	pblist_next = pblist + 1;
    	DEBUG("one root found");
    	return ROOT_01_FOUND;
    }

    double *palist_next1;
    double *pblist_next1;

    double midval = x_eval(deg, poly_coefs, 0.5);

    //cerr << "midval = " << midval << endl;

    if (midval == 0) {

    	//cerr << "doing midval = 0" << endl;

    	*palist = 0.5;
    	*pblist = 0.5;
    	palist_next1 = ++palist;
    	pblist_next1 = ++pblist;

    	//memcpy(poly_1, poly_coefs, copy_size);
    	//make_trafo_divide(deg, poly_1);

        palist = palist_next1;
        pblist = pblist_next1;

        // return core_alg_01(deg-1, poly_1, palist, pblist, palist_next, pblist_next);
    }

    memcpy(poly_1, poly_coefs, copy_size);

    x_make_trafo_2(deg, poly_1);

    memcpy(poly_2, poly_1, copy_size);

    DEBUG("trafo2 done");

    x_core_alg_01(deg, poly_1, palist, pblist, palist_next1, pblist_next1);
    for (double * pa = palist; pa < palist_next1; ++pa) {
    	*pa *= 0.5;
    }
    for (double * pb = pblist; pb < pblist_next1; ++pb) {
    	*pb *= 0.5;
    }
    palist = palist_next1;
    pblist = pblist_next1;

    // memcpy(poly_1, poly_coefs, copy_size);

    x_make_trafo_3(deg, poly_2);

    DEBUG("trafo3 done");

    x_core_alg_01(deg, poly_2, palist, pblist, palist_next1, pblist_next1);

    for (double * pa = palist; pa < palist_next1; ++pa) {
    	*pa += 1;
    	*pa *= 0.5;
    }
    for (double * pb = pblist; pb < pblist_next1; ++pb) {
    	*pb += 1;
    	*pb *= 0.5;
    }
    palist_next = palist_next1;
    pblist_next = pblist_next1;

    return 0;

}


__device__ void x_print_intervals(int len, double *alist, double *blist) {
	for(int i = 0; i < len; ++i){
		printf( "[ %f , %f ]\n", alist[i], blist[i]);
	}
}

__device__ void x_list_out(int len, double* list)
{
	printf( "( " );
	for(int i = 0; i < len; ++i) {
		printf("%f , ", list[i]);
	}
	printf( ")\n" );
}

/*

void test_root_list()
{
	double poly_coefs[POLY_MAX_DEG + 1];
	double alist[POLY_MAX_DEG + 1];
	double blist[POLY_MAX_DEG + 1];

	double x0_list[POLY_MAX_DEG + 1];

	int n_len;
	int deg;
	int x0_list_len;

	// (x-1)*(x-2)

	deg = 2;

	poly_coefs[0] = 1;
	poly_coefs[1] = -3;
	poly_coefs[2] = 2;

	int ret = root_interv_list(deg, poly_coefs, 10, 1, alist, blist, n_len);
	print_intervals(n_len, alist, blist);

	root_final_list(deg, poly_coefs, 10, x0_list, x0_list_len);
	list_out(x0_list_len, x0_list);
	cout << endl;

	// early an error from asxp.cpp, now works correct

	deg = 2;

	poly_coefs[0] = 312.5;
	poly_coefs[1] = 44921.9;
	poly_coefs[2] = -7791.31;

	ret = root_interv_list(deg, poly_coefs, 10, 1, alist, blist, n_len);
	print_intervals(n_len, alist, blist);

	root_final_list(deg, poly_coefs, 10, x0_list, x0_list_len);
	list_out(x0_list_len, x0_list);
	cout << endl;

	// currently active error from asxp.cpp

	deg = 2;

	// f = 21 * (x-2) * (x-5)
	// produced an interesting error:
	// intervals given back are [5,5] [0,10]
	// problem is: second interval is for f/(x-5) = 21 * (x-2)
	// for this [0, 10] is a perfectly sound interval
	// unfortunately it is entered in a common list of intervals
	// with no notification, that the reference polynomial has been
	// changed by the algorithm
	// solution: very simple: just add the midpoint to the original
	// interval list and do *not* divide out the factor, just
	// consider the intervals [0,1/2] and [1/2, 0] as if
	// at 1/2 there was no special case of the polynomial vanishing.

	poly_coefs[0] = 21;
	poly_coefs[1] = -147;
	poly_coefs[2] = 210;

	ret = root_interv_list(deg, poly_coefs, 10, 1, alist, blist, n_len);
	print_intervals(n_len, alist, blist);

	root_final_list(deg, poly_coefs, 10, x0_list, x0_list_len);
	list_out(x0_list_len, x0_list);
	cout << endl;
	cout << "(x) press key." << endl;
	getchar();


	// (x-1)*(x-3)*(x-5)
	deg = 3;

	poly_coefs[0] = 1;
	poly_coefs[1] = -9;
	poly_coefs[2] = 23;
	poly_coefs[3] = -15;

	ret = root_interv_list(deg, poly_coefs, 10, 1, alist, blist, n_len);
	print_intervals(n_len, alist, blist);

	root_final_list(deg, poly_coefs, 10, x0_list, x0_list_len);
	list_out(x0_list_len, x0_list);
	cout << endl;

	// (x+1)*(x+2)*(x+5)
	deg = 3;

	poly_coefs[0] = 1;
	poly_coefs[1] = 8;
	poly_coefs[2] = 17;
	poly_coefs[3] = 10;

	poly_scale(deg, poly_coefs, -1);
	ret = root_interv_list(deg, poly_coefs, 10, -1, alist, blist, n_len);
	print_intervals(n_len, alist, blist);

	poly_scale(deg, poly_coefs, -1);
	root_final_list(deg, poly_coefs, 10, x0_list, x0_list_len);
	list_out(x0_list_len, x0_list);
	cout << endl;

	// 2*x-3
	deg = 1;

	poly_coefs[0] = -3;
	poly_coefs[1] = 2;

	ret = root_interv_list(deg, poly_coefs, 10, -1, alist, blist, n_len);
	print_intervals(n_len, alist, blist);

	root_final_list(deg, poly_coefs, 10, x0_list, x0_list_len);
	list_out(x0_list_len, x0_list);
	cout << endl;

	// 2*x+3
	deg = 1;

	poly_coefs[0] = 3;
	poly_coefs[1] = 2;

	poly_scale(deg, poly_coefs, -1);
	ret = root_interv_list(deg, poly_coefs, 10, -1, alist, blist, n_len);
	print_intervals(n_len, alist, blist);

	poly_scale(deg, poly_coefs, -1);
	root_final_list(deg, poly_coefs, 10, x0_list, x0_list_len);
	list_out(x0_list_len, x0_list);
	cout << endl;

}
*/
