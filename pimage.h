
#ifndef __pimage_h
#define __pimage_h



class cudaPoly3 {

public:

	double coefs[40];

	int xexpo[40];
	int yexpo[40];
	int zexpo[40];

	int len;

#ifndef NOCUDA
	__device__ __host__ int degree();
	__device__ __host__ cudaPoly3 & diff(int iv);
	__device__ __host__ double eval(double* sl);
#else
	int degree();
	cudaPoly3 & diff(int iv);
	double eval(double* sl);
#endif



};


void gpu_compute_colmat(double a, double b, int xmax, int ymax, const cudaPoly3 & f3,
		double euler_phi, double euler_theta, double euler_psi,
		int *colmat_r, int *colmat_g, int *colmat_b);







#endif
