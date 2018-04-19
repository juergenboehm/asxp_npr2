
#ifndef __pimage_local_h
#define __pimage_local_h


__global__ void x_prepare_binom();

__global__ void x_print_binom();

__device__ int x_root_final_list(int deg, const double* poly_coefs, double b_bound, double * x0_list, int & x0_len);

__device__ double x_comp_disc(int deg, double* a);





#endif

