
#ifndef __roots_h
#define __roots_h

int prepare_binom();

int print_binom();

int root_final_list(int deg, const double* poly_coefs, double b_bound, double * x0_list, int & x0_len);

double comp_disc(int deg, double* a);


void test_root_list();





#endif // __sturm_h
