#ifndef __eigen_h
#define __eigen_h


double m_norm(double* vec, int n);

void m_skalmult(double cf, double * mat, double * mat_res, int n);

void m_skalmult_mat(double cf, double * mat, double * mat_res);


void m_set(double* lis, double val, int len);

double m_skalprod(double* v1, double* v2, int n);

void m_mult_vec_mat(double* vec, double *mat, double * vec_res);

void m_mult_mat_vec(double* vec, double *mat, double * vec_res);

void m_transpose_mat(double * mat, double* mat_res);

void m_mult(double * m1, double * m2, double * m3);

void m_print(double* m);



double cofactor_m3(double* mat, int i, int j);
void adjoint_m3(double* mat, double* adjmat);
double detmat_m3(double* mat);
double invmat_m3(double* mat, double* invmat);

int quadratic_solve(double a, double b, double c, double & x1, double & x2 );


void make_matP(double* nablaf, double norm_nabla_f, double* matP);

int restrict_hessian_to_tangentspace(double* vec_normal, double* hess33, double* shape_mat22,
		double* base_v1, double* base_v2);

int m_calculate_diagonal_form(double* shape_mat, double & x1, double & x2, double* base_v1, double* base_v2);


int eigen_test();


















#endif


