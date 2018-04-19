
#include <iostream>

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>


#include "eigen.h"

using namespace std;

#define MAT(m,i,j) (m[3 * (i) + (j)])

double m_norm(double* vec, int n) {
	double s = 0;
	for(int i = 0; i < n; ++i ) {
		s += vec[i] * vec[i];
	}
	return sqrt(s);
}


void m_skalmult(double cf, double * mat, double * mat_res, int n) {
	for(int i = 0; i < n; ++i) {
		mat_res[i] = cf * mat[i];
	}
}

void m_skalmult_mat(double cf, double * mat, double * mat_res) {
	m_skalmult(cf, mat, mat_res, 9);
}


void m_set(double* lis, double val, int len) {
	for(int i = 0; i < len; ++i) {
		lis[i] = val;
	}
}

double m_skalprod(double* v1, double* v2, int n) {

	double s = 0;

	for(int i = 0; i < n; ++i) {
		s += v1[i] * v2[i];
	}

	return s;

}


void m_mult_vec_mat(double* vec, double *mat, double * vec_res) {
	vec_res[0] = vec[0] * mat[0] + vec[1] * mat[3] + vec[2] * mat[6];
	vec_res[1] = vec[0] * mat[1] + vec[1] * mat[4] + vec[2] * mat[7];
	vec_res[2] = vec[0] * mat[2] + vec[1] * mat[5] + vec[2] * mat[8];
}

void m_mult_mat_vec(double* vec, double *mat, double * vec_res) {
	vec_res[0] = vec[0] * mat[0] + vec[1] * mat[1] + vec[2] * mat[2];
	vec_res[1] = vec[0] * mat[3] + vec[1] * mat[4] + vec[2] * mat[5];
	vec_res[2] = vec[0] * mat[6] + vec[1] * mat[7] + vec[2] * mat[8];
}

void m_transpose_mat(double * mat, double* mat_res) {

	for(int i = 0; i < 3; ++i) {
		for(int j = 0; j < 3; ++j) {
			mat_res[j * 3 + i] = mat[i * 3 + j];
		}
	}

}

void m_mult(double * m1, double * m2, double * m3)
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

void m_print(double* m)
{
	for(int i = 0; i < 3; ++i) {
		for(int j = 0; j < 3; ++j) {
			cout << m[i*3+j] << " , ";
		}
		cout << endl;
	}
}


void other_indices(int i, int & i1, int & i2) {

	if (i == 0) {
		i1 = 1;
		i2 = 2;
	} else if (i == 1) {
		i1 = 0;
		i2 = 2;
	} else if (i == 2) {
		i1 = 0;
		i2 = 1;
	}
}

double cofactor_m3(double* mat, int i, int j) {

	int i1, i2, j1, j2;
	i1 = i2 = j1 = j2 = -1;

	other_indices(i, i1, i2);
	other_indices(j, j1, j2);
	double cf = MAT(mat,i1,j1) * MAT(mat, i2, j2) - MAT(mat, i1, j2) * MAT(mat, i2, j1);
	return cf;
}


void adjoint_m3(double* mat, double* adjmat) {

	double sign = 1;
	for(int i = 0; i < 3; ++i) {
		for(int j = 0; j < 3; ++ j) {
			double cofactor_ij = cofactor_m3(mat, i, j);
			adjmat[3 * j + i] = sign * cofactor_ij;
			sign = -sign;
		}
	}

}

double detmat_m3(double* mat) {
	double adjmat[9];
	double prodmat[9];
	adjoint_m3(mat, adjmat);
	m_mult(mat, adjmat, prodmat);

	return prodmat[0];
}

double invmat_m3(double* mat, double* invmat) {

	double adjmat[9];
	adjoint_m3(mat, adjmat);

	double det = detmat_m3(mat);

	if (fabs(det) > 1e-12) {
		m_skalmult_mat(1/det, adjmat, invmat);
	}
	return det;
}

void make_matP(double* nablaf, double norm_nabla_f, double* matP) {

	for(int i = 0; i < 3; ++i) {
		for(int j = 0; j < 3; ++j ) {
			matP[3 * i + j] = ((i == j) * 1) - nablaf[i] * nablaf[j]/(norm_nabla_f * norm_nabla_f);
		}
	}
}

void m_gram_schmidt(double* v1, double* v2) {

	double norm1 = m_norm(v1, 3);
	double norm2 = m_norm(v2, 3);
	for(int i = 0; i < 3; ++i) {
		v1[i] /= norm1;
		v2[i] /= norm2;
	}

	double prod12 = m_skalprod(v1, v2, 3);

	for(int i = 0; i < 3; ++i) {
		v2[i] = v2[i] - prod12 * v1[i];
	}
	norm2 = m_norm(v2, 3);
	for(int i = 0; i < 3; ++i) {
		v2[i] /= norm2;
	}
}



int restrict_hessian_to_tangentspace(double* vec_normal, double* hess33, double* shape_mat22,
		double* base_v1, double* base_v2) {

	double nx = vec_normal[0];
	double ny = vec_normal[1];
	double nz = vec_normal[2];

	if (fabs(nz) < 1e-8) {
		return 2;
	}

	double Ymat[9];
	double Ymat_inv[9];

	for(int i=0; i < 3; ++i) {
		Ymat[i * 3 + 2] = vec_normal[i];
	}

	// two vectors v1 v2 perpendicular to (nx ny nz)
	// Ymat is Y = (v1 v2 n)

	Ymat[0] = 1; Ymat[3] = 0; Ymat[6] = -nx/nz;
	Ymat[1] = 0; Ymat[4] = 1; Ymat[7] = -ny/nz;

	for(int i = 0; i < 3; ++i) {
		base_v1[i] = Ymat[3 * i];
		base_v2[i] = Ymat[3 * i + 1];
	}

	m_gram_schmidt(base_v1, base_v2);

	for(int i = 0; i < 3; ++i) {
		Ymat[3 * i] = base_v1[i];
		Ymat[3 * i + 1] = base_v2[i];
	}


	// cout << "Ymat = " << endl;
	// m_print(Ymat);

	invmat_m3(Ymat, Ymat_inv);

	// cout << "Ymat_inv = " << endl;
	// m_print(Ymat_inv);

	double YinvH[9];

	m_mult(Ymat_inv, hess33, YinvH);

	double auxmat[9];

	m_mult(YinvH, Ymat, auxmat);

	//cout << "auxmat = " << endl;
	//m_print(auxmat);

	// auxmat is of the form
	//
	// lambda_1 mu_1  sigma_1
	// lambda_2 mu_2  sigma_2
	// 0        0     0
	//

	double a = MAT(auxmat,0,0);
	double b = MAT(auxmat,0,1);
	double c = MAT(auxmat,1,0);
	double d = MAT(auxmat,1,1);

	//cout << "a b c d = " << a << " , " << b << " , " << c << " , " << d << endl;

	shape_mat22[0] = a;
	shape_mat22[1] = b;
	shape_mat22[2] = c;
	shape_mat22[3] = d;

	//cin.get();

	return 0;

}


int quadratic_solve(double a, double b, double c, double & x1, double & x2 ) {

	if (a == 0) {
		if (b == 0) {
			x1 = x2 = DBL_MAX;
			return 3;
		}
		x1 = x2 = -c/b;
		return 1;
	}

	double p = b/a;
	double q = c/a;

	double disc = p*p/4 - q;

	if (fabs(disc) < 1e-11) {
		disc = 0;
	}

	if (disc < 0) {
		x1 = x2 = DBL_MAX;
		return 3;
	}

	double root = sqrt(disc);
	x1 = -p/2 + root;
	x2 = -p/2 - root;

	return 0;
}

void eigenvalues_22(double* mat22, double & x1, double & x2 ) {

	double a = mat22[0];
	double b = mat22[1];
	double c = mat22[2];
	double d = mat22[3];

	double p = -(a + d);
	double q = a * d - b * c;

	int res_code = quadratic_solve(1, p, q, x1, x2);

	if (x1 > x2) {
		::swap(x1, x2);
	}

	if (!(x1 <= x2)) {
		std::cout << "p = " << p << " q = " << q << " disc = " << p * p /4.0 - q << std::endl;
		std::cout << "err: x1 = " << x1 << " x2 = " << x2 << std::endl;
		assert(x1 <= x2);
	}
#if 0
	if (!(x1 <= x2)) {
		cout << "rescode = " << res_code << endl;
		cout << "x1 = " << x1 << "x2 = " << x2 << endl;
		assert(0);
	}
#endif

}

// x1 and x2 are the eigenvalues, x1 <= x2
// w1 and w2 are the eigenvectors to x1, x2
void eigenvectors_22(double* mat22, double x1, double x2, double* w1, double * w2) {

	double mat1[4];
	double mat2[4];

	for(int i = 0; i < 2; ++i) {
		for(int j = 0; j < 2; ++j) {
			mat1[2 * i + j] = x1 * (i == j) - mat22[2 * i + j];
			mat2[2 * i + j] = x2 * (i == j) - mat22[2 * i + j];
		}
	}

	double norm_max = 0;
	int j_ind_max = -1;
	double new_norm;

	for(int j = 0; j < 2; ++j) {
		for(int i = 0; i < 2; ++i) {
			w1[i] = mat2[2 * i + j];
		}
		new_norm = m_norm(w1, 2);
		if (new_norm > norm_max) {
			norm_max = new_norm;
			j_ind_max = j;
		}
	}

	if (j_ind_max >= 0) {
		for(int i = 0; i < 2; ++i) {
			w1[i] = mat2[2 * i + j_ind_max];
		}
	}

	norm_max = 0;
	j_ind_max = -1;

	for(int j = 0; j < 2; ++j) {
		for(int i = 0; i < 2; ++i) {
			w2[i] = mat1[2 * i + j];
		}
		new_norm = m_norm(w2, 2);
		if (new_norm > norm_max) {
			norm_max = new_norm;
			j_ind_max = j;
		}
	}

	if (j_ind_max >= 0) {
		for(int i = 0; i < 2; ++i) {
			w2[i] = mat1[2 * i + j_ind_max];
		}
	}


}


int m_calculate_diagonal_form(double* shape_mat, double & x1, double & x2, double* base_v1, double* base_v2) {

	double w1[2];
	double w2[2];

	double new_v1[3];
	double new_v2[3];

	eigenvalues_22(shape_mat, x1, x2);

	if (fabs(x1 - x2) > 0) {
		eigenvectors_22(shape_mat, x1, x2, w1, w2);
		for(int i = 0; i < 3; ++i) {
			new_v1[i] = w1[0] * base_v1[i] + w1[1] * base_v2[i];
			new_v2[i] = w2[0] * base_v1[i] + w2[1] * base_v2[i];
		}
	} else {
		for(int i = 0; i < 3; ++i) {
			new_v1[i] = base_v1[i];
			new_v2[i] = base_v2[i];
		}
	}

	memcpy(base_v1, new_v1, sizeof(new_v1));
	memcpy(base_v2, new_v2, sizeof(new_v2));

	return 0;

}



int eigen_test() {

	double testmat[9] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	double testmat1[9] = { 2, -1, 0, -1, 2, -1, 0, -1, 2 };

	m_print(testmat);

	double adjmat[9];

	double prodmat[9];

	adjoint_m3(testmat, adjmat);

	m_print(adjmat);

	m_mult(testmat, adjmat, prodmat);

	m_print(prodmat);

	double det = invmat_m3(testmat1, adjmat);

	cout << "det = " << det << endl;

	m_print(adjmat);

	cout << "Hessian restrict test" << endl;

	double vnormal[3] = {1, 1, 1};
	double hess33[9] = {1, 0, 0, 0, 2, 0, 0, 0, 3};
	double resmat22[4];

	double base_v1[3];
	double base_v2[3];

	restrict_hessian_to_tangentspace(vnormal, hess33, resmat22, base_v1, base_v2);

	cout << resmat22[0] << " " << resmat22[1] << endl;
	cout << resmat22[2] << " " << resmat22[2] << endl;

	return 0;

}



