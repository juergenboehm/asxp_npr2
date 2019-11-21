

#include <math.h>
#include <float.h>


#include "asxp_arrays.h"

#include "eigen.h"
#include "roots.h"

#include "poly.h"
#include "parser.h"



#include "screenwidget.h"


#include "procpoly.h"






#define UNUSED(expr) do { (void)(expr); } while (0)


#define EPS 1e-12
#define M_INF -1e38

#define POLY_DEG 3
#define POLY_DEG1 (POLY_DEG+1)



//Poly5 f5;
//Poly3 f3;

// stand alone
void m_rot_z(double a, double * mat)
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

// stand alone
void m_rot_x(double a, double * mat)
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



// private functions of ProcPoly

void ProcPoly::normalize_poly_coefs(double* poly_coefs, int akt_deg, int & deg_new)
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



void ProcPoly::lb_poly_mult(int deg, double * coefs, double a)
{
	// multiply coefs poly coefs[0] * t^deg + ... with (t + a)
	coefs[deg+1] = 0;
	for(int i = deg + 1; i >=  1; --i) {
		coefs[i] += a * coefs[i-1];
	}
}

void ProcPoly::lb_poly_scal_mult(int deg, double * coefs, double a)
{
	for(int i = 0; i <= deg; ++i) {
		coefs[i] *= a;
	}
}


void ProcPoly::lb_gen_lagrange_basis(int deg, double* xbase, double lagr_basis[][max_deg + 1])
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
			lb_poly_mult(curr_deg, aux_poly, -xbase[j]);
			++curr_deg;
		}
		lb_poly_scal_mult(deg, aux_poly, 1/pre_coef);

		memcpy(lagr_basis[i], aux_poly, sizeof(aux_poly));

	}
}


double ProcPoly::eval_f3_mat(Poly3 & f3, double *sl)
{
	double sl1[3];
	m_mult_mat_vec(sl, m_euler, sl1);

	double f = f3.eval(sl1);

	return f;

}


int ProcPoly::eval_poly_poly(Poly3 & f3,
		double x, double y, double z, double & f, double* fnormal, double* fhessian, bool only_normal)
{

	// hier muss
	// 1) sl mit m_euler substituiert werden vor dem Einsetzen
	// 2) (fx fy fz) * m_euler multipliziert werden

	double sl[3];
	double sl1[3];
	double nablaf_eul[3];

	double fhess[9];

	sl[0] = x;
	sl[1] = y;
	sl[2] = z;

	m_mult_mat_vec(sl, m_euler, sl1);

	f = f3.eval(sl1);

	sl[0] = f3x.eval(sl1);
	sl[1] = f3y.eval(sl1);
	sl[2] = f3z.eval(sl1);

	// sl2 is nabla f

	m_mult_vec_mat(sl, m_euler, nablaf_eul);

	memcpy(fnormal, nablaf_eul, 3 * sizeof(double));

	if (only_normal) {
		return 0;
	}


	double norm_nabla_f_eul = m_norm(nablaf_eul, 3);

	//
	// the hessian is A^t H A
	// where A is the euler matrix m_euler and H is
	// the matrix fhess of second derivatives
	//

	fhess[0] = f3xx.eval(sl1);
	fhess[1] = f3xy.eval(sl1);
	fhess[2] = f3xz.eval(sl1);

	fhess[3] = fhess[1];
	fhess[4] = f3yy.eval(sl1);
	fhess[5] = f3yz.eval(sl1);

	fhess[6] = fhess[2];
	fhess[7] = fhess[5];
	fhess[8] = f3zz.eval(sl1);

	double fhess1[9];

	//cout << "m_euler = " << endl;
	//m_print(m_euler);

	m_mult(fhess, m_euler, fhess1);
	m_mult(m_euler_transpose, fhess1, fhess);

	double matP[9];
	make_matP(nablaf_eul, norm_nabla_f_eul, matP);

	m_mult(fhess, matP, fhess1);
	m_skalmult_mat(1/norm_nabla_f_eul, fhess1, fhess );

	m_transpose_mat(fhess, fhessian);

	//memcpy(fhessian, fhess1, 9 * sizeof(double));
	//memcpy(fnormal, nablaf_eul, 3 * sizeof(double));

	return 0;
}



int ProcPoly::eval_coefs_poly(Poly3 & f3, double x, double y, double* coefs_lis)
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
		double fi = eval_f3_mat(f3, sl);
		for(int j = 0; j <= akt_deg_global; ++j) {
			coefs[j] += fi * lagrange_basis[i][j];
		}
	}

	memcpy(coefs_lis, coefs, sizeof(coefs));

	return 0;
}

int ProcPoly::eval_coefs_poly_point_normal(Poly3 & f3, double x, double y, double z,
			double nx, double ny, double nz, double* coefs_lis)
{
	// hier die Lagrangeinterpolation einfügen
	// z wird nacheinander z0, z1, z2,....zn gesetzt mit n = deg f3

	double sl[3];
	double coefs[max_deg + 1];
	for(int i = 0; i <= max_deg; ++i) {
		coefs[i] = 0;
	}
	for(int i = 0; i <= akt_deg_global; ++i) {
		double t = akt_xbase[i];
		sl[0] = x + t * nx;
		sl[1] = y + t * ny;
		sl[2] = z + t * nz;
		double fi = f3.eval(sl);
		for(int j = 0; j <= akt_deg_global; ++j) {
			coefs[j] += fi * lagrange_basis[i][j];
		}
	}

	memcpy(coefs_lis, coefs, sizeof(coefs));

	return 0;

}







const double clip_radius = 8; //20;

inline bool in_clip_radius(double x, double y, double z)
{
	return x * x + y * y + z * z <= clip_radius * clip_radius;
}


//static double min_disc_abs = DBL_MAX;


void m_fold(const double* m_image, const double* m_mask, int n, double &result) {
	double s = 0.0;
	for(int i = 0; i < n * n; ++i) {
		s += m_image[i] * m_mask[i];
	}
	result = s;
}




#define RAND_COL ((int)(255.0*((float)rand())/RAND_MAX))


const int max_roots = 10;

const int arr_size = 2 * gl_win_size;

ProcPoly pp(arr_size, max_roots, gl_win_size, xrast_to_x, yrast_to_y);




// public interface of ProcPoly


// can not stand alone uses m_euler
int ProcPoly::rotate_mat(double phi, double theta, double psi)
{
	const int ndim = 3;
	const int nsize = ndim * ndim;

	double m_z_phi[nsize];
	double m_x_theta[nsize];
	double m_z_psi[nsize];

	double m_aux[nsize];

	m_rot_z(phi, m_z_phi);
	m_rot_x(theta, m_x_theta);
	m_rot_z(psi, m_z_psi);

	// note that f(m_euler x) = 0
	// means f(x) = 0 is turned by m_euler^{-1}

	m_mult(m_x_theta, m_z_phi, m_aux);
	m_mult(m_z_psi, m_aux, m_euler);

	m_transpose_mat(m_euler, m_euler_transpose);

	return 0;

}



int ProcPoly::eval_parms_poly(Poly5 poly_in, Poly3 & poly_out, double *parm)
{
	Poly4 p4_aux;
	poly_in.eval_last(parm[1], p4_aux);
	p4_aux.eval_last(parm[0], poly_out);

	return 0;
}

void ProcPoly::eval_poly_f3(double x, double y, double z, double & resval)
{
	double sl[3];
	sl[0] = x;
	sl[1] = y;
	sl[2] = z;

	//m_mult_mat_vec(sl, m_euler, sl1);

	resval = f3.eval(sl);
}

int ProcPoly::eval_poly_poly_f3(double x, double y, double z, double & f, double* fnormal, double* fhessian)
{
	return eval_poly_poly(f3, x, y, z, f, fnormal, fhessian, false);
}



int ProcPoly::init_f3_diff(Poly3 & f3)
{
	f3x = f3;
	f3x.diff(0);
	f3y = f3;
	f3y.diff(1);
	f3z = f3;
	f3z.diff(2);

	f3xx = f3x;
	f3xx.diff(0);

	f3xy = f3x;
	f3xy.diff(1);

	f3xz = f3x;
	f3xz.diff(2);

	f3yy = f3y;
	f3yy.diff(1);

	f3yz = f3y;
	f3yz.diff(2);

	f3zz = f3z;
	f3zz.diff(2);


	akt_deg_global = f3.degree();

	lb_gen_lagrange_basis(akt_deg_global, akt_xbase, lagrange_basis);



	return 0;
}

int ProcPoly::init_f5(QString f5str)
{

	f5.read(f5str.toLatin1().data());

	return 0;

}

int ProcPoly::root_list_poly_point_normal(double x, double y, double z,
		double nx, double ny, double nz, double* root_list, int & root_list_len)
{
	double poly_coefs[max_deg + 1];

	eval_coefs_poly_point_normal(f3, x, y, z, nx, ny, nz, poly_coefs);

	// coefficient of leading monomial is in poly_coefs[0]
	// akt_deg is intended degree

	int deg_new;

	normalize_poly_coefs(poly_coefs, akt_deg_global, deg_new);

	root_final_list(deg_new, poly_coefs, 20, root_list, root_list_len);

	return 0;

}

// is private function
int ProcPoly::get_z_intersect_poly(double x, double y, double  & z, double  & n_z, bool & disc_zero, double* full_z_list,
		double* shape_mat, double* base_v1, double* base_v2, double & l1, double & l2, int & jsel, int & nsel,
		int shade_type)
{

	double poly_coefs[max_deg + 1];

	double z_erg;
	double z_erg_new = M_INF;
	double z_erg_list[max_deg + 1];

	eval_coefs_poly(f3, x, y, poly_coefs);

	// coefficient of leading monomial is in poly_coefs[0]
	// akt_deg is intended degree

	int deg;
	int num_z_erg;

	normalize_poly_coefs(poly_coefs, akt_deg_global, deg);

	double disc_poly = comp_disc(deg, poly_coefs);

	if ( 1 /*disc_poly != 0*/) {

		root_final_list(deg, poly_coefs, 20, z_erg_list, num_z_erg);

		nsel = num_z_erg;

		for(int i = 0; i < num_z_erg; ++i) {
			full_z_list[i] = z_erg_list[i];
		}

		int j;

		z_erg_new = M_INF;

		j = num_z_erg - 1;

		while (j >= 0) {
			if (in_clip_radius(x,y,z_erg_list[j])) {

				z_erg_new = z_erg_list[j];
				jsel = j;

				break;
			}
			--j;
		}

		disc_zero = false;

	} else {

		disc_zero = true;

		z_erg_new = M_INF;
	}

	z_erg = z_erg_new;


	if (! in_clip_radius(x, y, z_erg)) {
		z_erg = M_INF;
	}

	z = z_erg;

	m_set(shape_mat, 0, 4);
	m_set(base_v1, 0, 3);
	m_set(base_v2, 0, 3);


	if (z_erg > M_INF) {
		double f, fx, fy, fz;

		double fnormal[3];
		double fhessian[9];

		bool only_normal = (shade_type == 1) || (shade_type == 2);

		eval_poly_poly(f3, x, y, z_erg, f, fnormal, fhessian, only_normal);

		fx = fnormal[0];
		fy = fnormal[1];
		fz = fnormal[2];

		//cout << "before norm: fx = " << fx << " fy = " << fy << " fz = " << fz << endl;

		double normf = sqrt(fx * fx + fy * fy + fz * fz);

		if (normf > 0) {
			for(int i = 0; i < 3; ++i) {
				fnormal[i] *= (1/normf);
			}
		}

		fx = fnormal[0];
		fy = fnormal[1];
		fz = fnormal[2];

		if (shade_type == 1 || shade_type == 2) {
			n_z = fz;
			return 0;
		}

		//cout << "fx = " << fx << " fy = " << fy << " fz = " << fz << endl;

		int rescode = restrict_hessian_to_tangentspace(fnormal, fhessian, shape_mat, base_v1, base_v2);

		l1 = l2 = 0;

#if 0
		std::cout << "z_erg = " << z_erg << std::endl;
		std::cout << "rescode = " << rescode << std::endl;
		std::cout << "a = " << shape_mat[0] << " ";
		std::cout << "b = " << shape_mat[1] << " ";
		std::cout << "c = " << shape_mat[2] << " ";
		std::cout << "d = " << shape_mat[3] << " " << std::endl;
#endif

		if (rescode != 2) {
			m_calculate_diagonal_form(shape_mat, l1, l2, base_v1, base_v2);
		}

		assert(l1 <= l2);

		n_z = fz;
	} else {
		n_z = 0;
	}

	return 0;
}


void ProcPoly::fill_arrays(int xmax, int ymax, int shade_type) {

	int x;
	int y;

	bool disc_zero;


	for (x = 0; x < xmax; ++x) {
    	for (y = 0; y < ymax; ++y) {

    		//double y1 = (y - win_offset)/SCALE * local_scale;
    		//double x1 = (x - win_offset)/SCALE * local_scale;

    		double y1 = yrast_to_y(y);
    		double x1 = xrast_to_x(x);

    		double z;
    		double n_z;

    		double full_z_list[max_roots];
    		int jsel = -1;
    		int nsel = -1;

    		double shape_mat[4];
    		double vbase_1[3];
    		double vbase_2[3];

    		double l1 = 0;
    		double l2 = 0;

    		get_z_intersect_poly(x1, y1, z, n_z, disc_zero, full_z_list, shape_mat,
    				vbase_1, vbase_2, l1, l2, jsel, nsel, shade_type);


    		z_buf[x][y] = z;
    		jsel_buf[x][y] = jsel;
    		nsel_buf[x][y] = nsel;

    		if (z > M_INF) {
    			n_buf[x][y] = n_z;
    		} else {
    			n_buf[x][y] = 0;
    		}

    		if (shade_type == 1 || shade_type == 2) {
    			continue;
    		}


    		for(int i = 0; i < 4; ++i) {
    			shape_mat_buf[x][y][i] = shape_mat[i];
    		}

    		l1_buf[x][y] = l1;
    		l2_buf[x][y] = l2;

    		for(int i = 0; i < 3; ++i) {
    			vbase_1_buf[x][y][i] = vbase_1[i];
    			vbase_2_buf[x][y][i] = vbase_2[i];
    		}

    		double norm1 = m_norm(vbase_1, 2);
    		double norm2 = m_norm(vbase_2, 2);

    		for(int i = 0; i < 2; ++i) {
    			if (norm1 != 0)
    				vbase_1[i] /= norm1;
    			if (norm2 != 0)
    				vbase_2[i] /= norm2;
    		}

    		for(int i = 0; i < 2; ++i) {
    			v1_buf[x][y][i] = vbase_1[i];
    			v2_buf[x][y][i] = vbase_2[i];
    		}

    		assert(l1_buf[x][y] <= l2_buf[x][y]);

    		norm1 = m_norm(vbase_1, 2);
    		norm2 = m_norm(vbase_2, 2);

    		assert(fabs(norm1 - 1) < 1e-6 || norm1 < 1e-6);
    		assert(fabs(norm2 - 1) < 1e-6 || norm2 < 1e-6);

    		// memcpy(zfull_buf[x][y], full_z_list, max_roots * sizeof(double));

    		for(int i = 0; i < max_roots; ++i)
    			zfull_buf[x][y][i] = full_z_list[i];


#if CHECK
    		for(int ii = 0; ii < max_roots; ++ii) {
    			assert(zfull_buf[x][y][ii] == full_z_list[ii]);
    		}
#endif

    	}
    }

}

void ProcPoly::calc_silhouette(int xmax, int ymax, double & k_gauss_min, double & k_gauss_max) {

	double local_scale = akt_win_size/xmax;
	int win_offset = xmax/2;

	UNUSED(local_scale);
	UNUSED(win_offset);

	k_gauss_min = DBL_MAX;
    k_gauss_max = -DBL_MAX;

    int x;
    int y;

    silhouette_pointl.clear();

    for (x = 0; x < xmax; ++x) {
    	for (y = 0; y < ymax; ++y) {

    		//double y1 = (y - win_offset)/SCALE * local_scale;
    		//double x1 = (x - win_offset)/SCALE * local_scale;

			is_silhouette_mat[x][y] = false;

    		double a = shape_mat_buf[x][y][0];
    		double b = shape_mat_buf[x][y][1];
    		double c = shape_mat_buf[x][y][2];
    		double d = shape_mat_buf[x][y][3];

    		double k_gauss = a * d - b * c;
    		//double k_gauss_1 = l1_buf[x][y] * l2_buf[x][y];


    		if (z_buf[x][y] != M_INF) {

        		//assert(fabs(k_gauss - k_gauss_1) < 1e-6);


    			if (k_gauss > k_gauss_max) {
					k_gauss_max = k_gauss;
				}

				if (k_gauss < k_gauss_min) {
					k_gauss_min = k_gauss;
				}
    		}

    		double zmat[9];
    		double nmat[9];

    		int jnextmat[9];

    		UNUSED(zmat);
    		UNUSED(nmat);
    		UNUSED(jnextmat);

    		for(int i = 0; i < 9; ++i) {
    			zmat[i] = M_INF;
    			nmat[i] = 0.0;
    		}

    		for(int i = ::max(x-1,0); i <= ::min(x + 1, xmax - 1); ++i ) {
    			for(int j = ::max(y-1,0); j <= ::min(y + 1, ymax - 1); ++j) {

    				int r = i - x + 1;
    				int s = j - y + 1;

    				int jnext = -1;

    				double minval = DBL_MAX;
    				double testval = z_buf[x][y];
    				//double cnt_near = 0;

    				for(int ii = 0; ii < nsel_buf[i][j]; ++ii) {
    					double diffval = fabs(zfull_buf[i][j][ii] - testval);
    					//if (diffval < 0.1 && (i != x && j != y)) {
    					//	++cnt_near;
    					//}
    					if (diffval < minval) {
    						minval = diffval;
    						jnext = ii;
    					}
    				}

    				jnextmat[r + 3 * s] = jnext;

    				bool is_on_same_plane = jsel_buf[i][j] >= 0
    						&& fabs(zfull_buf[i][j][jnext] - zfull_buf[i][j][jsel_buf[i][j]]) < 1e-6;

    				if (!is_on_same_plane
    						/*jnext != jsel_buf[i][j] */|| (testval != M_INF && z_buf[i][j] == M_INF) /* || (cnt_near > 1) */) {

    					if (!is_silhouette_mat[x][y] && z_buf[x][y] != M_INF) {
    						silhouette_pointl.push_back(Point2D(x,y));
    					}

    					is_silhouette_mat[x][y] = true;

    				}

    			}
    		}

    	}
    }

}

