
#include <QtGui>
#include <QString>
#include <QLabel>

#include <list>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <timing.h>

#include "asxp_arrays.h"



#include "poly.h"
#include "parser.h"
#include "roots.h"

#include "screenwidget.h"
#include "asxp.h"

#include "eigen.h"

#include "streamline.h"
#include "pointlist.h"
#include "streamplot.h"

#define NOCUDA

#include "pimage.h"


int phong_exponent;
 

#define Q_SQ(x) ((x) * (x))
#define Q_NORM(x,y) (Q_SQ(x) + Q_SQ(y))
#define Q_SKALPROD(x1, y1, x2, y2) ((x1)*(x2) + (y1)*(y2))

#define TO_GRAD(x) ((x)/(2 * 3.14159265) * 360)
#define SIGN(x) ((x) > 0 ? 1 : (((x) == 0) ? 0 : -1))

#define UNUSED(expr) do { (void)(expr); } while (0)



inline int & mref(int* mat, int i, int j)
{
	return mat[i * 2 * gl_win_size + j];
}


PaintHelper::PaintHelper()
 {
    QLinearGradient gradient(QPointF(50, -20), QPointF(80, 20));
    gradient.setColorAt(0.0, Qt::white);
    gradient.setColorAt(1.0, QColor(0xa6, 0xce, 0x39));
 
    background = QBrush(QColor(200, 200, 200));
    circleBrush = QBrush(gradient);
    circlePen = QPen(Qt::black);
    circlePen.setWidth(1);
    textPen = QPen(Qt::white);
    textFont.setPixelSize(50);

    const int colmat_size = 4 * gl_win_size * gl_win_size * sizeof(int);

    colmat_r = (int*)malloc(colmat_size);
    colmat_g = (int*)malloc(colmat_size);
    colmat_b = (int*)malloc(colmat_size);

    init_colmat();

    full_plot1 = 0;
    full_plot2 = 0;

    dsep = 10;

    displayFlowLines = false;
    displayCrossField = true;

    streamLineColor = 0;

    streamgen_type = 0;

    colmat_valid = true;
}

void PaintHelper::init_colmat()
{
	for(int x = 0; x < gl_win_size; ++x) {
		for(int y = 0; y < gl_win_size; ++y) {
			mref(colmat_r, x, y) = 200; //64;
			mref(colmat_g, x, y) = 200;// 32;
			mref(colmat_b, x, y) = 200;  // 64;
		}
	}
}

PaintHelper::~PaintHelper()
{
	free(colmat_r);
	free(colmat_g);
	free(colmat_b);
}

#define EPS 1e-12
#define M_INF -1e38

#define POLY_DEG 3
#define POLY_DEG1 (POLY_DEG+1)

Poly5 f5;
Poly3 f3;
Poly3 f3x, f3y, f3z;

Poly3 f3xx, f3xy, f3xz, f3yy, f3yz, f3zz;

double m_euler[3 * 3];
double m_euler_transpose[3 * 3];

const int max_deg = 20;

int akt_deg_global;

double akt_xbase[max_deg+1] = { 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 };

double lagrange_basis[max_deg+1][max_deg+1];

void lb_poly_mult(int deg, double * coefs, double a)
{
	// multiply coefs poly coefs[0] * t^deg + ... with (t + a)
	coefs[deg+1] = 0;
	for(int i = deg + 1; i >=  1; --i) {
		coefs[i] += a * coefs[i-1];
	}
}

void lb_poly_scal_mult(int deg, double * coefs, double a)
{
	for(int i = 0; i <= deg; ++i) {
		coefs[i] *= a;
	}
}


void lb_gen_lagrange_basis(int deg, double* xbase, double lagr_basis[][max_deg + 1])
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


double eval_f3_mat(Poly3 & f3, double *sl)
{
	double sl1[3];
	m_mult_mat_vec(sl, m_euler, sl1);

	double f = f3.eval(sl1);

	return f;

}

void eval_poly_f3(double x, double y, double z, double & resval)
{
	double sl[3];
	sl[0] = x;
	sl[1] = y;
	sl[2] = z;

	//m_mult_mat_vec(sl, m_euler, sl1);

	resval = f3.eval(sl);
}


int eval_poly_poly(Poly3 & f3,
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

int eval_poly_poly_f3(double x, double y, double z, double & f, double* fnormal, double* fhessian)
{
	return eval_poly_poly(f3, x, y, z, f, fnormal, fhessian, false);
}


int eval_coefs_poly(Poly3 & f3, double x, double y, double* coefs_lis)
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

int eval_coefs_poly_point_normal(Poly3 & f3, double x, double y, double z,
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

void normalize_poly_coefs(double* poly_coefs, int akt_deg, int & deg_new)
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

int root_list_poly_point_normal(double x, double y, double z,
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

int eval_parms_poly(Poly5 poly_in, Poly3 & poly_out, double *parm)
{
	Poly4 p4_aux;
	poly_in.eval_last(parm[1], p4_aux);
	p4_aux.eval_last(parm[0], poly_out);

	return 0;
}


int rotate_mat(double phi, double theta, double psi)
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

int init_f3_diff(Poly3 & f3)
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

int init_f5(QString f5str)
{

	f5.read(f5str.toLatin1().data());

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

int get_z_intersect_poly(double x, double y, double  & z, double  & n_z, bool & disc_zero, double* full_z_list,
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

#define RAND_COL ((int)(255.0*((float)rand())/RAND_MAX))


const int win_size =  gl_win_size;
const int max_roots = 10;

const int arr_size = 2 * gl_win_size;

Array2d_double z_buf(boost::extents[arr_size][arr_size]);
Array2d_double n_buf(boost::extents[arr_size][arr_size]);


Array3d_double zfull_buf(boost::extents[arr_size][arr_size][max_roots]);
Array2d_int jsel_buf(boost::extents[arr_size][arr_size]);
Array2d_int nsel_buf(boost::extents[arr_size][arr_size]);

//Array2d_double zsobel_arr(boost::extents[arr_size][arr_size]);
//Array2d_double nsobel_arr(boost::extents[arr_size][arr_size]);

Array2d_bool is_silhouette_mat(boost::extents[arr_size][arr_size]);

Array3d_double shape_mat_buf(boost::extents[arr_size][arr_size][4]);

Array3d_double vbase_1_buf(boost::extents[arr_size][arr_size][3]);
Array3d_double vbase_2_buf(boost::extents[arr_size][arr_size][3]);

Array3d_double v1_buf(boost::extents[arr_size][arr_size][2]);
Array3d_double v2_buf(boost::extents[arr_size][arr_size][2]);

bool is_vfs_filled = false;

Array2d_double l1_buf(boost::extents[arr_size][arr_size]);
Array2d_double l2_buf(boost::extents[arr_size][arr_size]);

Point2DList silhouette_pointl;



Array2d_double* pz_buf = &z_buf;
Array2d_double* pn_buf = &n_buf;


Array3d_double* pzfull_buf = & zfull_buf;
Array2d_int* pjsel_buf = &jsel_buf;
Array2d_int* pnsel_buf = &nsel_buf;

Array2d_bool* pis_silhouette_mat = &is_silhouette_mat;

Array3d_double* pshape_mat_buf = & shape_mat_buf;

Array3d_double* pvbase_1_buf = & vbase_1_buf;
Array3d_double* pvbase_2_buf = & vbase_2_buf;

Array3d_double* pv1_buf = & v1_buf;
Array3d_double* pv2_buf = & v2_buf;

Array2d_double* pl1_buf = & l1_buf;
Array2d_double* pl2_buf = & l2_buf;



void PaintHelper::getImageInfo(int x, int y, double & a, double & b, double & c, double & d )
{

	a = l1_buf[x][y];
	b = l2_buf[x][y];

	c = a * b;
	d = a + b;

	a = nsel_buf[x][y];
	b = jsel_buf[x][y];

}



//#define SCALE 10.0
#define SCALE 50.0

void calc_brightness( double nz, bool & is_backside, double & brightness, double & highlight) {

	double nabs = fabs(nz);

	is_backside = (nz < 0);

	brightness = nabs/2;

	double phong_kernel = 2 * nz * nz - 1;
	double spec_coef = 0.3 * pow(phong_kernel, phong_exponent);

	highlight = spec_coef;

}


void shade_color(double nz, double & color_red, double & color_green, double & color_blue ) {

	double brightness;
	double highlight;
	bool is_backside;

	calc_brightness(nz, is_backside, brightness, highlight);

	if (is_backside) {
		color_red = brightness + highlight;
		color_green = 0;
		color_blue = 0;

	} else {
		color_red = 0;
		color_green = brightness + highlight;
		color_blue = 0;
	}

}

Array2d_double floyd_stein_buf(boost::extents[win_size + 1][win_size + 1]);

void shade_dots(double nz, double & color_red, double & color_green, double & color_blue ) {

	double brightness;
	double highlight;
	bool is_backside;

	calc_brightness(nz, is_backside, brightness, highlight);

	color_red = brightness + highlight;
	color_green = brightness + highlight;
	color_blue = brightness + highlight;

	double rand_thr = ((double (rand()))/(RAND_MAX));


	if (color_red > rand_thr) {
		color_red = 0.0;
		color_green = 0.0;
		color_blue = 0.0;
	} else {
		color_red = 1.0;
		color_green = 1.0;
		color_blue = 1.0;
	}


}

void shade_grey(double nz, double & c_grey)
{
	double brightness;
	double highlight;
	bool is_backside;

	calc_brightness(nz, is_backside, brightness, highlight);

	c_grey = brightness + highlight;

}

void shade_floyd_steinberg(int x, int y, double & color_red, double & color_green, double & color_blue)
{
	bool is_set = floyd_stein_buf[x][y] > 0.5;

	double error;

	if (is_set) {
		error = floyd_stein_buf[x][y] - 1;

		color_red = color_green = color_blue = 1;

	} else {

		error = floyd_stein_buf[x][y];

		color_red = color_green = color_blue = 0;
	}

	floyd_stein_buf[x + 1][y] = floyd_stein_buf[x + 1][y] + (3 * error) / 8;
	floyd_stein_buf[x][y + 1] = floyd_stein_buf[x][y + 1] + (3 * error) / 8;
	floyd_stein_buf[x + 1][y + 1] = floyd_stein_buf[x + 1][y + 1] + error / 4;

}

Scale xrast_to_x, yrast_to_y, x_to_xrast, y_to_yrast;


//  y
//  ^
//  |   a        b
//  |
//  |    ..sx *
//  |        sy
//  |         .
//  |   c     .  d
//  |
//  -----------------> x



static double bilin_interpolate(double a, double b, double c, double d, double sx, double sy) {

	double val_a = a * sy * (1 - sx);
	double val_b = b * sx * sy;
	double val_c = c * (1 - sx) * (1 - sy);
	double val_d = d * sx * (1 - sy);

	double val = val_a + val_b + val_c + val_d;

	return val;
}


void fill_arrays(int xmax, int ymax, int shade_type) {

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

void calc_silhouette(int xmax, int ymax, double & k_gauss_min, double & k_gauss_max) {

	double local_scale = gl_win_size/xmax;
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


void PaintHelper::compute_colmat(double a, double b, int xmax, int ymax)
{

	double local_scale = gl_win_size/xmax;

	double parm[2];
	int win_offset = xmax/2;
	bool disc_zero;

	UNUSED(win_offset);
	UNUSED(disc_zero);

	UNUSED(local_scale);

	int x;
	int y;

	cout << "xmax = " << xmax << "ymax = " << ymax << endl;


	parm[0] = a;
	parm[1] = b;

    eval_parms_poly(f5, f3, parm);
    rotate_mat(euler_phi, euler_theta, euler_psi);

	init_f3_diff(f3);

#ifdef RANDOM_COLMAT
	for (x = 0; x < gl_win_size; ++x) {
		for( y = 0; y < gl_win_size; ++y) {
			mref(colmat_r, x, y) = RAND_COL;
			mref(colmat_g, x, y) = RAND_COL;
			mref(colmat_b, x, y) = RAND_COL;
		}
	}

	return;
#endif

	fill_arrays(xmax, ymax, shade_type);

	double k_gauss_min, k_gauss_max;

	if (!(shade_type == 1 || shade_type == 2)) {
		calc_silhouette(xmax, ymax, k_gauss_min, k_gauss_max);
	}

	//cout << "k_gauss_min = " << k_gauss_min << " k_gauss_max = " << k_gauss_max << endl;

	if (shade_type == 2) {
		for(x = 0; x < xmax; ++x) {
			for(y = 0; y < ymax; ++y) {

				double n = n_buf[x][y];
				double z = z_buf[x][y];

				if (z > M_INF) {

					double c_grey; // between 0.0 and 1.0

					shade_grey(n, c_grey);

					floyd_stein_buf[x][y] = c_grey;
				} else
					floyd_stein_buf[x][y] = 1;
			}
		}
	}

    for (x = 0; x < xmax; ++x) {
    	for (y = 0; y < ymax; ++y) {

    		double z;
    		double n;

    		z = z_buf[x][y];
    		n = n_buf[x][y];

    		if (z > M_INF) {

    			double color_red = 1.0;
    			double color_green = 1.0;
    			double color_blue = 1.0;

    			bool is_silhouette = false;
    			bool is_true_silhouette = false;

    			if (!(shade_type == 1 || shade_type == 2)) {

    				is_silhouette = is_silhouette_mat[x][y];

					is_true_silhouette = false;

					if (fabs(n) <= nsobel_thresh * 0.2) {
						is_silhouette = true;
						is_true_silhouette = true;
					}

    			}

    			if (is_silhouette == true) {

    				v1_buf[x][y][0] = 0;
    				v1_buf[x][y][1] = 0;

    				v2_buf[x][y][0] = 0;
    				v2_buf[x][y][1] = 0;

    				color_red = 0.0;
    				color_green = 0.0;
    				color_blue = 0.0;

    				if (is_true_silhouette) {
    					color_red = 1.0;
    				}

    			} else {

    				if (shade_type == 1) {
    					shade_color(n, color_red, color_green, color_blue);
    				}

    				if (shade_type == 2) {
    					shade_floyd_steinberg(x, y, color_red, color_green, color_blue);
    				}

    			}

    			int col_red = (int) (250 * color_red);
    			int col_green = (int) (250 * color_green);
    			int col_blue = (int) (250 * color_blue);

    			mref(colmat_r, x, y) = col_red;
    			mref(colmat_g, x, y) = col_green;
    			mref(colmat_b, x, y) = col_blue;

    		} else {

    			// background color
    			mref(colmat_r, x, y) = 250; // 64
    			mref(colmat_g, x, y) = 250; // 32
    			mref(colmat_b, x, y) = 250; // 64

    		}
    	}
    }


}



int translate_poly(cudaPoly3 & f3_h, Poly3 & f3)
{

	Poly3::vecIt it;
	int i;

	memset(f3_h.xexpo, 0, 40 * sizeof(int));
	memset(f3_h.yexpo, 0, 40 * sizeof(int));
	memset(f3_h.zexpo, 0, 40 * sizeof(int));


	for(i = 0, it = f3.terms.begin(); it < f3.terms.end(); ++it, ++i) {
		f3_h.coefs[i] = it->coef;
		f3_h.xexpo[i] = it->expos[0];
		f3_h.yexpo[i] = it->expos[1];
		f3_h.zexpo[i] = it->expos[2];
	}

	f3_h.len = i;
}


void PaintHelper::compute_colmat_gpu(double a, double b, int xmax, int ymax)
{

	double local_scale = gl_win_size/xmax;

	double parm[2];
	int win_offset = xmax/2;
	bool disc_zero;

	int x;
	int y;

	cout << "xmax = " << xmax << "ymax = " << ymax << endl;


	parm[0] = a;
	parm[1] = b;

    eval_parms_poly(f5, f3, parm);

    cudaPoly3 f3_cuda_h;

    translate_poly(f3_cuda_h, f3);

	gpu_compute_colmat(a, b, xmax, ymax, f3_cuda_h, euler_phi, euler_theta, euler_psi,
			colmat_r, colmat_g, colmat_b);

}


int data[gl_win_size * gl_win_size * sizeof(int)];

// (x, y) raster coordinates

void FindSilhouette::set_xy_max(int xmaxa, int ymaxa) {
	xmax = xmaxa;
	ymax = ymaxa;
}

int FindSilhouette::operator()(double x, double y) {

	int cx = int(x);
	int cy = int(y);

	if ((cx < 0 || cx >= xmax) || (cy < 0 || cy >= ymax)) {
		return true;
	} else {
		return is_silhouette_mat[cx][cy];
	}
}

FindSilhouette is_silhouette_point;

void FindBackground::set_xy_max(int xmaxa, int ymaxa) {
	xmax = xmaxa;
	ymax = ymaxa;
}

int FindBackground::operator()(double x, double y) {

	int cx = int(x);
	int cy = int(y);

	if ((cx < 0 || cx >= xmax) || (cy < 0 || cy >= ymax)) {
		return true;
	} else {
		return (z_buf[cx][cy] == M_INF);
	}
}

FindBackground is_background_point;


void PaintHelper::compute_streamfield_CGAL(int xmax, int ymax)
{

	const double dsep_multiplier = 4 * 0.6;

	if (!full_plot1) {
		full_plot1 = new Streamplot(&v1_buf, &v2_buf, xmax, ymax, dsep * dsep_multiplier,
								xrast_to_x, yrast_to_y, x_to_xrast, y_to_yrast,
								&is_background_point);

		full_plot1->set_mode(1);

		full_plot1->compute_stream_field_CGAL();
	}

	if (!full_plot2) {

		full_plot2 = new Streamplot(&v1_buf, &v2_buf, xmax, ymax, dsep * dsep_multiplier,
								xrast_to_x, yrast_to_y, x_to_xrast, y_to_yrast,
								&is_background_point);

		full_plot2->set_mode(2);

		full_plot2->compute_stream_field_CGAL();
	}

}

void PaintHelper::compute_streamfield(int xmax, int ymax)
{

	bool init_streamfield1 = true;
	bool init_streamfield2 = true;

	if (streamgen_type == 1) {
		compute_streamfield_CGAL(xmax, ymax);
		return;
	}

	for(int tstx = 0; tstx < xmax; tstx += 10) {

		for(int tsty = 0; tsty < ymax; tsty += 10) {

			char buf[128];

			sprintf(buf, "tstx = %d, tsty = %d", tstx, tsty);

			displayLabel->setText(QString(buf));
			displayLabel->repaint();

			cout << "tstx = " << tstx << " tsty = " << tsty << endl;

			while (1) {

				if (z_buf[tstx][tsty] == M_INF) {
					break;
				}

				double xstart, ystart;

				xstart = xrast_to_x(tstx);
				ystart = yrast_to_y(tsty);

#if 0
				Streamline flow_line(&v1_buf, &v2_buf, xmax, ymax, x_to_xrast, y_to_yrast, &is_silhouette_point);

				flow_line.integrate_from(xstart, ystart);
#endif


				if (!full_plot1) {

					full_plot1 = new Streamplot(&v1_buf, &v2_buf, xmax, ymax, dsep,  xrast_to_x, yrast_to_y, x_to_xrast, y_to_yrast,
											&is_silhouette_point);

					full_plot1->set_mode(1);

					init_streamfield1 = true;

				}

				cout << "compute streamfield 1" << endl;

				full_plot1->compute_stream_field(xstart, ystart, init_streamfield1);

				init_streamfield1 = false;

				if (!full_plot2) {

					full_plot2 = new Streamplot(&v1_buf, &v2_buf, xmax, ymax, dsep,  xrast_to_x, yrast_to_y, x_to_xrast, y_to_yrast,
											&is_silhouette_point);

					full_plot2->set_mode(2);


					init_streamfield2 = true;
				}

				cout << "compute streamfield 2" << endl;

				full_plot2->compute_stream_field(xstart, ystart, init_streamfield2);

				init_streamfield2 = false;

				break;

			}
		}
	}



}

void PaintHelper::paint(QPainter *painter, QPaintEvent *event, bool mouse_moved,
		int mousex, int mousey, double scale_im)
 {

	UNUSED(event);

	const int win_offset = gl_win_size/2;

	double a = (mousex - win_offset)/10.0;
	double b = (mousey - win_offset)/10.0;

	cout << "colmat_valid = " << colmat_valid << endl;


	double divider = scale_im;

	int xmax = gl_win_size/divider;
	int ymax = gl_win_size/divider;

	cout << "xmax = " << xmax << " ymax = " << ymax << endl;

	double local_scale_x = gl_win_size/((double)xmax);
	double local_scale_y = gl_win_size/((double)ymax);
	int win_offset_x = xmax/2;
	int win_offset_y = ymax/2;

	is_silhouette_point.set_xy_max(xmax, ymax);

	// do
	// x1 = (xrast - win_offset_x)/SCALE * local_scale;
	// and analogously to y1 =

	xrast_to_x.set(1/SCALE * local_scale_x, -win_offset_x/SCALE * local_scale_x);
	yrast_to_y.set(1/SCALE * local_scale_y, -win_offset_y/SCALE * local_scale_y);

	// xrast = win_offset + x * SCALE/local_scale;
	// and analogously to yrast =

	x_to_xrast.set(SCALE/local_scale_x, win_offset_x);
	y_to_yrast.set(SCALE/local_scale_y, win_offset_y);

	cout << "xrast_to_x(x_to_xrast(1)) = " << xrast_to_x(x_to_xrast(1)) << endl;

	//cin.get();


	assert(xrast_to_x(x_to_xrast(1)) == 1);
	assert(yrast_to_y(y_to_yrast(1)) == 1);
	assert(xrast_to_x(x_to_xrast(0)) == 0);
	assert(yrast_to_y(y_to_yrast(0)) == 0);

	cout << "divider = " << divider << endl;

	if (!colmat_valid) {

		bool isPrint = (shade_type == -1);

		if (isPrint) {
			shade_type = 0;
		}

		Timing time0;

		time0.start();

		//compute_colmat(a, b, gl_win_size/divider, gl_win_size/divider);
		if (shade_type != 1) {
			compute_colmat(a, b, xmax, ymax);
		} else {
			compute_colmat_gpu(a, b, xmax, ymax);
		}

		time0.stop();
		double elapsed = time0.elapsed();

		cout << "time to render = " << elapsed << endl;
		colmat_valid = true;

		if (full_plot1 != 0) {
			delete full_plot1;
			full_plot1= 0;
		}

		if (full_plot2 != 0) {
			delete full_plot2;
			full_plot2 = 0;
		}
		if (isPrint) {

			compute_streamfield(xmax, ymax);

		}
	}

    //if (!mouse_moved)
    //	painter->fillRect(event->rect(), background);

	// painter->translate(win_size/2, win_size/2);

	painter->scale(scale_im, scale_im);


    //painter->save();
    //painter->setBrush(circleBrush);
    //painter->setPen(textPen);

    if (mouse_moved) {
		char buf[128];
		sprintf(buf, "%.3f %.3f", a, b);
		painter->setPen(QPen(Qt::blue));
		painter->fillRect(99,100 - 12, 150, 13, background);
		painter->drawText(100,100, QString(buf));
    }

    int x;
    int y;

    bool is_pen_image = displayFlowLines && (shade_type == 0);

    if (! is_pen_image) {

		Timing time0;

		time0.start();

		int stride = gl_win_size/divider;

		for (x = 0; x < xmax; ++x) {
			for (y = 0; y < xmax; ++y) {

				int col_red = mref(colmat_r, x, y);
				int col_green = mref(colmat_g, x, y);
				int col_blue = mref(colmat_b, x, y);

				data[y * stride + x] = (255 << 24) + (col_red << 16) + (col_green << 8) + col_blue;

			}
		}
		time0.stop();

		cout << "time to copy = " << time0.elapsed() << endl;


		int im_size = gl_win_size / divider;

		QImage qimage((uchar*)data, im_size, im_size, QImage::Format_ARGB32);

		painter->drawImage(0,0, qimage);

    }

    if (colmat_valid && displayFlowLines) {

    	// draw silhouette points

    	if (shade_type == 0) {

    		painter->fillRect(0, 0, xmax, ymax, QBrush(QColor(255,255,255)));

    		painter->setPen(QPen(QColor(0,0,0)));
    		painter->setBrush(QBrush(QColor(0,0,0)));

			for(Point2DListIterator it = silhouette_pointl.begin(); it != silhouette_pointl.end(); ++it) {
				painter->drawPoint(it->x, it->y);
			}
    	}

    	// draw cross-field

    	if (displayCrossField) {

			for(int x = 0; x < xmax; ++x) {
				for(int y = 0; y < ymax; ++y) {

					if (z_buf[x][y] == M_INF) {
						continue;
					}

					if (x % 16 != 0 || y % 16 != 0) {
						continue;
					}

					const double c = 7;

					double v1x = c * v1_buf[x][y][0];
					double v1y = c * v1_buf[x][y][1];

					double v2x = c * v2_buf[x][y][0];
					double v2y = c * v2_buf[x][y][1];

					painter->setPen(QPen(Qt::blue));

					painter->drawLine(x, y, x + v1x, y + v1y);
					painter->drawLine(x, y, x - v1x, y - v1y);


					painter->setPen(QPen(Qt::green));

					painter->drawLine(x, y, x + v2x, y + v2y);
					painter->drawLine(x, y, x - v2x, y - v2y);


				}
			}
    	}
    	// get mouse (x,y) as start-point for integrate flow-line

    	if (full_plot1 &&  full_plot2) {

    		switch (streamLineColor) {
    			case 0: painter->setPen(QPen(QColor(0,0,0)));
    					break;
    			case 1:	painter->setPen(QPen(Qt::red));
    					break;
    			default: painter->setPen(QPen(Qt::blue));
    					break;
    		}
    	    // dummy while, only for break

    		Streamplot* full_plot_arr[2];

    		full_plot_arr[0] = full_plot1;
    		full_plot_arr[1] = full_plot2;

    		for(int ifpa = 0; ifpa < 2; ++ifpa) {

    			Streamplot* full_plot = full_plot_arr[ifpa];

    			if (full_plot == 0) {
    				continue;
    			}

    			cout << "no of streamlines = " << full_plot->streaml_res.size() << endl;
    			// cin.get();

    			Point2DList fringe_d = full_plot->debug_points;

    			for(Point2DListIterator pit = fringe_d.begin(); pit != fringe_d.end(); ++pit) {
    				painter->drawPoint(pit->x, pit->y);
    			}

    			for(list<Streamline>::iterator psl = full_plot->streaml_res.begin();
    						psl != full_plot->streaml_res.end(); ++psl) {

    				Point2DList & line_samples(psl->line_samples);

    				Point2DListIterator pit = line_samples.begin();

    				cout << "Anfang size = " << line_samples.size() << endl;

#if 0

    				while (pit != line_samples.end()) {
    					double x = pit->x;
    					double y = pit->y;
    					++pit;
    					cout << "( " << x << " , " << y << " ) ";
    				}
#endif
    				cout << endl << "Ende" << endl;

    				// line_samples contains screen coordinates

    				pit = line_samples.begin();

    				double x0 = pit->x;
    				double y0 = pit->y;

    				++pit;

    				while (pit != line_samples.end()) {
    					double x1 = pit->x;
    					double y1 = pit->y;

    					++pit;

    					painter->drawLine(x0, y0, x1, y1);

    					x0 = x1;
    					y0 = y1;

    				}

    			}

    			cout << "no of streamlines = " << full_plot->streaml_res.size() << endl;
    		}

    	}


	}

    char buf[128];
    sprintf(buf, "%.3f %.3f", a, b);
    painter->setPen(QPen(Qt::blue));
    painter->drawText(100,100, QString(buf));

}
