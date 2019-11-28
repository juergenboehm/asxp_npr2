#ifndef __procpoly_h_
#define __procpoly_h_

#include <QString>

#include "asxp_arrays.h"

#include "poly.h"
#include "streamline.h"
#include "streamplot.h"
//#include "asxp_arrays.h"
#include "pointlist.h"


class ProcPoly {

public:

	ProcPoly(const int arr_size, const int max_rootsa, const int akt_win_sizea,
			Scale & xrast_to_xa, Scale & yrast_to_ya, Scale & x_to_xrasta, Scale & y_to_yrasta);

	int rotate_mat(double phi, double theta, double psi);

	int eval_parms_poly(Poly5 poly_in, Poly3 & poly_out, double *parm);


	int init_f5(QString f5str);

	int init_f3_diff(Poly3 & f3);


	void eval_poly_f3(double x, double y, double z, double & resval);

	int eval_poly_poly_f3(double x, double y, double z, double & f, double* fnormal, double* fhessian);


	int root_list_poly_point_normal(double x, double y, double z,
		double nx, double ny, double nz, double* root_list, int & root_list_len);


	void fill_arrays(int xmax, int ymax, int shade_type);

	void calc_silhouette(int xmax, int ymax, double & k_gauss_min, double & k_gauss_max);


	void reset_full_plots();



	void compute_streamfield_CGAL(int xmax, int ymax);

	void compute_streamfield(int xmax, int ymax);



	Poly5 f5;
	Poly3 f3;

	Array2d_double z_buf;
	Array2d_double n_buf;


	Array3d_double zfull_buf;
	Array2d_int jsel_buf;
	Array2d_int nsel_buf;

	//Array2d_double zsobel_arr(boost::extents[arr_size][arr_size]);
	//Array2d_double nsobel_arr(boost::extents[arr_size][arr_size]);

	Array2d_bool is_silhouette_mat;

	Array3d_double shape_mat_buf;

	Array3d_double vbase_1_buf;
	Array3d_double vbase_2_buf;

	Array3d_double v1_buf;
	Array3d_double v2_buf;

	Array2d_double l1_buf;
	Array2d_double l2_buf;


	Array2d_double* pz_buf;
	Array2d_double* pn_buf;


	Array3d_double* pzfull_buf;
	Array2d_int* pjsel_buf;
	Array2d_int* pnsel_buf;

	Array2d_bool* pis_silhouette_mat;

	Array3d_double* pshape_mat_buf;

	Array3d_double* pvbase_1_buf;
	Array3d_double* pvbase_2_buf;

	Array3d_double* pv1_buf;
	Array3d_double* pv2_buf;

	Array2d_double* pl1_buf;
	Array2d_double* pl2_buf;


	Point2DList silhouette_pointl;

	Streamplot* full_plot1;
	Streamplot* full_plot2;


	double dsep;

	int streamgen_type; // 0 for own 1 for CGAL




private:


	static const int max_deg = 20;

	int max_roots;


	Scale & xrast_to_x;
	Scale & yrast_to_y;

	Scale & x_to_xrast;
	Scale & y_to_yrast;


	int akt_win_size;


	int get_z_intersect_poly(double x, double y, double  & z, double  & n_z, bool & disc_zero, double* full_z_list,
			double* shape_mat, double* base_v1, double* base_v2, double & l1, double & l2, int & jsel, int & nsel,
			int shade_type);

	void lb_poly_mult(int deg, double * coefs, double a);

	void lb_poly_scal_mult(int deg, double * coefs, double a);

	void lb_gen_lagrange_basis(int deg, double* xbase, double lagr_basis[][max_deg + 1]);

	double eval_f3_mat(Poly3 & f3, double *sl);

	int eval_poly_poly(Poly3 & f3,
			double x, double y, double z, double & f, double* fnormal, double* fhessian, bool only_normal);

	int eval_coefs_poly(Poly3 & f3, double x, double y, double* coefs_lis);

	int eval_coefs_poly_point_normal(Poly3 & f3, double x, double y, double z,
				double nx, double ny, double nz, double* coefs_lis);


	void normalize_poly_coefs(double* poly_coefs, int akt_deg, int & deg_new);


	Poly3 f3x, f3y, f3z;

	Poly3 f3xx, f3xy, f3xz, f3yy, f3yz, f3zz;

	double m_euler[3 * 3];
	double m_euler_transpose[3 * 3];

	int akt_deg_global;

	double akt_xbase[max_deg+1] = { 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0 };

	double lagrange_basis[max_deg+1][max_deg+1];




};


class FindSilhouette : public PointClassifier {

public:

	virtual ~FindSilhouette() {};

	virtual int operator()(double x, double y);
	void set_xy_max(int xmaxa, int ymaxa);

	int xmax;
	int ymax;
};

class FindBackground : public PointClassifier {

public:

	virtual ~FindBackground() {};

	virtual int operator()(double x, double y);
	void set_xy_max(int xmaxa, int ymaxa);

	int xmax;
	int ymax;
};




extern ProcPoly pp;

extern Scale xrast_to_x, yrast_to_y, zrast_to_z, x_to_xrast, y_to_yrast, z_to_zrast;

extern FindSilhouette is_silhouette_point;




//#define EPS 1e-12
#define M_INF -1e38




#endif /* PROCPOLY_H_ */
