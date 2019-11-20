#ifndef __procpoly_h_
#define __procpoly_h_

#include <QString>

#include "asxp_arrays.h"

#include "poly.h"
#include "streamline.h"
#include "asxp_arrays.h"
#include "pointlist.h"


class ProcPoly {

public:

	ProcPoly(int arr_size, int max_roots):

		z_buf(boost::extents[arr_size][arr_size]),
		n_buf(boost::extents[arr_size][arr_size]),


		zfull_buf(boost::extents[arr_size][arr_size][max_roots]),
		jsel_buf(boost::extents[arr_size][arr_size]),
		nsel_buf(boost::extents[arr_size][arr_size]),


		is_silhouette_mat(boost::extents[arr_size][arr_size]),

		shape_mat_buf(boost::extents[arr_size][arr_size][4]),

		vbase_1_buf(boost::extents[arr_size][arr_size][3]),
		vbase_2_buf(boost::extents[arr_size][arr_size][3]),

		v1_buf(boost::extents[arr_size][arr_size][2]),
		v2_buf(boost::extents[arr_size][arr_size][2]),

		l1_buf(boost::extents[arr_size][arr_size]),
		l2_buf(boost::extents[arr_size][arr_size]),


		pz_buf(&z_buf),
		pn_buf(&n_buf),


		pzfull_buf(& zfull_buf),
		pjsel_buf(&jsel_buf),
		pnsel_buf(&nsel_buf),

		pis_silhouette_mat(&is_silhouette_mat),

		pshape_mat_buf(& shape_mat_buf),

		pvbase_1_buf(& vbase_1_buf),
		pvbase_2_buf(& vbase_2_buf),

		pv1_buf(& v1_buf),
		pv2_buf(& v2_buf),

		pl1_buf(& l1_buf),
		pl2_buf(& l2_buf)



	{};

	int eval_parms_poly(Poly5 poly_in, Poly3 & poly_out, double *parm);


	int init_f5(QString f5str);

	int init_f3_diff(Poly3 & f3);


	void eval_poly_f3(double x, double y, double z, double & resval);

	int eval_poly_poly_f3(double x, double y, double z, double & f, double* fnormal, double* fhessian);


	int root_list_poly_point_normal(double x, double y, double z,
		double nx, double ny, double nz, double* root_list, int & root_list_len);


	void fill_arrays(int xmax, int ymax, int shade_type);

	void calc_silhouette(int xmax, int ymax, double & k_gauss_min, double & k_gauss_max);


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



private:

	int get_z_intersect_poly(double x, double y, double  & z, double  & n_z, bool & disc_zero, double* full_z_list,
			double* shape_mat, double* base_v1, double* base_v2, double & l1, double & l2, int & jsel, int & nsel,
			int shade_type);



};

int rotate_mat(double phi, double theta, double psi);


extern ProcPoly pp;

extern Scale xrast_to_x, yrast_to_y, x_to_xrast, y_to_yrast;




extern const int win_size;

extern Point2DList silhouette_pointl;



//#define EPS 1e-12
#define M_INF -1e38




#endif /* PROCPOLY_H_ */
