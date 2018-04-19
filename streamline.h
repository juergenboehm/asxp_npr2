#ifndef __streamline_h
#define __streamline_h

#include <list>

#include "asxp_arrays.h"
#include "pointlist.h"


class Scale {

public:

	Scale();

	Scale(double aa, double ba);

	void set(double aa, double bb);

	double operator()(double x);

	double a, b;

};



class Streamline {

public:

	Streamline() {};

	Point2DList line_samples;

};


class Integrator {

public:

	Integrator(Array3d_double* vf1a, Array3d_double* vf2a, int xmaxa, int ymaxa, Scale x_to_xrasta, Scale y_to_yrasta,
				PointClassifier* is_stop_point_a, Streamline & streamlinea);

	void integrate_from(double x1, double y1);

	void set_aux(Array2d_bool* is_set_already_a);

	void set_h(double ha);

	void set_mode(int modea);

	void get_vector(Array3d_double* vf, double x, double y, int xmax, int ymax,
			double & vx, double & vy);


	void get_cross_field(Array3d_double* vf1, Array3d_double* vf2,
									double x, double y, int xmax, int ymax,
										double & v1x, double & v1y,
										double & v2x, double & v2y );

	void get_vector_from_crossfield( Array3d_double* vf1, Array3d_double* vf2,
						double x, double y, double vfx_old, double vfy_old, int xmax, int ymax,
							double & vx, double & vy);


	void integrate_flow(double xstart, double ystart, double h);


	Array3d_double* vf1;
	Array3d_double* vf2;

	double h;

	int xmax, ymax;

	Scale x_to_xrast;
	Scale y_to_yrast;


	PointClassifier* is_stop_point;


	Array2d_bool* is_set_already;

	int mode;

	Streamline & streamline;


};

























#endif
