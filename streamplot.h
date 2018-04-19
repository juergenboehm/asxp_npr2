#ifndef __streamplot_h
#define __streamplot_h

#include "asxp_arrays.h"
#include "pointlist.h"
#include "pointraster.h"
#include "streamline.h"

#include <list>

class OrPointClassifier: public PointClassifier {

public:

	OrPointClassifier(PointClassifier* pcf_aa, PointClassifier* pcf_ba);
	virtual ~OrPointClassifier() {};

	virtual int operator()(double x, double y);

	PointClassifier* pcf_a;
	PointClassifier* pcf_b;
};



class PrastPointClassifier: public PointClassifier {

public:

	PrastPointClassifier(Pointraster & prasta);
	virtual ~PrastPointClassifier() {};

	virtual int operator()(double x, double y);

	Pointraster & prast;

};

class Streamplot {

public:

	static const double d_default = 10;

	Streamplot(Array3d_double* vf1a, Array3d_double* vf2a, int xmaxa, int ymaxa, double da,
				Scale xrast_to_xa, Scale yrast_to_ya, Scale x_to_xrasta, Scale y_to_yrasta,
				PointClassifier* is_boundarya);

	void set_mode(int modea);

	void init_queue();

	void compute_stream_field(double xseed, double yseed, bool init_comp);

	void compute_d_fringe(Streamline & curr_streaml, double d, Pointraster & prast, Point2DList & curr_fringe);

	void prast_update(Point2DList & plist);

	void compute_stream_field_CGAL();


	Array3d_double* vf1;
	Array3d_double* vf2;

	int xmax, ymax;
	double d_main;

	Scale xrast_to_x;
	Scale yrast_to_y;
	Scale x_to_xrast;
	Scale y_to_yrast;

	double h;

	PointClassifier* is_boundary;

	Pointraster prast;

	Point2DList debug_points;

	int mode;

	std::list<Streamline> streaml_queue;

	std::list<Streamline> streaml_res;

};




















#endif
