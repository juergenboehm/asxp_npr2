
#include <iostream>

#include <QString>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Stream_lines_2.h>
#include <CGAL/Euler_integrator_2.h>
#include <CGAL/Runge_kutta_integrator_2.h>

//#define TRY 0

#ifndef TRY
#include <CGAL/Regular_grid_2.h>
#endif

#if TRY
#include "Regular_grid_xp_2.h"
#endif

#include "asxp_arrays.h"

#include "streamplot.h"
#include "pointraster.h"
#include "streamline.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

#ifndef TRY
typedef CGAL::Regular_grid_2<K> Regular_grid;
#endif

#if TRY
typedef CGAL::Regular_grid_xp_2<K> Regular_grid;
#endif

typedef CGAL::Runge_kutta_integrator_2<Regular_grid> Runge_kutta_integrator;
typedef CGAL::Euler_integrator_2<Regular_grid> Euler_integrator;
typedef CGAL::Stream_lines_2<Regular_grid, Runge_kutta_integrator> Stream_lines;
typedef CGAL::Stream_lines_2<Regular_grid, Runge_kutta_integrator>::Stream_line_iterator_2 Stream_line_iterator;
typedef CGAL::Stream_lines_2<Regular_grid, Runge_kutta_integrator>::Point_iterator_2 Point_iterator;
typedef CGAL::Stream_lines_2<Regular_grid, Runge_kutta_integrator>::Point_2 Point_2;
typedef CGAL::Stream_lines_2<Regular_grid, Runge_kutta_integrator>::Vector_2 Vector;


typedef K::Iso_rectangle_2 Iso_rectangle_2;

#define Q_SQ(x) ((x) * (x))
#define Q_NORM(x,y) (Q_SQ(x) + Q_SQ(y))
#define Q_SKALPROD(x1, y1, x2, y2) ((x1)*(x2) + (y1)*(y2))

#define TO_GRAD(x) ((x)/(2 * 3.14159265) * 360)
#define SIGN(x) ((x) > 0 ? 1 : (((x) == 0) ? 0 : -1))



using namespace std;

OrPointClassifier::OrPointClassifier(PointClassifier* pcf_aa, PointClassifier* pcf_ba):
					pcf_a(pcf_aa), pcf_b(pcf_ba)
{

}

int OrPointClassifier::operator()(double x, double y)
{
	return int(((bool) (*pcf_a)(x,y)) || ((bool) (*pcf_b)(x,y)));
}


PrastPointClassifier::PrastPointClassifier(Pointraster & prasta): prast(prasta)
{
}

int PrastPointClassifier::operator()(double x, double y)
{
	return prast.check_xy(x, y);
}



static void normalize_vec(double & x, double & y)
{
	double norm = sqrt(Q_NORM(x,y));

	if (norm > 0) {
		x /= norm;
		y /= norm;
	}
}

const double d_ratio = 0.5;

Streamplot::Streamplot(Array3d_double* vf1a, Array3d_double* vf2a, int xmaxa, int ymaxa, double da,
				Scale xrast_to_xa, Scale yrast_to_ya,
				Scale x_to_xrasta, Scale y_to_yrasta,
				PointClassifier *is_boundarya):
				vf1(vf1a), vf2(vf2a),
					xmax(xmaxa), ymax(ymaxa), d_main(da),
					xrast_to_x(xrast_to_xa), yrast_to_y(yrast_to_ya),
					x_to_xrast(x_to_xrasta), y_to_yrast(y_to_yrasta), h(1e-3),
					is_boundary(is_boundarya),
					prast(0, 0, xmaxa, ymaxa, d_ratio * da, d_ratio * da, d_ratio * da), mode(0)
{
}

void Streamplot::init_queue()
{
	streaml_queue.clear();
}

void Streamplot::compute_d_fringe(Streamline & curr_streaml, double d, Pointraster & prast, Point2DList & curr_fringe)
{
	Point2DList aux_list;
	PrastPointClassifier is_near_d(prast);

	distance_d_list(d,	curr_streaml.line_samples, aux_list);
	cout << "before remove points len = " << aux_list.size() << endl;
	remove_points(is_near_d, aux_list, curr_fringe);
	cout << "after remove points len = " << curr_fringe.size() << endl;

}

void Streamplot::prast_update(Point2DList & plist)
{
	cout << "prast_update: old_size = " << prast.count_points() << endl;
	for(Point2DListIterator pit = plist.begin(); pit != plist.end(); ++pit) {
		double xval = pit->x;
		double yval = pit->y;

		int xind, yind;

		prast.enter_xy(xval, yval, xind, yind);

	}
	cout << "prast_update: new_size = " << prast.count_points() << endl;

}

Array2d_bool is_set_already(boost::extents[1600 + 1][1600 + 1]);

const int minimal_line_len = 20;

void Streamplot::compute_stream_field(double xseed, double yseed, bool init_comp)
{

	for(int x = 0; x < 1601; ++x) {
		for(int y = 0; y < 1601; ++y) {
			is_set_already[x][y] = false;
		}
	}

	if (init_comp) {
		prast.init();
	}

	PrastPointClassifier is_near_d(prast);
	OrPointClassifier is_stop_point(is_boundary, &is_near_d);

	double xseed_rast = x_to_xrast(xseed);
	double yseed_rast = y_to_yrast(yseed);

	double dd = d_main * d_ratio;

	if (is_near_d(xseed_rast + dd, yseed_rast - dd) ||
		is_near_d(xseed_rast + dd, yseed_rast)||
		is_near_d(xseed_rast + dd, yseed_rast + dd)||

		is_near_d(xseed_rast, yseed_rast - dd) ||
		is_near_d(xseed_rast, yseed_rast)||
		is_near_d(xseed_rast, yseed_rast + dd)||

		is_near_d(xseed_rast - dd, yseed_rast - dd) ||
		is_near_d(xseed_rast - dd, yseed_rast)||
		is_near_d(xseed_rast - dd, yseed_rast + dd)) {

		return;
	}

	//Streamline curr_streaml(vf1, vf2, xmax, ymax, x_to_xrast, y_to_yrast, is_boundary);
	Streamline curr_streaml;

	Integrator curr_streaml_integr(vf1, vf2, xmax, ymax, x_to_xrast, y_to_yrast, is_boundary, curr_streaml);

	curr_streaml_integr.set_mode(mode);

	curr_streaml_integr.set_aux(&is_set_already);

	curr_streaml_integr.set_h(1e-4);

	curr_streaml_integr.integrate_from(xseed, yseed);

	//streaml_queue.push_back(curr_streaml);
	streaml_res.push_back(curr_streaml);
	prast_update(curr_streaml.line_samples);

	bool finished = false;

	do {
		Point2DList curr_fringe_d;
		Point2DList new_fringe_d;

		compute_d_fringe(curr_streaml, d_main, prast, curr_fringe_d);

		int cnt_seeds = 0;

		cout << "d_fringe_size start = " << curr_fringe_d.size() << endl;

#ifdef SHOW_FRINGE
		cin.get();

		debug_points = curr_fringe_d;

		break;
#endif

		while (curr_fringe_d.size() > 0) {
			Point2D new_sp = curr_fringe_d.front();
			curr_fringe_d.pop_front();

			if (is_stop_point(new_sp.x, new_sp.y)) {
				continue;
			}

			//Streamline csl(vf1, vf2, xmax, ymax, x_to_xrast, y_to_yrast, &is_stop_point);

			Streamline csl;
			Integrator csl_int(vf1, vf2, xmax, ymax, x_to_xrast, y_to_yrast, &is_stop_point, csl);

			csl_int.set_mode(mode);

			csl_int.set_aux(&is_set_already);

			csl_int.set_h(1e-4);

			csl_int.integrate_from(xrast_to_x(new_sp.x), yrast_to_y(new_sp.y));

			if (csl.line_samples.size() < minimal_line_len) {
				continue;
			}

			cout << "line_samples len = " << csl.line_samples.size() << endl;

			prast_update(csl.line_samples);

			streaml_queue.push_back(csl);

			streaml_res.push_back(csl);

			new_fringe_d.clear();
			remove_points(is_near_d, curr_fringe_d, new_fringe_d);
			curr_fringe_d = new_fringe_d;

			cout << "d_fringe_size iter = " << curr_fringe_d.size() << endl;

			++cnt_seeds;
		}

		// cin.get();

		if (streaml_queue.size() == 0) {
			finished = true;
		} else {
			curr_streaml = streaml_queue.front();
			streaml_queue.pop_front();
		}

	} while (!finished);

}

void Streamplot::set_mode(int modea) {

	mode = modea;

}

void Streamplot::compute_stream_field_CGAL() {

	double integrating = 1;

	Runge_kutta_integrator * runge_kutta_integrator = new Runge_kutta_integrator(integrating);

	//euler_integrator = new Euler_integrator(integrating);

	double xw_min = xrast_to_x(0);
	double xw_max = xrast_to_x(xmax);

	double yw_min = yrast_to_y(0);
	double yw_max = yrast_to_y(ymax);

	double ix_size, iy_size;

	ix_size = xmax;
	iy_size = ymax;

	double scale_i = 4;

	Regular_grid * regular_grid = new Regular_grid(xmax, ymax, ix_size * scale_i, iy_size * scale_i);

#if TRY
	regular_grid->set_point_classifier(is_boundary);
#endif

	/*fill the grid with the appropriate values*/

	for (int i = 0; i < xmax; i++) {
		for (int j = 0; j < ymax; j++) {

			double vec_x, vec_y;

			if (mode == 1) {

				vec_x = (*vf1)[i][j][0];
				vec_y = (*vf1)[i][j][1];

			} else if (mode == 2) {

				vec_x = (*vf2)[i][j][0];
				vec_y = (*vf2)[i][j][1];

			}

			normalize_vec(vec_x, vec_y);

			regular_grid->set_field(i, j, Vector(vec_x, vec_y));
		}
	}

	std::cout << "enter_generate..." << std::endl;

	double density = d_main;
	double ratio = 1.2;
	int sampling = 1;

	Stream_lines * stream_lines = new Stream_lines(*regular_grid,
			*runge_kutta_integrator, density, ratio, sampling);

	Stream_line_iterator sit = stream_lines->begin();
	Stream_line_iterator sit_end = stream_lines->end();
	int cnt = 0;

	QString fname;

	if (mode == 1) {
		fname = "intpoints1.dat";
	} else if (mode == 2) {
		fname = "intpoints2.dat";
	}

	ofstream intpointss(fname.toLatin1().data());

	while (sit != sit_end) {

		intpointss << "Streamline cnt = " << cnt << endl;

		// only dummy, no integration, just fill csl.line_samples
		Streamline csl;

		Point_iterator pit = sit->first;

		double xnew, ynew;
		double xold, yold;

		int xnew_i, ynew_i;
		int xold_i, yold_i;

		xnew = pit->x();
		ynew = pit->y();

		xnew_i = int(xnew/scale_i);
		ynew_i = int(ynew/scale_i);

		intpointss << xnew << " " << ynew << "\n";

		csl.line_samples.push_back(Point2D(xnew_i, ynew_i));

		++pit;

		while (pit != sit->second) {

			xold = xnew;
			yold = ynew;

			xold_i = xnew_i;
			yold_i = ynew_i;

			xnew = pit->x();
			ynew = pit->y();

			xnew_i = int(xnew/scale_i);
			ynew_i = int(ynew/scale_i);

			intpointss << xnew << " " << ynew << "\n";

			if (xnew_i != xold_i || ynew_i != yold_i) {

				csl.line_samples.push_back(Point2D(xnew_i, ynew_i));
			}

			++pit;

		}

		streaml_res.push_back(csl);

		++cnt;
		++sit;
	}


	delete stream_lines;
	delete regular_grid;
	delete runge_kutta_integrator;

}
