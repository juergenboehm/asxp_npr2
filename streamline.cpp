
#include <iostream>

#include <assert.h>

#include <math.h>
#include <float.h>

#include <list>

#include <stdlib.h>
#include <stdio.h>

#include "streamline.h"

#define Q_SQ(x) ((x) * (x))
#define Q_NORM(x,y) (Q_SQ(x) + Q_SQ(y))
#define Q_SKALPROD(x1, y1, x2, y2) ((x1)*(x2) + (y1)*(y2))

#define TO_GRAD(x) ((x)/(2 * 3.14159265) * 360)
#define SIGN(x) ((x) > 0 ? 1 : (((x) == 0) ? 0 : -1))


using namespace std;

Scale::Scale(): a(1), b(0)
{
}

Scale::Scale(double aa, double ba): a(aa), b(ba)
{
}

void Scale::set(double aa, double ba)
{
	a = aa;
	b = ba;
}

double Scale::operator()(double x)
{
	return a * x + b;
}


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



inline double bilin_interpolate(double a, double b, double c, double d, double sx, double sy) {

	double val_a = a * sy * (1 - sx);
	double val_b = b * sx * sy;
	double val_c = c * (1 - sx) * (1 - sy);
	double val_d = d * sx * (1 - sy);

	double val = val_a + val_b + val_c + val_d;

	return val;
}



static void normalize_vec(double & x, double & y)
{
	double norm = sqrt(Q_NORM(x,y));

	if (norm > 0) {
		x /= norm;
		y /= norm;
	}
}

void full_normalize_vec(double & x, double & y)
{
	const double eps = 1e-6;

	double x_old = x;

	if (fabs(x_old) > eps) {
		x /= x_old;
		y /= x_old;

		normalize_vec(x, y);

	} else {
		x = 0;
		y = (fabs(y) > eps) ? 1 : 0;
	}
}



// Integrator

Integrator::Integrator(Array3d_double* vf1a, Array3d_double* vf2a, int xmaxa, int ymaxa,
		Scale x_to_xrasta, Scale y_to_yrasta, PointClassifier* is_stop_point_a, Streamline & streamlinea):
		vf1(vf1a), vf2(vf2a), h(1e-3), xmax(xmaxa), ymax(ymaxa), x_to_xrast(x_to_xrasta), y_to_yrast(y_to_yrasta),
		is_stop_point(is_stop_point_a), streamline(streamlinea)
{

}

void Integrator::set_aux(Array2d_bool* is_set_already_a)
{
	is_set_already = is_set_already_a;
}


void Integrator::integrate_from(double x1, double y1)
{
	Point2DList forward_part;

	Point2DList & line_samples = streamline.line_samples;


	integrate_flow(x1, y1, h);
	forward_part = line_samples;

	line_samples.clear();

	integrate_flow(x1, y1, -h);

	Point2DListIterator pit = line_samples.begin();
	++pit;

	while (pit != line_samples.end()) {
		forward_part.push_front(*pit);
		++pit;
	}

	line_samples = forward_part;

	cout << "integrate_from: total points of integration = " << line_samples.size() << endl;
}

//
// (x, y) is in the model coordinate system
//
void Integrator::get_vector(Array3d_double* vf, double x, double y, int xmax, int ymax,
		double & vx, double & vy)
{

	double xrast = x_to_xrast(x);
	double yrast = y_to_yrast(y);

	int cx = floor(xrast);
	int cy = floor(yrast);

#if 0
	double sx = xrast - cx;
	double sy = yrast - cy;
#endif

	bool outofx = (cx < 0) || (cx > xmax - 1);
	bool outofy = (cy < 0) || (cy > ymax - 1);

	if (outofx || outofy) {
		vx = vy = 0;
	} else {

#if 0
		vx = bilin_interpolate((*vf)[cx][cy + 1][0], (*vf)[cx+1][cy+1][0], (*vf)[cx][cy][0], (*vf)[cx+1][cy][0], sx, sy);
		vy = bilin_interpolate((*vf)[cx][cy + 1][1], (*vf)[cx+1][cy+1][1], (*vf)[cx][cy][1], (*vf)[cx+1][cy][1], sx, sy);
#endif

#if 1
		vx = (*vf)[cx][cy][0];
		vy = (*vf)[cx][cy][1];
#endif

	}

}



void Integrator::get_cross_field(Array3d_double* vf1, Array3d_double* vf2,
								double x, double y, int xmax, int ymax,
									double & v1x, double & v1y,
									double & v2x, double & v2y )
{

	get_vector(vf1, x, y, xmax, ymax, v1x, v1y);
	get_vector(vf2, x, y, xmax, ymax, v2x, v2y);

}

void Integrator::get_vector_from_crossfield( Array3d_double * vf1, Array3d_double * vf2,
					double x, double y, double vfx_old, double vfy_old, int xmax, int ymax,
						double & vx, double & vy) {


	double v1x_new, v1y_new, v2x_new, v2y_new;

	if (mode == 0) {
		get_cross_field(vf1, vf2, x, y, xmax, ymax, v1x_new, v1y_new, v2x_new, v2y_new );

		double p11 = v1x_new * vfx_old + v1y_new * vfy_old;
		double p21 = v2x_new * vfx_old + v2y_new * vfy_old;

		if (fabs(p11) < 1e-3 && fabs(p21) < 1e-3) {
			vx = 0;
			vy = 0;
		} else if (fabs(p11) >= fabs(p21)) {
			vx = v1x_new * SIGN(p11);
			vy = v1y_new * SIGN(p11);
		} else {
			vx = v2x_new * SIGN(p21);
			vy = v2y_new * SIGN(p21);
		}
	} else if (mode == 1) {
		get_vector(vf1, x, y, xmax, ymax, vx, vy);
		normalize_vec(vx , vy);
	} else if (mode == 2) {
		get_vector(vf2, x, y, xmax, ymax, vx, vy);
		normalize_vec(vx, vy);
	}


}

void Integrator::set_mode(int modea) {

	mode = modea;

}

// (xstart, ystart) is in the model coordinate system

void Integrator::integrate_flow(double xstart, double ystart, double h) {

	double x0 = xstart, y0 = ystart;

	double x0rast = x_to_xrast(x0);
	double y0rast = y_to_yrast(y0);

	double x1, y1;

	double v1x;
	double v1y;

	//double v2x;
	//double v2y;

	cout << "enter integrate_flow..." << endl;

	get_vector(vf1, x0, y0, xmax, ymax, v1x, v1y);
	//get_vector(vf2, x0, y0, xmax, ymax, v2x, v2y);

	double vxc, vyc;

	vxc = v1x;
	vyc = v1y;

	normalize_vec( vxc, vyc );

	int len = 0;
	int cnt_points = 0;

	Point2DList & line_samples = streamline.line_samples;

	while (len < 10 * 1000 * 1000 && cnt_points <= 10000) {

		double vxc_new, vyc_new;

		get_vector_from_crossfield(vf1, vf2, x0, y0, vxc, vyc, xmax, ymax, vxc_new, vyc_new);

#if 1
		if (Q_NORM(vxc, vyc) > 0 && Q_NORM(vxc_new, vyc_new) > 0) {
			assert(fabs(Q_NORM(vxc, vyc) - 1) < 1e-6 && fabs(Q_NORM(vxc_new, vyc_new) - 1) < 1e-6);
		}
#endif

		double ptst = Q_SKALPROD(vxc, vyc, vxc_new, vyc_new);

		double angle = TO_GRAD(acos(ptst));

		if (fabs(angle) > 10) {
			//printf("ptst problem: phi = %f (%f, %f), (%f, %f) \n", angle, vxc, vyc, vxc_new, vyc_new);
		}

		vxc = vxc_new;
		vyc = vyc_new;

		x1 = x0 + vxc * h;
		y1 = y0 + vyc * h;


#if 0
		if (len % 1000 == 0 ) {
			cout << "step = " << len << endl;
		}
#endif

		double x1rast = x_to_xrast(x1);
		double y1rast = y_to_yrast(y1);

		if (int(x1rast) != int(x0rast) || int(y1rast) != int(y0rast)) {

			double x0rast1 = ::min(x0rast, (double(xmax) - 1.0));
			double y0rast1 = ::min(y0rast, (double(ymax) - 1.0));

			x0rast1 = ::max(x0rast1, 0.0);
			y0rast1 = ::max(y0rast1, 0.0);

			if ((*is_set_already)[int(x0rast1)][int(y0rast1)]) {
				break;
			}

			if ((*is_stop_point)(int(x0rast1), int(y0rast1))) {
				break;
			}

			++cnt_points;

			(*is_set_already)[int(x0rast1)][int(y0rast1)] = true;

			line_samples.push_back(Point2D(x0rast1, y0rast1));

			//xtestlist.push_back(x0rast1);
			//ytestlist.push_back(y0rast1);
		}

		const double eps = 1e-3 * fabs(h);

		if (fabs(x1 - x0) < eps && fabs(y1 - y0) < eps) {
			break;
		}

		x0 = x1;
		y0 = y1;


		x0rast = x1rast;
		y0rast = y1rast;

		++len;

	}
	for(Point2DListIterator pit = line_samples.begin(); pit != line_samples.end(); ++pit) {
		double x = pit->x;
		double y = pit->y;

		(*is_set_already)[int(x)][int(y)] = false;
	}

}


void Integrator::set_h(double ha)
{
	h = ha;
}
















