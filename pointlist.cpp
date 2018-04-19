
#include <iostream>

#include <math.h>

#include "clipper.hpp"

#include "pointlist.h"

using namespace std;


Point2D::Point2D(double xa, double ya): x(xa), y(ya)
{
}


int linpol_pointlist(Point2DList & list_in, Point2DList & list_out)
{

	return 0;
}

double distance(double x1, double y1, double x2, double y2)
{
	x1 -= x2;
	y1 -= y2;

	return sqrt(x1 * x1 + y1 * y1);
}

static void line(int x0, int y0, int x1, int y1, Point2DList & list_add, double scale) {

  int dx = abs(x1-x0), sx = x0<x1 ? 1 : -1;
  int dy = abs(y1-y0), sy = y0<y1 ? 1 : -1;
  int err = (dx>dy ? dx : -dy)/2, e2;

  for(;;){

    list_add.push_back(Point2D(x0/scale,y0/scale));

    if (x0==x1 && y0==y1) break;
    e2 = err;
    if (e2 >-dx) { err -= dy; x0 += sx; }
    if (e2 < dy) { err += dx; y0 += sy; }
  }
}

int distance_d_list(double d, Point2DList & list_in, Point2DList & list_out)
{

	const double scale = 10;

	ClipperLib::ClipperOffset clip_offset;

	ClipperLib::Path list_in_path;

	if (list_in.size() == 0) {
		return 0;
	}

	Point2DListIterator pit = list_in.begin();

	double x0 = pit->x;
	double y0 = pit->y;

	for(pit = ++pit; pit != list_in.end(); ++pit) {
		double x1 = pit->x;
		double y1 = pit->y;

		Point2DList aux_list;

		line(x0 * scale, y0 * scale, x1 * scale, y1 * scale, aux_list, 1);

		for(Point2DListIterator qit = aux_list.begin(); qit != aux_list.end(); ++qit) {

			list_in_path.push_back(ClipperLib::IntPoint(qit->x, qit->y));
		}

		x0 = x1;
		y0 = y1;

	}

	cout << "list_in_path.size = " << list_in_path.size() << endl;

	ClipperLib::Paths list_out_paths;

	clip_offset.AddPath(list_in_path, ClipperLib::jtRound, ClipperLib::etOpenRound);

	clip_offset.Execute(list_out_paths, d * scale);

	cout << "list_out_paths.size = " << list_out_paths.size() << endl;

	for(ClipperLib::Paths::iterator plit = list_out_paths.begin();
			plit != list_out_paths.end(); ++plit ){

		ClipperLib::Path::iterator pit = plit->begin();

		int X0, Y0;
		int X1, Y1;

		X0 = pit->X;
		Y0 = pit->Y;

		for(pit = ++pit; pit != plit->end(); ++pit) {
			X1 = pit->X;
			Y1 = pit->Y;

			line(X0, Y0, X1, Y1, list_out, scale);

			X0 = X1;
			Y0 = Y1;
		}
	}

#if SLOW_ROUTINE
	for(Point2DListIterator pit = list_in.begin(); pit != list_in.end(); ++pit){

		double d_used = d + 1;

		double xc = pit->x;
		double yc = pit->y;

		double xmin = xc - d_used;
		double xmax = xc + d_used;

		int xmini = int(xmin);
		int xmaxi = int(xmax);

		for(int xi = xmini; xi <= xmaxi; ++xi){

			double xval = xi;
			double xs = xi - xc;
			double yval1 = yc + sqrt(d_used * d_used - xs * xs);
			double yval2 = yc -sqrt(d_used * d_used - xs * xs);

			bool ok1 = true;
			bool ok2 = true;

			for(Point2DListIterator qit = list_in.begin(); qit != list_in.end(); ++qit) {
				if (distance(xval, yval1, qit->x, qit->y) <= d) {
					ok1 = false;
					break;
				}
			}
			if (ok1) {
				list_out.push_back((Point2D(xval, yval1)));
			}
			for(Point2DListIterator qit = list_in.begin(); qit != list_in.end(); ++qit) {
				if (distance(xval, yval2, qit->x, qit->y) <= d) {
					ok2 = false;
					break;
				}
			}
			if (ok2) {
				list_out.push_back((Point2D(xval, yval2)));
			}
		}
	}

#endif

	return 0;
}

int select_points(PointClassifier & pc, Point2DList & list_in, Point2DList & list_out)
{
	for(Point2DListIterator pit = list_in.begin(); pit != list_in.end(); ++pit) {
		if (pc(pit->x, pit->y)) {
			list_out.push_back(*pit);
		}
	}
	return 0;
}


int remove_points(PointClassifier & pc, Point2DList & list_in, Point2DList & list_out)
{
	for(Point2DListIterator pit = list_in.begin(); pit != list_in.end(); ++pit) {
		if (!pc(pit->x, pit->y)) {
			list_out.push_back(*pit);
		}
	}
	return 0;
}
