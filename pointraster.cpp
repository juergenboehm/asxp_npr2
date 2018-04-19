
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "pointraster.h"

using namespace std;

Pointraster::Pointraster(double x1, double y1, double x2, double y2, double hxa, double hya, double da):
	xul(x1), yul(y1), xdr(x2), ydr(y2), hx(hxa), hy(hya), d_main(da)
{
}

void Pointraster::set_h(double hxa, double hya)
{
	hx = hxa;
	hy = hya;

}

void Pointraster::set_d(double da)
{
	d_main = da;
}


void Pointraster::set_xy_size()
{
	double xmax1 = (xdr - xul) / hx;
	double ymax1 = (ydr - yul) / hy;

	int maxx = floor(xmax1);
	int maxy = floor(ymax1);

	xsize = maxx + 1;
	ysize = maxy + 1;

}

void Pointraster::init()
{
	set_xy_size();

	plist_arr_size = xsize * ysize;

	plist_arr.clear();

	plist_arr.resize(plist_arr_size);

}


Point2DList & Pointraster::get_arr(int x, int y)
{
	return plist_arr[x * xsize + y];
}

void Pointraster::get_xy_ind(double x, double y, int & xind, int & yind)
{
	int xi = (int)(floor((x - xul) / hx));
	int yi = (int)(floor((y - yul) / hy));

	//assert(xi >= 0 && xi < xsize);
	//assert(yi >= 0 && yi < ysize);

	xind = xi;
	yind = yi;
}

void Pointraster::enter_xy(double x, double y, int & xind, int & yind)
{
	get_xy_ind(x, y, xind, yind);

	if (xind >= 0 && xind < xsize && yind >= 0 && yind < ysize) {
		Point2DList & xy_list = get_arr(xind, yind);

		xy_list.push_back(Point2D(x, y));
	}


}

int Pointraster::count_points()
{
	int cnt = 0;
	for(int xi = 0; xi < xsize; ++xi ) {
		for(int yi = 0; yi < ysize; ++yi) {
			cnt += get_arr(xi, yi).size();
		}
	}
	return cnt;
}

bool Pointraster::check_d(double x, double y, double d, Point2DList & xy_list)
{
	for(Point2DListIterator it = xy_list.begin(); it != xy_list.end(); ++it) {
		double xit = it->x;
		double yit = it->y;

		xit -= x;
		yit -= y;

		double dit = sqrt(xit * xit + yit * yit);

		if (dit <= d) {
			return true;
		}

	}
	return false;
}

// returns true if distance of (x,y) to a point in raster is <= d
bool Pointraster::check_xy_d(double x, double y, double d)
{
	int xi, yi;

	get_xy_ind(x, y, xi, yi);

	int ximin = ::max(xi - 1, 0);
	int ximax = ::min(xi + 1, xsize - 1);

	int yimin = ::max(yi - 1, 0);
	int yimax = ::min(yi + 1, ysize - 1);


	for(int xi1 = ximin; xi1 <= ximax; ++xi1) {
		for(int yi1 = yimin; yi1 <= yimax; ++yi1) {

			Point2DList & xy_list = get_arr(xi1, yi1);
			bool is_near = check_d(x, y, d, xy_list);
			if (is_near) {
				return true;
			}

		}
	}
	return false;

}

bool Pointraster::check_xy(double x, double y)
{
	return check_xy_d(x, y, d_main);
}













