#ifndef __pointraster_h
#define __pointraster_h

#include <vector>

#include "pointlist.h"


class Pointraster {

public:

	Pointraster(double x1, double y1, double x2, double y2, double hxa, double hya, double d);

	void init();

	void set_xy_size();

	void set_h(double hxa, double hya);

	void set_d(double da);

	Point2DList & get_arr(int x, int y);

	int count_points();

	void get_xy_ind(double x, double y, int & xind, int & yind);

	void enter_xy(double x, double y, int & xind, int & yind);

	bool check_d(double x, double y, double d, Point2DList & xy_list);

	bool check_xy_d(double x, double y, double d);

	bool check_xy(double x, double y);


	double xul, yul;
	double xdr, ydr;

	double hx, hy;

	double d_main;

	int xsize, ysize;

	int plist_arr_size;


	std::vector<Point2DList> plist_arr;

};

#endif


