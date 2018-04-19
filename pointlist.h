#ifndef __pointlist_h
#define __pointlist_h

#include <list>


class Point2D {

public:

	Point2D(double xa, double ya);

	double x, y;

};

class PointClassifier {

public:

	virtual ~PointClassifier() {};
	virtual int operator()(double x, double y) = 0;
};


typedef std::list<Point2D> Point2DList;
typedef Point2DList::iterator Point2DListIterator;
typedef Point2DList::const_iterator Point2DListConstIterator;

int linpol_pointlist(Point2DList & list_in, Point2DList & list_out);

int distance_d_list(double d, Point2DList & list_in, Point2DList & list_out);

int select_points(PointClassifier & pc, Point2DList & list_in, Point2DList & list_out);

int remove_points(PointClassifier & pc, Point2DList & list_in, Point2DList & list_out);
















#endif

