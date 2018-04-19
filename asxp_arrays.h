#ifndef __asxp_arrays_h
#define __asxp_arrays_h

//#ifndef Q_MOC_RUN
#include <boost/multi_array.hpp>
//#endif

// Create a 3D array that is 3 x 4 x 2
//typedef boost::multi_array<double, 3> array_type;
//typedef array_type::index Aindex;
//array_type A(boost::extents[3][4][2]);

typedef boost::multi_array<double, 2> Array2d_double;

typedef boost::multi_array<double, 3> Array3d_double;

typedef boost::multi_array<int, 2> Array2d_int;

typedef boost::multi_array<int, 3> Array3d_int;

typedef boost::multi_array<bool, 2> Array2d_bool;













#endif
