
#include <sstream>

#include <stdlib.h>
#include <string.h>

#include "poly.h"

Poly5 calcpoly(5);


template<>
int Poly<double>::Term::sprintcoef(char* buf, bool with_sign) const
{
	double cf = with_sign ? coef : fabs(coef);
	return sprintf(buf, "%.2f", cf);
}


template<>
std::string Poly<double>::Term::to_str_coef(bool with_sign) const
{
	double cf = with_sign ? coef: fabs(coef);

	std::stringstream ss;

	ss << cf;

	return ss.str();
}


template<>
double Poly<double>::Term::powx(double val, int x) const
{
	return ::pow(val, x);
}
