

#include "timing.h"


void Timing::start()
{
	::gettimeofday(&begin, 0);
}

void Timing::stop()
{
    ::gettimeofday(&current, (struct timezone*) 0);

}

double Timing::elapsed()
{
	double used = (current.tv_sec - begin.tv_sec) + ((current.tv_usec
	            - begin.tv_usec) / 1000000.0F);
	return used;
}
