
#ifndef __timing_h
#define __timing_h

#include <sys/time.h>

class Timing {

public:
	Timing() {};

	void start();
	void stop();

	double elapsed();

private:

	timeval begin;
    timeval current;

};






#endif
