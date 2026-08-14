#ifndef UTILS_H
#define UTILS_H

#include "variables.h"

inline int idx(int x, int t, int mu) {
    //x ranges from 0 to width_x+1
    //t ranges from 0 to width_t+1
    //The physical volume runs from 1 to width_x (or width_t)
    //mu = 0, 1
    return ((x*(mpi::width_t+2) + t)*2 + mu);
}


/*
Generate a random U(1) variable
*/
inline c_double RandomU1() {
	//Random angle in (0,2*pi) with uniform distribution 
	double cociente = ((double) rand() / (RAND_MAX));
    double theta = 2.0*pi * cociente;
	c_double z(cos(theta), sin(theta));
	return z;
}


#endif