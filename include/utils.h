#ifndef UTILS_H
#define UTILS_H

#include "variables.h"
#include <algorithm>
#include <cmath>
#include <random>

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


/*
    dot product between two spinors of the form psi[ntot][2]
    A.B = sum_i A_i conj(B_i) 
*/
inline c_double dot(const spinor_v2& x, const spinor_v2& y) {
    c_double local_z = 0;
    //reduction over all lattice points and spin components
    for (int n = 0; n < mpi::maxSize; n++) {
        local_z += x.val[2*n]   * std::conj(y.val[2*n]);
        local_z += x.val[2*n+1] * std::conj(y.val[2*n+1]);
    }
    c_double z;
    MPI_Allreduce(&local_z, &z, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi::cart_comm);
    return z;
}


//mean of a vector
template <typename T>
double mean(std::vector<T> x){ 
    double prom = 0;
    for (T i : x) {
        prom += i*1.0;
    }   
    prom = prom / x.size();
    return prom;
}

//random double number in the inteval [a,b] a = min, b = max
inline double rand_range(double a, double b){
    double cociente = ((double) rand() / (RAND_MAX));
    double x = (b-a) * cociente + a;
    return x;
}

//----------Jackknife---------//
std::vector<double> samples_mean(std::vector<double> dat, int bin); 
double Jackknife_error(std::vector<double> dat, int bin); 
double Jackknife(std::vector<double> dat, std::vector<int> bins); 

//---------------Linspace (similar to python)----------------------//
template <typename T>
std::vector<double> linspace(T min, T max, int n) {
    std::vector<double> linspace;
    double h = (1.0*max - 1.0*min) / (n - 1);
    for (int i = 0; i < n; ++i) {
        linspace.insert(linspace.begin() + i, min*1.0 + i * h); 
    }
    return linspace;
}
 
#endif