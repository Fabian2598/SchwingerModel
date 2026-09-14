#ifndef UTILS_H
#define UTILS_H

#include "variables.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <fstream>
#include <sstream>

/*
    Linearized index for a spinor
*/
inline int idx(int x, int t, int mu) {
    //x ranges from 0 to width_x+1
    //t ranges from 0 to width_t+1
    //The physical volume runs from 1 to width_x (or width_t)
    //mu = 0, 1
    return ((x*(mpi::width_t+2) + t)*2 + mu);
}
/*
    Linearized index for a vector
*/
inline int idx_vec(int x, int t){
    return (x*(mpi::width_t+2) + t);
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

// Converts a double like -0.4568 → "-04568"
inline std::string format(double val) {
    char buf[8];
    int sign    = (val < 0) ? -1 : 1;
    int digits  = static_cast<int>(std::round(std::abs(val) * 10000));
    std::snprintf(buf, sizeof(buf), "%s%05d",
                  (sign < 0 ? "-" : ""),
                  digits);
    return buf;
}


/*
    dot product between two spinors
    A.B = sum_i A_i conj(B_i) 
*/
inline c_double dot(const spinor& X, const spinor& Y) {
    c_double local_z = 0;
    //reduction over all lattice points and spin components
    int n;
    for(int x = 1; x<=mpi::width_x; x++){
        for(int t = 1; t<=mpi::width_t; t++){
            n = x*(mpi::width_t+2)+t;
            local_z += X.val[2*n]   * std::conj(Y.val[2*n]);
            local_z += X.val[2*n+1] * std::conj(Y.val[2*n+1]);
        }
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

void print_parameters();
 


#endif