#ifndef UTILS_H
#define UTILS_H

#include "variables.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <fstream>
#include <sstream>
#include <type_traits>

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
    Executes f(n) for every site n = x*(width_t+2)+t of the local lattice (halo excluded).
    With CUDA this is a GPU kernel with one thread per lattice site, check backend.h for the rules
    the lambda f has to follow. Usage:
        for_each_site([=] HOST_DEVICE (int n){ ... });
*/
template <typename F>
inline void for_each_site(F f){
    backend::for_sites(mpi::width_x, mpi::width_t, f);
}

/*
    Sum of f(n) over all the lattice sites of all the ranks, T = double or c_double.
    f can also update fields at site n, which fuses vector updates and reductions in one kernel.
        double s = lattice_sum<double>([=] HOST_DEVICE (int n){ return ...; });
*/
template <typename T, typename F>
inline T lattice_sum(F f){
    T local_sum = backend::sum_sites<T>(mpi::width_x, mpi::width_t, f);
    T sum;
    MPI_Allreduce(&local_sum, &sum, 1, std::is_same_v<T, double> ? MPI_DOUBLE : MPI_DOUBLE_COMPLEX,
                  MPI_SUM, mpi::cart_comm);
    return sum;
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
c_double dot(const spinor& X, const spinor& Y);


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