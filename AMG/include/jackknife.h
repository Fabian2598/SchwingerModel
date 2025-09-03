#ifndef JACKKNIFE_H_INCLUDED
#define JACKKNIFE_H_INCLUDED
#include <vector>
#include <algorithm>
#include <cmath>

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

//----------Jackknife---------//
std::vector<double> samples_mean(std::vector<double> dat, int bin); 
double Jackknife_error(std::vector<double> dat, int bin); 
double Jackknife(std::vector<double> dat, std::vector<int> bins); 


#endif