#ifndef JACKKNIFE_H_INCLUDED
#define JACKKNIFE_H_INCLUDED
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

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

template <typename T>
double standard_deviation(const std::vector<T>& data) {
    if (data.empty()) return 0.0;
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    double sq_sum = 0.0;
    for (const auto& val : data) {
        sq_sum += (static_cast<double>(val) - mean) * (static_cast<double>(val) - mean);
    }
    return std::sqrt(sq_sum / data.size());
}

//----------Jackknife---------//
std::vector<double> samples_mean(std::vector<double> dat, int bin); 
double Jackknife_error(std::vector<double> dat, int bin); 
double Jackknife(std::vector<double> dat, std::vector<int> bins); 


#endif