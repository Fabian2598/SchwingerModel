#ifndef VARIABLES_H_INCLUDED
#define VARIABLES_H_INCLUDED
#include "config.h"
#include <vector>
#include <complex>

extern double pi;
typedef std::complex<double> c_double;


//------------Lattice parameters--------------//
namespace LV {
    //Lattice dimensions//
    constexpr int Nx= NS; //We extract this value from config.h
    constexpr int Nt = NT; //We extract this value from config.h
    constexpr int Ntot = Nx*Nt; //Total number of lattice points
}


void Coordinates(); //Vectorized coordinates. Coords[x][t]. Computed only once
extern std::vector<std::vector<int>>Coords; //coordinates x * Nt + t

//We store the neighbor coordinates to prevent recomputing them each time D is called.
extern std::vector<std::vector<int>>RightPB; //Right periodic boundary
extern std::vector<std::vector<int>>LeftPB; //Left periodic boundary
extern std::vector<std::vector<c_double>>SignR; //Right fermionic boundary
extern std::vector<std::vector<c_double>>SignL; //Left fermionic boundary
extern std::vector<int>x_1_t1; //x-\hat{x}+\hat{t}
extern std::vector<int>x1_t_1; //x+\hat{x}-\hat{t}
extern std::vector<std::vector<c_double>>Y;


#endif 