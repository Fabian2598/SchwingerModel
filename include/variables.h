#ifndef VARIABLES_H_INCLUDED
#define VARIABLES_H_INCLUDED
#include "config.h"
#include <vector>
#include <complex>

extern double pi;
typedef std::complex<double> c_double;

struct spinor{
    c_double mu0[LV::Ntot];
    c_double mu1[LV::Ntot];
};


//------------Lattice parameters--------------//
namespace LV {
    //Lattice dimensions//
    constexpr int Nx= NS; //We extract this value from config.h
    constexpr int Nt = NT; //We extract this value from config.h
    constexpr int Ntot = Nx*Nt; //Total number of lattice points
}

namespace CG{
    extern int max_iter; //Maximum number of iterations for the conjugate gradient method
    extern double tol; //Tolerance for convergence
}

int Coords(const int& x, const int& t);
extern int LeftPB[LV::Ntot*2]; 
extern int RightPB[LV::Ntot*2]; 
extern c_double SignL[LV::Ntot*2]; 
extern c_double SignR[LV::Ntot*2]; 
extern int x_1_t1[LV::Ntot]; 
extern int x1_t_1[LV::Ntot]; 


//Memory preallocation
extern std::vector<std::vector<c_double>> TEMP;
extern std::vector<std::vector<c_double>> DTEMP;



#endif 