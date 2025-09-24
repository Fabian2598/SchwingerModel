#include "variables.h"

double pi=3.14159265359;

/*
	Vectorized lattice coords.*/
int Coords(const int& x, const int& t){
	return x*LV::Nt + t;
}

namespace CG{
	int max_iter = 1000; //Maximum number of iterations for the conjugate gradient method
	double tol = 1e-10; //Tolerance for convergence
}

int* LeftPB = nullptr;
int* RightPB = nullptr;
c_double* SignL = nullptr;
c_double* SignR = nullptr;
int* x_1_t1 = nullptr;
int* x1_t_1 = nullptr;

void allocate_lattice_arrays() {
    LeftPB  = new int[LV::Ntot * 2];
    RightPB = new int[LV::Ntot * 2];
    SignL   = new c_double[LV::Ntot * 2];
    SignR   = new c_double[LV::Ntot * 2];
    x_1_t1  = new int[LV::Ntot];
    x1_t_1  = new int[LV::Ntot];
}

void free_lattice_arrays() {
    delete[] LeftPB;
    delete[] RightPB;
    delete[] SignL;
    delete[] SignR;
    delete[] x_1_t1;
    delete[] x1_t_1;
}

//Memory preallocation
spinor DTEMP;
spinor TEMP; 