#include "variables.h"

double pi=3.14159265359;

/*
	Vectorized lattice coords.*/
int Coords(const int& x, const int& t){
	return x*mpi::width + t;
}

namespace mpi{
    int rank = 0;
    int size = 1; 
    int maxSize = LV::Ntot; //Default value, will be updated in main
    int width = LV::Nx;
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
    using namespace mpi;
    LeftPB  = new int[maxSize * 2];
    RightPB = new int[maxSize * 2];
    SignL   = new c_double[maxSize * 2];
    SignR   = new c_double[maxSize * 2];
    x_1_t1  = new int[maxSize];
    x1_t_1  = new int[maxSize];
}

void free_lattice_arrays() {
    delete[] LeftPB;
    delete[] RightPB;
    delete[] SignL;
    delete[] SignR;
    delete[] x_1_t1;
    delete[] x1_t_1;
}

//Memory preallocation I need to give them the right dimension mpi::maxSize, but I only know it in main after MPI_Init ...
spinor DTEMP;
spinor TEMP; 


//Buffers for MPI communication
spinor TopRow(LV::Nt);
spinor BottomRow(LV::Nt);