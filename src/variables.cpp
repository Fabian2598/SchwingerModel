#include "variables.h"

double pi=3.14159265359;

MPI_Datatype sub_block_type;
MPI_Datatype sub_block_resized;

/*
	Vectorized lattice coords.*/
int Coords(const int& x, const int& t){
	return x*mpi::width_t + t;
}

namespace mpi{
    int rank = 0;
    int size = 1; 
    int maxSize = LV::Ntot; //Default value, will be updated in main
    int ranks_x = 1;
    int ranks_t = 1;
    int width_x = LV::Nx;
    int width_t = LV::Nt;
    int rank2d = 0; //linearize rank
    int coords[2] = {0,0}; //rank coords
    int top = 0; 
    int bot = 0; 
    int right = 0; 
    int left = 0;
    MPI_Comm cart_comm;
}

namespace CG{
	int max_iter = 10000; //Maximum number of iterations for the conjugate gradient method
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

spinor TopRow(LV::Nt); //Should be width_t
spinor BottomRow(LV::Nt);
spinor RightCol(LV::Nx);
spinor LeftCol(LV::Nx); //Shoul be width_x
//Buffers for MPI communication
