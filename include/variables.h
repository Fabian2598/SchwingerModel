#ifndef VARIABLES_H_INCLUDED
#define VARIABLES_H_INCLUDED
#include "config.h"
#include <vector>
#include <complex>
#include "mpi.h"

extern double pi;
typedef std::complex<double> c_double;

namespace mpi{
    extern int rank;
    extern int size; 
    extern int maxSize;
    extern int ranks_x;
    extern int ranks_t;
    extern int width_x;
    extern int width_t;
    extern int rank2d; //Rank id in the 2D communicator
    extern int coords[2];
    extern int top; 
    extern int bot; 
    extern int right; 
    extern int left;
    extern MPI_Comm cart_comm;
}


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


struct spinor {
    c_double* mu0;
    c_double* mu1;
    int size;
    //Constructor
    spinor(int N = LV::Ntot) : size(N) {
        mu0 = new c_double[N]();
        mu1 = new c_double[N]();
    }

    // Copy constructor (deep copy)
    spinor(const spinor& other) : size(other.size) {
        mu0 = new c_double[size];
        mu1 = new c_double[size];
        std::copy(other.mu0, other.mu0 + size, mu0);
        std::copy(other.mu1, other.mu1 + size, mu1);
    }

    // Assignment operator (deep copy)
    spinor& operator=(const spinor& other) {
        if (this != &other) {
            if (size != other.size) {
                delete[] mu0;
                delete[] mu1;
                size = other.size;
                mu0 = new c_double[size];
                mu1 = new c_double[size];
            }
            std::copy(other.mu0, other.mu0 + size, mu0);
            std::copy(other.mu1, other.mu1 + size, mu1);
        }
        return *this;
    }

    // Destructor
    ~spinor() {
        delete[] mu0;
        delete[] mu1;
    }

    inline void clearBuffer(){
        for(int n = 0; n<size; n++){
            mu0[n] = 0;
            mu1[n] = 0;
        }
    }
};

struct re_field {
    double* mu0;
    double* mu1;
    int size;
    //Constructor
    re_field(int N = LV::Ntot) : size(N) {
        mu0 = new double[N]();
        mu1 = new double[N]();
    }

    // Copy constructor (deep copy)
    re_field(const re_field& other) : size(other.size) {
        mu0 = new double[size];
        mu1 = new double[size];
        std::copy(other.mu0, other.mu0 + size, mu0);
        std::copy(other.mu1, other.mu1 + size, mu1);
    }

    // Assignment operator (deep copy)
    re_field& operator=(const re_field& other) {
        if (this != &other) {
            if (size != other.size) {
                delete[] mu0;
                delete[] mu1;
                size = other.size;
                mu0 = new double[size];
                mu1 = new double[size];
            }
            std::copy(other.mu0, other.mu0 + size, mu0);
            std::copy(other.mu1, other.mu1 + size, mu1);
        }
        return *this;
    }

    // Destructor
    ~re_field() {
        delete[] mu0;
        delete[] mu1;
    }
};

typedef spinor c_matrix;

int Coords(const int& x, const int& t);
extern int* LeftPB;
extern int* RightPB;
extern c_double* SignL;
extern c_double* SignR;
extern int* x_1_t1;
extern int* x1_t_1;

void allocate_lattice_arrays();
void free_lattice_arrays();


//Memory preallocation
extern spinor DTEMP;
extern spinor TEMP; 

//Buffers for MPI communication
extern spinor TopRow;
extern spinor BottomRow;
extern spinor RightCol;
extern spinor LeftCol; //Shoul be width_x


/*
	Modulo operation
*/
inline int mod(int a, int b) {
	int r = a % b;
	return r < 0 ? r + b : r;
}






#endif 