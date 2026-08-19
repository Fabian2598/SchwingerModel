#ifndef VARIABLES_H_INCLUDED
#define VARIABLES_H_INCLUDED
#include "config.h"
#include <iostream>
#include <vector>
#include <complex>
#include <iomanip>
#include <string>
#include "mpi.h"

constexpr double pi = 3.14159265359;
typedef std::complex<double> c_double;
const c_double I_number(0,1); //imaginary number

namespace mpi{
    extern int rank;     //rank ID
    extern int size;     //number of ranks
    extern int maxSizeH; //number of variables (including halo and spin)
    extern int sitesH;   //number of lattice sites (including halo)
    extern int ranks_x;  //number of ranks on x 
    extern int ranks_t;  //number of ranks on t
    extern int width_x;  //number of lattice sites on x (no halo)
    extern int width_t;  //number of lattice sites on t (no halo)
    extern int rank2d;   //rank ID for the cartesian communicator
    extern int coords[2]; //rank 2D coordinates
    //Neighboring ranks of rank2d
    extern int top;
    extern int bot;
    extern int right;
    extern int left;
    extern int bot_left;
    extern int bot_right;
    extern int top_left;
    extern int top_right;

    extern MPI_Comm cart_comm; //cartesian communicator
    //Datatypes for reading/writing gauge confs and rhs
    extern MPI_Datatype column_type;
    extern MPI_Datatype global_conf_type;
    extern MPI_Datatype global_conf_resized;
    extern MPI_Datatype local_conf_type;
    extern MPI_Datatype local_conf_resized;

}

//Lattice dimensions
namespace LV {
    //Lattice dimensions//
    constexpr int Nx= NS; //We extract this value from config.h
    constexpr int Nt = NT; //We extract this value from config.h
    constexpr int Ntot = Nx*Nt; //Total number of lattice points
}

//CG parameters
namespace CG{
    extern int max_iter; //Maximum number of iterations for the conjugate gradient method
    extern double tol; //Tolerance for convergence
    extern bool print_convergence_message;
}

//BiCGstab parameters
namespace BiCG{
    extern int max_iter; //Maximum number of iterations
    extern double tol; //Tolerance for convergence
    extern bool print_convergence_message;
}

//Simulation parameters
namespace sim_params {
    extern double beta;
    extern double m0;
    extern int MD_steps;
    extern double trajectory_length;
    extern int Ntherm;
    extern int Nmeas;
    extern int Nsteps;
    extern std::string start_time_str;
}

//Flattened spinor
template <typename T>
struct spinor_template{
    T* val;
    int size;

    explicit spinor_template(int N = LV::Ntot) : size(N), val(new T[N]()) {}

    spinor_template(const spinor_template& other) : size(other.size), val(new T[size]) {
        std::copy(other.val, other.val + size, val);
    }

    spinor_template& operator=(const spinor_template& other) {
        if (this != &other) {
            if (size != other.size) {
                delete[] val;
                size = other.size;
                val = new T[size];
            }
            std::copy(other.val, other.val + size, val);
        }
        return *this;
    }

    ~spinor_template() {
        delete[] val;
    }

    void clearBuffer() {
        std::fill(val, val + size, T{});
    }
};

using spinor = spinor_template<c_double>;
using re_field = spinor_template<double>;

/*
	Modulo operation
*/
inline int mod(int a, int b) {
	int r = a % b;
	return r < 0 ? r + b : r;
}

#endif 