#ifndef VARIABLES_H_INCLUDED
#define VARIABLES_H_INCLUDED
#include "config.h"
#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include <utility>
#include "mpi.h"
#include "backend.h"
#include "complex_double.h" //c_double (complex number usable on the GPU) and I_number

constexpr double pi = 3.14159265359;
constexpr double l_8 = 1.0/8.0;
constexpr double l_16 = 1.0/16.0;

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
    extern double tm; //twisted mass
    extern double csw; //Sheikholeslami-Wohlert constant
}

/*
    Flattened field (spinor, gauge links, momenta, ...), halo included.
    When compiled with CUDA the data is stored in GPU memory: val can only be dereferenced inside the
    for_each_site / lattice_sum lambdas (see backend.h). toHost() and fromHost() copy the whole field
    between the GPU and the CPU.
*/
template <typename T>
struct spinor_template{
    T* val;
    int size;

    explicit spinor_template(int N = LV::Ntot) : val(static_cast<T*>(backend::allocate_bytes(sizeof(T)*N))), size(N) {
        clearBuffer();
    }

    spinor_template(const spinor_template& other) :
        val(static_cast<T*>(backend::allocate_bytes(sizeof(T)*other.size))), size(other.size) {
        backend::copy_bytes(val, other.val, sizeof(T)*size);
    }

    spinor_template(spinor_template&& other) noexcept : val(other.val), size(other.size) {
        other.val = nullptr;
        other.size = 0;
    }

    spinor_template& operator=(const spinor_template& other) {
        if (this != &other) {
            if (size != other.size) {
                backend::free_bytes(val);
                size = other.size;
                val = static_cast<T*>(backend::allocate_bytes(sizeof(T)*size));
            }
            backend::copy_bytes(val, other.val, sizeof(T)*size);
        }
        return *this;
    }

    spinor_template& operator=(spinor_template&& other) noexcept {
        std::swap(val, other.val);
        std::swap(size, other.size);
        return *this;
    }

    ~spinor_template() {
        backend::free_bytes(val);
    }

    void clearBuffer() {
        backend::zero_bytes(val, sizeof(T)*size);
    }

    //Copy of the field (halo included) in host memory
    std::vector<T> toHost() const {
        std::vector<T> host(size);
        backend::copy_bytes_to_host(host.data(), val, sizeof(T)*size);
        return host;
    }

    //Overwrites the field with data from host memory (host has to have the same size as the field)
    void fromHost(const std::vector<T>& host) {
        if (static_cast<int>(host.size()) != size) {
            std::cerr << "fromHost: the host buffer has " << host.size() << " elements, the field " << size << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        backend::copy_bytes_from_host(val, host.data(), sizeof(T)*size);
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
