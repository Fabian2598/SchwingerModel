#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include <chrono>
#include <sstream>
#include <iomanip>
#include "mpi_setup.h"
#include "hmc.h"
#include "tests.h"


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi::size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi::rank);
        
    srand((mpi::rank+1)*time(0));
    
    int Ntherm, Nmeas, Nsteps, Nm0; //Simulation parameters
    double beta; //Beta range
    double trajectory_length; //HMC parameters
    int MD_steps;
    double m0; //bare mass
	int saveconf = 0; //Save configurations
    

    //To call the sequential program one has to choose ranks_x = ranks_t = 1
    if (mpi::rank == 0){
         //---Input data---//
        std::cerr << "  -----------------------------" << std::endl;
        std::cerr << "|  Two-flavor Schwinger model   |" << std::endl;
        std::cerr << "| Hybrid Monte Carlo simulation |" << std::endl;
        std::cerr << "  -----------------------------" << std::endl;
        std::cerr << "Nx " << LV::Nx << " Nt " << LV::Nt << std::endl;
        std::cerr << "ranks_x: " << std::endl;
        std::cin >> mpi::ranks_x;
        std::cerr << "ranks_t: " << std::endl;
        std::cin >> mpi::ranks_t;
    }
    
    MPI_Bcast(&mpi::ranks_x, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(&mpi::ranks_t, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(&m0, 1, MPI_DOUBLE,  0, MPI_COMM_WORLD);
    MPI_Bcast(&MD_steps, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(&trajectory_length, 1, MPI_DOUBLE,  0, MPI_COMM_WORLD);
    MPI_Bcast(&beta, 1, MPI_DOUBLE,  0, MPI_COMM_WORLD);
    MPI_Bcast(&Ntherm, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(&Nmeas, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(&Nsteps, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(&saveconf, 1, MPI_INT,  0, MPI_COMM_WORLD);

    m0 = 0;
    sim_params::m0 = m0;
    sim_params::beta = beta;

    
    initializeMPI(); //2D rank topology
    allocate_lattice_arrays(); //Allocates memory for arrays of coordinates
    periodic_boundary(); //Stores neighbors

    
    Tests tests(sim_params::m0);
    tests.test_D_operator();
    tests.test_D_dagger_operator();
    tests.test_phi_dag_partialD_phi();
    tests.test_D_D_dagger_phi();
    tests.test_CG();

    //Free coordinate arrays
    free_lattice_arrays();
    MPI_Finalize();

	return 0;
}

