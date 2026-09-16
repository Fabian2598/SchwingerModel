#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include "hmc.h"


TEST_CASE("HMC", "[HMC]") {

    sim_params::beta = 4;
    sim_params::trajectory_length = 1;
    sim_params::MD_steps = 5;
    sim_params::Ntherm = 1;
    sim_params::Nmeas = 1;
    sim_params::Nsteps = 1;
    sim_params::m0 = 1;
    sim_params::tm = 0.5;
    sim_params::csw = 1;
    using namespace sim_params;
    int saveconf=0;

    GaugeConf GConf = GaugeConf();  //Initial gauge configuration         
    HMC hmc = HMC(GConf,MD_steps, trajectory_length, Ntherm, Nmeas, Nsteps, beta, LV::Nx, LV::Nt, m0,saveconf);   
    double begin = MPI_Wtime();
    hmc.HMC_algorithm();
    double end = MPI_Wtime();

    REQUIRE(true);

}


int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi::rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi::size);

    srand((mpi::rank + 1) * time(0));
    mpi::ranks_t = 2;
    mpi::ranks_x = 2;

    if (mpi::size!=4 && mpi::rank == 0){
        std::cerr << "ERROR: This test is meant to be run with 4 ranks" << std::endl;
        MPI_Finalize();
        return 0;
    }
    initializeMPI();
    
    int result = Catch::Session().run(argc, argv);
    MPI_Finalize();

    return result;
}