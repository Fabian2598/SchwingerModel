#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include "gauge_conf.h"

TEST_CASE("GaugeConf plaquette", "[Plaquette]") {
    // Create a simple gauge configuration (identity links)
    GaugeConf g;
    g.initialization();
    exchange_halo(g.Conf.val);
    g.Compute_Plaquette01();

    bool u1_vars = true; //we check that the plaquettes are U1 vars
    for(int x = 1; x<=mpi::width_x; x++){
		for(int t = 1; t<=mpi::width_t; t++){
			int n = x*(mpi::width_t+2)+t;
            if (std::abs(g.Plaquette01[n] * std::conj(g.Plaquette01[n]) - 1.0) > 1e-10){
                u1_vars = false;
                break;
            }
        }
    }
    REQUIRE(u1_vars);
}

TEST_CASE("GaugeConf Q01", "[Q01]") {
    // Create a simple gauge configuration (identity links)
    GaugeConf g;
    g.initialization();
    exchange_halo(g.Conf.val);
    g.Compute_Q();

    bool notzero = true; //check that Q01 are not all zero
    for(int x = 1; x<=mpi::width_x; x++){
		for(int t = 1; t<=mpi::width_t; t++){
			int n = x*(mpi::width_t+2)+t;
            if (std::abs(g.Q01[n]) < 1e-10){
                notzero = false;
                break;
            }
        }
    }
    REQUIRE(notzero);
}

TEST_CASE("GaugeConf Q10", "[Q10]") {
    // Create a simple gauge configuration (identity links)
    GaugeConf g;
    g.initialization();
    exchange_halo(g.Conf.val);
    g.Compute_Q();

    bool notzero = true; //check that Q01 are not all zero
    for(int x = 1; x<=mpi::width_x; x++){
		for(int t = 1; t<=mpi::width_t; t++){
			int n = x*(mpi::width_t+2)+t;
            if (std::abs(g.Q10[n]) < 1e-10){
                notzero = false;
                break;
            }
        }
    }
    REQUIRE(notzero);
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