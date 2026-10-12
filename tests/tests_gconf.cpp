#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include "gauge_conf.h"

TEST_CASE("GaugeConf plaquette", "[Plaquette]") {
    // Create a simple gauge configuration (identity links)
    GaugeConf g;
    g.initialization();
    exchange_halo_gauge(g.Conf.val);
    g.Compute_Plaquette01();

    std::vector<c_double> plaquette = g.Plaquette01.toHost();
    int u1_vars = 1; //we check that the plaquettes are U1 vars
    for(int x = 1; x<=mpi::width_x; x++){
		for(int t = 1; t<=mpi::width_t; t++){
			int n = x*(mpi::width_t+2)+t;
            if (abs(plaquette[n] * conj(plaquette[n]) - 1.0) > 1e-10){
                u1_vars = 0;
            }
        }
    }
    int global_u1_vars;
    MPI_Allreduce(&u1_vars, &global_u1_vars, 1, MPI_INT, MPI_LAND, mpi::cart_comm);
    REQUIRE(global_u1_vars == 1);
}

TEST_CASE("GaugeConf Q01 = Q10^+", "[Q01 and Q10]") {
    // Create a simple gauge configuration (identity links)
    GaugeConf g;
    g.initialization();
    exchange_halo_gauge(g.Conf.val);
    g.Compute_Q();

    std::vector<c_double> Q01 = g.Q01.toHost();
    std::vector<c_double> Q10 = g.Q10.toHost();
    int notzero = 1; //check that Q01 = Q10^+ for all n
    for(int x = 1; x<=mpi::width_x; x++){
		for(int t = 1; t<=mpi::width_t; t++){
			int n = x*(mpi::width_t+2)+t;
            if (abs( Q01[n] - conj(Q10[n]) ) > 1e-10){
                notzero = 0;
            }
        }
    }
    int global_notzero;
    MPI_Allreduce(&notzero, &global_notzero, 1, MPI_INT, MPI_LAND, mpi::cart_comm);
    REQUIRE(global_notzero == 1);
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
