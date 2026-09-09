#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include "conjugate_gradient.h"

TEST_CASE("CG test convergence", "[CG]") {
    // Create a simple gauge configuration (identity links)
    GaugeConf g;
    g.initialization();
    exchange_halo(g.Conf.val);
    
    // Create input spinor with simple values
    spinor phi(mpi::maxSizeH);
    int n;
	for (int x = 1; x <= mpi::width_x; x++) {
        for (int t = 1; t <= mpi::width_t; t++) {
            n = x*(mpi::width_t+2)+t;
            phi.val[2*n] = RandomU1();
            phi.val[2*n+1] = RandomU1();
        }
    }
        
    spinor D_phi_output(mpi::maxSizeH);    
    sim_params::m0 = 0.5;
    sim_params::tm = 0.1;
    sim_params::csw = 1;

    CG::print_convergence_message = true;
    int it = conjugate_gradient(g,phi,D_phi_output);
    //When it == 0 it means the algorithm didn't converge
    REQUIRE(it != 0);
}

TEST_CASE("CG converges to the right solution", "[CG sol]") {
    // Create a simple gauge configuration (identity links)
    GaugeConf g;
    g.initialization();
    exchange_halo(g.Conf.val);
    
    // Create input spinor with simple values
    spinor phi(mpi::maxSizeH);
    int n;
	for (int x = 1; x <= mpi::width_x; x++) {
        for (int t = 1; t <= mpi::width_t; t++) {
            n = x*(mpi::width_t+2)+t;
            phi.val[2*n] = RandomU1();
            phi.val[2*n+1] = RandomU1();
        }
    }
        
    spinor sol(mpi::maxSizeH);    
    spinor D_D_dagg_sol(mpi::maxSizeH);
    sim_params::m0 = 0.5;
    sim_params::tm = 0.1;
    sim_params::csw = 1;

    CG::print_convergence_message = true;
    CG::tol = 1e-15;
    int it = conjugate_gradient(g,phi,sol);
    
    D_D_dagger_phi(g,sol,D_D_dagg_sol);

    bool equal = true;
    for (int x = 1; x <= mpi::width_x; x++) {
        for (int t = 1; t <= mpi::width_t; t++) {
            n = x*(mpi::width_t+2)+t;
            for(int mu = 0; mu<2; mu++){
                if (std::abs(phi.val[2*n+mu]-D_D_dagg_sol.val[2*n+mu]) > 1e-10) {
                    equal = false;
                    break;
                }
            }
        }
    }
    
    REQUIRE(equal);
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
