#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include "dirac_operator.h"

// Test: Check that D_phi produces non-zero output for non-zero input
TEST_CASE("D_phi produces non-zero output", "[dirac_operator]") {
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
    D_phi(g, phi, D_phi_output);
    
    // Check that output is not all zeros
    bool has_nonzero = false;
    for (int x = 1; x <= mpi::width_x; x++) {
        for (int t = 1; t <= mpi::width_t; t++) {
            n = x*(mpi::width_t+2)+t;
            for(int mu = 0; mu<2; mu++){
                if (std::abs(D_phi_output.val[2*n+mu]) > 1e-10) {
                    has_nonzero = true;
                    break;
                }
            }
        }
    }
    
    REQUIRE(has_nonzero);
}

// Test: Check that D_phi produces non-zero output for non-zero input
TEST_CASE("D_Ddagg_phi produces non-zero output", "[dirac_operator]") {
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
    D_dagger_phi(g, phi, D_phi_output);
    
    // Check that output is not all zeros
    bool has_nonzero = false;
    for (int x = 1; x <= mpi::width_x; x++) {
        for (int t = 1; t <= mpi::width_t; t++) {
            n = x*(mpi::width_t+2)+t;
            for(int mu = 0; mu<2; mu++){
                if (std::abs(D_phi_output.val[2*n+mu]) > 1e-10) {
                    has_nonzero = true;
                    break;
                }
            }
        }
    }
    
    REQUIRE(has_nonzero);
}


 // Test: Check that phi_dag_partialD_phi produces non-zero output for non-zero input
TEST_CASE("phi_dag_partialD_phi produces non-zero output", "[dirac_operator]") {
    // Create a simple gauge configuration (identity links)
    GaugeConf g;
    g.initialization();
    exchange_halo(g.Conf.val);
    
    // Create input spinor with simple values
    spinor right_term(mpi::maxSizeH);
    spinor left_term(mpi::maxSizeH);
    re_field D_phi_output(mpi::maxSizeH);
    int n;
	for (int x = 1; x <= mpi::width_x; x++) {
        for (int t = 1; t <= mpi::width_t; t++) {
            n = x*(mpi::width_t+2)+t;
            right_term.val[2*n] = RandomU1();
            right_term.val[2*n+1] = RandomU1();
            left_term.val[2*n] = RandomU1();
            left_term.val[2*n+1] = RandomU1();
        }
    }
  
    sim_params::m0 = 0.5;
    phi_dag_partialD_phi(g, left_term,right_term,D_phi_output);
    
    // Check that output is not all zeros
    bool has_nonzero = false;
    for (int x = 1; x <= mpi::width_x; x++) {
        for (int t = 1; t <= mpi::width_t; t++) {
            n = x*(mpi::width_t+2)+t;
            for(int mu = 0; mu<2; mu++){
                if (std::abs(D_phi_output.val[2*n+mu]) > 1e-10) {
                    has_nonzero = true;
                    break;
                }
            }
        }
    }
    
    REQUIRE(has_nonzero);
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi::rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi::size);

    srand((mpi::rank + 1) * time(0));

    mpi::ranks_t = 2;
    mpi::ranks_x = 2;

    initializeMPI();

    int result = Catch::Session().run(argc, argv);

    MPI_Finalize();

    return result;
}
