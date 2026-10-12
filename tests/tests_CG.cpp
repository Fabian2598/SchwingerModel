#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include "conjugate_gradient.h"

//Spinor with random U(1) entries (generated on the CPU and copied to the GPU)
spinor random_spinor(){
    spinor phi(mpi::maxSizeH);
    std::vector<c_double> phi_host(phi.size);
    int n;
	for (int x = 1; x <= mpi::width_x; x++) {
        for (int t = 1; t <= mpi::width_t; t++) {
            n = x*(mpi::width_t+2)+t;
            phi_host[2*n] = RandomU1();
            phi_host[2*n+1] = RandomU1();
        }
    }
    phi.fromHost(phi_host);
    return phi;
}

//Checks on all ranks that two spinors are equal up to tol
bool equal_spinors(const spinor& a, const spinor& b, const double tol){
    std::vector<c_double> a_host = a.toHost();
    std::vector<c_double> b_host = b.toHost();
    int equal = 1;
    int n;
    for (int x = 1; x <= mpi::width_x; x++) {
        for (int t = 1; t <= mpi::width_t; t++) {
            n = x*(mpi::width_t+2)+t;
            for(int mu = 0; mu<2; mu++){
                if (abs(a_host[2*n+mu]-b_host[2*n+mu]) > tol) {
                    equal = 0;
                }
            }
        }
    }
    int global_equal;
    MPI_Allreduce(&equal, &global_equal, 1, MPI_INT, MPI_LAND, mpi::cart_comm);
    return global_equal == 1;
}

TEST_CASE("CG test convergence", "[CG]") {
    // Create a simple gauge configuration (identity links)
    GaugeConf g;
    g.initialization();
    exchange_halo_gauge(g.Conf.val);
    #ifdef CLOVER
        g.Compute_Q();
    #endif
    sim_params::m0 = 0.5;
    sim_params::tm = 0.01;
    sim_params::csw = 1;
    // Create input spinor with simple values
    spinor phi = random_spinor();
    spinor D_phi_output(mpi::maxSizeH);

    CG::print_convergence_message = false;
    int it = conjugate_gradient(g,phi,D_phi_output);
    //When it == 0 it means the algorithm didn't converge
    REQUIRE(it != 0);
}

TEST_CASE("CG converges to the right solution", "[CG sol]") {
    // Create a simple gauge configuration (identity links)
    GaugeConf g;
    g.initialization();
    exchange_halo_gauge(g.Conf.val);

    #ifdef CLOVER
        g.Compute_Q();
    #endif
    // Create input spinor with simple values
    spinor phi = random_spinor();
    spinor sol(mpi::maxSizeH);
    spinor D_D_dagg_sol(mpi::maxSizeH);

    sim_params::m0 = 0.5;
    sim_params::tm = 0.01;
    sim_params::csw = 1;

    CG::print_convergence_message = false;
    CG::tol = 1e-12;
    int it = conjugate_gradient(g,phi,sol);


    D_D_dagger_phi(g,sol,D_D_dagg_sol);

    REQUIRE(equal_spinors(phi, D_D_dagg_sol, 1e-10));
}

TEST_CASE("Clover vs No Clover", "[CG clover and no clover]") {
    // Create a simple gauge configuration (identity links)
    GaugeConf g;
    g.initialization();
    exchange_halo_gauge(g.Conf.val);

    #ifdef CLOVER
        g.Compute_Q();
    #endif
    // Create input spinor with simple values
    spinor phi = random_spinor();
    spinor sol(mpi::maxSizeH);
    spinor D_D_dagg_sol(mpi::maxSizeH);

    sim_params::m0 = 0.5;
    sim_params::tm = 0.01;
    sim_params::csw = 1;
    CG::print_convergence_message = true;
    CG::tol = 1e-15;

    int it = conjugate_gradient(g,phi,sol);
    if (mpi::rank == 0)
        std::cout << "Test with clover converged in " << it << " iterations\n";

    D_D_dagger_phi(g,sol,D_D_dagg_sol);

    REQUIRE(equal_spinors(phi, D_D_dagg_sol, 1e-10));

    sim_params::csw = 0;

    spinor sol_no_clover(mpi::maxSizeH);
    spinor D_D_dagg_sol_no_clover(mpi::maxSizeH);
    it = conjugate_gradient(g,phi,sol_no_clover);
    if (mpi::rank == 0)
        std::cout << "Test with no clover converged in " << it << " iterations\n";
    D_D_dagger_phi(g,sol_no_clover,D_D_dagg_sol_no_clover);

    REQUIRE(equal_spinors(phi, D_D_dagg_sol_no_clover, 1e-10));
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
