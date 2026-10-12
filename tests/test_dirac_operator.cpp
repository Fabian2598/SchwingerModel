#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include "dirac_operator.h"

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

double magnitude(const double v){ return std::fabs(v); }
double magnitude(const c_double v){ return abs(v); }

//Checks on all ranks that the field is not zero
template <typename T>
bool has_nonzero(const spinor_template<T>& field){
    std::vector<T> field_host = field.toHost();
    int nonzero = 0;
    int n;
    for (int x = 1; x <= mpi::width_x; x++) {
        for (int t = 1; t <= mpi::width_t; t++) {
            n = x*(mpi::width_t+2)+t;
            for(int mu = 0; mu<2; mu++){
                if (magnitude(field_host[2*n+mu]) > 1e-10) {
                    nonzero = 1;
                }
            }
        }
    }
    int global_nonzero;
    MPI_Allreduce(&nonzero, &global_nonzero, 1, MPI_INT, MPI_LAND, mpi::cart_comm);
    return global_nonzero == 1;
}

// Test: Check that D_phi produces non-zero output for non-zero input
TEST_CASE("D_phi produces non-zero output", "[dirac_operator]") {
    // Create a simple gauge configuration (identity links)
    GaugeConf g;
    g.initialization();
    exchange_halo_gauge(g.Conf.val);
    #ifdef CLOVER
        g.Compute_Q();
    #endif

    // Create input spinor with simple values
    spinor phi = random_spinor();
    spinor D_phi_output(mpi::maxSizeH);

    D_phi(g, phi, D_phi_output);

    // Check that output is not all zeros
    REQUIRE(has_nonzero(D_phi_output));
}

// Test: Check that D_phi produces non-zero output for non-zero input
TEST_CASE("D_Ddagg_phi produces non-zero output", "[dirac_operator]") {
    // Create a simple gauge configuration (identity links)
    GaugeConf g;
    g.initialization();
    exchange_halo_gauge(g.Conf.val);
    #ifdef CLOVER
        g.Compute_Q();
    #endif

    // Create input spinor with simple values
    spinor phi = random_spinor();
    spinor D_phi_output(mpi::maxSizeH);
    D_dagger_phi(g, phi, D_phi_output);

    // Check that output is not all zeros
    REQUIRE(has_nonzero(D_phi_output));
}

// Test: D^+ is the adjoint of D, (chi, D phi) = (D^+ chi, phi). It checks the hopping terms, the boundary
// signs, the clover term and the halo exchange of D and D^+ against each other.
TEST_CASE("D_dagger_phi is the adjoint of D_phi", "[dirac_operator]") {
    GaugeConf g;
    g.initialization();
    exchange_halo_gauge(g.Conf.val);
    #ifdef CLOVER
        g.Compute_Q();
    #endif

    spinor phi = random_spinor();
    spinor chi = random_spinor();
    spinor D_phi_output(mpi::maxSizeH);
    spinor D_dagger_chi(mpi::maxSizeH);
    D_phi(g, phi, D_phi_output);
    D_dagger_phi(g, chi, D_dagger_chi);

    //dot(A,B) = sum_i A_i conj(B_i) = B^+ A
    c_double chi_D_phi = dot(D_phi_output, chi);      //chi^+ D phi
    c_double Ddagg_chi_phi = dot(phi, D_dagger_chi);  //(D^+ chi)^+ phi
    REQUIRE(abs(chi_D_phi - Ddagg_chi_phi) < 1e-10 * abs(chi_D_phi));
}


 // Test: Check that phi_dag_partialD_phi produces non-zero output for non-zero input
TEST_CASE("phi_dag_partialD_phi produces non-zero output", "[dirac_operator]") {
    // Create a simple gauge configuration (identity links)
    GaugeConf g;
    g.initialization();
    exchange_halo_gauge(g.Conf.val);

    // Create input spinor with simple values
    spinor right_term = random_spinor();
    spinor left_term = random_spinor();
    re_field D_phi_output(mpi::maxSizeH);

    phi_dag_partialD_phi(g, left_term,right_term,D_phi_output);

    // Check that output is not all zeros
    REQUIRE(has_nonzero(D_phi_output));
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi::rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi::size);

    srand((mpi::rank + 1) * time(0));

    mpi::ranks_t = 2;
    mpi::ranks_x = 2;
    sim_params::m0 = 0.5;
    sim_params::tm = 0.1;
    sim_params::csw = 1;

    initializeMPI();

    int result = Catch::Session().run(argc, argv);

    MPI_Finalize();

    return result;
}
