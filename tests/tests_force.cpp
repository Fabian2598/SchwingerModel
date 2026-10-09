#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <array>
#include "hmc.h"

/*
    Tests of the pseudofermion heatbath and of the HMC force. They do not depend on the fermion
    discretization, so they apply to every build (Wilson, twisted mass, clover and four flavors):
        - phi_k = D chi_k has to satisfy phi_k^dagger (DD^dagger)^-1 phi_k = chi_k^dagger chi_k
        - The force has to be minus the derivative of the action, F_mu(n) = -dS/d omega_mu(n). This is
          checked against a central finite difference of HMC::Action at fixed pseudofermion fields, for
          links on the time boundary and on the boundaries between ranks.
*/

//Access to the private members of HMC (declared as friend in hmc.h)
struct HMC_test_access {
    static void RandomCHI(HMC& hmc) { hmc.RandomCHI(); }
    static const spinor& chi(const HMC& hmc) { return hmc.chi; }
    static GaugeConf& GConf(HMC& hmc) { return hmc.GConf; }
    static double Action(HMC& hmc, GaugeConf& GConfig, const std::vector<spinor>& phi) { return hmc.Action(GConfig, phi); }
    static void Force(HMC& hmc, GaugeConf& GConfig, const std::vector<spinor>& phi) { hmc.Force(GConfig, phi); }
    static const re_field& Forces(const HMC& hmc) { return hmc.Forces; }
};
using priv = HMC_test_access; //short name (access would clash with access() from unistd.h)

//Parameters for both tests. A heavy mass keeps the solver fast on a random gauge configuration.
void set_parameters() {
    sim_params::beta = 2;
    sim_params::m0 = 0.5;
    sim_params::tm = 0.5;  //only used with TWISTED_MASS
    sim_params::csw = 1;   //only used with CLOVER
    CG::tol = 1e-14;       //the finite differences need an accurate action
    CG::max_iter = 100000;
}

TEST_CASE("Pseudofermion heatbath", "[HMC heatbath]") {
    using namespace sim_params;
    set_parameters();
    GaugeConf GConf = GaugeConf();
    HMC hmc = HMC(GConf, MD_steps, trajectory_length, Ntherm, Nmeas, Nsteps, beta, LV::Nx, LV::Nt, m0, 0);
    GaugeConf& U = priv::GConf(hmc);
    U.initialization();

    std::vector<spinor> phi(N_pf, spinor(mpi::maxSizeH));
    spinor sol(mpi::maxSizeH);
    //Pseudofermions generated as in HMC::HMC_Update, checking each field before chi is overwritten
    exchange_halo(U.Conf.val);
    #ifdef CLOVER
        U.Compute_Q();
    #endif
    for (int k = 0; k < N_pf; k++) {
        priv::RandomCHI(hmc);
        const spinor& chi = priv::chi(hmc);
        D_phi(U, chi, phi[k]);
        conjugate_gradient(U, phi[k], sol); //(DD^dagger)^-1 phi_k
        double action = std::real(dot(sol, phi[k]));
        double chi_norm2 = std::real(dot(chi, chi));
        INFO("Pseudofermion field " << k << ": phi^dagger (DD^dagger)^-1 phi = " << action << ", chi^dagger chi = " << chi_norm2);
        REQUIRE(std::abs(action - chi_norm2) < 1e-10 * chi_norm2);
    }
}

TEST_CASE("Force is minus the derivative of the action", "[HMC force]") {
    using namespace sim_params;
    set_parameters();
    GaugeConf GConf = GaugeConf();
    HMC hmc = HMC(GConf, MD_steps, trajectory_length, Ntherm, Nmeas, Nsteps, beta, LV::Nx, LV::Nt, m0, 0);
    GaugeConf& U = priv::GConf(hmc);
    U.initialization();

    //Pseudofermions generated as in HMC::HMC_Update
    std::vector<spinor> phi(N_pf, spinor(mpi::maxSizeH));
    exchange_halo(U.Conf.val);
    #ifdef CLOVER
        U.Compute_Q();
    #endif
    for (int k = 0; k < N_pf; k++) {
        priv::RandomCHI(hmc);
        D_phi(U, priv::chi(hmc), phi[k]);
    }
    priv::Force(hmc, U, phi);
    const re_field& F = priv::Forces(hmc);

    //Global coordinates (x, t, mu) of the links. They include the last time slice (antiperiodic boundary)
    //and, for 2x2 ranks, links on both sides of the boundaries between ranks.
    const int Nx = LV::Nx, Nt = LV::Nt;
    std::vector<std::array<int, 3>> links;
    for (int mu = 0; mu < 2; mu++) {
        links.push_back({0, 0, mu});
        links.push_back({0, Nt-1, mu});
        links.push_back({Nx-1, Nt-1, mu});
        links.push_back({Nx/2-1, Nt/2-1, mu});
        links.push_back({Nx/2, Nt/2-1, mu});
        links.push_back({Nx/2-1, Nt/2, mu});
    }

    const double eps = 1e-5;
    for (const auto& [x, t, mu] : links) {
        //Only the rank that owns the link perturbs it
        bool owner = (mpi::coords[0] == x / mpi::width_x) && (mpi::coords[1] == t / mpi::width_t);
        int i = idx(x % mpi::width_x + 1, t % mpi::width_t + 1, mu);
        double F_local = owner ? F.val[i] : 0.0, F_link;
        MPI_Allreduce(&F_local, &F_link, 1, MPI_DOUBLE, MPI_SUM, mpi::cart_comm);

        //S[U_mu(x) exp(+i eps)] and S[U_mu(x) exp(-i eps)]
        double S[2];
        for (int s = 0; s < 2; s++) {
            GaugeConf U_eps = U;
            if (owner) U_eps.Conf.val[i] *= std::exp((s == 0 ? 1.0 : -1.0) * eps * I_number);
            S[s] = priv::Action(hmc, U_eps, phi);
        }
        double dS = (S[0] - S[1]) / (2 * eps);

        //Every rank has the same values, so all of them reach the same verdict
        INFO("Link x = " << x << ", t = " << t << ", mu = " << mu << ": F = " << F_link << ", -dS/domega = " << -dS);
        REQUIRE(std::abs(F_link + dS) < 1e-4);
    }
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
