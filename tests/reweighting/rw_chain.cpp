#include "hmc.h"

/*
    HMC chain for the reweighting test (see run_reweighting_test.sh). It runs on a single rank and is only
    meant for very small lattices, since det(D) is computed from the dense Dirac matrix.
    For every trajectory after thermalization it writes one line with
        plaquette    dH    Re log det(D)    Im log det(D) (mod 2 pi)
    Usage:
        rw_chain beta m0 Ntherm Ntraj MD_steps trajectory_length output_file
*/

//Access to the private members of HMC (declared as friend in hmc.h)
struct HMC_test_access {
    static void HMC_Update(HMC& hmc) { hmc.HMC_Update(); }
    static GaugeConf& GConf(HMC& hmc) { return hmc.GConf; }
};
using priv = HMC_test_access; //short name (access would clash with access() from unistd.h)

//log det of a dense complex n x n matrix (column-major), using LU decomposition with partial pivoting
c_double log_det(std::vector<c_double> A, int n) {
    c_double logdet = 0;
    for (int k = 0; k < n; k++) {
        int p = k;
        for (int i = k + 1; i < n; i++)
            if (std::abs(A[i + k*n]) > std::abs(A[p + k*n])) p = i;
        if (p != k) {
            for (int j = 0; j < n; j++) std::swap(A[k + j*n], A[p + j*n]);
            logdet += c_double(0, pi); //row swap -> factor -1
        }
        c_double pivot = A[k + k*n];
        logdet += std::log(pivot);
        for (int i = k + 1; i < n; i++) {
            c_double f = A[i + k*n] / pivot;
            for (int j = k + 1; j < n; j++) A[i + j*n] -= f * A[k + j*n];
        }
    }
    return logdet;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi::size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi::rank);
    if (argc != 8) {
        std::cerr << "Usage: " << argv[0] << " beta m0 Ntherm Ntraj MD_steps trajectory_length output_file" << std::endl;
        MPI_Finalize();
        return 1;
    }
    //Chains are started in parallel, so the seed cannot depend on the time only
    std::random_device rd;
    srand(rd());

    mpi::ranks_x = 1;
    mpi::ranks_t = 1;
    initializeMPI();

    double beta = std::atof(argv[1]), m0 = std::atof(argv[2]);
    int Ntherm = std::atoi(argv[3]), Ntraj = std::atoi(argv[4]), MD_steps = std::atoi(argv[5]);
    double trajectory_length = std::atof(argv[6]);
    sim_params::beta = beta;
    sim_params::m0 = m0;

    GaugeConf GConf = GaugeConf();
    HMC hmc = HMC(GConf, MD_steps, trajectory_length, Ntherm, Ntraj, 1, beta, LV::Nx, LV::Nt, m0, 0);
    GaugeConf& U = priv::GConf(hmc);
    U.initialization();
    for (int i = 0; i < Ntherm; i++) priv::HMC_Update(hmc);

    const int dim = 2 * LV::Ntot;
    std::vector<c_double> D(dim * dim);
    spinor e(mpi::maxSizeH), De(mpi::maxSizeH);
    std::ofstream out(argv[7]);
    out << std::setprecision(15);
    for (int i = 0; i < Ntraj; i++) {
        priv::HMC_Update(hmc);
        exchange_halo(U.Conf.val);
        U.Compute_Plaquette01();
        double plaquette = U.MeasureSp_HMC() / LV::Ntot;

        //Dense Dirac matrix, one column per unit vector
        int col = 0;
        for (int x = 1; x <= LV::Nx; x++)
            for (int t = 1; t <= LV::Nt; t++)
                for (int mu = 0; mu < 2; mu++, col++) {
                    e.clearBuffer();
                    e.val[idx(x, t, mu)] = 1;
                    D_phi(U, e, De);
                    int row = 0;
                    for (int y = 1; y <= LV::Nx; y++)
                        for (int s = 1; s <= LV::Nt; s++)
                            for (int nu = 0; nu < 2; nu++, row++)
                                D[row + col * dim] = De.val[idx(y, s, nu)];
                }
        c_double logdet = log_det(D, dim);
        out << plaquette << " " << hmc.getDeltaH() << " " << logdet.real() << " "
            << std::remainder(logdet.imag(), 2 * pi) << "\n";
    }
    out.close();
    MPI_Finalize();
    return 0;
}
