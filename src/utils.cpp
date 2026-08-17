#include "utils.h"

//----------Jackknife---------//
std::vector<double> samples_mean(std::vector<double> dat, int bin) {
    std::vector<double> samples_mean(bin);
    int dat_bin = dat.size() / bin;
    double prom = 0;
    for (int i = 0; i < bin; i++) {
        for (int k = 0; k < bin; k++) {
            for (int j = k * dat_bin; j < k * dat_bin + dat_bin; j++) {
                if (k != i) {
                    prom += dat[j];
                }
            }
        }
        prom = prom / (dat.size() - dat_bin);
        samples_mean[i] = prom;
        prom = 0;
    }
    return samples_mean;
}

double Jackknife_error(std::vector<double> dat, int bin) {
    double error = 0;
    std::vector<double> sm = samples_mean(dat, bin);
    double normal_mean = mean(dat);
    for (int m = 0; m < bin; m++) {
        error += pow((sm[m] - normal_mean), 2);
    }
    error = sqrt(error * (bin - 1) / bin);
    return error;
}


double Jackknife(std::vector<double> dat, std::vector<int> bins) {
    std::vector<double> errores(bins.size());
    double error;
    for (int i = 0; i < bins.size(); i++) {
        errores[i] = Jackknife_error(dat, bins[i]);
    }
    error = *std::max_element(errores.begin(), errores.end());
    return error;
}
//-------------End of Jackknife--------------//

void print_parameters(){
    using namespace sim_params;
    if (mpi::rank == 0){
        std::cout << "**********************************************************************" << std::endl;
        std::cout << "*                              PARAMETERS" << std::endl;
        std::cout << "* Nx = " << LV::Nx << ", Nt = " << LV::Nt << std::endl;
        std::cout << "* m0 = " << m0 << ", kappa = " << 1/(2*(m0+2)) << std::endl;
        std::cout << "* beta = " << beta << std::endl;
        std::cout << "* Thermalization confs = " << Ntherm << std::endl;
        std::cout << "* Measurement confs = " << Nmeas << std::endl;
        std::cout << "* Decorrelation steps (confs dropped between measurements) = " << Nsteps << std::endl;
        std::cout << "* Trajectory length = " << trajectory_length << ", Leapfrog steps = " << MD_steps << 
        ", Integration step = " << trajectory_length/MD_steps << std::endl;
        std::cout << "* CG max iterations = " << CG::max_iter << ", CG tolerance = " << CG::tol << std::endl;
        std::cout << "* Number of ranks on x = " << mpi::ranks_x << ", Number of ranks on t = "  << mpi::ranks_t << std::endl;
        std::cout << "* Total number of MPI ranks = " << mpi::size << std::endl;
        std::cout << "* Each rank has " << mpi::width_x*mpi::width_t << " lattice sites" << std::endl;
        std::cout << "* Host: " << std::getenv("HOSTNAME") << std::endl;
        std::cout << "* Start time: " << start_time_str << std::endl;
        std::cout << "**********************************************************************" << std::endl;
    }
        
}