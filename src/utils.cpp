#include "utils.h"

 //----------Jackknife---------//
 std::vector<double> samples_mean(std::vector<double> dat, int bin) {
    if (bin < 2 || static_cast<int>(dat.size()) < 2 * bin) {
        throw std::runtime_error("samples_mean: need bin >= 2 and at least "
                                 "2 samples per block");
    }
    const int blk  = dat.size() / bin; //block size
    const int used = bin * blk;        //trailing dat.size()%bin samples dropped
    double total = 0.0;
    for (int j = 0; j < used; j++) total += dat[j];

     std::vector<double> samples_mean(bin);

     for (int i = 0; i < bin; i++) {
        double block = 0.0;
        for (int j = i * blk; j < (i + 1) * blk; j++) block += dat[j];
        samples_mean[i] = (total - block) / (used - blk);
    }
    return samples_mean;
 }



 double Jackknife_error(std::vector<double> dat, int bin) {
    double error = 0;
    std::vector<double> sm = samples_mean(dat, bin);
    double sm_mean = mean(sm);
     for (int m = 0; m < bin; m++) {
        error += pow((sm[m] - sm_mean), 2);
    }
    error = sqrt(error * (bin - 1) / bin);
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
        std::cout << "* mu (twisted mass) = " << sim_params::tm << std::endl;
        std::cout << "* csw (clover term constant) = " << sim_params::csw << std::endl;
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
        const char* hostname = std::getenv("HOSTNAME");
        std::cout << "* Host: " << (hostname ? hostname : "unknown") << std::endl;
        std::cout << "* Start time: " << start_time_str << std::endl;
        std::cout << "**********************************************************************" << std::endl;
    }
        
}