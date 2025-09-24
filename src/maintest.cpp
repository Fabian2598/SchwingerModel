#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include "conjugate_gradient.h"
#include <format>

//Formats decimal numbers
//For opening file with confs 
static std::string format(const double& number) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << number;
    std::string str = oss.str();
    str.erase(str.find('.'), 1); //Removes decimal dot 
    return str;
}

int main() {
    allocate_lattice_arrays();

    srand(19);
    
    
    double beta = 2; 
    double m0 = -0.18840579710144945;
    int nconf;
    if (LV::Nx == 64)
        nconf = 0;
    else if (LV::Nx == 128)
        nconf = 3;
    else if (LV::Nx == 256)
        nconf = 20;  

    CG::max_iter = 10000; //Maximum number of iterations for the conjugate gradient method
    CG::tol = 1e-10; //Tolerance for convergence
    //---Input data---//
    std::cout << "  -----------------------------" << std::endl;
    std::cout << "        CG tests        " << std::endl;
    std::cout << "  -----------------------------" << std::endl;
    std::cout << "Nx " << LV::Nx << " Nt " << LV::Nt << std::endl;
    std::cout << "m0 " << m0 << std::endl;
    std::cout << "beta " << beta << std::endl;
    std::cout << " " << std::endl;
    
    periodic_boundary(); //Compute right and left periodic boundary
    
	GaugeConf GConf = GaugeConf(LV::Nx, LV::Nt);  //Gauge configuration
    GConf.initialization(); //Random initialization of the gauge configuration 
    

    spinor sol,rhs;
    for(int n = 0; n < LV::Ntot; n++) {
        rhs.mu0[n] = RandomU1(); //spin up
        rhs.mu1[n] = RandomU1(); //spin down
    }
    
    

    std::ostringstream NameData;
    NameData << "../confs/b" << beta << "_" << LV::Nx << "x" << LV::Nt << "/m-018/2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt << "_b" << 
    format(beta).c_str() << "_m" << format(m0).c_str() << "_" << nconf << ".ctxt";
    GConf.read_conf(NameData.str());

    //Checking that the operatons have no problem 
    GConf.Compute_Staple(); //Compute staples 
    GConf.Compute_Plaquette01(); //Compute plaquettes
    double Sp = GConf.MeasureSp_HMC(); //Measure average plaquette
    double Sg = GConf.Compute_gaugeAction(beta); //Compute gauge action
    std::cout << std::format("Average plaquette: {:.6f}", Sp) << std::endl;
    std::cout << std::format("Gauge action: {:.6f}", Sg) << std::endl;
    ///----------------------------------///
    
    double start, end;
    start = omp_get_wtime();
    int iter = conjugate_gradient(GConf.Conf, rhs, sol, m0);
    end = omp_get_wtime();
    double time_taken = end - start; //in seconds
    std::cout << std::format("Time taken by CG: {:.6f} seconds", time_taken) << std::endl;
    


    

    //Free coordinate arrays
    free_lattice_arrays();

	return 0;
}

