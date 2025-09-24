#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include "hmc.h"
#include <format>


int main() {
    allocate_lattice_arrays();
    srand(time(0));
    
    int Ntherm, Nmeas, Nsteps, Nm0; //Simulation parameters
    double beta; //Beta range
    double trajectory_length; //HMC parameters
    int MD_steps;
    double m0_min, m0_max; //bare mass
	int saveconf = 0; //Save configurations

    CG::max_iter = 10000; //Maximum number of iterations for the conjugate gradient method
    CG::tol = 1e-10; //Tolerance for convergence
    //---Input data---//
    std::cout << "  -----------------------------" << std::endl;
    std::cout << "|  Two-flavor Schwinger model   |" << std::endl;
    std::cout << "| Hybrid Monte Carlo simulation |" << std::endl;
    std::cout << "  -----------------------------" << std::endl;
    std::cout << "Nx " << LV::Nx << " Nt " << LV::Nt << std::endl;
    std::cout << "m0 min: ";
    std::cin >> m0_min;
    std::cout << "m0 max: ";
    std::cin >> m0_max;
    std::cout << "Number of masses in [m0_min, m0_max]: ";
    std::cin >> Nm0;
    std::cout << "Molecular dynamics steps: ";
    std::cin >> MD_steps;
    std::cout << "Trajectory length: ";
    std::cin >> trajectory_length; 
    std::cout << "beta: ";
    std::cin >> beta;
    std::cout << "Thermalization: ";
    std::cin >> Ntherm;
    std::cout << "Measurements: ";
    std::cin >> Nmeas;
    std::cout << "Step (sweeps between measurements): ";
    std::cin >> Nsteps;
    std::cout << "Save configurations yes/no (1 or 0): ";
    std::cin >> saveconf;
    std::cout << " " << std::endl;
    
    std::vector<double> Masses(Nm0);

    periodic_boundary(); //Compute right and left periodic boundary
    
	GaugeConf GConf = GaugeConf(LV::Nx, LV::Nt);  //Gauge configuration
 
    if (Nm0 == 1) {
        Masses = { m0_min };
    }
    else {
        Masses = linspace(m0_min, m0_max, Nm0);
    }

    std::ostringstream NameData;
    NameData << "2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt << "_Meas" << Nmeas << ".txt";
    std::ofstream Datfile;
    Datfile.open(NameData.str());
    Datfile << std::format("{:<30.17g}{:<30d}{:<30d}{:<30d}\n", beta, Ntherm, Nmeas, Nsteps);
    Datfile << std::format("{:<30.17g}{:<30d}\n", trajectory_length, MD_steps);
    for (double m0 : Masses) {
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
        std::cout << "**********************************************************************" << std::endl;


        HMC hmc = HMC(GConf,MD_steps, trajectory_length, Ntherm, Nmeas, Nsteps, beta, LV::Nx, LV::Nt, m0,saveconf);   
        double begin = omp_get_wtime();
        hmc.HMC_algorithm();
        double end = omp_get_wtime();
        std::cout << "Average plaquette value / volume: Ep = " << hmc.getEp() << " dEp = " << hmc.getdEp() << std::endl;
        std::cout << "Average gauge action / volume: gS = " << hmc.getgS() << " dgS = " << hmc.getdgS() << std::endl;
        std::cout << "Acceptance rate: " << hmc.getacceptance_rate() << std::endl;
        double elapsed_secs = end - begin;

        Datfile << std::format("{:<30.17g}\n", m0);
        Datfile << std::format("{:<30.17g}{:<30.17g}\n", hmc.getEp(), hmc.getdEp());
        Datfile << std::format("{:<30.17g}{:<30.17g}\n", hmc.getgS(), hmc.getdgS());
        Datfile << std::format("{:<30.17g}\n", hmc.getacceptance_rate());
        Datfile << std::format("{:<30.17g}", elapsed_secs);
        
 
        std::cout << "Time = " << elapsed_secs << " s" << std::endl;
        std::cout << "-------------------------------" << std::endl;

    }
    Datfile.close();
    
    //Free coordinate arrays
    free_lattice_arrays();

	return 0;
}

