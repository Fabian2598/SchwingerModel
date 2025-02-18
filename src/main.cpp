#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include "hmc.h"


int main() {
     srand(time(0));
    int Ntherm, Nmeas, Nsteps, Nbeta; //Simulation parameters
    double beta_min, beta_max; //Beta range
    double trajectory_length; //HMC parameters
    int MD_steps;
    double m0; //bare mass
    //---Input data---//
    std::cout << "----------------------------" << std::endl;
    std::cout << "|  Two-flavor Schwinger model   |" << std::endl;
    std::cout << "| Hybrid Monte Carlo simulation |" << std::endl;
    std::cout << "----------------------------" << std::endl;
    std::cout << "Ns " << Ns << " Nt " << Nt << std::endl;
    std::cout << "m0: ";
    std::cin >> m0;
    std::cout << "Molecular dynamics steps: ";
    std::cin >> MD_steps;
    std::cout << "Trajectory length: ";
    std::cin >> trajectory_length; 
    std::cout << "beta min: ";
    std::cin >> beta_min;
    std::cout << "beta max: ";
    std::cin >> beta_max;
    std::cout << "Number of betas: ";
    std::cin >> Nbeta;
    std::cout << "Thermalization: ";
    std::cin >> Ntherm;
    std::cout << "Measurements: ";
    std::cin >> Nmeas;
    std::cout << "Step (sweeps between measurements): ";
    std::cin >> Nsteps;
    std::cout << " " << std::endl;

    std::vector<double> Betas(Nbeta);
    initialize_matrices(); //Intialize gamma matrices, identity and unit vectors
	Coordinates(); //Compute vectorized coordinates
    periodic_boundary(); //Compute right and left periodic boundary
    int Ntot = Ns * Nt;
	GaugeConf GConf = GaugeConf(Ns, Nt);  //Gauge configuration
    if (Nbeta == 1) {
        Betas = { beta_min };
    }
    else {
        Betas = linspace(beta_min, beta_max, Nbeta);
    }
    char NameData[500], Data_str[500];
    sprintf(NameData, "2D_U1_Ns%d_Nt%d_Meas%d.txt", Ns, Nt, Nmeas);

    std::ofstream Datfile;
    Datfile.open(NameData);
    for (double beta : Betas) {
        std::cout << "beta = " << beta << std::endl;
        HMC hmc = HMC(GConf,MD_steps, trajectory_length, Ntherm, Nmeas, Nsteps, beta, Ns, Nt, Ntot, m0);   
        clock_t begin = clock();
        hmc.HMC_algorithm();
        clock_t end = clock();
        sprintf(Data_str, "%-30.17g%-30.17g%-30.17g\n", beta, hmc.getEp(), hmc.getdEp());
        std::cout << "Ep = " << hmc.getEp() << " dEp = " << hmc.getdEp() << std::endl;
        Datfile << Data_str;
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        std::cout << "Time = " << elapsed_secs << " s" << std::endl;
        std::cout << "------------------------------" << std::endl;

    }
    Datfile.close();

	return 0;
}
