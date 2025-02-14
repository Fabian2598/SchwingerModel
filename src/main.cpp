#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include "gauge_conf.h"
#include "matrix_operations.h"
#include "conjugate_gradient.h"
#include "hmc.h"

int main() {
    srand(0);//time(0);
	initialize_matrices(); //Intialize gamma matrices, identity and unit vectors
	Coordinates(); //Compute vectorized coordinates
    int Ntot = Ns * Nt;
	GaugeConf GConf = GaugeConf(Ns, Nt);  //Gauge configuration
	GConf.initialization(); //Random initialization of the gauge configuration
    double m0 = 1, beta = 1; //bare mass and beta
    
    int MD_steps = 10;
    double trajectory_length = 1;
    HMC hmc = HMC(GConf,MD_steps, trajectory_length, 1, 1, 1, beta, Ns, Nt, Ntot, m0);   
    hmc.HMC_Update(); 
   
    return 0;
}

/*
int main() {
    srand(time(0));
    int Ntherm, Nmeas, Nsteps, Nbeta;
    double beta_min, beta_max;
    double trajectory_length;
    int MD_steps;
    //---Input data---//
    std::cout << "----------------------------" << std::endl;
    std::cout << "|  Two-flavor Schwinger model   |" << std::endl;
    std::cout << "| Hybrid Monte Carlo simulation |" << std::endl;
    std::cout << "----------------------------" << std::endl;
    std::cout << "Ns " << Ns << " Nt " << Nt << std::endl;
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
    GaugeConf Configuration = GaugeConf(Ns, Nt);

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
	Configuration.Coordinates(); //Compute vectorized coordinates
    for (double beta : Betas) {
        clock_t begin = clock();
        std::cout << "beta = " << beta << std::endl;
        Configuration.setBeta(beta);
        Configuration.HMC(MD_steps, trajectory_length, Ntherm, Nmeas, Nsteps);
        sprintf(Data_str, "%-30.17g%-30.17g%-30.17g\n", beta, Configuration.getEp(), Configuration.getdEp());
        std::cout << "Ep = " << Configuration.getEp() << " dEp = " << Configuration.getdEp() << std::endl;
        Datfile << Data_str;
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        std::cout << "Time = " << elapsed_secs << " s" << std::endl;
        std::cout << "------------------------------" << std::endl;

    }
    Datfile.close();

	return 0;
}

*/
