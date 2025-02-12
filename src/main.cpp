#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include "gauge_conf.h"
#include "matrix_operations.h"
#include "conjugate_gradient.h"
#include "hmc.h"

int main() {
    srand(time(0));
	initialize_matrices(); //Intialize gamma matrices, identity and unit vectors
	Coordinates(); //Compute vectorized coordinates
    int Ntot = Ns * Nt;
    std::vector<std::vector<std::complex<double>>>Conf(Ntot, std::vector<std::complex<double>>(2, 0));
    std::vector<std::vector<std::complex<double>>> chi = RandomChi();
	//Test configuration//
    for (int i = 0; i < Ntot; i++) {
         for (int mu = 0; mu < 2; mu++) {
             Conf[i][mu] = RandomU1(); //Conf[Ns x Nt][mu in {0,1}]
         }
    }
    
    char Name[500];
    sprintf(Name, "Conf.txt");
    SaveConf(Conf, Name);
    sprintf(Name, "Chi.txt");
    SaveConf(chi, Name);
    double m0 = 1;
    std::vector<std::vector<std::complex<double>>> phi = D_phi(Conf, chi, m0);
    std::cout << "---------psi = Dchi---------" << std::endl;
    for (int x = 0; x < Ns; x++) {
        for (int t = 0; t < Nt; t++) {
            int n = Coords[x][t];
            std::cout << x << " " << t << " " << phi[n][0]
                << " " << phi[n][1] << std::endl;
        }
    }
    std::cout  << std::endl;

    std::vector<std::vector<std::complex<double>>> psi = conjugate_gradient(Conf, phi, m0);  //(DD^dagger)^-1 phi
    std::cout << "---------(DD^dagger)^-1 psi---------" << std::endl;
    for (int x = 0; x < Ns; x++) {
        for (int t = 0; t < Nt; t++) {
            int n = Coords[x][t];
            std::cout << x << " " << t << " " << psi[n][0]
                << " " << psi[n][1] << std::endl;
        }
    }
    std::cout  << std::endl;
    std::vector<std::vector<std::complex<double>>> D_dagger_psi = D_dagger_phi(Conf, psi, m0);
    std::cout << "---------D^dagger (DD^dagger)^-1 psi---------" << std::endl;
    for (int x = 0; x < Ns; x++) {
        for (int t = 0; t < Nt; t++) {
            int n = Coords[x][t];
            std::cout << x << " " << t << " " << D_dagger_psi[n][0]
                << " " << D_dagger_psi[n][1] << std::endl;
        }
    }
    std::cout  << std::endl;
    std::vector<std::vector<double>> Forces = phi_dag_partialD_phi(Conf,psi,D_dagger_psi); //psi^dagger partial D / partial omega(n) D psi
     
    for (int x = 0; x < Ns; x++) {
        for (int t = 0; t < Nt; t++) {
            int n = Coords[x][t];
            std::cout << x << " " << t << " " << Forces[n][0]
                << " " << Forces[n][1] << std::endl;
        }
    }
    
 
  
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
    std::cout << "-----------------------" << std::endl;
    std::cout << "|Pure Gauge U(1) theory|" << std::endl;
    std::cout << "-----------------------" << std::endl;
    std::cout << "Ns " << Ns << " Nt " << Nt << std::endl;
    std::cout << "----HMC----" << std::endl;
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
