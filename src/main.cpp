#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include "hmc.h"
#include <format>


int main() {
    srand(time(0));
    
    int Ntherm, Nmeas, Nsteps, Nm0; //Simulation parameters
    double beta; //Beta range
    double trajectory_length; //HMC parameters
    int MD_steps;
    double m0_min, m0_max; //bare mass
	int saveconf = 0; //Save configurations
    //---Input data---//
    
    std::cout << "----------------------------" << std::endl;
    std::cout << "|  Two-flavor Schwinger model   |" << std::endl;
    std::cout << "| Hybrid Monte Carlo simulation |" << std::endl;
    std::cout << "----------------------------" << std::endl;
    std::cout << "Nx " << LV::Nx << " Nt " << LV::Nt << std::endl;
    std::cout << "m0 min: ";
    std::cin >> m0_min;
    std::cout << "m0 max: ";
    std::cin >> m0_max;
    std::cout << "Number of masses in [m0_min, m0_max] ";
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
    
	Coordinates(); //Compute vectorized coordinates
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
        std::cout << "**********************************************************************" << std::endl;


        HMC hmc = HMC(GConf,MD_steps, trajectory_length, Ntherm, Nmeas, Nsteps, beta, LV::Nx, LV::Nt, m0,saveconf);   
        clock_t begin = clock();
        hmc.HMC_algorithm();
        clock_t end = clock();
        std::cout << "Average plaquette value / volume: Ep = " << hmc.getEp() << " dEp = " << hmc.getdEp() << std::endl;
        std::cout << "Average gauge action / volume: gS = " << hmc.getgS() << " dgS = " << hmc.getdgS() << std::endl;
        std::cout << "Acceptance rate: " << hmc.getacceptance_rate() << std::endl;
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

        Datfile << std::format("{:<30.17g}\n", m0);
        Datfile << std::format("{:<30.17g}{:<30.17g}\n", hmc.getEp(), hmc.getdEp());
        Datfile << std::format("{:<30.17g}{:<30.17g}\n", hmc.getgS(), hmc.getdgS());
        Datfile << std::format("{:<30.17g}", elapsed_secs);
        
 
        std::cout << "Time = " << elapsed_secs << " s" << std::endl;
        std::cout << "------------------------------" << std::endl;

    }
    Datfile.close();
    

	return 0;
}



/*

for(int at = 0; at<10000; at++){
    spinor right(LV::Ntot,c_vector(2,0));
    spinor left(LV::Ntot,c_vector(2,0));

    for(int i = 0; i< LV::Ntot; i++){
        right[i][0] = RandomU1();
        left[i][0] = RandomU1();
        right[i][1] = RandomU1();
        left[i][1] = RandomU1();
    }

    re_field X = phi_dag_partialD_phi(GConf.Conf,left,right);
    re_field Xold = phi_dag_partialD_phi_old(GConf.Conf,left,right);

    for(int i = 0; i < LV::Ntot; i++) {
        if(std::abs(X[i][0] - Xold[i][0]) > 1e-12 || std::abs(X[i][1] - Xold[i][1]) > 1e-12) {
            std::cout << "Error in the operator phi_dag_partialD_phi " << std::endl;
            std::cout << "n " << i << std::endl;
            return 1;
        }
    }
    
    }

    std::cout << "derivative working fine\n";
*/