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
    double m0 = 0.1; //bare mass
	int saveconf = 0; //Save configurations
    
    initialize_matrices(); //Intialize gamma matrices, identity and unit vectors
	Coordinates(); //Compute vectorized coordinates
    periodic_boundary(); //Compute right and left periodic boundary
    int Ntot = Ns * Nt;
	GaugeConf GConf = GaugeConf(Ns, Nt);  //Gauge configuration
    GConf.initialization();

    double beta = 1;
    std::cout << "beta = " << beta << std::endl;


    double start;
    double end;
    c_matrix chi = RandomChi();
    c_matrix phi;
    const int N = 100;
    double serialTime[N];
	double parallelTime[N];
    for (int i = 0; i < N; i++) {
        start = omp_get_wtime();
        phi = D_phi(GConf.Conf, chi, m0);
        end = omp_get_wtime();
		serialTime[i] = end - start;
        if (i == 10){
            std::cout << "phi[0][0]= " << phi[0][0] << "  phi[0][1]= " << phi[0][1] << std::endl;
        }
       

        start = omp_get_wtime();
        phi = D_phi_parallel(GConf.Conf, chi, m0);
        end = omp_get_wtime();
		parallelTime[i] = end - start;
        if (i == 10) {
            std::cout << "phi[0][0]= " << phi[0][0] << "  phi[0][1]= " << phi[0][1] << std::endl;
        }
    }

	double serialSum = 0;
	double parallelSum = 0;
	for (int i = 0; i < N; i++) {
		serialSum += serialTime[i];
		parallelSum += parallelTime[i];
	}
	double serialAvg = serialSum / N;
	double parallelAvg = parallelSum / N;
	std::cout << "Serial average time: " << serialAvg << " seconds" << std::endl;
	std::cout << "Parallel average time: " << parallelAvg << " seconds" << std::endl;
    


    


	return 0;
}
