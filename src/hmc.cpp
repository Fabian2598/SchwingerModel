#include "hmc.h"
#include <iomanip>

//Compute staple at coordinate (x,t) in the mu-direction
void HMC::StapleHMC(const std::vector<std::vector<std::complex<double>>>& U) {
    //x and t --> lattice sites. x: rows, t: columns
    // WARNING: Some references define the staple as the conjugate of this:
    //U_v(x) U_m(x+v) U*_v(x+m) + U*_v(x-v) U_m(x-v) U_v(x+m-v)
    //mu = 0 time direction, mu = 1 space direction
    for (int x = 0; x < Ns; x++) {
        for (int t = 0; t < Nt; t++) {
            //These coordinates could change depending on the conventions 
            int x1 = modulo(x - 1, Ns); 
            int x_1 = modulo(x + 1, Ns);
            int t1 = modulo(t - 1, Nt);
            int t_1 = modulo(t + 1, Nt);
            int i = Coords[x][t];
            for (int mu = 0; mu < 2; mu++) {
                if (mu == 0) {
                    const std::complex<double>& conf1 = U[i][1];
                    const std::complex<double>& conf2 = U[Coords[x1][t]][0];
                    const std::complex<double>& conf3 = U[Coords[x][t1]][1];
                    const std::complex<double>& conf4 = U[Coords[x_1][t]][1];
                    const std::complex<double>& conf5 = U[Coords[x_1][t]][0];
                    const std::complex<double>& conf6 = U[Coords[x_1][t1]][1];
                    staples[i][mu] = conf1 * conf2 * std::conj(conf3) +
                        std::conj(conf4) * conf5 * conf6;
                }
                else {
                    const std::complex<double>& conf1 = U[i][0];
                    const std::complex<double>& conf2 = U[Coords[x][t1]][1];
                    const std::complex<double>& conf3 = U[Coords[x1][t]][0];
                    const std::complex<double>& conf4 = U[Coords[x][t_1]][0];
                    const std::complex<double>& conf5 = U[Coords[x][t_1]][1];
                    const std::complex<double>& conf6 = U[Coords[x1][t_1]][0];
                    staples[i][mu] = conf1 * conf2 * std::conj(conf3) +
                        std::conj(conf4) * conf5 * conf6;
                }
            }
        }
    }
}

//Pure gauge force
void HMC::Force_G(const std::vector<std::vector<std::complex<double>>>& U) {
    StapleHMC(U); //Computes staples
    //NOTE: I HAVE TO CALL phi_dag_partialD_phi FIRST
	for (int x = 0; x < Ns; x++) {
		for (int t = 0; t < Nt; t++) {
			int i = Coords[x][t];
			for (int mu = 0; mu < 2; mu++) {
				Forces[i][mu] += -beta * std::imag(U[i][mu] * std::conj(staples[i][mu]));
			}
		}
	}
}

//Fermions force
//2* Re[ Psi^dagger partial D / partial omega(n) D Psi], where Psi = (DD^dagger)^(-1)phi, phi = D chi
void HMC::Force(const std::vector<std::vector<std::complex<double>>>& U,const std::vector<std::vector<std::complex<double>>>& phi) {
    std::vector<std::vector<std::complex<double>>> psi;
    psi = conjugate_gradient(U, phi, m0);  //(DD^dagger)^-1 phi
    Forces = phi_dag_partialD_phi(U,psi,D_dagger_phi(U, psi, m0)); //psi^dagger partial D / partial omega(n) D psi
    for (int x = 0; x < Ns; x++) {
		for (int t = 0; t < Nt; t++) {
			int i = Coords[x][t];
			for (int mu = 0; mu < 2; mu++) {
				Forces[i][mu] *= beta;
			}
		}
    }
    Force_G(U); //Gauge force 
}

//Generates new configuration [U,Pi]
//LATER REMOVE PHI AND CONF AS ARGUMENTS ...
void HMC::Leapfrog( 
const std::vector<std::vector<std::complex<double>>>& Conf,const std::vector<std::vector<std::complex<double>>>& phi ){
        double StepSize = trajectory_length / (MD_steps * 1.0);
        PConf = RandomMomentum(); //this will go in another function later
        char Name[500];
        sprintf(Name, "PConf.txt");
        SavePIConf(PConf, Name);
        PConf_copy = PConf;
        Conf_copy = Conf;
        std::complex<double> inumber(0.0, 1.0); //imaginary number

        //Conf_copy = Conf*exp(0.5i * StepSize * PConf_copy)
        for (int x = 0; x < Ns; x++) {
            for (int t = 0; t < Nt; t++) {
                int i = Coords[x][t];
                for (int mu = 0; mu < 2; mu++) {
                    Conf_copy[i][mu] = Conf_copy[i][mu] * exp(0.5 * inumber * StepSize * PConf_copy[i][mu]);
                }
            }
        }
		Force(Conf_copy,phi);
        for (int step = 1; step < MD_steps - 1; step++) {
            //PConf_copy += StepSize*force
            //Conf_copy *= exp(i * StepSize * PConf_copy)
            for (int x = 0; x < Ns; x++) {
                for (int t = 0; t < Nt; t++) {
                    int i = Coords[x][t];
                    for (int mu = 0; mu < 2; mu++) {
                        PConf_copy[i][mu] += StepSize *  Forces[i][mu];
                        Conf_copy[i][mu] *= exp(inumber * StepSize * PConf_copy[i][mu]);        
                    }
                }
            }
            Force(Conf_copy,phi);
        }
        //PConf_copy += StepSize*force
        //Conf_copy = Conf*exp(0.5i * StepSize* PConf_copy)
        for (int x = 0; x < Ns; x++) {
            for (int t = 0; t < Nt; t++) {
                int i = Coords[x][t];
                for (int mu = 0; mu < 2; mu++) {
                    PConf_copy[i][mu] += StepSize * Forces[i][mu];
                    Conf_copy[i][mu] *= exp(0.5 * inumber * StepSize * PConf_copy[i][mu]);
                    
                }
            }
        }
       // SaveConf(Conf_copy,  "testing_confCopy.txt");

}

/*
double HMC::DeltaH() {
    Compute_Plaquette01();
    double deltaH = 0.0;
    double dS = 0.0;
    for (int x = 0; x < Ns; x++) {
        for (int t = 0; t < Nt; t++) {
            dS += beta * std::real(Plaquette01[Coords[x][t]] - Plaquette01_prime[Coords[x][t]]);
            for (int mu = 0; mu < 2; mu++) {
                deltaH += 0.5 * (PConf_copy[Coords[x][t]][mu] * PConf_copy[Coords[x][t]][mu] - PConf[Coords[x][t]][mu] * PConf[Coords[x][t]][mu]);
            }  
        }
    }
    deltaH += dS;
    return deltaH;
}


void HMC::HMC_Update(const int& MD_steps, const double& trajectory_length){
    PConf = RandomMomentum();//random momentum conf sampled from a normal distribution
    Leapfrog(MD_steps, trajectory_length); //Evolve [Pi] and [U] (updates PConf_copy and Conf_copy)
    double deltaH = DeltaH(); //deltaH = Hamiltonian[U'][Pi'] - [U][Pi]
    double r = rand_range(0, 1);
    if (r<= exp(-deltaH)){
        //Accept the new configuration
        Conf = Conf_copy;
    }
    //Else configuration is not modified.
}

//HMC algorithm
void HMC::HMC_algorithm(const int& MD_steps, 
	const double& trajectory_length, const int& Ntherm, 
	const int& Nmeas, const int& Nsteps){

    std::vector<double> SpVector(Nmeas);
    initialization();
    for(int i = 0; i < Ntherm; i++) {HMC_Update(MD_steps,trajectory_length);} //Thermalization
    for(int i = 0; i < Nmeas; i++) {
        HMC_Update(MD_steps,trajectory_length);
        SpVector[i] = MeasureSp_HMC(); 
        for(int j = 0; j < Nsteps; j++) {HMC_Update(MD_steps,trajectory_length);} //Decorrelation
    }
    Ep = mean(SpVector) / (Ntot * 1.0); dEp = Jackknife_error(SpVector, 20) / (Ntot * 1.0);
} 

*/