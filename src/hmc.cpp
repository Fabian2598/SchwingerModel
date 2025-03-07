#include "hmc.h"
#include <iomanip>
#include <string>
#include <sstream>



//Pure gauge force
//NOTE: phi_dag_partialD_phi HAS TO BE CALLED FIRST
void HMC::Force_G(GaugeConf& GConfig) {
    GConfig.Compute_Staple(); //Computes staples
	for (int x = 0; x < Ns; x++) {
		for (int t = 0; t < Nt; t++) {
			int i = Coords[x][t];
			for (int mu = 0; mu < 2; mu++) {
				Forces[i][mu] += -beta * std::imag(GConfig.Conf[i][mu] * std::conj(GConfig.Staples[i][mu]));
			}
		}
	}
}

//Fermions force
//2* Re[ Psi^dagger partial D / partial omega(n) D Psi], where Psi = (DD^dagger)^(-1)phi, phi = D chi
void HMC::Force(GaugeConf& GConfig,const c_matrix& phi) {
    c_matrix psi;
    psi = conjugate_gradient(GConfig.Conf, phi, m0);  //(DD^dagger)^-1 phi
    Forces = phi_dag_partialD_phi(GConfig.Conf,psi,D_dagger_phi(GConfig.Conf, psi, m0)); //psi^dagger partial D / partial omega(n) D psi
    for (int x = 0; x < Ns; x++) {
		for (int t = 0; t < Nt; t++) {
			int i = Coords[x][t];
			for (int mu = 0; mu < 2; mu++) {
				Forces[i][mu] *= beta; 
			}
		}
    }
    Force_G(GConfig); //Gauge force 
}

//Generates new configuration [U,Pi]
void HMC::Leapfrog(const c_matrix& phi){
        double StepSize = trajectory_length / (MD_steps * 1.0);
        PConf_copy = PConf;
        c_double inumber(0.0, 1.0); //imaginary number
		GConf_copy = GConf; //Copy of the gauge configuration
        //Conf_copy = Conf*exp(0.5i * StepSize * PConf_copy)
        for (int x = 0; x < Ns; x++) {
            for (int t = 0; t < Nt; t++) {
                int i = Coords[x][t];
                for (int mu = 0; mu < 2; mu++) {
                    GConf_copy.Conf[i][mu] = GConf_copy.Conf[i][mu] * exp(0.5 * inumber * StepSize * PConf_copy[i][mu]);
                }
            }
        }
		Force(GConf_copy,phi); 
        for (int step = 1; step < MD_steps - 1; step++) {
            //PConf_copy += StepSize*force
            //Conf_copy *= exp(i * StepSize * PConf_copy)
            for (int x = 0; x < Ns; x++) {
                for (int t = 0; t < Nt; t++) {
                    int i = Coords[x][t];
                    for (int mu = 0; mu < 2; mu++) {
                        PConf_copy[i][mu] += StepSize *  Forces[i][mu];
                        GConf_copy.Conf[i][mu] *= exp(inumber * StepSize * PConf_copy[i][mu]);
                    }
                }
            }
            Force(GConf_copy,phi);
        }
        //PConf_copy += StepSize*force
        //Conf_copy = Conf*exp(0.5i * StepSize* PConf_copy)
        for (int x = 0; x < Ns; x++) {
            for (int t = 0; t < Nt; t++) {
                int i = Coords[x][t];
                for (int mu = 0; mu < 2; mu++) {
                    PConf_copy[i][mu] += StepSize * Forces[i][mu];
                    GConf_copy.Conf[i][mu] *= exp(0.5 * inumber * StepSize * PConf_copy[i][mu]);
                    
                }
            }
        }

}

double HMC::Action(GaugeConf& GConfig, const c_matrix& phi) {
    double action = 0;
    c_matrix phi_dagg(Ntot, c_vector(2, 0)); //phi^dagger
    GConfig.Compute_Plaquette01();
    //Gauge contribution
	for (int i = 0; i < Ntot; i++) {
        action += beta * std::real(1.0-GConfig.Plaquette01[i]);
		phi_dagg[i][0] = std::conj(phi[i][0]);
        phi_dagg[i][1] = std::conj(phi[i][1]);
	}
    //Fermions contribution
    //Phi^dagger (DD^dagger)^-1 Phi = dot(Phi,(DD^dagger)^-1 Phi) (the dot function takes into account the dagger)
	action += beta*std::real( dot(phi, conjugate_gradient(GConfig.Conf, phi, m0))); 
    return action;
}

double HMC::Hamiltonian(GaugeConf& GConfig, const std::vector<std::vector<double>>& Pi,const c_matrix& phi) {
    double H = 0;
    //Momentum contribution
    for (int i = 0; i < Ntot; i++) {
        for (int mu = 0; mu < 2; mu++) {
			H += 0.5 * Pi[i][mu] * Pi[i][mu];
        }
    }
    //Action contribution
    H += Action(GConfig,phi);

    return H;
}

void HMC::HMC_Update() {
	PConf = RandomMomentum(); //random momentum conf sampled from a normal distribution
    //pseudofermions phi = D chi, where chi is normaly sampled
    c_matrix chi = RandomChi();
    c_matrix phi = D_phi(GConf.Conf, chi, m0);
    Leapfrog(phi); //Evolve [Pi] and [U] 
    double deltaH = Hamiltonian(GConf_copy, PConf_copy, phi) - Hamiltonian(GConf, PConf, phi); //deltaH = Hamiltonian[U'][Pi'] - [U][Pi]
    double r = rand_range(0, 1);
    if (r <= exp(-deltaH)) {
        //Accept the new configuration
        GConf = GConf_copy;
    }
    //Else configuration is not modified.
}

//Formats decimal numbers
static std::string format(const double& number) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << number;
    std::string str = oss.str();
    str.erase(str.find('.'), 1); //Removes decimal dot
    return str;
}

void HMC::HMC_algorithm(){
    std::vector<double> SpVector(Nmeas);
	GConf.initialization(); //Initialize the gauge configuration
    for(int i = 0; i < Ntherm; i++) {HMC_Update();} //Thermalization
    for(int i = 0; i < Nmeas; i++) {
        HMC_Update();
        SpVector[i] = GConf.MeasureSp_HMC(); //Plaquettes are computed when the action is called
		if (saveconf == 1) {
			char NameData[500];
            sprintf(NameData, "2D_U1_Ns%d_Nt%d_b%s_m%s_%d.txt", Ns, Nt, format(beta).c_str(), format(m0).c_str(), i); 
            SaveConf(GConf.Conf, NameData);
		}
		for (int j = 0; j < Nsteps; j++) { HMC_Update(); } //Decorrelation
    }
    Ep = mean(SpVector) / (Ntot * 1.0); dEp = Jackknife_error(SpVector, 20) / (Ntot * 1.0); //Average Plaquette Value
} 




