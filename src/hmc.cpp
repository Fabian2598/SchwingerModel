#include "hmc.h"
#include <iomanip>
#include <string>
#include <sstream>

 void HMC::RandomPI() {

	static std::random_device rd;
	static std::default_random_engine generator(rd());
	static std::normal_distribution<double> distribution(0.0, 1.0); //mu, std
	
	for (int n = 0; n < Ntot; n++) {
		PConf[n][0] = distribution(generator);
		PConf[n][1] = distribution(generator);
	}

}

//Random Chi vector 
void HMC::RandomCHI() {
	static std::random_device rd;
	static std::default_random_engine generator(rd());
	static std::normal_distribution<double> distribution(0.0, 1/sqrt(2)); //mu, standard deviation

	for (int n = 0; n < Ntot; n++) {
		chi[n][0] = 1.0 * distribution(generator) + I_number * distribution(generator);
		chi[n][1] = 1.0 * distribution(generator) + I_number * distribution(generator);
	}
}

//Pure gauge force
//NOTE: phi_dag_partialD_phi HAS TO BE CALLED FIRST
void HMC::Force_G(GaugeConf& GConfig) {
    GConfig.Compute_Staple(); //Computes staples
	for (int n = 0; n < Ntot; n++) {
		Forces[n][0] += -beta * std::imag(GConfig.Conf[n][0] * std::conj(GConfig.Staples[n][0]));
        Forces[n][1] += -beta * std::imag(GConfig.Conf[n][1] * std::conj(GConfig.Staples[n][1]));
	}
		
}

//Fermions force
//2* Re[ Psi^dagger partial D / partial omega(n) D Psi], where Psi = (DD^dagger)^(-1)phi, phi = D chi
void HMC::Force(GaugeConf& GConfig,const spinor& phi) {
    spinor psi;
    psi = conjugate_gradient(GConfig.Conf, phi, m0);  //(DD^dagger)^-1 phi
    Forces = phi_dag_partialD_phi(GConfig.Conf,psi,D_dagger_phi(GConfig.Conf, psi, m0)); //psi^dagger partial D / partial omega(n) D psi
    Force_G(GConfig); //Gauge force 
}

//Generates new configuration [U,Pi]
void HMC::Leapfrog(const spinor& phi){
    double StepSize = trajectory_length / (MD_steps * 1.0);
    PConf_copy = PConf;
    c_double inumber(0.0, 1.0); //imaginary number
	GConf_copy = GConf; //Copy of the gauge configuration

    //Conf_copy = Conf*exp(0.5i * StepSize * PConf_copy)
    for (int n = 0; n < Ntot; n++) {
        GConf_copy.Conf[n][0] = GConf_copy.Conf[n][0] * exp(0.5 * inumber * StepSize * PConf_copy[n][0]);
        GConf_copy.Conf[n][1] = GConf_copy.Conf[n][1] * exp(0.5 * inumber * StepSize * PConf_copy[n][1]);   
    }

	Force(GConf_copy,phi); 

    for (int step = 1; step < MD_steps - 1; step++) {
        //PConf_copy += StepSize*force
        //Conf_copy *= exp(i * StepSize * PConf_copy)
        for (int n = 0; n < Ntot; n++) {
            //mu = 0
            PConf_copy[n][0] += StepSize *  Forces[n][0];
            GConf_copy.Conf[n][0] *= exp(inumber * StepSize * PConf_copy[n][0]);

            //mu = 1
            PConf_copy[n][1] += StepSize *  Forces[n][1];
            GConf_copy.Conf[n][1] *= exp(inumber * StepSize * PConf_copy[n][1]);
        }
        Force(GConf_copy,phi);
    }

    //PConf_copy += StepSize*force
    //Conf_copy = Conf*exp(0.5i * StepSize* PConf_copy)
    for (int n = 0; n < Ntot; n++) {
        //mu = 0
        PConf_copy[n][0] += StepSize * Forces[n][0];
        GConf_copy.Conf[n][0] *= exp(0.5 * inumber * StepSize * PConf_copy[n][0]);

        //mu = 1
        PConf_copy[n][1] += StepSize * Forces[n][1];
        GConf_copy.Conf[n][1] *= exp(0.5 * inumber * StepSize * PConf_copy[n][1]);
    }

}

double HMC::Action(GaugeConf& GConfig, const spinor& phi) {
    double action = 0;
    GConfig.Compute_Plaquette01();
    //Gauge contribution
	for (int i = 0; i < Ntot; i++) {
        action += beta * std::real(1.0-GConfig.Plaquette01[i]);
	}
    //Fermions contribution
    //Phi^dagger (DD^dagger)^-1 Phi = dot(Phi,(DD^dagger)^-1 Phi) (the dot function takes into account the dagger)
	action += std::real( dot( conjugate_gradient(GConfig.Conf, phi, m0), phi)); 
    return action;
}

double HMC::Hamiltonian(GaugeConf& GConfig, const re_field& Pi,const spinor& phi) {
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
   
	//PConf = RandomMomentum(); //random momentum conf sampled from a normal distribution
    RandomPI(); 
    //pseudofermions phi = D chi, where chi is normaly sampled
    //spinor chi = RandomChi();
    RandomCHI();


    spinor phi = D_phi(GConf.Conf, chi, m0);
    Leapfrog(phi); //Evolve [Pi] and [U] 
    double deltaH = Hamiltonian(GConf_copy, PConf_copy, phi) - Hamiltonian(GConf, PConf, phi); //deltaH = Hamiltonian[U'][Pi'] - [U][Pi]
    double r = rand_range(0, 1);
    if (r <= exp(-deltaH)) {
        //Accept the new configuration
        GConf = GConf_copy;
        if (therm == true) {
            acceptance_rate += 1.0;
        }
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
    re_vector SpVector(Nmeas);
    re_vector gAction(Nmeas);
	GConf.initialization(); //Initialize the gauge configuration
    for(int i = 0; i < Ntherm; i++) {HMC_Update();} //Thermalization
    therm = true; //Set the flag to true
    for(int i = 0; i < Nmeas; i++) {
        HMC_Update();
        SpVector[i] = GConf.MeasureSp_HMC(); //Plaquettes are computed when the action is called
        gAction[i] = GConf.Compute_gaugeAction(beta); //Gauge action
		if (saveconf == 1) {
			std::ostringstream NameData;
                NameData << "2D_U1_Ns" << Nx << "_Nt" << Nt
                << "_b" << format(beta)
                << "_m" << format(m0)
                << "_" << i << ".ctxt";
            SaveConf(GConf.Conf, NameData.str());
		}
		for (int j = 0; j < Nsteps; j++) { HMC_Update(); } //Decorrelation
    }
    Ep = mean(SpVector) / (Ntot * 1.0); dEp = Jackknife_error(SpVector, 20) / (Ntot * 1.0); //Average Plaquette Value
    gS = mean(gAction) / (Ntot * 1.0); dgS = Jackknife_error(gAction, 20) / (Ntot * 1.0);
} 




