#include "hmc.h"
#include <iomanip>
#include <string>
#include <sstream>

 void HMC::RandomPI() {

	static std::random_device rd;
	static std::default_random_engine generator(rd());
	std::normal_distribution<double> distribution(0.0, 1.0); //mu, std

	#pragma omp parallel for
	for (int n = 0; n < mpi::maxSize; n++) {
		PConf.mu0[n] = distribution(generator);
		PConf.mu1[n] = distribution(generator);
	}

}

//Random Chi vector 
void HMC::RandomCHI() {
	static std::random_device rd;
	static std::default_random_engine generator(rd());
	std::normal_distribution<double> distribution(0.0, 1/sqrt(2)); //mu, standard deviation

    #pragma omp parallel for
	for (int n = 0; n < mpi::maxSize; n++) {
		chi.mu0[n] = 1.0 * distribution(generator) + I_number * distribution(generator);
		chi.mu1[n] = 1.0 * distribution(generator) + I_number * distribution(generator);
	}
}

//Pure gauge force
//NOTE: phi_dag_partialD_phi HAS TO BE CALLED FIRST
void HMC::Force_G(GaugeConf& GConfig) {
    GConfig.Compute_Staple(); //Computes staples

    #pragma omp parallel for
	for (int n = 0; n < mpi::maxSize; n++) {
		Forces.mu0[n] += -beta * std::imag(GConfig.Conf.mu0[n] * std::conj(GConfig.Staples.mu0[n]));
        Forces.mu1[n] += -beta * std::imag(GConfig.Conf.mu1[n] * std::conj(GConfig.Staples.mu1[n]));
	}
		
}

//Fermions force
//2* Re[ Psi^dagger partial D / partial omega(n) D Psi], where Psi = (DD^dagger)^(-1)phi, phi = D chi
void HMC::Force(GaugeConf& GConfig,const spinor& phi) {
    spinor psi(mpi::maxSize); 
    conjugate_gradient(GConfig.Conf, phi,psi, m0);  //(DD^dagger)^-1 phi
    D_dagger_phi(GConfig.Conf, psi,TEMP, m0);
    Forces = phi_dag_partialD_phi(GConfig.Conf,psi,TEMP); //psi^dagger partial D / partial omega(n) D psi
    Force_G(GConfig); //Gauge force 
}

//Generates new configuration [U,Pi]
void HMC::Leapfrog(const spinor& phi){
    double StepSize = trajectory_length / (MD_steps * 1.0);
    PConf_copy = PConf;
    c_double inumber(0.0, 1.0); //imaginary number
	GConf_copy = GConf; //Copy of the gauge configuration
    //Conf_copy = Conf*exp(0.5i * StepSize * PConf_copy)
    #pragma omp parallel for
    for (int n = 0; n < mpi::maxSize; n++) {
        GConf_copy.Conf.mu0[n] = GConf_copy.Conf.mu0[n] * exp(0.5 * inumber * StepSize * PConf_copy.mu0[n]);
        GConf_copy.Conf.mu1[n] = GConf_copy.Conf.mu1[n] * exp(0.5 * inumber * StepSize * PConf_copy.mu1[n]);   
    }

	Force(GConf_copy,phi); 

    for (int step = 1; step < MD_steps - 1; step++) {
        //PConf_copy += StepSize*force
        //Conf_copy *= exp(i * StepSize * PConf_copy)
        #pragma omp parallel for
        for (int n = 0; n < mpi::maxSize; n++) {
            //mu = 0
            PConf_copy.mu0[n] += StepSize *  Forces.mu0[n];
            GConf_copy.Conf.mu0[n] *= exp(inumber * StepSize * PConf_copy.mu0[n]);

            //mu = 1
            PConf_copy.mu1[n] += StepSize *  Forces.mu1[n];
            GConf_copy.Conf.mu1[n] *= exp(inumber * StepSize * PConf_copy.mu1[n]);
        }
        Force(GConf_copy,phi);
    }

    //PConf_copy += StepSize*force
    //Conf_copy = Conf*exp(0.5i * StepSize* PConf_copy)
    #pragma omp parallel for
    for (int n = 0; n < mpi::maxSize; n++) {
        //mu = 0
        PConf_copy.mu0[n] += StepSize * Forces.mu0[n];
        GConf_copy.Conf.mu0[n] *= exp(0.5 * inumber * StepSize * PConf_copy.mu0[n]);

        //mu = 1
        PConf_copy.mu1[n] += StepSize * Forces.mu1[n];
        GConf_copy.Conf.mu1[n] *= exp(0.5 * inumber * StepSize * PConf_copy.mu1[n]);
    }

}

double HMC::Action(GaugeConf& GConfig, const spinor& phi) {
    double local_action = 0.0;
    double action;
    GConfig.Compute_Plaquette01();
    //Gauge contribution
    #pragma omp parallel for reduction(+:local_action)
	for (int n = 0; n < mpi::maxSize; n++) {
        local_action += beta * std::real(1.0-GConfig.Plaquette01[n]);
	}
    MPI_Allreduce(&local_action, &action, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //Fermions contribution
    //Phi^dagger (DD^dagger)^-1 Phi = dot(Phi,(DD^dagger)^-1 Phi) (the dot function takes into account the dagger)
    conjugate_gradient(GConfig.Conf, phi,TEMP, m0);
    action += std::real( dot( TEMP, phi)); 
    return action;
}

double HMC::Hamiltonian(GaugeConf& GConfig, const re_field& Pi,const spinor& phi) {
    double local_H = 0;
    //Momentum contribution
    #pragma omp parallel for reduction(+:local_H)
    for (int n = 0; n < mpi::maxSize; n++) {
        local_H += 0.5 * Pi.mu0[n] * Pi.mu0[n];
		local_H += 0.5 * Pi.mu1[n] * Pi.mu1[n];
        
    }
    double H;
    MPI_Allreduce(&local_H, &H, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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

    spinor phi(mpi::maxSize);
    D_phi(GConf.Conf, chi,phi, m0);
    Leapfrog(phi); //Evolve [Pi] and [U] 
    double deltaH = Hamiltonian(GConf_copy, PConf_copy, phi) - Hamiltonian(GConf, PConf, phi); //deltaH = Hamiltonian[U'][Pi'] - [U][Pi]
    double r;

    if (mpi::rank == 0)
        r = rand_range(0, 1); 

    MPI_Bcast(&r, 1, MPI_DOUBLE,  0, MPI_COMM_WORLD);

    if (r <= exp(-deltaH)) {
        //Accept the new configuration
        GConf = GConf_copy;
        if (therm == true) {
            acceptance_rate += 1.0;    
        }
    }
    //Else configuration is not modified.
    //if (therm == true && mpi::rank==0)
        //std::cout << "Conf number " << conf_i << " acceptance rate " << getacceptance_rate(conf_i) << std::endl;
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
    std::vector<double> gAction(Nmeas);
	GConf.initialization(); //Initialize the gauge configuration randomly
    for(int i = 0; i < Ntherm; i++) {
        HMC_Update();
        if (i%100 == 0 && mpi::rank == 0)
            std::cout << "Conf " << i << " out of " << Ntherm << " for thermalization" << std::endl;
    } //Thermalization
    therm = true; //Set the flag to true
    if (mpi::rank == 0)
        std::cout << "Thermalization done" <<std::endl; 
    conf_i = 0;
    for(int i = 0; i < Nmeas; i++) {
        conf_i += 1;
        HMC_Update();
        SpVector[i] = GConf.MeasureSp_HMC(); //Plaquettes are computed when the action is called
        gAction[i] = GConf.Compute_gaugeAction(beta); //Gauge action
		if (saveconf == 1) {
			std::ostringstream NameData;
            NameData << "2D_U1_Ns" << Nx << "_Nt" << Nt
                << "_b" << format(beta)
                << "_m" << format(m0)
                << "_" << i << ".ctxt";
            SaveConf(GConf, NameData.str());
		}
		for (int j = 0; j < Nsteps; j++) { conf_i+=1; HMC_Update(); } //Decorrelation
    }
    Ep = mean(SpVector) / (Ntot * 1.0); dEp = Jackknife_error(SpVector, 20) / (Ntot * 1.0); //Average Plaquette Value
    gS = mean(gAction) / (Ntot * 1.0); dgS = Jackknife_error(gAction, 20) / (Ntot * 1.0);
} 




