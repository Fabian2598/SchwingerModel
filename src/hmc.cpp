#include "hmc.h"
#include <string>
#include <sstream>

 void HMC::RandomPI() {

	static std::random_device rd;
	static std::default_random_engine generator(rd());
	std::normal_distribution<double> distribution(0.0, 1.0); //mu, std

    int n;
    for(int x = 1; x<=mpi::width_x; x++){
		for(int t = 1; t<=mpi::width_t; t++){
			n = x*(mpi::width_t+2)+t;
            PConf.val[2*n]   = distribution(generator);
		    PConf.val[2*n+1] = distribution(generator);
        }
    }
}

//Random Chi vector 
void HMC::RandomCHI() {
	static std::random_device rd;
	static std::default_random_engine generator(rd());
	std::normal_distribution<double> distribution(0.0, 1/sqrt(2)); //mu, standard deviation

    int n;
    for(int x = 1; x<=mpi::width_x; x++){
		for(int t = 1; t<=mpi::width_t; t++){
            n = x*(mpi::width_t+2)+t;
            chi.val[2*n]   = 1.0 * distribution(generator) + I_number * distribution(generator);
		    chi.val[2*n+1] = 1.0 * distribution(generator) + I_number * distribution(generator);
        }
    }
}

//Pure gauge force
//NOTE: phi_dag_partialD_phi HAS TO BE CALLED FIRST
void HMC::Force_G(GaugeConf& GConfig) {
    GConfig.Compute_Staple(); //Computes staples

	int n;
    for(int x = 1; x<=mpi::width_x; x++){
		for(int t = 1; t<=mpi::width_t; t++){
            n = x*(mpi::width_t+2)+t;
		    Forces.val[2*n]   += -beta * std::imag(GConfig.Conf.val[2*n]   * std::conj(GConfig.Staples.val[2*n]));
            Forces.val[2*n+1] += -beta * std::imag(GConfig.Conf.val[2*n+1] * std::conj(GConfig.Staples.val[2*n+1]));
        }
	}
		
}

//Fermions force
//2* Re[ Psi^dagger partial D / partial omega(n) D Psi], where Psi = (DD^dagger)^(-1)phi, phi = D chi
void HMC::Force(GaugeConf& GConfig,const spinor& phi) {
    spinor psi(mpi::maxSizeH); 
    CG_convergence = conjugate_gradient(GConfig.Conf, phi,psi);  //(DD^dagger)^-1 phi
    //Save gauge configuration if CG does not converge
    if (CG_convergence == 0){
        std::ostringstream NameData;
        NameData << "2D_U1_" << Nx << "x" << Nt
                 << "_b" << format(beta)
                 << "_m" << format(m0)
                 << "_illConf" << illConfId << ".ctxt";
        GConf.SaveConf(NameData.str());
        illConfId += 1;
    } 
    D_dagger_phi(GConfig.Conf, psi,TEMP);
    phi_dag_partialD_phi(GConfig.Conf,psi,TEMP,Forces); //psi^dagger partial D / partial omega(n) D psi
    Force_G(GConfig); //Gauge force 
}

//Generates new configuration [U,Pi]
void HMC::Leapfrog(const spinor& phi){
    double StepSize = trajectory_length / (MD_steps * 1.0);
    PConf_copy = PConf;
	GConf_copy = GConf; //Copy of the gauge configuration
    //Conf_copy = Conf*exp(0.5i * StepSize * PConf_copy)
    int n;
    for(int x = 1; x<=mpi::width_x; x++){
		for(int t = 1; t<=mpi::width_t; t++){
            n = x*(mpi::width_t+2)+t;
            GConf_copy.Conf.val[2*n]   = GConf_copy.Conf.val[2*n]   * exp(0.5 * I_number * StepSize * PConf_copy.val[2*n]);
            GConf_copy.Conf.val[2*n+1] = GConf_copy.Conf.val[2*n+1] * exp(0.5 * I_number * StepSize * PConf_copy.val[2*n+1]);   
        }
    }

	Force(GConf_copy,phi); 

    for (int step = 1; step < MD_steps - 1; step++) {
        //PConf_copy += StepSize*force
        //Conf_copy *= exp(i * StepSize * PConf_copy)
        for(int x = 1; x<=mpi::width_x; x++){
		    for(int t = 1; t<=mpi::width_t; t++){
                n = x*(mpi::width_t+2)+t;
                //mu = 0
                PConf_copy.val[2*n] += StepSize *  Forces.val[2*n];
                GConf_copy.Conf.val[2*n] *= exp(I_number * StepSize * PConf_copy.val[2*n]);

                //mu = 1
                PConf_copy.val[2*n+1] += StepSize *  Forces.val[2*n+1];
                GConf_copy.Conf.val[2*n+1] *= exp(I_number * StepSize * PConf_copy.val[2*n+1]);
            }
        }
        Force(GConf_copy,phi);
    }

    //PConf_copy += StepSize*force
    //Conf_copy = Conf*exp(0.5i * StepSize* PConf_copy)
    for(int x = 1; x<=mpi::width_x; x++){
		for(int t = 1; t<=mpi::width_t; t++){
            n = x*(mpi::width_t+2)+t;
            //mu = 0
            PConf_copy.val[2*n] += StepSize * Forces.val[2*n];
            GConf_copy.Conf.val[2*n] *= exp(0.5 * I_number * StepSize * PConf_copy.val[2*n]);

            //mu = 1
            PConf_copy.val[2*n+1] += StepSize * Forces.val[2*n+1];
            GConf_copy.Conf.val[2*n+1] *= exp(0.5 * I_number * StepSize * PConf_copy.val[2*n+1]);
        }
    }
}

double HMC::Action(GaugeConf& GConfig, const spinor& phi) {
    double local_action = 0.0;
    double action;
    GConfig.Compute_Plaquette01();
    //Gauge contribution
    int n;
	for(int x = 1; x<=mpi::width_x; x++){
	    for(int t = 1; t<=mpi::width_t; t++){
            n = x*(mpi::width_t+2)+t;
            local_action += beta * std::real(1.0-GConfig.Plaquette01[n]);
        }
	}

    MPI_Allreduce(&local_action, &action, 1, MPI_DOUBLE, MPI_SUM, mpi::cart_comm);
    //Fermions contribution
    //Phi^dagger (DD^dagger)^-1 Phi = dot(Phi,(DD^dagger)^-1 Phi) (the dot function takes into account the dagger)
    CG_convergence = conjugate_gradient(GConfig.Conf, phi,TEMP);
    action += std::real( dot( TEMP, phi)); 

    //Save gauge configuration if CG does not converge
    /*
    if (CG_convergence == 0){
        std::ostringstream NameData;
        NameData << "2D_U1_" << Nx << "x" << Nt
                 << "_b" << format(beta)
                 << "_m" << format(m0)
                 << "_illConf" << illConfId << ".ctxt";
        GConf.SaveConf(NameData.str());
        illConfId += 1;
    } 
    */
   
    return action;
}

double HMC::Hamiltonian(GaugeConf& GConfig, const re_field& Pi,const spinor& phi) {
    double local_H = 0;
    //Momentum contribution
    int n;
    for(int x = 1; x<=mpi::width_x; x++){
	    for(int t = 1; t<=mpi::width_t; t++){
            n = x*(mpi::width_t+2)+t;
            local_H += 0.5 * Pi.val[2*n] * Pi.val[2*n];
		    local_H += 0.5 * Pi.val[2*n+1] * Pi.val[2*n+1];
        }
    }
    double H;
    MPI_Allreduce(&local_H, &H, 1, MPI_DOUBLE, MPI_SUM, mpi::cart_comm);
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

    spinor phi(mpi::maxSizeH);
    D_phi(GConf.Conf, chi,phi);
    Leapfrog(phi); //Evolve [Pi] and [U] 
    double deltaH = Hamiltonian(GConf_copy, PConf_copy, phi) - Hamiltonian(GConf, PConf, phi); //deltaH = Hamiltonian[U'][Pi'] - [U][Pi]
    double r;
    
    //Same random number for all ranks
    if (mpi::rank == 0)
        r = rand_range(0, 1); 

    MPI_Bcast(&r, 1, MPI_DOUBLE,  0, mpi::cart_comm);

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
            GConf.SaveConf(NameData.str());
		}
		if (i != Nmeas-1){
            for (int j = 0; j < Nsteps; j++) { conf_i+=1; HMC_Update(); } //Decorrelation
        }
    }
    Ep = mean(SpVector) / (Ntot * 1.0); dEp = Jackknife_error(SpVector, 20) / (Ntot * 1.0); //Average Plaquette Value
    gS = mean(gAction) / (Ntot * 1.0); dgS = Jackknife_error(gAction, 20) / (Ntot * 1.0);
} 




