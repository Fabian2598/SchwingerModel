#include "hmc.h"
#include <string>
#include <sstream>

/*
    The loops over the lattice are executed by for_each_site / lattice_sum (GPU kernels when compiled with CUDA).
    The kernels cannot access class members, so the pointers to the data of the fields (field.val) and the
    parameters they need are copied to local variables first (check backend.h).
*/

 void HMC::RandomPI() {

	static std::random_device rd;
	static std::default_random_engine generator(rd());
	std::normal_distribution<double> distribution(0.0, 1.0); //mu, std

    //The momenta are generated on the CPU and then copied to the GPU
    std::vector<double> P(PConf.size);
    int n;
    for(int x = 1; x<=mpi::width_x; x++){
		for(int t = 1; t<=mpi::width_t; t++){
			n = x*(mpi::width_t+2)+t;
            P[2*n]   = distribution(generator);
		    P[2*n+1] = distribution(generator);
        }
    }
    PConf.fromHost(P);
}

//Random Chi vector
void HMC::RandomCHI() {
	static std::random_device rd;
	static std::default_random_engine generator(rd());
	std::normal_distribution<double> distribution(0.0, 1/sqrt(2)); //mu, standard deviation

    //chi is generated on the CPU and then copied to the GPU
    std::vector<c_double> chi_host(chi.size);
    int n;
    for(int x = 1; x<=mpi::width_x; x++){
		for(int t = 1; t<=mpi::width_t; t++){
            n = x*(mpi::width_t+2)+t;
            chi_host[2*n]   = 1.0 * distribution(generator) + I_number * distribution(generator);
		    chi_host[2*n+1] = 1.0 * distribution(generator) + I_number * distribution(generator);
        }
    }
    chi.fromHost(chi_host);
}

//Pure gauge force
//NOTE: phi_dag_partialD_phi HAS TO BE CALLED FIRST
void HMC::Force_G(GaugeConf& GConfig) {
    GConfig.Compute_Staple(); //Computes staples

    double* F = Forces.val;
    const c_double* U = GConfig.Conf.val;
    const c_double* staple = GConfig.Staples.val;
    const double b = beta;
    for_each_site([=] HOST_DEVICE (int n){
        F[2*n]   += -b * imag(U[2*n]   * conj(staple[2*n]));
        F[2*n+1] += -b * imag(U[2*n+1] * conj(staple[2*n+1]));
	});

}

//Fermions force
//2* Re[ Psi^dagger partial D / partial omega(n) D Psi], where Psi = (DD^dagger)^(-1)phi, phi = D chi
void HMC::Force(GaugeConf& GConfig,const spinor& phi) {
    spinor psi(mpi::maxSizeH);
    exchange_halo_gauge(GConfig.Conf.val); //the corners are needed for the staples and the clover term
    #ifdef CLOVER
        GConfig.Compute_Q();
    #endif
    CG_convergence = conjugate_gradient(GConfig, phi,psi);  //(DD^dagger)^-1 phi
    //Save gauge configuration if CG does not converge
    if (CG_convergence == 0){
        std::ostringstream NameData;
        NameData << "2D_U1_" << Nx << "x" << Nt
                 << "_b" << format(beta)
                 << "_m" << format(m0)
                 << "_illConf" << illConfId << ".ctxt";
        //GConf.SaveConf(NameData.str());
        illConfId += 1;
    }

    D_dagger_phi(GConfig, psi,TEMP);
    phi_dag_partialD_phi(GConfig,psi,TEMP,Forces); //psi^dagger partial D / partial omega(n) D psi
    Force_G(GConfig); //Gauge force
    #ifdef CLOVER
        Force_Clover(GConfig,psi,TEMP); //update clover force
    #endif
}

void HMC::Force_Clover(const GaugeConf& GConf,const spinor& left_in, const spinor& right_in){
    const c_double factor = I_number * l_8 * sim_params::csw;
    const LocalLattice lat = local_lattice();
    const c_double* left_term = left_in.val;
    const c_double* right_term = right_in.val;
    const c_double* P1 = GConf.P1.val;
    const c_double* P2 = GConf.P2.val;
    const c_double* P3 = GConf.P3.val;
    const c_double* P4 = GConf.P4.val;
    c_double* j1 = J1.val;
    c_double* j2 = J2.val;
    c_double* j3 = J3.val;
    c_double* j4 = J4.val;
    double* F = Forces.val;

    for_each_site([=] HOST_DEVICE (int n){
        const c_double f = factor * (-conj(left_term[2*n]) * right_term[2*n]
                            + conj(left_term[2*n+1]) * right_term[2*n+1]);
        j1[n] = f * (P1[n] + conj(P1[n]));
        j2[n] = f * (P2[n] + conj(P2[n]));
        j3[n] = f * (P3[n] + conj(P3[n]));
        j4[n] = f * (P4[n] + conj(P4[n]));
    });

    //Halo exchange including corners
    exchange_halo_vec(j1);
    exchange_halo_vec(j2);
    exchange_halo_vec(j3);
    exchange_halo_vec(j4);

    for_each_site([=] HOST_DEVICE (int n){
        int right, down, left, up;
        double lsign, rsign;
        int x1_t_1, x_1_t_1, x_1_t1, x1_t1; //n-0+1, n-0-1, n+0-1, n+0+1
        get_neighbors(lat, n, right, down, left, up, rsign, lsign);
        get_corners(lat, n, x1_t_1, x_1_t_1, x_1_t1, x1_t1);
        F[2*n]   -= imag( j1[n]-j1[up]
                            - j2[x_1_t1] + j2[right]
                            - j3[right] + j3[x1_t1]
                            + j4[down] - j4[n]
                            );
        F[2*n+1] -= imag(j1[left] - j1[n]
                            + j2[n] - j2[right]
                            - j3[x1_t1] + j3[down]
                            - j4[down] + j4[x1_t_1]
                            );
	});


}

//Generates new configuration [U,Pi]
void HMC::Leapfrog(const spinor& phi){
    const double StepSize = trajectory_length / (MD_steps * 1.0);
    PConf_copy = PConf;
	GConf_copy = GConf; //Copy of the gauge configuration

    c_double* U = GConf_copy.Conf.val;
    double* P = PConf_copy.val;
    const double* F = Forces.val;

    //Conf_copy = Conf*exp(0.5i * StepSize * PConf_copy)
    for_each_site([=] HOST_DEVICE (int n){
        U[2*n]   = U[2*n]   * exp(0.5 * I_number * StepSize * P[2*n]);
        U[2*n+1] = U[2*n+1] * exp(0.5 * I_number * StepSize * P[2*n+1]);
    });

	Force(GConf_copy,phi);

    for (int step = 1; step < MD_steps; step++) {
        //PConf_copy += StepSize*force
        //Conf_copy *= exp(i * StepSize * PConf_copy)
        for_each_site([=] HOST_DEVICE (int n){
            //mu = 0
            P[2*n] += StepSize *  F[2*n];
            U[2*n] *= exp(I_number * StepSize * P[2*n]);

            //mu = 1
            P[2*n+1] += StepSize *  F[2*n+1];
            U[2*n+1] *= exp(I_number * StepSize * P[2*n+1]);
        });
        Force(GConf_copy,phi);
    }

    //PConf_copy += StepSize*force
    //Conf_copy = Conf*exp(0.5i * StepSize* PConf_copy)
    for_each_site([=] HOST_DEVICE (int n){
        //mu = 0
        P[2*n] += StepSize * F[2*n];
        U[2*n] *= exp(0.5 * I_number * StepSize * P[2*n]);

        //mu = 1
        P[2*n+1] += StepSize * F[2*n+1];
        U[2*n+1] *= exp(0.5 * I_number * StepSize * P[2*n+1]);
    });
}

double HMC::Action(GaugeConf& GConfig, const spinor& phi) {
    double action;
    exchange_halo_gauge(GConfig.Conf.val);
    #ifdef CLOVER
        GConfig.Compute_Q(); //The clover term has to be computed with the links of GConfig, not with the ones of the last force evaluation
    #endif
    GConfig.Compute_Plaquette01();
    //Gauge contribution
    const c_double* plaquette = GConfig.Plaquette01.val;
    const double b = beta;
    action = lattice_sum<double>([=] HOST_DEVICE (int n){
        return b * real(1.0-plaquette[n]);
    });

    //Fermions contribution
    //Phi^dagger (DD^dagger)^-1 Phi = dot(Phi,(DD^dagger)^-1 Phi) (the dot function takes into account the dagger)
    CG_convergence = conjugate_gradient(GConfig, phi,TEMP);
    action += real( dot( TEMP, phi));

    return action;
}

double HMC::Hamiltonian(GaugeConf& GConfig, const re_field& Pi,const spinor& phi) {
    //Momentum contribution
    const double* P = Pi.val;
    double H = lattice_sum<double>([=] HOST_DEVICE (int n){
        return 0.5 * P[2*n] * P[2*n] + 0.5 * P[2*n+1] * P[2*n+1];
    });
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

    exchange_halo_gauge(GConf.Conf.val);
    #ifdef CLOVER
        GConf.Compute_Q(); //Terms needed for clover
    #endif
    D_phi(GConf, chi,phi);
    #ifdef TWISTED_MASS
        //phi = (D + i mu gamma_5) chi, with gamma_5 = diag(-1,1), since (D + i mu gamma_5)(D + i mu gamma_5)^dagger = DD^dagger + mu^2.
        //Then phi is distributed according to exp(-phi^dagger (DD^dagger + mu^2)^-1 phi), as required by the action.
        c_double* phi_v = phi.val;
        const c_double* chi_v = chi.val;
        const double tm = sim_params::tm;
        for_each_site([=] HOST_DEVICE (int n){
            phi_v[2*n]   -= I_number * tm * chi_v[2*n];
            phi_v[2*n+1] += I_number * tm * chi_v[2*n+1];
        });
    #endif
    Leapfrog(phi); //Evolve [Pi] and [U] GConf is copied to GConf_copy inside leapfrog
    deltaH = Hamiltonian(GConf_copy, PConf_copy, phi) - Hamiltonian(GConf, PConf, phi); //deltaH = Hamiltonian[U'][Pi'] - [U][Pi]
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

    const int adapt_begin = static_cast<int>(0.2 * Ntherm);
    const int adapt_end   = static_cast<int>(0.8 * Ntherm);
    //Adapt from the very first trajectory so the chain can actually move, then
    //throw away the hot-start statistics at adapt_begin and average only over
    //equilibrated trajectories. Freeze at adapt_end so the remaining
    //thermalization trajectories run the exact kernel used for measurement.
    hmc::HMCTuner tuner(trajectory_length, MD_steps, p_target, adapt_end);
    tuner.set_verbose(mpi::rank == 0);
    for(int i = 0; i < Ntherm; i++) {
        const int    md_used  = MD_steps;
        const double eps_used = trajectory_length / MD_steps;

        HMC_Update();
        if (tune_MD && i < adapt_end) {
             tuner.record(deltaH);
             MD_steps = tuner.n_steps();
             MPI_Bcast(&MD_steps, 1, MPI_INT, 0, mpi::cart_comm);
         }

        if (tune_MD && i == adapt_begin) tuner.reset(adapt_end - adapt_begin);
        if (tune_MD && i == adapt_end) {
            tuner.freeze();
            MD_steps = tuner.n_steps();
            MPI_Bcast(&MD_steps, 1, MPI_INT, 0, mpi::cart_comm);
            sim_params::MD_steps = MD_steps;
        }

        if ((i+1) % 100 == 0 && mpi::rank == 0) {
            const char* phase = !tune_MD          ? "fixed"
                              : (i <  adapt_begin) ? "burn "
                              : (i <  adapt_end)   ? "adapt"
                                                   : "froze";
            std::cout << "therm " << std::setw(6) << i+1 << "/" << Ntherm
                      << "  [" << phase << "]"
                      << "  MD_steps = " << std::setw(4) << md_used
                      << "  eps = "      << std::setw(9) << eps_used
                      << "  dH = "     << std::setw(9) << deltaH;
            std::cout << std::endl;
        }
     } //Thermalization

    therm = true;
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




