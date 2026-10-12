#ifndef HMC_INCLUDED
#define HMC_INCLUDED

#include "gauge_conf.h"
#include "conjugate_gradient.h"
#include "hmc_tuner.h"


class HMC {

public:
	HMC(GaugeConf& GConf, const int& MD_steps, const double& trajectory_length, const int& Ntherm, const int& Nmeas,
		const int& Nsteps, const double& beta, const int& Nspace, const int& Ntime, const double& m0, const int& saveconf) :
		Nx(Nspace), Nt(Ntime), Ntot(Nspace*Ntime), MD_steps(MD_steps), Ntherm(Ntherm), Nmeas(Nmeas), Nsteps(Nsteps),
		saveconf(saveconf), trajectory_length(trajectory_length), beta(beta), m0(m0),
		PConf(mpi::maxSizeH), 		//Momenta PI
		PConf_copy(mpi::maxSizeH), 	//Momenta PI copy
		Forces(mpi::maxSizeH), 		//Forces
		GConf(GConf),
		chi(mpi::maxSizeH),
		TEMP(mpi::maxSizeH), 		//buffer
		//Needed to compute the fermion force of the clover term
		J1(mpi::sitesH), J2(mpi::sitesH), J3(mpi::sitesH), J4(mpi::sitesH) {

		Ep = 0; dEp = 0;
		acceptance_rate = 0;
		illConfId = 0;

	}

	void HMC_algorithm();
	//Enable automatic tuning of MD_steps during thermalization.
	//p_acc_target: 0.65 is cost-optimal in theory, 0.75-0.80 is safer with clover.
	void enableTuning(double p_acc_target = 0.78) { tune_MD = true; p_target = p_acc_target; }
	double getDeltaH() const { return deltaH; }
	int getMDsteps() const { return MD_steps; }

	double getEp() { return Ep; }
	double getdEp() { return dEp; }
	double getgS() { return gS; }
	double getdgS() { return dgS; }
	double getacceptance_rate(int conf_number) { return acceptance_rate/((conf_number)*1.0); }

	//The following functions are only used internally. They are public because CUDA does not allow
	//GPU kernels (device lambdas) inside private or protected member functions.
	double Action(GaugeConf& GConfig, const spinor& phi);
	void Force_G(GaugeConf& GConfig); //force for gauge part
	void Force_Clover(const GaugeConf& GConfig,const spinor& left_term, const spinor& right_term); //Clover term contribution
	void Force(GaugeConf& GConfig, const spinor& phi); //force_G + fermions + Clover
	void Leapfrog(const spinor& phi );
	double Hamiltonian(GaugeConf& GConfig, const re_field& Pi, const spinor& phi);
	void HMC_Update();

	void RandomPI();
	void RandomCHI();

private:
	int Nx, Nt, Ntot;
	int MD_steps, Ntherm, Nmeas, Nsteps;
	int saveconf;
	int conf_i;
	double trajectory_length;
	double beta;
	double m0;
	double Ep, dEp;
	double gS, dgS;
	double acceptance_rate;
	int CG_convergence; //1->converge, 0->not converged
	int illConfId;
	bool therm = false;
	double deltaH = 0.0;     //dH of the most recent trajectory
	bool   tune_MD = false;  //adapt MD_steps during thermalization
	double p_target = 0.78;  //target acceptance rate
	re_field PConf; //Momenta PI
	re_field PConf_copy; //Momenta PI copy
	re_field Forces; //Forces
	GaugeConf GConf; //Gauge configuration
	GaugeConf GConf_copy; //Copy of the gauge configuration
	spinor chi;
	spinor TEMP;  //buffer

	//For the force of the clover term (one complex number per site)
	spinor J1;
	spinor J2;
	spinor J3;
	spinor J4;

};


#endif
