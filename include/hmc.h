#ifndef HMC_INCLUDED
#define HMC_INCLUDED

#include "gauge_conf.h"
#include "conjugate_gradient.h"


class HMC {

public:
	HMC(GaugeConf& GConf, const int& MD_steps, const double& trajectory_length, const int& Ntherm, const int& Nmeas, 
		const int& Nsteps, const double& beta, const int& Nspace, const int& Ntime, const double& m0, const int& saveconf) : 
		MD_steps(MD_steps), trajectory_length(trajectory_length), Ntherm(Ntherm), Nmeas(Nmeas), Nsteps(Nsteps), 
		beta(beta), Nx(Nspace), Nt(Ntime), Ntot(Nspace*Ntime), m0(m0), saveconf(saveconf), GConf(GConf) {

		Ep = 0; dEp = 0;
		acceptance_rate = 0;
		illConfId = 0;

		PConf = re_field(mpi::maxSizeH); //Momenta PI
		PConf_copy = re_field(mpi::maxSizeH); //Momenta PI copy
		Forces = re_field(mpi::maxSizeH); //Forces
		chi = spinor(mpi::maxSizeH);
		TEMP = spinor(mpi::maxSizeH); //buffer

		//Needed to compute the fermion force of the clover term
		J1 = new c_double[mpi::sitesH];	 
		J2 = new c_double[mpi::sitesH];	
		J3 = new c_double[mpi::sitesH];	
		J4 = new c_double[mpi::sitesH];

	}
	~HMC() {
		delete[] J1;
		delete[] J2;
		delete[] J3;
		delete[] J4;
	} 
	
	void HMC_algorithm();
	double getEp() { return Ep; }
	double getdEp() { return dEp; }
	double getgS() { return gS; }
	double getdgS() { return dgS; }
	double getacceptance_rate(int conf_number) { return acceptance_rate/((conf_number)*1.0); }

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
	re_field PConf; //Momenta PI
	re_field PConf_copy; //Momenta PI copy
	re_field Forces; //Forces
	GaugeConf GConf; //Gauge configuration
	GaugeConf GConf_copy; //Copy of the gauge configuration
	spinor chi;
	spinor TEMP;  //buffer

	//For the force of the clover term
	c_double* J1;
	c_double* J2;
	c_double* J3;
	c_double* J4;

	double Action(GaugeConf& GConfig, const spinor& phi);
	void Force_G(GaugeConf& GConfig); //force for gauge part
	void Force_Clover(GaugeConf& GConfig,const spinor& left_term, const spinor& right_term); //Clover term contribution
	void Force(GaugeConf& GConfig, const spinor& phi); //force_G + fermions + Clover
	void Leapfrog(const spinor& phi );
	double Hamiltonian(GaugeConf& GConfig, const re_field& Pi, const spinor& phi);
	void HMC_Update();

	void RandomPI();
	void RandomCHI();
	
};


#endif