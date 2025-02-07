#ifndef HMC_INCLUDED
#define HMC_INCLUDED

#include "gauge_conf.h"

class HMC {

public:
	HMC(const int& MD_steps, const double& trajectory_length, const int& Ntherm, const int& Nmeas, 
		const int& Nsteps, const double& beta, const int& Nspace, const int& Ntime) : 
		MD_steps(MD_steps), trajectory_length(trajectory_length), Ntherm(Ntherm), Nmeas(Nmeas), Nsteps(Nsteps), beta(beta), Ns(Nspace), Nt(Ntime) {		
		PConf = std::vector<std::vector<double>>(Ntot, std::vector<double>(2, 0));//Momenta PI
		PConf_copy = std::vector<std::vector<double>>(Ntot, std::vector<double>(2, 0));//Momenta PI copy
		Forces = std::vector<std::vector<double>>(Ntot, std::vector<double>(2, 0)); //Forces
		//std::cout << "HMC constructor" << std::endl;
	}
	~HMC() {
		//std::cout << "HMC destructor" << std::endl;
	}
	void HMC_Update(const int& MD_steps, const double& trajectory_length);
	void Leapfrog(const int& MD_steps, const double& trajectory_length);
	void Force_G(const std::vector<std::vector<std::complex<double>>>& U); //force for gauge part
	double DeltaH();;
	void HMC_algorithm(const int& MD_steps, const double& trajectory_length, const int& Ntherm, const int& Nmeas, const int& Nsteps);

private:
	int Ns, Nt;
	int MD_steps, Ntherm, Nmeas, Nsteps;
	double trajectory_length;
	double beta;
	std::vector<std::vector<std::complex<double>>> PConf;
	std::vector<std::vector<std::complex<double>>> PConf_copy;
	std::vector<std::vector<std::complex<double>>> Forces;
};
#endif