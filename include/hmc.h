#ifndef HMC_INCLUDED
#define HMC_INCLUDED

#include "gauge_conf.h"
#include "matrix_operations.h"
#include "conjugate_gradient.h"

class HMC {

public:
	HMC(const int& MD_steps, const double& trajectory_length, const int& Ntherm, const int& Nmeas, 
		const int& Nsteps, const double& beta, const int& Nspace, const int& Ntime, const int& Ntot, const double& m0) : 
		MD_steps(MD_steps), trajectory_length(trajectory_length), Ntherm(Ntherm), Nmeas(Nmeas), Nsteps(Nsteps), 
		beta(beta), Ns(Nspace), Nt(Ntime), Ntot(Ntot), m0(m0) {		
		PConf = std::vector<std::vector<double>>(Ntot, std::vector<double>(2, 0));//Momenta PI
		PConf_copy = std::vector<std::vector<double>>(Ntot, std::vector<double>(2, 0));//Momenta PI copy
		Forces = std::vector<std::vector<double>>(Ntot, std::vector<double>(2, 0)); //Forces
		staples = std::vector<std::vector<std::complex<double>>>(Ntot, std::vector<std::complex<double>>(2, 0)); //staples
		Conf_copy = std::vector<std::vector<std::complex<double>>>(Ntot, std::vector<std::complex<double>>(2, 0)); 
		//std::cout << "HMC constructor" << std::endl;
	}
	~HMC() {
		//std::cout << "HMC destructor" << std::endl;
	}
	//Most of these functions will be private in the future, for now I leave them this way for testing
	void HMC_Update(const int& MD_steps, const double& trajectory_length);
	void Leapfrog(const std::vector<std::vector<std::complex<double>>>& Conf,
	const std::vector<std::vector<std::complex<double>>>& phi );
	void StapleHMC(const std::vector<std::vector<std::complex<double>>>& U);
	void Force_G(const std::vector<std::vector<std::complex<double>>>& U); //force for gauge part
	void Force(const std::vector<std::vector<std::complex<double>>>& U,const std::vector<std::vector<std::complex<double>>>& phi); //force_G + fermions
	//double DeltaH();
	//void HMC_algorithm(const int& MD_steps, const double& trajectory_length, const int& Ntherm, const int& Nmeas, const int& Nsteps);
	std::vector<std::vector<double>> getForce() { return Forces; }
	std::vector<std::vector<double>> getMomentum() { return PConf_copy; }
	std::vector<std::vector<std::complex<double>>> getConfcopy() { return Conf_copy; }
private:
	int Ns, Nt, Ntot;
	int MD_steps, Ntherm, Nmeas, Nsteps;
	double trajectory_length;
	double beta;
	double m0;
	std::vector<std::vector<double>> PConf;
	std::vector<std::vector<double>> PConf_copy;
	std::vector<std::vector<double>> Forces;
	std::vector<std::vector<std::complex<double>>> staples;
	std::vector<std::vector<std::complex<double>>> Conf_copy;//for testing, this should be in gauge_conf.h
};
#endif