#ifndef HMC_INCLUDED
#define HMC_INCLUDED

#include "gauge_conf.h"
#include "matrix_operations.h"
#include "conjugate_gradient.h"

class HMC {

public:
	HMC(GaugeConf& GConf, const int& MD_steps, const double& trajectory_length, const int& Ntherm, const int& Nmeas, 
		const int& Nsteps, const double& beta, const int& Nspace, const int& Ntime, const int& Ntot, const double& m0) : 
		MD_steps(MD_steps), trajectory_length(trajectory_length), Ntherm(Ntherm), Nmeas(Nmeas), Nsteps(Nsteps), 
		beta(beta), Ns(Nspace), Nt(Ntime), Ntot(Ntot), m0(m0), GConf(GConf) {

		PConf = std::vector<std::vector<double>>(Ntot, std::vector<double>(2, 0));//Momenta PI
		PConf_copy = std::vector<std::vector<double>>(Ntot, std::vector<double>(2, 0));//Momenta PI copy
		Forces = std::vector<std::vector<double>>(Ntot, std::vector<double>(2, 0)); //Forces
		//std::cout << "HMC constructor" << std::endl;
	}
	~HMC() {} //std::cout << "HMC destructor" << std::endl;	
	//---All these functions will be private in the future, except for the action---//
	void Force_G(GaugeConf& GConfig); //force for gauge part
	void Force(GaugeConf& GConfig, const std::vector<std::vector<std::complex<double>>>& phi); //force_G + fermions
	void Leapfrog(const std::vector<std::vector<std::complex<double>>>& phi );
	double Hamiltonian(GaugeConf& GConfig, const std::vector<std::vector<double>>& Pi, const std::vector<std::vector<std::complex<double>>>& phi);
	//---------------------------------------------------------------------------------//
	
	double Action(GaugeConf& GConfig, const std::vector<std::vector<std::complex<double>>>& phi);
	void HMC_Update();
	void HMC_algorithm();
	std::vector<std::vector<double>> getForce() { return Forces; }
	double getEp() { return Ep; }
	double getdEp() { return dEp; }

private:
	int Ns, Nt, Ntot;
	int MD_steps, Ntherm, Nmeas, Nsteps;
	double trajectory_length;
	double beta;
	double m0;
	double Ep, dEp;
	std::vector<std::vector<double>> PConf;
	std::vector<std::vector<double>> PConf_copy;
	std::vector<std::vector<double>> Forces;
	GaugeConf GConf; //Gauge configuration
	GaugeConf GConf_copy; //Copy of the gauge configuration
};


//Save Momentum configuration (useful for testing)
inline void SavePIConf(const std::vector<std::vector<double>>& PConf, char* Name) {
	char NameData[500], Data_str[500];
	sprintf(NameData, Name);
	std::ofstream Datfile;
	Datfile.open(NameData);
	for (int x = 0; x < Ns; x++) {
		for (int t = 0; t < Nt; t++) {
			int i = x * Ns + t;
			for (int mu = 0; mu < 2; mu++) {
				sprintf(Data_str, "%-30d%-30d%-30d%-30.17g\n", x, t, mu, PConf[i][mu]);
				Datfile << Data_str;
			}
		}
	}
	Datfile.close();
}

inline std::vector<std::vector<double>> RandomMomentum() {
	//We use a Gaussian distribution to sample the momenta
	std::random_device rd;
	std::default_random_engine generator;
	generator.seed(rd()); //This generator has to be seeded differently. srand does't work here
	std::normal_distribution<double> distribution(0.0, 1.0); //mu, std
	std::vector<std::vector<double>> RandPI(Ns * Nt, std::vector<double>(2, 0));
	for (int i = 0; i < Ns * Nt; i++) {
		for (int mu = 0; mu < 2; mu++) {
			RandPI[i][mu] = distribution(generator);
		}
	}
	return RandPI;
}

//Random Chi vector 
inline std::vector<std::vector<std::complex<double>>> RandomChi() {
	//We use a Gaussian distribution to sample the momenta
	std::random_device rd;
	std::default_random_engine generator;
	generator.seed(rd()); //This generator has to be seeded differently. srand doesn't work here
	std::normal_distribution<double> distribution(0.0, 1.0); //mu, standard deviation
	std::vector<std::vector<std::complex<double>>> RandPI(Ns * Nt, std::vector<std::complex<double>>(2, 0));
	for (int i = 0; i < Ns * Nt; i++) {
		for (int mu = 0; mu < 2; mu++) {
			RandPI[i][mu] = 1.0 * distribution(generator) + 0 * (0, 1);

		}
	}

	return RandPI;
}

#endif