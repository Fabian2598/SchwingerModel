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

		PConf = re_field(Ntot, re_vector(2, 0));//Momenta PI
		PConf_copy = re_field(Ntot, re_vector(2, 0));//Momenta PI copy
		Forces = re_field(Ntot, re_vector(2, 0)); //Forces
		Ep = 0; dEp = 0;
		acceptance_rate = 0;

	}
	~HMC() {} 
	
	void HMC_algorithm();
	double getEp() { return Ep; }
	double getdEp() { return dEp; }
	double getgS() { return gS; }
	double getdgS() { return dgS; }
	double getacceptance_rate() { return acceptance_rate/((Nmeas+Nsteps*Nmeas)*1.0); }

private:
	int Nx, Nt, Ntot;
	int MD_steps, Ntherm, Nmeas, Nsteps;
	int saveconf;
	double trajectory_length;
	double beta;
	double m0;
	double Ep, dEp;
	double gS, dgS;
	double acceptance_rate;
	bool therm = false;
	re_field PConf;
	re_field PConf_copy;
	re_field Forces;
	GaugeConf GConf; //Gauge configuration
	GaugeConf GConf_copy; //Copy of the gauge configuration

	double Action(GaugeConf& GConfig, const spinor& phi);
	void Force_G(GaugeConf& GConfig); //force for gauge part
	void Force(GaugeConf& GConfig, const spinor& phi); //force_G + fermions
	void Leapfrog(const spinor& phi );
	double Hamiltonian(GaugeConf& GConfig, const re_field& Pi, const spinor& phi);
	void HMC_Update();
	
};


//Save Momentum configuration (useful for testing)
inline void SavePIConf(const spinor& PConf, char* Name) {

	using namespace LV;
	char NameData[500], Data_str[500];
	sprintf(NameData, Name);
	std::ofstream Datfile;
	Datfile.open(NameData);
	for (int x = 0; x < Nx; x++) {
		for (int t = 0; t < Nt; t++) {
			int i = x * Nx + t;
			for (int mu = 0; mu < 2; mu++) {
				sprintf(Data_str, "%-30d%-30d%-30d%-30.17g\n", x, t, mu, PConf[i][mu]);
				Datfile << Data_str;
			}
		}
	}
	Datfile.close();
}

inline re_field RandomMomentum() {
	using namespace LV;
	//We use a Gaussian distribution to sample the momenta
	std::random_device rd;
	std::default_random_engine generator;
	generator.seed(rd()); //This generator has to be seeded differently. srand does't work here
	std::normal_distribution<double> distribution(0.0, 1.0); //mu, std
	re_field RandPI(Nx * Nt, re_vector(2, 0));
	for (int i = 0; i < Nx * Nt; i++) {
		for (int mu = 0; mu < 2; mu++) {
			RandPI[i][mu] = distribution(generator);
		}
	}
	return RandPI;
}

//Random Chi vector 
inline spinor RandomChi() {
	//We use a Gaussian distribution to sample the momenta
	std::random_device rd;
	std::default_random_engine generator;
	generator.seed(rd()); //This generator has to be seeded differently. srand doesn't work here
	// The result is very much affected by the standard deviation
	std::normal_distribution<double> distribution(0.0, 1/sqrt(2)); //mu, standard deviation
	using namespace LV;
	spinor RandPI(Nx * Nt, c_vector(2, 0));
	for (int i = 0; i < Nx * Nt; i++) {
		for (int mu = 0; mu < 2; mu++) {
			RandPI[i][mu] = 1.0 * distribution(generator) + I_number * distribution(generator);//distribution(generator);

		}
	}
	return RandPI;
}

#endif