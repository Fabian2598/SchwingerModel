#ifndef GAUGECONF_H_INCLUDED
#define GAUGECONF_H_INCLUDED
#include <iostream>
#include <cmath> 
#include <complex>  
#include "variables.h"
#include "statistics.h"
#include <fstream>

std::complex<double> RandomU1(); //defined on gauge_conf.cpp
void Coordinates();

inline std::vector<std::vector<std::complex<double>>> RandomChi() {
	//We use a Gaussian distribution to sample the momenta
	std::random_device rd;
	std::default_random_engine generator;
	generator.seed(rd()); //This generator has to be seeded differently. srand does't work here
	std::normal_distribution<double> distribution(0.0, 1.0); //mu, std
	std::vector<std::vector<std::complex<double>>> RandPI(Ns * Nt, std::vector<std::complex<double>>(2, 0));
	for (int i = 0; i < Ns * Nt; i++) {
		for (int mu = 0; mu < 2; mu++) {
			RandPI[i][mu] = 1.0 * distribution(generator) + 0 * (0, 1);

		}
	}

	return RandPI;
}

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


//Save Gauge configuration
inline void SaveConf(const std::vector<std::vector<std::complex<double>>>& Conf, char* Name) {
	char NameData[500], Data_str[500];
	sprintf(NameData, Name);
	std::ofstream Datfile;
	Datfile.open(NameData);
	for (int x = 0; x < Ns; x++) {
		for (int t = 0; t < Nt; t++) {
			int i = x * Ns + t;
			for (int mu = 0; mu < 2; mu++) {
				sprintf(Data_str, "%-30d%-30d%-30d%-30.17g%-30.17g\n", x, t, mu, std::real(Conf[i][mu]), std::imag(Conf[i][mu]));
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


class GaugeConf {
public:
	//GaugeConf[Ntot][mu]; //The first entry is the number of lattice points, the second entry 
	//indicates the link variable in x-direction or y-direction.
	GaugeConf(const int& Nspace, const int& Ntime) : Ns(Nspace), Nt(Ntime), Ntot(Nspace*Ntime) {
		Conf = std::vector<std::vector<std::complex<double>>>(Ntot, std::vector<std::complex<double>>(2, 0)); //Gauge configurationion copy
		Plaquette01 = std::vector<std::complex<double>> (Ntot,0);
		Plaquette01_prime = std::vector<std::complex<double>>(Ntot, 0);
	}
	~GaugeConf() {};
	void initialization(); 
	int getNs() const { return Ns; }
	int getNt() const { return Nt; }
	std::vector<std::vector<std::complex<double>>> getConf() const { return Conf; }
	double getEp() const { return Ep; }
	double getdEp() const { return dEp; }

	void printPlaquette(); //Prints plaquettes
		
private:
	int Ns, Nt, Ntot;
	std::vector<std::vector<std::complex<double>>> Conf; //Gauge Conf
	//std::vector<std::vector<int>> Coords;
	//double beta;
	double Ep, dEp;
	std::vector<std::vector<std::complex<double>>> Conf_copy;//Gauge Conf copy
	std::vector<std::complex<double>> Plaquette01;
	std::vector<std::complex<double>> Plaquette01_prime;	
};

std::ostream& operator<<(std::ostream& out, const GaugeConf& GConf);

template <typename T>
void printConf(std::vector<std::vector<T>> V) {
	for (int i = 0; i < V.size(); i++) {
		for (int j = 0; j < V[i].size(); j++) {
			std::cout << V[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

#endif
