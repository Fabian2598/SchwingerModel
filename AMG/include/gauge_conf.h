#ifndef GAUGECONF_H_INCLUDED
#define GAUGECONF_H_INCLUDED
#include "variables.h"
#include <complex>  
#include <iostream>
#include <fstream>

std::complex<double> RandomU1(); //defined on gauge_conf.cpp
void Coordinates(); //Vectorized coordinates. Coords[x][t]. Computed only once
void Aggregates(); //Compute aggregates for a lattice blocking. 
void SaveConf(std::vector<std::vector<std::complex<double>>>& Conf, char* Name); //Save Gauge configuration
void PrintAggregates();

//class GaugeConf;

class GaugeConf {
public:
	GaugeConf() {} //Default constructor
	//Constructor
	GaugeConf(const int& Nspace, const int& Ntime) : Ns(Nspace), Nt(Ntime), Ntot(Nspace* Ntime) {
		Conf = std::vector<std::vector<std::complex<double>>>(Ntot, std::vector<std::complex<double>>(2, 0)); //Gauge configurationion copy
		Plaquette01 = std::vector<std::complex<double>>(Ntot, 0); //Plaquettes
		Staples = std::vector<std::vector<std::complex<double>>>(Ntot, std::vector<std::complex<double>>(2, 0)); //Staples
	}
	//Copy constructor
	GaugeConf(const GaugeConf& GConfig) : Ns(GConfig.getNs()), Nt(GConfig.getNt()), Ntot(Ns*Nt) {
		Conf = GConfig.Conf; 
		Plaquette01 = GConfig.Plaquette01; 
		Staples = GConfig.Staples; 
	}
	~GaugeConf() {}; //Destructor

	void initialization(); //Random initialization of the gauge configuration
	int getNs() const { return Ns; }
	int getNt() const { return Nt; }

	//HMC needs access to these variables
	std::vector<std::vector<std::complex<double>>> Conf; //Conf[Ntot][mu]; //The first entry is the number of lattice points, the second entry is mu
	std::vector<std::vector<std::complex<double>>> Staples; //Staples
	std::vector<std::complex<double>> Plaquette01; //Plaquette U_01(x)


	
	void Compute_Staple(); //Computes staples
	void Compute_Plaquette01(); //Computes plaquettes
	double MeasureSp_HMC();
private:
	int Ns, Nt, Ntot;
};




#endif
