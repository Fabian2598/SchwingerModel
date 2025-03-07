#ifndef GAUGECONF_H_INCLUDED
#define GAUGECONF_H_INCLUDED
#include "operator_overloads.h"
#include <complex>  
#include "variables.h"
#include "statistics.h"
#include <fstream>

c_double RandomU1(); //defined on gauge_conf.cpp
void Coordinates(); //Vectorized coordinates. Coords[x][t]. Computed only once
void SaveConf(c_matrix& Conf, char* Name); //Save Gauge configuration

//class GaugeConf;

class GaugeConf {
public:
	GaugeConf() {} //Default constructor
	//Constructor
	GaugeConf(const int& Nspace, const int& Ntime) : Ns(Nspace), Nt(Ntime), Ntot(Nspace* Ntime) {
		Conf = c_matrix(Ntot, c_vector(2, 0)); //Gauge configurationion copy
		Plaquette01 = c_vector(Ntot, 0); //Plaquettes
		Staples = c_matrix(Ntot, c_vector(2, 0)); //Staples
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
	c_matrix Conf; //Conf[Ntot][mu]; //The first entry is the number of lattice points, the second entry is mu
	c_matrix Staples; //Staples
	c_vector Plaquette01; //Plaquette U_01(x)

	void Compute_Staple(); //Computes staples
	void Compute_Plaquette01(); //Computes plaquettes
	double MeasureSp_HMC();
private:
	int Ns, Nt, Ntot;
};




#endif
