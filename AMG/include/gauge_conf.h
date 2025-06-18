#ifndef GAUGECONF_H_INCLUDED
#define GAUGECONF_H_INCLUDED
#include "variables.h"
#include "operator_overloads.h"
#include <complex>  
#include <iostream>
#include <fstream>

std::complex<double> RandomU1(); //defined on gauge_conf.cpp



//class GaugeConf;
class GaugeConf {
public:
	GaugeConf() {} //Default constructor
	//Constructor
	GaugeConf(const int& Nspace, const int& Ntime) : Nx(Nspace), Nt(Ntime), Ntot(Nspace* Ntime) {
		Conf = std::vector<std::vector<std::complex<double>>>(Ntot, std::vector<std::complex<double>>(2, 0)); //Gauge configurationion copy
		
	}
	//Copy constructor
	GaugeConf(const GaugeConf& GConfig) : Nx(GConfig.getNx()), Nt(GConfig.getNt()), Ntot(Nx*Nt) {
		Conf = GConfig.Conf; 
	}
	~GaugeConf() {}; //Destructor

	void initialization(); //Random initialization of the gauge configuration
	void set_gconf(const std::vector<std::vector<std::complex<double>>>& CONF) {Conf = CONF;}
	void SaveConf(char* Name); //Save Gauge configuration
	int getNx() const { return Nx; }
	int getNt() const { return Nt; }

	std::vector<std::vector<std::complex<double>>> Conf; //Conf[Ntot][mu]; 
	//The first entry is the number of lattice points, the second entry is mu
	void PrintConf();
private:
	int Nx, Nt, Ntot;
};




#endif
