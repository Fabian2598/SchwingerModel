#ifndef GAUGECONF_H_INCLUDED
#define GAUGECONF_H_INCLUDED
#include "variables.h"
#include "statistics.h"
#include <iomanip>
#include <fstream>
#include <iostream>
#include <string>
#include "mpi.h"


/*
Generate a random U(1) variable
*/
c_double RandomU1(); 


class GaugeConf {
public:

	GaugeConf() {
		Plaquette01 = new c_double[mpi::maxSize];
		Conf = spinor(mpi::maxSize); //Gauge configuration
		Staples = spinor(mpi::maxSize); //Staples
	}

	/*
	Copy constructor
	*/
	GaugeConf(const GaugeConf& GConfig) {
		Conf = GConfig.Conf; 
		Staples = GConfig.Staples; 
		Plaquette01 = new c_double[mpi::maxSize];
        std::copy(GConfig.Plaquette01, GConfig.Plaquette01 + mpi::maxSize, Plaquette01);
	}

	/*
	Assignment operator
	*/
	GaugeConf& operator=(const GaugeConf& GConfig) {
		if (this != &GConfig) {
			Conf = GConfig.Conf;
			Staples = GConfig.Staples;
			delete[] Plaquette01;
			Plaquette01 = new c_double[mpi::maxSize];
			std::copy(GConfig.Plaquette01, GConfig.Plaquette01 + mpi::maxSize, Plaquette01);
		}
		return *this;
	}

	/*
	Destructor
	*/
	~GaugeConf() {
		delete[] Plaquette01;
	}; 

	/*
		Random initialization of the gauge configuration
		It calls RandomU(1) for every site
	*/
	void initialization(); 

	/*
		std::vector with the gauge configuration. Conf[Nx*Nt][2]
		The first entry is the number of lattice points, the second entry is mu
	*/
	spinor Conf; 
	spinor Staples; //Staples
	c_double* Plaquette01; //Plaquette U_01(x)

	/*
		Computes staple
		U_v(x) U_m(x+v) U*_v(x+m) + U*_v(x-v) U_m(x-v) U_v(x+m-v)
    	mu = 0 time direction, mu = 1 space direction
		WARNING: Some references define the staple as the conjugate of the term I just wrote above. 
		PERIODIC BOUNDARIES MUST BE PRECOMPUTED FIRST
	*/
	void Compute_Staple();

	/*
		Compute the plaquette
		U_mv(x) = U_m(x) U_v(x+m) U*_m(x+v) U*_v(x)
		with m = 0, nu = 1
		PERIODIC BOUNDARIES MUST BE PRECOMPUTED FIRST
	*/
	void Compute_Plaquette01(); 

	/*
		Measures average plaquette real value
		Sp = < U_01(x) >
		Plaquettes have to be measured before
	*/
	double MeasureSp_HMC();

	/*
		Gauge action
		S_G = beta * sum_x (1 - U_01(x))
		Plaquettes have to be measured before
	*/
	double Compute_gaugeAction(const double& beta); //Computes the gauge action

	//Reads a gauge configuration from a file
	void read_conf(const std::string& name);

	void readBinary(const std::string& name);
	
};


/*
	Save Gauge configuration
*/
void SaveConf(const GaugeConf& GConf, const std::string& Name); 



#endif
