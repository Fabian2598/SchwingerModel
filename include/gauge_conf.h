#ifndef GAUGECONF_H_INCLUDED
#define GAUGECONF_H_INCLUDED
#include "utils.h"
#include <fstream>
#include "mpi.h"
#include "halo_exchange.h"
#include "boundary.h"

class GaugeConf {
public:

	GaugeConf() {
		Plaquette01 = new c_double[mpi::sitesH];
		Conf = spinor(mpi::maxSizeH); //Gauge configuration
		Staples = spinor(mpi::maxSizeH); //Staples
	}

	/*
	Copy constructor
	*/
	GaugeConf(const GaugeConf& GConfig) {
		Conf = GConfig.Conf; 
		Staples = GConfig.Staples; 
		Plaquette01 = new c_double[mpi::sitesH];
        std::copy(GConfig.Plaquette01, GConfig.Plaquette01 + mpi::sitesH, Plaquette01);
	}

	/*
	Assignment operator
	*/
	GaugeConf& operator=(const GaugeConf& GConfig) {
		if (this != &GConfig) {
			Conf = GConfig.Conf;
			Staples = GConfig.Staples;
			delete[] Plaquette01;
			Plaquette01 = new c_double[mpi::sitesH];
			std::copy(GConfig.Plaquette01, GConfig.Plaquette01 + mpi::sitesH, Plaquette01);
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

	spinor Conf; 
	spinor Staples; //Staples
	c_double* Plaquette01; //Plaquette U_01(x)

	/*
		Computes staple
		U_v(x) U_m(x+v) U*_v(x+m) + U*_v(x-v) U_m(x-v) U_v(x+m-v)
    	mu = 0 time direction, mu = 1 space direction
		WARNING: Some references define the staple as the conjugate of the term I just wrote above. 
	*/
	void Compute_Staple();

	/*
		Compute the plaquette
		U_mv(x) = U_m(x) U_v(x+m) U*_m(x+v) U*_v(x)
		with m = 0, nu = 1
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

	//

	/*
		Read/Write Gauge configuration
	*/
	void ReadConf(const std::string& name);
	void SaveConf(const std::string& Name); 
	
	
};






#endif
