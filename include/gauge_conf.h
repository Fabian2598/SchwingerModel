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
		Q01 = new c_double[mpi::sitesH];	//Clover
		Q10 = new c_double[mpi::sitesH];	
		P1 = new c_double[mpi::sitesH];	//U_01 
		P2 = new c_double[mpi::sitesH];	//U_{1,-0} 
		P3 = new c_double[mpi::sitesH];	//U_{-0,-1} 
		P4 = new c_double[mpi::sitesH];	//U_{-10}
	}

	/*
	Copy constructor
	*/
	GaugeConf(const GaugeConf& GConfig) {
		Conf = GConfig.Conf; 
		Staples = GConfig.Staples; 
		Plaquette01 = new c_double[mpi::sitesH];
		Q01 = new c_double[mpi::sitesH];	//Clover
		Q10 = new c_double[mpi::sitesH];	
		P1 = new c_double[mpi::sitesH];	
		P2 = new c_double[mpi::sitesH];	
		P3 = new c_double[mpi::sitesH];	
		P4 = new c_double[mpi::sitesH];	
        std::copy(GConfig.Plaquette01, GConfig.Plaquette01 + mpi::sitesH, Plaquette01);
		std::copy(GConfig.Q01, GConfig.Q01 + mpi::sitesH, Q01);
		std::copy(GConfig.Q10, GConfig.Q10 + mpi::sitesH, Q10);
		std::copy(GConfig.P1, GConfig.P1 + mpi::sitesH, P1);
		std::copy(GConfig.P2, GConfig.P2 + mpi::sitesH, P2);
		std::copy(GConfig.P3, GConfig.P3 + mpi::sitesH, P3);
		std::copy(GConfig.P4, GConfig.P4 + mpi::sitesH, P4);
	}

	/*
	Assignment operator
	*/
	GaugeConf& operator=(const GaugeConf& GConfig) {
		if (this != &GConfig) {
			Conf = GConfig.Conf;
			Staples = GConfig.Staples;
			delete[] Plaquette01;
			delete[] Q01;
			delete[] Q10;
			delete[] P1;
			delete[] P2;
			delete[] P3;
			delete[] P4;
			Plaquette01 = new c_double[mpi::sitesH];
			Q01 = new c_double[mpi::sitesH];	
			Q10 = new c_double[mpi::sitesH];	
			P1 = new c_double[mpi::sitesH];	 
			P2 = new c_double[mpi::sitesH];	
			P3 = new c_double[mpi::sitesH];	
			P4 = new c_double[mpi::sitesH];
			std::copy(GConfig.Plaquette01, GConfig.Plaquette01 + mpi::sitesH, Plaquette01);
			std::copy(GConfig.Q01, GConfig.Q01 + mpi::sitesH, Q01);
			std::copy(GConfig.Q10, GConfig.Q10 + mpi::sitesH, Q10);
			std::copy(GConfig.P1, GConfig.P1 + mpi::sitesH, P1);
			std::copy(GConfig.P2, GConfig.P2 + mpi::sitesH, P2);
			std::copy(GConfig.P3, GConfig.P3 + mpi::sitesH, P3);
			std::copy(GConfig.P4, GConfig.P4 + mpi::sitesH, P4);
		}
		return *this;
	}

	/*
	Destructor
	*/
	~GaugeConf() {
		delete[] Plaquette01;
		delete[] Q01;
		delete[] Q10;
		delete[] P1;
		delete[] P2;
		delete[] P3;
		delete[] P4;
	}; 

	/*
		Random initialization of the gauge configuration
		It calls RandomU(1) for every site
	*/
	void initialization(); 

	spinor Conf; 
	spinor Staples; //Staples
	c_double* Plaquette01; //Plaquette U_01(x)
	c_double* Q01; 
	c_double* Q10; 
	c_double* P1;
	c_double* P2;
	c_double* P3;
	c_double* P4;

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
		Compute Q_01(x) and Q_10(x) for the clover term and 4-plaquettes on the mv plane at each 
		Q_mv(x) = U_{m,v}(x) + U_{v,-m}(x) + U_{-m,-v}(x) + U_{-v,m}(x)
	*/
	void Compute_Q();

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
