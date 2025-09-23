#ifndef GAUGECONF_H_INCLUDED
#define GAUGECONF_H_INCLUDED
#include "operator_overloads.h"
#include <complex>  
#include "variables.h"
#include "statistics.h"
#include <iomanip>
#include <fstream>

/*
Generate a random U(1) variable
*/
c_double RandomU1(); 

/*
	Save Gauge configuration
*/
void SaveConf(const c_double (&Conf)[2*LV::Ntot], const std::string& Name); 


class GaugeConf {
public:
	/*
	Nspace: number of lattice points in the space direction
	Ntime: number of lattice points in the time direction
	*/
	GaugeConf() {} 

	GaugeConf(const int& Nspace, const int& Ntime) : Nx(Nspace), Nt(Ntime), Ntot(Nspace* Ntime) {
	}

	/*
	Copy constructor
	*/
	GaugeConf(const GaugeConf& GConfig) : Nx(GConfig.getNx()), Nt(GConfig.getNt()), Ntot(Nx*Nt) {
		Conf = GConfig.Conf; 
		Plaquette01 = GConfig.Plaquette01; 
		Staples = GConfig.Staples; 
	}
	/*
	Destructor
	*/
	~GaugeConf() {
	}; 

	/*
		Random initialization of the gauge configuration
		It calls RandomU(1) for every site
	*/
	void initialization(); 
	int getNx() const { return Nx; }
	int getNt() const { return Nt; }

	/*
		std::vector with the gauge configuration. Conf[Nx*Nt][2]
		The first entry is the number of lattice points, the second entry is mu
	*/
	spinor Conf; 
	spinor Staples; //Staples
	c_double Plaquette01[2*LV::Ntot]; //Plaquette U_01(x)

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
private:
	int Nx, Nt, Ntot;
};




#endif
