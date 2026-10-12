#ifndef GAUGECONF_H_INCLUDED
#define GAUGECONF_H_INCLUDED
#include "utils.h"
#include <fstream>
#include "mpi.h"
#include "halo_exchange.h"
#include "boundary.h"

class GaugeConf {
public:

	GaugeConf() :
		Conf(mpi::maxSizeH),		//Gauge configuration
		Staples(mpi::maxSizeH),		//Staples
		Plaquette01(mpi::sitesH),
		Q01(mpi::sitesH),			//Clover
		Q10(mpi::sitesH),
		P1(mpi::sitesH),			//U_01
		P2(mpi::sitesH),			//U_{1,-0}
		P3(mpi::sitesH),			//U_{-0,-1}
		P4(mpi::sitesH)				//U_{-10}
	{}

	//Copy constructor, assignment operator and destructor are the default ones: every field is copied
	//(GPU to GPU when compiled with CUDA)

	/*
		Random initialization of the gauge configuration
		It calls RandomU(1) for every site
	*/
	void initialization();

	spinor Conf; //U_0(n) = Conf.val[2n], U_1(n) = Conf.val[2n+1]
	spinor Staples; //Staples
	//The following fields have one complex number per site
	spinor Plaquette01; //Plaquette U_01(x)
	spinor Q01;
	spinor Q10;
	spinor P1;
	spinor P2;
	spinor P3;
	spinor P4;

	/*
		Computes staple
		U_v(x) U_m(x+v) U*_v(x+m) + U*_v(x-v) U_m(x-v) U_v(x+m-v)
    	mu = 0 time direction, mu = 1 space direction
		WARNING: Some references define the staple as the conjugate of the term I just wrote above.
		The halo of Conf (corners included) has to be communicated before with exchange_halo_gauge.
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
		The halo of Conf (corners included) has to be communicated before with exchange_halo_gauge.
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
