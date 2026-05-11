#ifndef GAUGECONF_H_INCLUDED
#define GAUGECONF_H_INCLUDED
#include "variables.h"
#include <iomanip>
#include <fstream>
#include <iostream>
#include <string>


/*
Generate a random U(1) variable
*/
c_double RandomU1(); 


class GaugeConf {
public:

	GaugeConf() {
		Conf = spinor(mpi::maxSize); //Gauge configuration
	}

	/*
	Copy constructor
	*/
	GaugeConf(const GaugeConf& GConfig) {
		Conf = GConfig.Conf; 
	}

	/*
	Assignment operator
	*/
	GaugeConf& operator=(const GaugeConf& GConfig) {
		if (this != &GConfig) {
			Conf = GConfig.Conf;
		}
		return *this;
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

	/*
		std::vector with the gauge configuration. Conf[Nx*Nt][2]
		The first entry is the number of lattice points, the second entry is mu
	*/
	spinor Conf; 

	//Reads a gauge configuration from a file
	void read_conf(const std::string& name);

	void readBinary(const std::string& name);
	
};




#endif
