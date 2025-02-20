#ifndef AMG_H_INCLUDED
#define AMG_H_INCLUDED
#include "gauge_conf.h"
#include <algorithm>


class AMG() {
public:
	AMG(const GaugeConf & Gconf,const int& Ns, const int& Nt, const int& Ntest, const double& m0) : GConf(GConf), Ns(Ns), Nt(Nt), Ntot(Ntot), Ntest(Ntest), m0(m0) {
		test_vectors_init = std::vector<std::vector<std::vector<std::complex<double>>>>(Ntest,
		std::vector<std::vector<std::complex<double>>>( Ns * Nt, std::vector<std::complex<double>> (2,0))); //test_vectors[Number of test vectors][Ns Nt][2], two spin directions, no color
	
	};
	~AMG();
	//--CHECK LATER WHICH MEMBER FUNCTIONS WILL BE PRIVATE--//

	void tv_init(const double& eps);
	void tv_update(); //update test vectors
	std::vector<std::vector<std::complex<double>>> P_v(); //Prolongation operator times vector


	std::vector<std::vector<std::vector<std::complex<double>>>> test_vectors; //test vectors[Ntest][Ntot][spin components]


private:
	GaugeConf GConf;
	double m0; 
	int Ns, Nt, Ntest;
	int Ntot;
};

#endif