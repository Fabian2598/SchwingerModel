#ifndef AMG_H_INCLUDED
#define AMG_H_INCLUDED
#include "gauge_conf.h"
#include "conjugate_gradient.h"
#include <algorithm>
#include <random>

void PrintVector(const std::vector<std::vector<std::complex<double>>>& v );
std::vector<std::vector<std::complex<double>>> canonical_vector(const int& i, const int& Ntot);

class AMG {
public:
	AMG(const GaugeConf & GConf,const int& Ns, const int& Nt, const int& Ntest, const double& m0) : GConf(GConf), 
	Ns(Ns), Nt(Nt), Ntot(Ns*Nt), Ntest(Ntest), m0(m0) {	
		test_vectors = std::vector<std::vector<std::vector<std::complex<double>>>>(Ntest,
		std::vector<std::vector<std::complex<double>>>( Ntot, std::vector<std::complex<double>> (2,0))); //test_vectors[Number of test vectors][Ns Nt][2], two spin directions, no color
	}
	~AMG() { };
	//--CHECK LATER WHICH MEMBER FUNCTIONS WILL BE PRIVATE--//

	void tv_init(const double& eps);
	//void tv_update(); //update test vectors
	std::vector<std::vector<std::complex<double>>> P_v(const std::vector<std::vector<std::complex<double>>>& v); //P v
	std::vector<std::vector<std::complex<double>>> Pt_v(const std::vector<std::vector<std::complex<double>>>& v); // Pt v
	std::vector<std::vector<std::vector<std::complex<double>>>> test_vectors; //test vectors[Ntest][Ntot][spin components]
private:
	GaugeConf GConf;
	double m0; 
	int Ns, Nt, Ntest;
	int Ntot;
};


#endif