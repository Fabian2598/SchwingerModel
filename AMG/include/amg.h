#ifndef AMG_H_INCLUDED
#define AMG_H_INCLUDED
#include "gauge_conf.h"
#include "conjugate_gradient.h"
#include <algorithm>
#include <random>

void PrintVector(const std::vector<std::vector<std::complex<double>>>& v );
std::vector<std::vector<std::complex<double>>> canonical_vector(const int& i, const int& N1, const int& N2);

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
	std::vector<std::vector<std::complex<double>>> TwoGrid(const int& nu1, const int& nu2, const std::vector<std::vector<std::complex<double>>>& x0, const std::vector<std::vector<std::complex<double>>>& phi);
	//void tv_update(); //update test vectors
	std::vector<std::vector<std::complex<double>>> P_v(const std::vector<std::vector<std::complex<double>>>& v); //P v
	std::vector<std::vector<std::complex<double>>> Pt_v(const std::vector<std::vector<std::complex<double>>>& v); // Pt v
	std::vector<std::vector<std::complex<double>>> Pt_D_P(const std::vector<std::vector<std::complex<double>>>& v); //Dc v = P^T D P v 

	std::vector<std::vector<std::complex<double>>> bi_cgstab_Dc(const c_matrix& U, const c_matrix& phi, const c_matrix& x0,
		 const double& m0, const int& max_iter, const double& tol, const bool& print_message); //Dc^-1 phi

	std::vector<std::vector<std::vector<std::complex<double>>>> test_vectors; //test vectors[Ntest][Ntot][spin components]
private:
	GaugeConf GConf;
	double m0; 
	int Ns, Nt, Ntest;
	int Ntot;
};


#endif