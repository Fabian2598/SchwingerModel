#ifndef AMG_H_INCLUDED
#define AMG_H_INCLUDED
#include "gauge_conf.h"
#include "conjugate_gradient.h"
#include <algorithm>
#include <random>

void PrintVector(const std::vector<std::vector<std::complex<double>>>& v );
void PPrintVector(const std::vector<std::complex<double>>& v );
std::vector<std::vector<std::complex<double>>> canonical_vector(const int& i, const int& N1, const int& N2);
void normalize(std::vector<std::vector<std::complex<double>>>& v);

template <typename T>
c_vector operator*(const T& lambda, const c_vector& A);
c_vector operator+(const c_vector& A, const c_vector& B);
c_vector operator-(const c_vector& A, const c_vector& B);
c_vector operator*(const c_matrix& A, const c_vector& v); //Av matrix times a vector
c_double dot(const c_vector& x, const c_vector& y);
c_double dot_v2(const c_matrix& x, const c_matrix& y);

class AMG {
public:
	AMG(const GaugeConf & GConf,const int& Ns, const int& Nt, const int& Ntest, const double& m0) : GConf(GConf), 
	Ns(Ns), Nt(Nt), Ntot(Ns*Nt), Ntest(Ntest), m0(m0) {	
		test_vectors = std::vector<std::vector<std::vector<std::complex<double>>>>(Ntest,
		std::vector<std::vector<std::complex<double>>>( Ntot, std::vector<std::complex<double>> (2,0))); //test_vectors[Number of test vectors][Ns Nt][2], two spin directions, no color
	}
	~AMG() { };
	//--CHECK LATER WHICH MEMBER FUNCTIONS WILL BE PRIVATE--//

	void tv_init(const double& eps, const int& Nit);
	std::vector<std::vector<std::complex<double>>> TwoGrid(const int& nu1, const int& nu2, 
		const std::vector<std::vector<std::complex<double>>>& x0, const std::vector<std::vector<std::complex<double>>>& phi,const bool& print_message);
	//void tv_update(); //update test vectors
	std::vector<std::vector<std::complex<double>>> P_v(const std::vector<std::vector<std::complex<double>>>& v); //P v
	std::vector<std::vector<std::complex<double>>> Pt_v(const std::vector<std::vector<std::complex<double>>>& v); // Pt v
	std::vector<std::vector<std::complex<double>>> Pt_D_P(const std::vector<std::vector<std::complex<double>>>& v); //Dc v = P^T D P v 

	std::vector<std::vector<std::complex<double>>> bi_cgstab_Dc(const c_matrix& U, const c_matrix& phi, const c_matrix& x0,
		 const double& m0, const int& max_iter, const double& tol, const bool& print_message); //Dc^-1 phi
	
	
	c_vector bi_cg_test(const c_matrix& U, const c_matrix& phi, const c_matrix& x0, const double& m0, const int& max_iter, const double& tol);
	
	std::vector<std::vector<std::vector<std::complex<double>>>> test_vectors; //test vectors[Ntest][Ntot][spin components]
private:
	GaugeConf GConf;
	double m0; 
	int Ns, Nt, Ntest;
	int Ntot;
};

//These function are provitioanl
void save_matrix(std::vector<std::vector<std::complex<double>>>& Matrix,char* Name);

void save_vector(std::vector<std::complex<double>>& Vector,char* Name);

#endif