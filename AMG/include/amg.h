#ifndef AMG_H_INCLUDED
#define AMG_H_INCLUDED
#include "gauge_conf.h"
#include "conjugate_gradient.h"
#include "operator_overloads.h"
#include "gmres.h"
#include <algorithm>
#include <random>

void PrintVector(const c_matrix& v );
void PrintVector(const c_vector& v );

void normalize(c_matrix& v);

class AMG {
public:
	AMG(const GaugeConf & GConf,const int& Ns, const int& Nt, const int& Ntest, const double& m0, const int& nu1, const int& nu2, const int& rpc) : GConf(GConf), 
	Ns(Ns), Nt(Nt), Ntot(Ns*Nt), Ntest(Ntest), m0(m0), nu1(nu1), nu2(nu2) {	
		test_vectors = std::vector<c_matrix>(Ntest,
		c_matrix( Ntot, c_vector (2,0))); //test_vectors[Number of test vectors][Ns Nt][2], two spin directions, no color
		std::vector<c_matrix> test_vectors_copy = std::vector<c_matrix>(Ntest,
			c_matrix( Ntot, c_vector (2,0))); 
	}
	~AMG() { };

	void tv_init(const double& eps, const int& Nit);
	c_matrix TwoGrid(const int& max_iter, const int& rpc,const double& tol,
		const c_matrix& x0, const c_matrix& phi,const bool& print_message);

	
private:
	GaugeConf GConf;
	double m0; 
	int Ns, Nt, Ntest;
	int Ntot;
	int nu1, nu2; //pre and post smoothing iterations

	c_matrix P_v(const c_matrix& v); //P v
	c_matrix Pt_v(const c_matrix& v); // Pt v
	c_matrix Pt_D_P(const c_matrix& v); //Dc v = P^T D P v 

	void orthonormalize(); //Orthonormalize the test vectors

	c_matrix bi_cgstab_Dc(const c_matrix& U, const c_matrix& phi, const c_matrix& x0,
		 const double& m0, const int& max_iter, const double& tol, const bool& print_message); //Dc^-1 phi
	
	std::vector<c_matrix> test_vectors; //test vectors[Ntest][Ntot][spin components]
	std::vector<c_matrix> test_vectors_copy; 
};

//These functions are provitioanl
void save_matrix(c_matrix& Matrix,char* Name);
void save_vector(c_vector& Vector,char* Name);

#endif