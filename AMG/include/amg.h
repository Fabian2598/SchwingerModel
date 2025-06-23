#ifndef AMG_H_INCLUDED
#define AMG_H_INCLUDED
#include "gauge_conf.h"
#include "operator_overloads.h"
#include "gmres.h"
#include "sap.h"
#include <algorithm>
#include <random>

void Aggregates(); //Aggregates A_j_0 = L_j x {0}, A_j_1 = L_j x {1}
void PrintAggregates();
void normalize(spinor& v); //Normalize a spinor  


class AMG {
public:
	AMG(const GaugeConf & GConf, const double& m0, const int& nu1, const int& nu2) : GConf(GConf), m0(m0), nu1(nu1), nu2(nu2) {	
		test_vectors = std::vector<spinor>(AMGV::Ntest,
		spinor( LV::Ntot, c_vector (2,0))); //test_vectors[Number of test vectors][Nx Nt][2], two spin directions, no color
		std::vector<spinor> test_vectors_copy = std::vector<spinor>(AMGV::Ntest,
			spinor( LV::Ntot, c_vector (2,0))); 
	}

	
	~AMG() { };

	void tv_init(const double& eps, const int& Nit);
	spinor TwoGrid(const int& max_iter, const double& tol,
		const spinor& x0, const spinor& phi,const bool& print_message);

private:
	GaugeConf GConf;
	double m0; 
	int nu1, nu2; //pre and post smoothing iterations

	spinor P_v(const spinor& v); //P v
	spinor Pt_v(const spinor& v); // Pt v
	spinor Pt_D_P(const spinor& v); //Dc v = P^T D P v 

	void orthonormalize(); //Orthonormalize the test vectors

	//These are the same iterative methods I defined for D^-1 phi. In principle I should be a able to pass
	//the matrix-vector operation as a function pointer. However, this is not so straightforward for 
	//member functions
	//Coarse grid solvers
	spinor bi_cgstab_Dc(const c_matrix& U, const spinor& phi, const spinor& x0,
		 const double& m0, const int& max_iter, const double& tol, const bool& print_message); //Dc^-1 phi
	spinor gmres(const int& dim1, const int& dim2,const c_matrix& U, const spinor& phi, const spinor& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message);

		
	
	std::vector<spinor> test_vectors; //test vectors[Ntest][Ntot][spin components]
	std::vector<spinor> test_vectors_copy; 
	
};

//Save spinor
void save_spinor(spinor& phi,char* Name);


#endif