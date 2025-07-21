#ifndef AMG_H_INCLUDED
#define AMG_H_INCLUDED
#include "gauge_conf.h"
#include "operator_overloads.h"
#include "gmres.h"
#include "sap.h"
#include <algorithm>
#include <random>

/*
	Build the aggregates A_j_0 = L_j x {0}, A_j_1 = L_j x {1} based on the lattice blocking
	indicated in the CMakeLists.txt. It is only called once at the beginning of the program.
*/
void Aggregates();

/*
	Print the aggregates. Useful for debugging.
*/
void PrintAggregates();

/*
	Normalize a spinor.
*/
void normalize(spinor& v); 

/*
	Canonical vector written in spinor notation
*/
spinor canonical_vector(const int& i, const int& N1, const int& N2);

/*
	Two-grid method for solving the linear system D x = phi
*/
class AMG {
public:
	/*
	GaugeConf GConf: Gauge configuration
    m0: Mass parameter for the Dirac matrix
    nu1: Number of pre-smoothing iterations
    nu2: Number of post-smoothing iterations
	test_vectors: Test vectors for the AMG method
	interpolator_columns: locally orthonormalized columns of the interpolator
	*/
	AMG(const GaugeConf & GConf, const double& m0, const int& nu1, const int& nu2) : GConf(GConf), m0(m0), nu1(nu1), nu2(nu2) {	
		test_vectors = std::vector<spinor>(AMGV::Ntest,
		spinor( LV::Ntot, c_vector (2,0))); 

		interpolator_columns = std::vector<spinor>(AMGV::Ntest,
			spinor( LV::Ntot, c_vector (2,0))); 

		valuesDc = c_vector(AMGV::Ntest * AMGV::Nagg * AMGV::Ntest * AMGV::Nagg, 0.0); //CSR matrix for the coarse grid operator
		rowsDc = std::vector<int>(AMGV::Ntest * AMGV::Nagg * AMGV::Ntest * AMGV::Nagg, 0); //Row indices of the coarse grid operator
		colsDc = std::vector<int>(AMGV::Ntest * AMGV::Nagg * AMGV::Ntest * AMGV::Nagg,0);
		v_chopped = std::vector<spinor>(AMGV::Ntest*AMGV::Nagg, spinor(LV::Ntot, c_vector(2,0)));
		//c_matrix DcMatrix = c_matrix(AMGV::Ntest*AMGV::Nagg, c_vector(AMGV::Ntest*AMGV::Nagg,0));
		P_TEMP = spinor(LV::Ntot, c_vector(2,0)); //Temporary spinor for the coarse grid operator
		Pt_TEMP = spinor(AMGV::Ntest, c_vector(AMGV::Nagg, 0)); //Temporary spinor for the coarse grid operator
	}
	~AMG() { };

	/*
	Set up phase of the AMG method
		Intializes the test vectors, v_i, with random values and updates them with the SAP method by performing 
		AMGV::SAP_test_vectors_iterations iterations to approximately solve the linear system D x = v_i.
	    The test vectors are then locally orthonormalized according to the aggregation.abort

	eps: Norm of the test vectors during random initialization
	Nit: Number of iterations for improving the interpolator
	*/
	void setUpPhase(const double& eps, const int& Nit);

	/*
	Two-grid method for solving the linear system D x = phi
		It can also be used as preconditioner for FGMRES

	max_iter: Maximum number of two-grid iterations
	tol: Tolerance for the solver
	x0: Initial guess
	phi: Right hand side
	print_message: If true, print the convergence message

	The convergence criterion is ||D x - phi|| < ||phi|| * tol
	*/
	spinor TwoGrid(const int& max_iter, const double& tol,
		const spinor& x0, const spinor& phi,const bool& print_message);

private:
	GaugeConf GConf;
	double m0; 
	int nu1, nu2; 
	int nonzero = 0; //Count the number of non-zero elements in the coarse grid operator
	c_vector valuesDc; //CSR matrix for the coarse grid operator
	std::vector<int> rowsDc;
	std::vector<int> colsDc;
	 
	//c_matrix DcMatrix;

	/*
	Interpolator times a spinor
	*/
	void P_v(const spinor& v,spinor& out); //P v

	/*
	Restriction operator times a spinor
	*/
	void Pt_v(const spinor& v,spinor& out); // P^T v

	/*
	Assemble the coarse grid operator Dc = P^H D P 
	dim(Dc) = Ntest Nagg x Ntest Nagg
	*/
	void assembleDc();

	/*
	Coarse grid operator Dc = P^H D P times a spinor
    */
	void Pt_D_P(const spinor& v,spinor& out); //Dc v = P^H D P v 

	/*
	Local orthonormalization of the test vectors
	*/
	void orthonormalize(); 

	//-----Coarse grid solvers-----// 

	/*
	bi_cgstab for the coarse grid operator

	U: gauge configuration
 	phi: right hand side
 	x0: initial guess
 	m0: mass parameter
	max_iter: maximum number of iterations
 	tol: tolerance for the solver 
 	print_message: if true, print the convergence message

	The convergence criterion is ||D x - phi|| < ||phi|| * tol
	*/
	spinor bi_cgstab(const c_matrix& U, const spinor& phi, const spinor& x0,
		 const double& m0, const int& max_iter, const double& tol, const bool& print_message); //Dc^-1 phi
	
	/*
	gmres for the coarse grid operator.
		It can also be used for smoothing at the fine level. 

	dim1: dimension of the first index of the spinor 
  	dim2: dimension of the second index of the spinor 	
		-> For instance, at the finest level dim1 = Nx*Nt, dim2 = 2 (spin components).
 		   For the coarse level, dim1 = Ntest, dim2 = Nagg (test vectors and aggregates).

	U: gauge configuration
 	phi: right hand side
 	x0: initial guess
 	m0: mass parameter
	m: restart length
  	restarts: number of restarts of length m
 	tol: tolerance for the solver 
 	print_message: if true, print the convergence message

	The convergence criterion is ||D x - phi|| < ||phi|| * tol
	*/
	spinor gmres(const int& dim1, const int& dim2,const c_matrix& U, const spinor& phi, const spinor& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message);

	std::vector<spinor> test_vectors; //test vectors[Ntest][Nx Nt][spin components], no color
	std::vector<spinor> interpolator_columns; 
	std::vector<spinor> v_chopped;
	spinor Pt_TEMP;
	spinor P_TEMP;
	
};

/*
	This function is for writing a spinor to a file. Useful for testing.
*/
void save_spinor(spinor& phi,char* Name);


#endif