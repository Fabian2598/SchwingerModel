#ifndef DIRAC_OPERATOR_H
#define DIRAC_OPERATOR_H
#include "variables.h"
#include "operator_overloads.h"
#include "omp.h"

extern std::vector<c_matrix> gamma_mat;  //Pauli matrices
extern c_double I_number; //imaginary number
extern c_matrix Identity; // 2 x 2 identity matrix

//unit vectors in the "mu" direction
//mu = 0 -> time, mu = 1 -> space
extern std::vector<std::vector<int>> hat_mu; //hat_mu[mu][2] = {hat_mu_t, hat_mu_x}

/*
	Initialize gamma matrices, identity and unit vectors
*/
void initialize_matrices();

/*
	Modulo operation
*/
inline int mod(int a, int b) {
	int r = a % b;
	return r < 0 ? r + b : r;
}

/*
	Periodic boundary conditions used for the link variables U_mu(n).
	This function builds the arrays x_1_t1, x1_t_1, RightPB and LeftPB, which
	store the neighbor coordinates for the periodic boundary conditions.
	This prevents recalculation every time we call the operator D.
	The function is only called once at the beginning of the program.

	right periodic boundary x+hat{mu}
	left periodic boundary x-hat{mu}
	hat_mu[0] = { 1, 0 } --> hat_t
	hat_mu[1] = { 0, 1 } --> hat_x
*/
inline void periodic_boundary() {
	using namespace LV; //Lattice parameters namespace
	for (int x = 0; x < Nx; x++) {
		for (int t = 0; t < Nt; t++) {
			for (int mu = 0; mu < 2; mu++) {
				RightPB[x][t][mu] = Coords[mod(x + hat_mu[mu][1], Nx)][mod(t + hat_mu[mu][0], Nt)]; 
				LeftPB[x][t][mu] = Coords[mod(x - hat_mu[mu][1], Nx)][mod(t - hat_mu[mu][0], Nt)];	
			}
		}
	}
}


/*
	Right boundary phi(n+hat{mu}) used for fermions (antiperiodic in time, periodic in space) 
	phi: spinor
	x: coordinate in the x direction
 	t: coordinate in the t direction
	mu: neighbor direction (0 for time, 1 for space)
*/
inline std::complex<double> rfb(const spinor& phi, const int& x, const int& t, const int& mu, const int& bet) {
	//time
	using namespace LV;	
	if (mu == 0) {
		if (t == Nt - 1) {
			return -phi[Coords[x][0]][bet];
		}
		else {
			return phi[Coords[x][t + 1]][bet];
		}
	}
	else {
	//space
		return phi[ Coords[mod(x + 1, Nx)][t] ][bet];
	}
}

/*
	Left boundary phi(n-hat{mu}) used for fermions (antiperiodic in time, periodic in space) 
	phi: spinor
	x: coordinate in the x direction
 	t: coordinate in the t direction
	mu: neighbor direction (0 for time, 1 for space)
*/
inline std::complex<double> lfb(const spinor& phi, const int& x, const int& t, const int& mu, const int& bet) {
	//time
	using namespace LV;	
	if (mu == 0) {
		if (t == 0) {
			return -phi[Coords[x][Nt-1]][bet];
		}
		else {
			return phi[Coords[x][t-1]][bet];
		}
	}
	else {
	//space
		return phi[ Coords[mod(x - 1, Nx)][t] ][bet];
	}
}

/*
	Dirac operator application D phi
	U: gauge configuration
	phi: spinor to apply the operator to
	m0: mass parameter
*/
spinor D_phi(const c_matrix& U, const spinor& phi, const double& m0);

#endif