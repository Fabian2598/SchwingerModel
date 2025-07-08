#ifndef DIRAC_OPERATOR_INCLUDED
#define DIRAC_OPERATOR_INCLUDED
#include <complex>
#include "variables.h"
#include "operator_overloads.h"

extern std::vector<c_matrix> gamma_mat;  //Pauli matrices
extern c_double I_number; //imaginary number
extern c_matrix Identity; //2 x 2 identity matrix

//unit vectors in the "mu" direction
//mu = 0 -> time, mu = 1 -> space
extern std::vector<std::vector<int>> hat_mu; //hat_mu[mu][2] = {hat_mu_x, hat_mu_t}

/*
	Intialize gamma matrices, identity and unit vectors
	Has to be called at the beginning of the program once
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
			x_1_t1[x][t] = Coords[mod(x - 1, Nx)][mod(t + 1, Nt)];
			x1_t_1[x][t] = Coords[mod(x + 1, Nx)][mod(t - 1, Nt)];
			for (int mu = 0; mu < 2; mu++) {
				RightPB[x][t][mu] = Coords[mod(x + hat_mu[mu][1], Nx)][mod(t + hat_mu[mu][0], Nt)]; 
				LeftPB[x][t][mu] = Coords[mod(x - hat_mu[mu][1], Nx)][mod(t - hat_mu[mu][0], Nt)];


				RightPBT[Coords[x][t]][mu] = RightPB[x][t][mu];
				LeftPBT[Coords[x][t]][mu] = LeftPB[x][t][mu];
				SignRT[Coords[x][t]][mu] = (mu == 0 && t == Nt - 1) ? -1 : 1; //sign for the right boundary in time
				SignLT[Coords[x][t]][mu] = (mu == 0 && t == 0) ? -1 : 1; //sign for the left boundary in time
				
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
inline c_double rfb(const spinor& phi, const int& x, const int& t, const int& mu, const int& bet) {
	//time
	if (mu == 0) {
		if (t == LV::Nt - 1) {
			return -phi[Coords[x][0]][bet];
		}
		else {
			return phi[Coords[x][t + 1]][bet];
		}
	}
	else {
	//space
		return phi[ Coords[mod(x + 1, LV::Nx)][t] ][bet];
	}
}

/*
	Left boundary phi(n-hat{mu}) used for fermions (antiperiodic in time, periodic in space) 
	phi: spinor
	x: coordinate in the x direction
 	t: coordinate in the t direction
	mu: neighbor direction (0 for time, 1 for space)
*/
inline c_double lfb(const spinor& phi, const int& x, const int& t, const int& mu, const int& bet) {
	//time
	if (mu == 0) {
		if (t == 0) {
			return -phi[Coords[x][LV::Nt-1]][bet];
		}
		else {
			return phi[Coords[x][t-1]][bet];
		}
	}
	else {
	//space
		return phi[ Coords[mod(x - 1, LV::Nx)][t] ][bet];
	}
}

/*
	Dirac operator application D phi
	U: gauge configuration
	phi: spinor to apply the operator to
	m0: mass parameter
*/
spinor D_phi(const c_matrix& U, const spinor& phi, const double& m0);

spinor D_phi_old(const c_matrix& U, const spinor& phi, const double& m0);

/*
	Dirac dagger operator application D^+ phi
	U: gauge configuration
	phi: spinor to apply the operator to
	m0: mass parameter
*/
spinor D_dagger_phi(const c_matrix& U, const spinor& phi, const double& m0);

spinor D_dagger_phi_old(const c_matrix& U, const spinor& phi, const double& m0);

/*
	Application of D D^+
	It just calls the previous functions
*/
spinor D_D_dagger_phi(const c_matrix& U, const spinor& phi, const double& m0);

/*
	2* Re ( left^+ d D / d omega(z) right )
	This derivative is needed for the fermion force
*/
re_field phi_dag_partialD_phi(const c_matrix& U, const spinor& left, const spinor& right);

#endif