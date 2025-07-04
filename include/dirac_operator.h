#ifndef DIRAC_OPERATOR_INCLUDED
#define DIRAC_OPERATOR_INCLUDED
#include <complex>
#include "variables.h"
#include "operator_overloads.h"

extern std::vector<c_matrix> gamma_mat;  //Pauli matrices
extern c_double I_number; //imaginary number
extern c_matrix Identity;
extern std::vector<std::vector<int>> hat_mu; //hat_mu[mu][2] = {hat_mu_x, hat_mu_t}

//Intialize gamma matrices, identity and unit vectors
void initialize_matrices();

inline int mod(int a, int b) {
	int r = a % b;
	return r < 0 ? r + b : r;
}

//right periodic boundary x+hat{mu}
//left periodic boundary x-hat{mu}
//hat_mu[0] = { 1, 0 }; //hat_t
//hat_mu[1] = { 0, 1 }; //hat_x
inline void periodic_boundary() {
	using namespace LV;
	for (int x = 0; x < Nx; x++) {
		for (int t = 0; t < Nt; t++) {
			x_1_t1[x][t] = Coords[mod(x - 1, Nx)][mod(t + 1, Nt)];
			x1_t_1[x][t] = Coords[mod(x + 1, Nx)][mod(t - 1, Nt)];
			for (int mu = 0; mu < 2; mu++) {
				RightPB[x][t][mu] = Coords[mod(x + hat_mu[mu][1], Nx)][mod(t + hat_mu[mu][0], Nt)]; 
				LeftPB[x][t][mu] = Coords[mod(x - hat_mu[mu][1], Nx)][mod(t - hat_mu[mu][0], Nt)];
				
			}
		}
	}
}


//right fermionic boundary (antiperiodic in time) x+hat{mu}
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
	//periodic
		return phi[ Coords[mod(x + 1, LV::Nx)][t] ][bet];
	}
}

//left fermionic boundary (antiperiodic in time) x-hat{mu}
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
	//periodic
		return phi[ Coords[mod(x - 1, LV::Nx)][t] ][bet];
	}
}

//D phi
spinor D_phi(const c_matrix& U, const spinor& phi, const double& m0);
//D^dagger phi
spinor D_dagger_phi(const c_matrix& U, const spinor& phi, const double& m0);
//D D^dagger phi
spinor D_D_dagger_phi(const c_matrix& U, const spinor& phi, const double& m0);


//psi^dag \partial D / \partial omega(z) psi
re_field phi_dag_partialD_phi(const c_matrix& U,
 const spinor& left, const spinor& right);

#endif