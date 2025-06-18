#ifndef MATRIX_OPERATIONS_INCLUDED
#define MATRIX_OPERATIONS_INCLUDED
#include "variables.h"
#include "operator_overloads.h"

extern std::vector<c_matrix> gamma_mat;  //Pauli matrices
extern c_double I_number; //imaginary number
extern c_matrix Identity;
extern std::vector<std::vector<int>> hat_mu; //hat_mu[mu][2] = {hat_mu_x, hat_mu_t}

spinor canonical_vector(const int& i, const int& N1, const int& N2);

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
			x_1_t1[x][t] = Coords[mod(x + 1, Nx)][mod(t - 1, Nt)];
			x1_t_1[x][t] = Coords[mod(x - 1, Nx)][mod(t + 1, Nt)];
			for (int mu = 0; mu < 2; mu++) {
				RightPB[x][t][mu] = Coords[mod(x + hat_mu[mu][1], Nx)][mod(t + hat_mu[mu][0], Nt)]; 
				LeftPB[x][t][mu] = Coords[mod(x - hat_mu[mu][1], Nx)][mod(t - hat_mu[mu][0], Nt)];
				
			}
		}
	}
}


//right fermionic boundary (antiperiodic in time) x+hat{mu}
inline std::complex<double> rfb(const std::vector<std::vector<std::complex<double>>>& phi, const int& x, const int& t, const int& mu, const int& bet) {
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
	//periodic
		return phi[ Coords[mod(x + 1, Nx)][t] ][bet];
	}
}

//left fermionic boundary (antiperiodic in time) x-hat{mu}
inline std::complex<double> lfb(const std::vector<std::vector<std::complex<double>>>& phi, const int& x, const int& t, const int& mu, const int& bet) {
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
	//periodic
		return phi[ Coords[mod(x - 1, Nx)][t] ][bet];
	}
}

//D phi
spinor D_phi(const std::vector<std::vector<std::complex<double>>>& U, const std::vector<std::vector<std::complex<double>>>& phi, const double& m0);

#endif