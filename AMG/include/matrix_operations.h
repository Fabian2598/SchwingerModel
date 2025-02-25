#ifndef MATRIX_OPERATIONS_INCLUDED
#define MATRIX_OPERATIONS_INCLUDED
#include "variables.h"

extern std::vector<std::vector<std::vector<std::complex<double>>>> gamma_mat;  //Pauli matrices
extern std::complex<double> I_number; //imaginary number
extern std::vector<std::vector<std::complex<double>>> Identity;
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
	for (int x = 0; x < Ns; x++) {
		for (int t = 0; t < Nt; t++) {
			x_1_t1[x][t] = Coords[mod(x + 1, Ns)][mod(t - 1, Nt)];
			x1_t_1[x][t] = Coords[mod(x - 1, Ns)][mod(t + 1, Nt)];
			for (int mu = 0; mu < 2; mu++) {
				RightPB[x][t][mu] = Coords[mod(x + hat_mu[mu][1], Ns)][mod(t + hat_mu[mu][0], Nt)]; 
				LeftPB[x][t][mu] = Coords[mod(x - hat_mu[mu][1], Ns)][mod(t - hat_mu[mu][0], Nt)];
				
			}
		}
	}
}


//right fermionic boundary (antiperiodic in time) x+hat{mu}
inline std::complex<double> rfb(const std::vector<std::vector<std::complex<double>>>& phi, const int& x, const int& t, const int& mu, const int& bet) {
	//time
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
		return phi[ Coords[mod(x + 1, Ns)][t] ][bet];
	}
}

//left fermionic boundary (antiperiodic in time) x-hat{mu}
inline std::complex<double> lfb(const std::vector<std::vector<std::complex<double>>>& phi, const int& x, const int& t, const int& mu, const int& bet) {
	//time
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
		return phi[ Coords[mod(x - 1, Ns)][t] ][bet];
	}
}

//D phi
std::vector<std::vector<std::complex<double>>> D_phi(const std::vector<std::vector<std::complex<double>>>& U, const std::vector<std::vector<std::complex<double>>>& phi, const double& m0);
//D^dagger phi
std::vector<std::vector<std::complex<double>>> D_dagger_phi(const std::vector<std::vector<std::complex<double>>>& U, const std::vector<std::vector<std::complex<double>>>& phi, const double& m0);
//D D^dagger phi
std::vector<std::vector<std::complex<double>>> D_D_dagger_phi(const std::vector<std::vector<std::complex<double>>>& U, const std::vector<std::vector<std::complex<double>>>& phi, const double& m0);


//psi^dag \partial D / \partial omega(z) psi
std::vector<std::vector<double>> phi_dag_partialD_phi(const std::vector<std::vector<std::complex<double>>>& U,
 const std::vector<std::vector<std::complex<double>>>& left, const std::vector<std::vector<std::complex<double>>>& right);

#endif