#ifndef MATRIX_OPERATIONS_INCLUDED
#define MATRIX_OPERATIONS_INCLUDED
#include <complex>
#include "variables.h"

constexpr int Ntot = Ns * Nt;
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
inline int rpb(const int& x, const int& t, const int& mu) {
	return Coords[mod(x + hat_mu[mu][1], Ns)][mod(t + hat_mu[mu][0], Nt)];
}
//left periodic boundary x-hat{mu}
inline int lpb(const int& x, const int& t, const int& mu) {
	return Coords[mod(x - hat_mu[mu][1], Ns)][mod(t - hat_mu[mu][0], Nt)];
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