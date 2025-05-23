#include "matrix_operations.h"

std::vector<std::vector<std::vector<std::complex<double>>>> gamma_mat(2,
	std::vector<std::vector<std::complex<double>>>(2, std::vector<std::complex<double>>(2, 0))
);  //Pauli matrices
std::complex<double> I_number(0, 1); //imaginary number
std::vector<std::vector<std::complex<double>>>Identity(Ntot, std::vector<std::complex<double>>(2, 0));
std::vector<std::vector<int>>hat_mu(Ntot, std::vector<int>(2, 0)); //hat_mu[mu][2] = {hat_mu_t, hat_mu_x}

//Intialize gamma matrices, identity and unit vectors
void initialize_matrices() {
	//sigma_0 --> time
	gamma_mat[0][0][0] = 0; gamma_mat[0][0][1] = 1;
	gamma_mat[0][1][0] = 1; gamma_mat[0][1][1] = 0;
	//simgma_1 --> space
	gamma_mat[1][0][0] = 0; gamma_mat[1][0][1] = -I_number;
	gamma_mat[1][1][0] = I_number; gamma_mat[1][1][1] = 0;
	//2d identity
	Identity[0][0] = 1; Identity[0][1] = 0;
	Identity[1][0] = 0; Identity[1][1] = 1;
	hat_mu[0] = { 1, 0 }; //hat_t
	hat_mu[1] = { 0, 1 }; //hat_x
}

//D phi
std::vector<std::vector<std::complex<double>>> D_phi(const std::vector<std::vector<std::complex<double>>>& U, const std::vector<std::vector<std::complex<double>>>& phi, const double& m0) {
	int Ntot = Ns * Nt;
	std::vector<std::vector<std::complex<double>>> Dphi(Ntot, std::vector<std::complex<double>>(2, 0)); //Dphi[Ntot][2]
	for (int x = 0; x < Ns; x++) {
		for (int t = 0; t < Nt; t++) {
			int n = Coords[x][t];
			for (int alf = 0; alf < 2; alf++) {
				Dphi[n][alf] = (m0 + 2) * phi[n][alf];
				for (int bet = 0; bet < 2; bet++) {
					for (int mu = 0; mu < 2; mu++) {
						Dphi[n][alf] += -0.5 * (
							(Identity[alf][bet] - gamma_mat[mu][alf][bet]) * U[n][mu] * rfb(phi, x, t, mu, bet)
							+ (Identity[alf][bet] + gamma_mat[mu][alf][bet]) * std::conj(U[LeftPB[x][t][mu]][mu]) * lfb(phi, x, t, mu, bet)
							);
					}
				}
			}
		}
	}
	return Dphi;
}

//D^dagger phi
std::vector<std::vector<std::complex<double>>> D_dagger_phi(const std::vector<std::vector<std::complex<double>>>& U, const std::vector<std::vector<std::complex<double>>>& phi, const double& m0) {
	int Ntot = Ns * Nt;
	std::vector<std::vector<std::complex<double>>> Dphi(Ntot, std::vector<std::complex<double>>(2, 0)); //Dphi[Ntot][2]
	for (int x = 0; x < Ns; x++) {
		for (int t = 0; t < Nt; t++) {
			int n = Coords[x][t];
			for (int alf = 0; alf < 2; alf++) {
				Dphi[n][alf] = (m0 + 2) * phi[n][alf];
				for (int bet = 0; bet < 2; bet++) {
					for (int mu = 0; mu < 2; mu++) {
						Dphi[n][alf] += -0.5 * (
							(Identity[alf][bet] - gamma_mat[mu][alf][bet]) * std::conj(U[LeftPB[x][t][mu]][mu]) * lfb(phi, x, t, mu, bet)
							+ (Identity[alf][bet] + gamma_mat[mu][alf][bet]) * U[n][mu] * rfb(phi, x, t, mu, bet)
							);
					}
				}
			}
		}
	}
	return Dphi;
}

//D D^dagger phi
std::vector<std::vector<std::complex<double>>> D_D_dagger_phi(const std::vector<std::vector<std::complex<double>>>& U, const std::vector<std::vector<std::complex<double>>>& phi, const double& m0) {
	return D_phi(U, D_dagger_phi(U, phi, m0), m0);
}


	