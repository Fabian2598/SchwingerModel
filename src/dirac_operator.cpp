#include "dirac_operator.h"


std::vector<c_matrix> gamma_mat(2,
	c_matrix(2, c_vector(2, 0))
);  //Pauli matrices
c_double I_number(0, 1); //imaginary number
c_matrix Identity(LV::Ntot, c_vector(2, 0));
std::vector<std::vector<int>>hat_mu(LV::Ntot, std::vector<int>(2, 0)); //hat_mu[mu][2] = {hat_mu_t, hat_mu_x}

//Intialize gamma matrices, identity and unit vectors
void initialize_matrices() {
	/*
	sigma_0 --> time
		0 1 
		1 0
	*/
	gamma_mat[0][0][0] = 0; gamma_mat[0][0][1] = 1;
	gamma_mat[0][1][0] = 1; gamma_mat[0][1][1] = 0;
	/*
	simgma_1 --> space
	  	0 -i
  		i  0
	*/
	gamma_mat[1][0][0] = 0; gamma_mat[1][0][1] = -I_number;
	gamma_mat[1][1][0] = I_number; gamma_mat[1][1][1] = 0;
	//2d identity
	Identity[0][0] = 1; Identity[0][1] = 0;
	Identity[1][0] = 0; Identity[1][1] = 1;
	hat_mu[0] = { 1, 0 }; //hat_t
	hat_mu[1] = { 0, 1 }; //hat_x
}

//D phi
spinor D_phi(const c_matrix& U, const spinor& phi, const double& m0) {
	using namespace LV;
	spinor Dphi(Ntot, c_vector(2, 0)); //Dphi[Ntot][2]
	for (int x = 0; x < Nx; x++) {
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
spinor D_dagger_phi(const c_matrix& U, const spinor& phi, const double& m0) {
	using namespace LV;
	spinor Dphi(Ntot, c_vector(2, 0)); //Dphi[Ntot][2]
	for (int x = 0; x < Nx; x++) {
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
spinor D_D_dagger_phi(const c_matrix& U, const spinor& phi, const double& m0) {
	return D_phi(U, D_dagger_phi(U, phi, m0), m0);
}


//2* Re ( left^dag \partial D / \partial omega(z) right )
re_field phi_dag_partialD_phi(const c_matrix& U,
 const spinor& left,const spinor& right){
	using namespace LV;
	re_field Dphi(Ntot, std::vector<double>(2, 0)); //Dphi[Ntot][2]
	for (int x = 0; x < Nx; x++) {
		for (int t = 0; t < Nt; t++) {
			int n = Coords[x][t];
			for (int mu = 0; mu < 2; mu++) {
				Dphi[n][mu] = 0;
				for (int alf = 0; alf < 2; alf++) {
					for (int bet = 0; bet < 2; bet++) {
						Dphi[n][mu] +=  std::imag(
						std::conj(left[n][alf]) * (Identity[alf][bet] - gamma_mat[mu][alf][bet]) * U[n][mu] * rfb(right, x, t, mu, bet)
						- std::conj(rfb(left, x, t, mu, alf)) * (Identity[alf][bet] + gamma_mat[mu][alf][bet]) * std::conj(U[n][mu]) * right[n][bet]
						);
					}
				}
				
			}
		}
	}

	return Dphi;
 }
	