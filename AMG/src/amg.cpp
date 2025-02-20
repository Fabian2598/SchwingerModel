#include "amg.h"

//test vectors initialization
void AMG::tv_init(const double& eps) {
	for (int i = 0; i < Ntest; i++) {
		for (int j = 0; j < Ntot; j++) {
			for (int k = 0; k < 2; k++) {
				test_vectors[i][j][k] = eps * RandomU1(); //epsilon fixes the norm
			}
		}
	}
	//Apply three steps of the smoother to approximately solve D x = v
	for (int n = 0; n < 3; n++) {
		for (int i = 0; i < Ntest; i++) {
			test_vectors[i] = conjugate_gradient_D(GConf.Conf, test_vectors[i], m0);
		}
	}
}

//x_i = P_ij v_j. dim(P) = 2 Ntot x Ntest Na, Na = block_x * block_t
std::vector<std::vector<std::complex<double>>> AMG::P_v(const std::vector<std::vector<std::complex<double>>>& v) {
	//Prolongation operator times vector
	std::vector<std::vector<std::complex<double>>> x(Ntot, std::vector<std::complex<double>>(2, 0));
	for (int i = 0; i < 2 * Ntot; i++) {
		//These two for loops run over the number of columns of P
		for (int k = 0; k < Ntest; k++) {
			for (int a = 0; a < block_x * block_t; a++) {
				//Checking if i belongs to the aggregate
				if (std::find(Agg[a].begin(), Agg[a].end(), i) != Agg[a].end()){
					//both spins belong to the same aggregate
					for (int alf = 0; alf < 2; alf++) {
						x[i][alf] += test_vectors[k][i][alf] * v[i][alf];
					}
				}
			}
		}
	}
	return x;
}

void AMG::tv_update() {
	//Build interpolator
	//Apply D^-1 to the test vectors, using MULTIGRID
}
	

