#include "amg.h"

void PrintVector(const std::vector<std::vector<std::complex<double>>>& v ){
    for(int i = 0; i < v.size(); i++){
        for(int j = 0; j < v[i].size(); j++){
            std::cout << v[i][j] << " ";
        }
        //std::cout << std::endl;
    }
    std::cout << std::endl;
}

//CHECK THIS FUNCTION
std::vector<std::vector<std::complex<double>>> canonical_vector(const int& i, const int& N1, const int& N2) {
	std::vector<std::vector<std::complex<double>>> e_i(N1, std::vector<std::complex<double> >(N2,0.0));
	e_i[i][0] = 1.0;
	e_i[i][1] = 0.0;
	return e_i;
}

//Random Chi vector 
inline std::vector<std::vector<std::complex<double>>> RandomChi() {
	//We use a Gaussian distribution to sample the momenta
	std::random_device rd;
	std::default_random_engine generator;
	generator.seed(rd()); //This generator has to be seeded differently. srand doesn't work here
	// The result is very much affected by the standard deviation
	std::normal_distribution<double> distribution(0.0, 0.5); //mu, standard deviation
	std::vector<std::vector<std::complex<double>>> RandPI(Ns * Nt, std::vector<std::complex<double>>(2, 0));
	for (int i = 0; i < Ns * Nt; i++) {
		for (int mu = 0; mu < 2; mu++) {
			RandPI[i][mu] = 1.0 * distribution(generator) + 0 * (0, 1);

		}
	}

	return RandPI;
}


//test vectors initialization
void AMG::tv_init(const double& eps) {
	for (int i = 0; i < Ntest; i++) {
		//std::vector<std::vector<std::complex<double>>> chi = RandomChi();
		//std::vector<std::vector<std::complex<double>>> phi = D_phi(GConf.Conf, chi, m0);
		//test_vectors[i] =  D_phi(GConf.Conf, chi, m0);
		
		for (int j = 0; j < Ntot; j++) {
			for (int k = 0; k < 2; k++) {
				test_vectors[i][j][k] = eps * RandomU1(); //epsilon fixes the norm
			}
		}
		
	}
	
	//Apply three steps of the smoother to approximately solve D x = v
	for (int n = 0; n < 3; n++) {
		for (int i = 0; i < Ntest; i++) {
			test_vectors[i] = conjugate_gradient(GConf.Conf, test_vectors[i], m0); //I have to implement Bi-GCR for this ... later SAP
		}
	}
}

//x_i = P_ij v_j. dim(P) = 2 Ntot x Ntest Na, Na = block_x * block_t
//dim(v) = Ntest Na, dim(x) = 2 Ntot
std::vector<std::vector<std::complex<double>>> AMG::P_v(const std::vector<std::vector<std::complex<double>>>& v) {
	//Prolongation operator times vector
	std::vector<std::vector<std::complex<double>>> x(Ntot, std::vector<std::complex<double>>(2, 0));
	for (int i = 0; i < Ntot; i++) {
		//These two for loops run over the number of columns of P
		for (int k = 0; k < Ntest; k++) {
			for (int a = 0; a < block_x * block_t; a++) {
				//Checking if i belongs to the aggregate	
				if (std::find(Agg[a].begin(), Agg[a].end(), i) != Agg[a].end()){
					//both spins belong to the same aggregate
					for (int alf = 0; alf < 2; alf++) {
						x[i][alf] += test_vectors[k][i][alf] * v[k][a];
					}
				}
			}
		}
	}
	return x;
}

//x_i = P^T_ij v_j. dim(P) = 2 Ntot x Ntest Na, Na = block_x * block_t
//dim(v) = Ntest Na, dim(x) = Ntest Nagg
std::vector<std::vector<std::complex<double>>> AMG::Pt_v(const std::vector<std::vector<std::complex<double>>>& v) {
	//Prolongation operator times vector
	std::vector<std::vector<std::complex<double>>> x(Ntot, std::vector<std::complex<double>>(2, 0));
	((Ntot, std::vector<std::complex<double>>(2, 0)));
	for (int i = 0; i < Ntest*block_x*block_t; i++) {
		std::vector<std::vector<std::complex<double>>> e_i = canonical_vector(i, Ntest,block_x*block_t);
		std::vector<std::vector<std::complex<double>>> P_row = P_v();
		//These two for loops run over the number of columns of Pt
		for (int k = 0; k < Ntest; k++) {
			for (int a = 0; a < block_x * block_t; a++) {
				//Checking if i belongs to the aggregate	
				//Maybe later I could do some kind of dictionary to store which rows belong to which aggregate
				//This would make the search faster after repeated calls
				if (std::find(Agg[a].begin(), Agg[a].end(), i) != Agg[a].end()){
					//both spins belong to the same aggregate
					for (int alf = 0; alf < 2; alf++) {
						x[i][alf] += test_vectors[k][i][alf] * v[k][a];
					}
				}
			}
		}
	}
	return x;
}


//void AMG::tv_update() {
	//Build interpolator
	//Apply D^-1 to the test vectors, using MULTIGRID
//}
	

