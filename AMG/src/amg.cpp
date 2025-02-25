#include "amg.h"

void PrintVector(const std::vector<std::vector<std::complex<double>>>& v ){
    for(int i = 0; i < v.size(); i++){
        for(int j = 0; j < v[i].size(); j++){
            std::cout << v[i][j] << " ";
        }
    }
    std::cout << std::endl;
}

std::vector<std::vector<std::complex<double>>> canonical_vector(const int& i, const int& N1, const int& N2) {
	std::vector<std::vector<std::complex<double>>> e_i(N1, std::vector<std::complex<double> >(N2,0.0));
	int j = i / N2;
	int k = i % N2;
	e_i[j][k] = 1.0;
	return e_i;
}

//test vectors initialization
void AMG::tv_init(const double& eps) {
	//Random initialization
	for (int i = 0; i < Ntest; i++) {
		for (int j = 0; j < Ntot; j++) {
			for (int k = 0; k < 2; k++) {
				test_vectors[i][j][k] = eps * RandomU1(); //epsilon fixes the norm
			}
		}
		
	}
	//Apply three steps of the smoother to approximately solve D x = v
	for (int i = 0; i < Ntest; i++) {
		test_vectors[i] = bi_cgstab(GConf.Conf, test_vectors[i], m0,3,1e-10,false); //The tolerance is not really relevant here
	}
	
}

//x_i = P_ij v_j. dim(P) = 2 Ntot x Ntest Na, Na = block_x * block_t
//dim(v) = Ntest Na, dim(x) = 2 Ntot
std::vector<std::vector<std::complex<double>>> AMG::P_v(const std::vector<std::vector<std::complex<double>>>& v) {
	//Prolongation operator times vector
	std::vector<std::vector<std::complex<double>>> x(Ntot, std::vector<std::complex<double>>(2, 0));
	//Loop over rows of P
	for (int i = 0; i < Ntot; i++) {
		//Loops over columns
		for (int j = 0; j < Ntest*Nagg; j++) {
			int k = j / Nagg; //Number of test vector
			int a = j % Nagg; //Number of aggregate
			//Checking if i belongs to the aggregate. We do it using the sets, complexitiy is log(N).
			if (Agg_sets[a].find(i) != Agg_sets[a].end()) {
				for (int alf = 0; alf < 2; alf++) {
					x[i][alf] += test_vectors[k][i][alf] * v[k][a];
				}
			}
			
		}
	}
	return x;
}

//x_i = P^T_ij v_j. dim(P) = 2 Ntot x Ntest Na, Na = block_x * block_t
//dim(v) = 2 NTot, dim(x) = Ntest Nagg
std::vector<std::vector<std::complex<double>>> AMG::Pt_v(const std::vector<std::vector<std::complex<double>>>& v) {
	std::vector<std::vector<std::complex<double>>> x(Ntest, std::vector<std::complex<double>>(Nagg, 0));
	//Loop over rows of Pt
	for (int i = 0; i < Ntest*Nagg; i++) {
		std::vector<std::vector<std::complex<double>>> e_i = canonical_vector(i, Ntest,Nagg);
		std::vector<std::vector<std::complex<double>>> Pt_row = P_v(e_i); 
		int j = i / Nagg;
		int k = i % Nagg;
		x[j][k] = dot(Pt_row,v);
	}
	return x;
}


//void AMG::tv_update() {
	//Build interpolator
	//Apply D^-1 to the test vectors, using MULTIGRID
//}
	

