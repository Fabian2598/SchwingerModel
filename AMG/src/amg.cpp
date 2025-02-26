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
		test_vectors[i] = bi_cgstab(GConf.Conf, test_vectors[i], test_vectors[i], m0,3,1e-10,false); //The tolerance is not really relevant here
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

//Dc = P^T D P 
std::vector<std::vector<std::complex<double>>> AMG::Pt_D_P(const std::vector<std::vector<std::complex<double>>>& v){
	return Pt_v(D_phi(GConf.Conf,P_v(v),m0));
}

//x = D^-1 phi
std::vector<std::vector<std::complex<double>>> AMG::TwoGrid(const int& nu1, const int& nu2, const std::vector<std::vector<std::complex<double>>>& x0, const std::vector<std::vector<std::complex<double>>>& phi){
	//nu1 --> pre-smoothing steps
	//nu2 --> post-smoothing steps
	//x0 --> initial guess
	//phi --> right hand side
	std::vector<std::vector<std::complex<double>>> x = bi_cgstab(GConf.Conf, phi,x0,m0,nu1,1e-10,false);
	//x = x + P*Dc^-1 * P^T * (phi-D*x); 
	std::vector<std::vector<std::complex<double>>> Pt_r = Pt_v(phi - D_phi(GConf.Conf,x,m0)); //TP^T (phi - D x)
	x = x + P_v(    bi_cgstab_Dc(GConf.Conf, Pt_r, Pt_r, m0,1000,1e-10,false)); //The bi_cgstab here has to be for DC^-1 not D^-1
	x = bi_cgstab(GConf.Conf, phi,x,m0,nu2,1e-10,false);

	return x;
}



c_matrix AMG::bi_cgstab_Dc(const c_matrix& U, const c_matrix& phi, const c_matrix& x0, const double& m0, const int& max_iter, const double& tol, const bool& print_message) {
    int k = 0; //Iteration number
    double err = 1;

    //D_D_dagger_phi(U, phi, m0); //DD^dagger  
    c_matrix r(Ntest, c_vector(Nagg, 0));  //r[coordinate][spin] residual
    c_matrix r_tilde(Ntest, c_vector(Nagg, 0));  //r[coordinate][spin] residual
    c_matrix d(Ntest, c_vector(Nagg, 0)); //search direction
    c_matrix s(Ntest, c_vector(Nagg, 0));
    c_matrix t(Ntest, c_vector(Nagg, 0));
    c_matrix Ad(Ntest, c_vector(Nagg, 0)); //DD^dagger*d
    c_matrix x(Ntest, c_vector(Nagg, 0)); //solution
    c_double alpha, beta, rho_i, omega, rho_i_2;
    x = x0; //initial solution
    r = phi - Pt_D_P(x); //r = b - A*x
    r_tilde = r;
    while (k<max_iter && err>tol) {
        rho_i = dot(r, r_tilde); //r . r_dagger
        if (k == 0) {
            d = r; //d_1 = r_0
        }
        else {
            beta = alpha * rho_i / (omega * rho_i_2); //beta_{i-1} = alpha_{i-1} * rho_{i-1} / (omega_{i-1} * rho_{i-2})
            d = r + beta * (d - omega * Ad); //d_i = r_{i-1} + beta_{i-1} * (d_{i-1} - omega_{i-1} * Ad_{i-1})
        }
        Ad = Pt_D_P(x);  //A d_i 
        alpha = rho_i / dot(Ad, r_tilde); //alpha_i = rho_{i-1} / (Ad_i, r_tilde)
        s = r - alpha * Ad; //s = r_{i-1} - alpha_i * Ad_i
        err = std::real(dot(s, s));
        if (err < tol) {
            x = x + alpha * d;
            return x;
        }
        t = Pt_D_P(x);   //A s
        omega = dot(s, t) / dot(t, t); //omega_i = t^dagg . s / t^dagg . t
        r = s - omega * t; //r_i = s - omega_i * t
        x = x + alpha * d + omega * s; //x_i = x_{i-1} + alpha_i * d_i + omega_i * s
        rho_i_2 = rho_i; //rho_{i-2} = rho_{i-1}
        k++;
    }
    if (print_message == true) {
        std::cout << "Did not converge in " << max_iter << " iterations" << " Error " << err << std::endl;
    }
    return x;
}
