#include "amg.h"

void PrintVector(const std::vector<std::vector<std::complex<double>>>& v ){
	//for c_matrix
    for(int i = 0; i < v.size(); i++){
        for(int j = 0; j < v[i].size(); j++){
            std::cout << v[i][j] << " ";
        }
    }
    std::cout << std::endl;
}

void PPrintVector(const std::vector<std::complex<double>>& v ){
	//for c_vector
    for(int i = 0; i < v.size(); i++){
        std::cout << v[i] << " ";
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

void normalize(std::vector<std::vector<std::complex<double>>>& v){
	std::complex<double> norm = sqrt(std::real(dot(v,v))) + 0.0*I_number; 
	v =  1.0/norm * v; 
}

//test vectors initialization
void AMG::tv_init(const double& eps,const int& Nit) {
	//eps --> the norm of the test vectors during random initialization
	//Nit --> number of iterations for improving the interpolator
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

	//Iterating over the multigrid to improve the test vectors
	for(int n=0; n<Nit; n++){
		std::vector<std::vector<std::vector<std::complex<double>>>> test_vectors_copy = test_vectors;
		for(int i = 0; i < Ntest; i++){
			//std::vector<std::vector<std::complex<double>>> x0(Ntot, std::vector<std::complex<double>>(2, 0));
			//I need to make a copy of the test vectors ...
			test_vectors_copy[i] = TwoGrid(1,1,test_vectors[i],test_vectors[i],false); 
			normalize(test_vectors_copy[i]); //normalizing the test vectors
		}
		test_vectors = test_vectors_copy;
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

//x_i = P^T_ij v_j. dim(P^T) =  Ntest Na x 2 Ntot, Nagg = block_x * block_t
//dim(v) = 2 NTot, dim(x) = Ntest Nagg
std::vector<std::vector<std::complex<double>>> AMG::Pt_v(const std::vector<std::vector<std::complex<double>>>& v) {
	std::vector<std::vector<std::complex<double>>> x(Ntest, std::vector<std::complex<double>>(Nagg, 0));
	//Loop over rows of Pt
	for (int i = 0; i < Ntest*Nagg; i++) {
		std::vector<std::vector<std::complex<double>>> e_i = canonical_vector(i, Ntest,Nagg);
		std::vector<std::vector<std::complex<double>>> Pt_row = P_v(e_i); //P^T row = P column, dimension: 2 Ntot
		int j = i / Nagg; //number of test vector
		int k = i % Nagg; //number of aggregate
		
		if (i == 0){
			std::cout << "------Pt_row------" << std::endl;
			for(int i = 0; i < Pt_row.size(); i++){
				for(int j = 0; j < Pt_row[i].size(); j++){
					std::cout << Pt_row[i][j] << " ";
				}
				std::cout << std::endl;
			}
			std::cout << "------v------" << std::endl;
			for(int i = 0; i < v.size(); i++){
				for(int j = 0; j < v[i].size(); j++){
					std::cout << v[i][j] << " ";
				}
				std::cout << std::endl;
			}
		}
		
		x[j][k] = dot_v2(Pt_row,v); //This not a complex dot product ...
	}
	return x;
}

//Dc = P^T D P, dim(Dc) = Ntest Nagg x Ntest Nagg, dim(v) = Ntest Nagg, 
std::vector<std::vector<std::complex<double>>> AMG::Pt_D_P(const std::vector<std::vector<std::complex<double>>>& v){
	return Pt_v(D_phi(GConf.Conf,P_v(v),m0));
}


//x = D^-1 phi
std::vector<std::vector<std::complex<double>>> AMG::TwoGrid(const int& nu1, const int& nu2, const std::vector<std::vector<std::complex<double>>>& x0, 
	const std::vector<std::vector<std::complex<double>>>& phi, const bool& print_message) {
	//nu1 --> pre-smoothing steps
	//nu2 --> post-smoothing steps
	//x0 --> initial guess
	//phi --> right hand side
	std::vector<std::vector<std::complex<double>>> x = bi_cgstab(GConf.Conf, phi,x0,m0,nu1,1e-10,false);
	//x = x + P*Dc^-1 * P^T * (phi-D*x); 
	std::vector<std::vector<std::complex<double>>> Pt_r = Pt_v(phi - D_phi(GConf.Conf,x,m0)); //TP^T (phi - D x)
	x = x + P_v(    bi_cgstab_Dc(GConf.Conf, Pt_r, Pt_r, m0,1000,1e-10,true)); //The bi_cgstab here has to be for DC^-1 not D^-1
	x = bi_cgstab(GConf.Conf, phi,x,m0,nu2,1e-10,false);
	double err = std::real(dot(phi - D_phi(GConf.Conf,x,m0),phi - D_phi(GConf.Conf,x,m0)));
	if (print_message == true){
		std::cout << "Error = " << err << std::endl;
	}
	return x;
}

//This version works properly (I compared with python())
c_vector AMG::bi_cg_test(const c_matrix& U, const c_matrix& phi, const c_matrix& x0, const double& m0, const int& max_iter, const double& tol){
	//Dc^-1 phi = x
	//Dc = P^T D P 
	int k = 0; //Iteration number
    double err = 1;
	c_matrix Dc(Ntest*Nagg, c_vector(Ntest*Nagg, 0));

	for(int col = 0; col < Ntest*Nagg; col++){
		c_matrix v = Pt_D_P(canonical_vector(col,Ntest,Nagg)); //column
		int count = 0;
		for(int i = 0; i < v.size(); i++){
			for(int j = 0; j < v[i].size(); j++){
				Dc[count][col] = v[i][j];
				count += 1;
			}
		}	
	} 
	
	std::cout << "------Dc------" << std::endl;
	for(int i = 0; i < Dc.size(); i++){
		for(int j = 0; j < Dc[i].size(); j++){
			std::cout << Dc[i][j] << " ";
		}
		std::cout << std::endl;
	}
	//Flatten phi
	c_vector phi_flat(Ntest*Nagg, 0);
	c_vector x0_flat(Ntest*Nagg, 0);
	int count = 0;
	for(int i = 0; i < phi.size(); i++){
		for(int j = 0; j < phi[i].size(); j++){
			x0_flat[count] = x0[i][j];
			phi_flat[count] = phi[i][j];
			count += 1;
		}
	}
	save_matrix(Dc,"Dc.dat");
	save_vector(phi_flat,"phi.dat");
	std::cout << "------phi_flat------" << std::endl;
	for(int i = 0; i < phi_flat.size(); i++){
		std::cout << phi_flat[i] << " ";
	}
	std::cout << std::endl;

	c_vector r(Ntest*Nagg, 0);
	c_vector r_tilde(Ntest*Nagg, 0);
	c_vector d(Ntest*Nagg, 0);
	c_vector s(Ntest*Nagg, 0);
	c_vector t(Ntest*Nagg, 0);	
	c_vector Ad(Ntest*Nagg, 0);
	c_vector x(Ntest*Nagg, 0);
    c_double alpha, beta, rho_i, omega, rho_i_2;
	
    x = x0_flat; //initial solution
	std::cout << "Dc phi " << std::endl;
	PPrintVector(Dc*x);
    r = phi_flat -   Dc*x; //r = b - A*x
	
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
        Ad = Dc*d;  //A d_i 
        alpha = rho_i / dot(Ad, r_tilde); //alpha_i = rho_{i-1} / (Ad_i, r_tilde)
        s = r - alpha * Ad; //s = r_{i-1} - alpha_i * Ad_i
        err = std::real(dot(s, s));
        if (err < tol) {
            x = x + alpha * d;
            return x;
        }
        t = Dc*s;   //A s
        omega = dot(s, t) / dot(t, t); //omega_i = t^dagg . s / t^dagg . t
        r = s - omega * t; //r_i = s - omega_i * t
        x = x + alpha * d + omega * s; //x_i = x_{i-1} + alpha_i * d_i + omega_i * s
        rho_i_2 = rho_i; //rho_{i-2} = rho_{i-1}
        k++;
    }
	std::cout << "Bi-CG-stab for Dc did not converge in " << max_iter << " iterations" << " Error " << err << std::endl;
    
	return x;

}

//This function is not working properly yet
c_matrix AMG::bi_cgstab_Dc(const c_matrix& U, const c_matrix& phi, const c_matrix& x0, const double& m0, const int& max_iter, const double& tol, const bool& print_message) {
    //Dc^-1 phi = x
	//Dc = P^T D P 
	int k = 0; //Iteration number
    double err = 1;


    c_matrix r(Ntest, c_vector(Nagg, 0));  //r[coordinate][spin] residual
    c_matrix r_tilde(Ntest, c_vector(Nagg, 0));  //r[coordinate][spin] residual
    c_matrix d(Ntest, c_vector(Nagg, 0)); //search direction
    c_matrix s(Ntest, c_vector(Nagg, 0));
    c_matrix t(Ntest, c_vector(Nagg, 0));
    c_matrix Ad(Ntest, c_vector(Nagg, 0)); //DD^dagger*d
    c_matrix x(Ntest, c_vector(Nagg, 0)); //solution
    c_double alpha, beta, rho_i, omega, rho_i_2;
    x = x0; //initial solution
	std::cout << "Dc phi" << std::endl;
	PrintVector(Pt_D_P(x));
    r = phi - Pt_D_P(x); //r = b - A*x //Already this operation gives something different
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
        Ad = Pt_D_P(d);  //A d_i 
        alpha = rho_i / dot(Ad, r_tilde); //alpha_i = rho_{i-1} / (Ad_i, r_tilde)
        s = r - alpha * Ad; //s = r_{i-1} - alpha_i * Ad_i
        err = std::real(dot(s, s));
        if (err < tol) {
            x = x + alpha * d;
			std::cout << "Converged in " << k << " iterations" << " Error " << err << std::endl;
            return x;
        }
        t = Pt_D_P(s);   //A s
        omega = dot(s, t) / dot(t, t); //omega_i = t^dagg . s / t^dagg . t
        r = s - omega * t; //r_i = s - omega_i * t
        x = x + alpha * d + omega * s; //x_i = x_{i-1} + alpha_i * d_i + omega_i * s
        rho_i_2 = rho_i; //rho_{i-2} = rho_{i-1}
        k++;
    }
    if (print_message == true) {
        std::cout << "Bi-CG-stab for Dc did not converge in " << max_iter << " iterations" << " Error " << err << std::endl;
    }
    return x;
}


//Complex dot product (not matrix multiplication)
c_double dot(const c_vector& x, const c_vector& y) {
    //Dot product of two vectors
    c_double z = 0;
    for (int i = 0; i < x.size(); i++) {
        z += x[i] * std::conj(y[i]);
    }
    return z;
}


c_double dot_v2(const c_matrix& x, const c_matrix& y) {
    //Dot product of two vectors
    c_double z = 0;
    for (int i = 0; i < x.size(); i++) {
        for (int j = 0; j < x[i].size(); j++) {
            z += x[i][j] * y[i][j];
        }
    }
    return z;
}

template <typename T>
c_vector operator*(const T& lambda, const c_vector& A){
    c_vector B(A.size(),0);
    for(int i = 0; i < A.size(); i++){
        B[i] = lambda*A[i];
    }
    return B;
}

c_vector operator+(const c_vector& A, const c_vector& B){
    c_vector C(A.size(),0);
    for(int i = 0; i < A.size(); i++){
        C[i] = A[i] + B[i];
    }
    return C;
}

c_vector operator-(const c_vector& A, const c_vector& B){
    c_vector C(A.size(),0);
    for(int i = 0; i < A.size(); i++){
        C[i] = A[i] - B[i];
    }
    return C;
}

c_vector operator*(const c_matrix& A, const c_vector& v) {
    c_vector w(A.size(), 0);
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            w[i] += A[i][j] * v[j];
        }
    }
    return w;
}

//These function are provitioanl
void save_matrix(std::vector<std::vector<std::complex<double>>>& Matrix,char* Name){
    char NameData[500], Data_str[500];
	sprintf(NameData, Name);
	std::ofstream Datfile;
	Datfile.open(NameData);
	for (int i = 0; i < Matrix.size(); i++) {
		for (int j = 0; j < Matrix[i].size(); j++) {
			sprintf(Data_str, "%-30d%-30d%-30.17g%-30.17g\n", i, j, std::real(Matrix[i][j]), std::imag(Matrix[i][j]));
			Datfile << Data_str;
		}
	}
	Datfile.close();
}

void save_vector(std::vector<std::complex<double>>& Vector,char* Name){
    char NameData[500], Data_str[500];
    sprintf(NameData, Name);
    std::ofstream Datfile;
    Datfile.open(NameData);
    for (int i = 0; i < Vector.size(); i++) {
        sprintf(Data_str, "%-30d%-30.17g%-30.17g\n", i, std::real(Vector[i]), std::imag(Vector[i]));
        Datfile << Data_str;
    }
    Datfile.close();
}