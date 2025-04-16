#include "amg.h"

void PrintVector(const c_matrix& v ){
	//for c_matrix
    for(int i = 0; i < v.size(); i++){
        for(int j = 0; j < v[i].size(); j++){
            std::cout << v[i][j] << " ";
        }
    }
    std::cout << std::endl;
}

void PrintVector(const c_vector& v ){
	//for c_vector
    for(int i = 0; i < v.size(); i++){
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;
}

void normalize(c_matrix& v){
	c_double norm = sqrt(std::real(dot(v,v))) + 0.0*I_number; 
	v =  1.0/norm * v; 
}

//Orthonormalize the P columns.
void AMG::orthonormalize(){
	//Gram-Schmidt orthonormalization
	std::vector<c_matrix> temp(Ntest, c_matrix( Ntot, c_vector (2,0)));
	std::vector<c_matrix> v_chopped(Ntest*Nagg, c_matrix(Ntot, c_vector(2,0)));
	for(int i = 0; i < Ntest*Nagg; i++){
		c_matrix e_i = canonical_vector(i, Ntest, Nagg);
		v_chopped[i] = P_v(e_i); //Columns of the interpolator
	}

	
	//clock_t begin = clock();
	for (int i = 0; i < Nagg; i++) {
		for (int nt = 0; nt < Ntest; nt++) {
			for (int j = 0; j < nt; j++) {
				c_double proj = dot(v_chopped[nt*Nagg+i], v_chopped[j*Nagg+i]);
				v_chopped[nt*Nagg+i] = v_chopped[nt*Nagg+i] - proj * v_chopped[j*Nagg+i];
			}
			normalize(v_chopped[nt*Nagg+i]);
		}
	}
	//clock_t end = clock();
	//double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//std::cout << "Time for orthonormalizing V_chopped with overload " << elapsed_secs << " s" << std::endl;
	
	for(int i = 0; i < Ntest; i++){
		for(int j = 0; j < Nagg; j++){
			for(int n = 0; n < Ntot; n++){
				for(int alf=0; alf<2; alf++){
					temp[i][n][alf] += v_chopped[i*Nagg + j][n][alf]; 
				}
			}
			
		}
	}
	test_vectors_copy = temp;
}; 

//test vectors initialization
void AMG::tv_init(const double& eps,const int& Nit) {
	//eps --> the norm of the test vectors during random initialization
	//Nit --> number of iterations for improving the interpolator
	//Random initialization
	for (int i = 0; i < Ntest; i++) {
		for (int j = 0; j < Ntot; j++) {
			for (int k = 0; k < 2; k++) {
				test_vectors[i][j][k] = eps * RandomU1();//i * Ntot * 2 + j * 2 + k + 1; 
			}
		}
	}

	for (int i = 0; i < Ntest; i++) {
		//20 restarts, 20 iterations per restart
		test_vectors[i] = gmres(GConf.Conf, test_vectors[i], test_vectors[i], m0, 20, 20, 1e-10, false);
	}

	test_vectors_copy = test_vectors; //This step essentially assembles the inerpolator
	orthonormalize(); //This modifies test_vectors_copy which is used in the interpolator NOT test_vectors
	std::cout << "Improving interpolator" << std::endl;
	for (int n = 0; n < Nit; n++) {
		std::cout << "****** Nit " << n << " ******" << std::endl;
		for (int i = 0; i < Ntest; i++) {
			test_vectors[i] = TwoGrid(1,20, 1e-10, test_vectors[i], test_vectors[i], false); //Updates the test vectors with the current P
		}
		test_vectors_copy = test_vectors; //"Assemble" interpolator
		orthonormalize();	
	}

}

//x_i = P_ij v_j. dim(P) = 2 Ntot x Ntest Na, Na = block_x * block_t
//dim(v) = Ntest Na, dim(x) = 2 Ntot
c_matrix AMG::P_v(const c_matrix& v) {
	//Prolongation operator times vector
	c_matrix x(Ntot, c_vector(2, 0));
	//Loop over columns
	for (int j = 0; j < Ntest * Nagg; j++) {
		int k = j / Nagg; //Number of test vector
		int a = j % Nagg; //Number of aggregate
		for (int i = 0; i < Agg[a].size(); i++) {
			if (Nagg ==  block_x * block_t){
				//Aggregation scheme 1
				for (int alf = 0; alf < 2; alf++) {
					x[Agg[a][i]][alf] += test_vectors_copy[k][Agg[a][i]][alf] * v[k][a];
				}
			}
			else{
				//Aggregation scheme 2 
				int x_coord = XCoord[Agg[a][i]], t_coord = TCoord[Agg[a][i]], s_coord = SCoord[Agg[a][i]];
				x[Coords[x_coord][t_coord]][s_coord] += test_vectors_copy[k][Coords[x_coord][t_coord]][s_coord] * v[k][a];		
			}
		}
	}
	return x;
}

//x_i = P^H_ij v_j. dim(P^H) =  Ntest Na x 2 Ntot, Nagg = block_x * block_t
//dim(v) = 2 NTot, dim(x) = Ntest Nagg
c_matrix AMG::Pt_v(const c_matrix& v) {
	c_matrix x(Ntest, c_vector(Nagg, 0));
	for (int i = 0; i < Ntest*Nagg; i++) {
		int k = i / Nagg; //number of test vector
		int a = i % Nagg; //number of aggregate
		for (int j = 0; j < Agg[a].size(); j++) {
			if (Nagg == block_x * block_t){
				//Aggregation scheme 1
				for (int alf = 0; alf < 2; alf++) {
					x[k][a] += std::conj(test_vectors_copy[k][Agg[a][j]][alf]) * v[Agg[a][j]][alf];
				}
			}
			else{
				//Aggregation scheme 2
				int x_coord = XCoord[Agg[a][j]], t_coord = TCoord[Agg[a][j]], s_coord = SCoord[Agg[a][j]];
				x[k][a] += std::conj(test_vectors_copy[k][Coords[x_coord][t_coord]][s_coord]) * v[Coords[x_coord][t_coord]][s_coord];
			}
		}
	}
	return x;
}


//Dc = P^H D P, dim(Dc) = Ntest Nagg x Ntest Nagg, dim(v) = Ntest Nagg, 
c_matrix AMG::Pt_D_P(const c_matrix& v){
	return Pt_v(D_phi(GConf.Conf,P_v(v),m0));
}

//x = D^-1 phi
c_matrix AMG::TwoGrid(const int& max_iter, const int& rpc, const double& tol, const c_matrix& x0, 
	const c_matrix& phi, const bool& print_message) {
	//nu1 --> pre-smoothing steps
	//nu2 --> post-smoothing steps
	//x0 --> initial guess
	//phi --> right hand side
	//max_iter --> maximum number of iterations
	//tol --> tolerance
	c_matrix x = x0;
	c_matrix r;
	double err = 1;
	int k = 0;
	double norm = sqrt(std::real(dot(phi,phi)));
	while(k < max_iter && err > tol*norm){
		//Pre-smoothing
		if (nu1>0){
			//x = kaczmarz(GConf.Conf, phi,x,m0,nu1,1e-10,false);
			x = gmres(GConf.Conf, phi, x, m0, rpc, nu1, 1e-10, false);

		} 
		//x = x + P*Dc^-1 * P^H * (phi-D*x);  Coarse grid correction
		c_matrix Pt_r = Pt_v(phi - D_phi(GConf.Conf,x,m0)); //P^H (phi - D x)

		x = x + P_v(bi_cgstab_Dc(GConf.Conf, Pt_r, Pt_r, m0,1000,1e-10,true)); //The bi_cgstab here is for Dc. This matrix is 
		//highly ill-conditioned, but we don't need a very precise solution for the coarse system. 
		
		//Post-smoothing
		if (nu2>0){
			x = gmres(GConf.Conf, phi, x, m0, rpc, nu2, 1e-10, false);
		}
		r = phi - D_phi(GConf.Conf, x, m0);
		err = sqrt(std::real(dot(r,r)));
		if (print_message == true){
			std::cout << "Two-grid method " << k+1 << " iterations" << " Error " << err << std::endl;
		}
		if (err < tol*norm){
			it_count = k + 1;
			if (print_message == true){
				std::cout << "Two-grid method converged in " << k+1 << " iterations" << " Error " << err << std::endl;
			}
		return x;
		}
		k++;
	}
	it_count = max_iter;
	if (print_message == true){
		std::cout << "Two-grid did not converge in " << max_iter << " iterations" << " Error " << err << std::endl;
	}
	return x;
}

//Bi-cgstab for Dc^-1 phi = x
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
    r = phi - Pt_D_P(x); //r = b - A*x 
    r_tilde = r;
	double norm_phi = sqrt(std::real(dot(phi, phi))); //norm of the right hand side
    while (k<max_iter && err>tol*norm_phi) {
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
        err = sqrt(std::real(dot(s, s)));
        if (err < tol*norm_phi) {
            x = x + alpha * d;
			if (print_message == true) {
				std::cout << "Bi-CG-stab for Dc converged in " << k << " iterations" << " Error " << err << std::endl;
			}
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


//These functions are for saving the matrices. Useful for testing.
void save_matrix(c_matrix& Matrix,char* Name){
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

void save_vector(c_vector& Vector,char* Name){
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