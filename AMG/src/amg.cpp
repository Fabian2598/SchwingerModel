#include "amg.h"

//Aggregates A_j_0 = L_j x {0}, A_j_1 = L_j x {1} (Volume times Spin component)
void Aggregates() {
	using namespace LV;
	for (int x = 0; x < block_x; x++) {
		for (int t = 0; t < block_t; t++) {
			for (int s = 0; s < 2; s++) {
				//vectorize the coordinates of these three loops
				int x0 = x * x_elements, t0 = t * t_elements;
				int x1 = (x + 1) * x_elements, t1 = (t + 1) * t_elements;
				int aggregate = x * block_t * 2 + t * 2 + s;
				int count = 0;
				//x and t are redefined in the following loop
				for (int x = x0; x < x1; x++) {
					for (int t = t0; t < t1; t++) {
						int i = x * Nt * 2  + t * 2 + s;
						XCoord[i] = x; TCoord[i] = t; SCoord[i] = s;
						Agg[aggregate][count] = i;
						count++;
					}
				}
				if (count != x_elements * t_elements) {
					std::cout << "Aggregate " << aggregate << " has " << count << " elements" << std::endl;
				}
				//Once the loops are finished count should be x_elements*t_elements
			}	
		}
	}
	AMGV::aggregates_initialized = true; //Set aggregates as initialized
}


void normalize(spinor& v){
	c_double norm = sqrt(std::real(dot(v,v))) + 0.0*I_number; 
	v =  1.0/norm * v; 
}


//Print aggregates. Useful for debugging
void PrintAggregates() {
	for (int i = 0; i < AMGV::Nagg; i++) {
		std::cout << "-------Aggregate-----" << i << std::endl;
		for (int j = 0; j < LV::x_elements * LV::t_elements; j++) {
			std::cout << Agg[i][j] << " ";
		}
		std::cout << std::endl;
	}
}


//Orthonormalize the P columns.
void AMG::orthonormalize(){
	//Gram-Schmidt orthonormalization
	using namespace AMGV;
	std::vector<spinor> temp(Ntest, spinor( LV::Ntot, c_vector (2,0)));
	std::vector<spinor> v_chopped(Ntest*Nagg, spinor(LV::Ntot, c_vector(2,0)));
	for(int i = 0; i < Ntest*Nagg; i++){
		spinor e_i = canonical_vector(i, Ntest, Nagg);
		v_chopped[i] = P_v(e_i); //Columns of the interpolator
	}


	for (int i = 0; i < Nagg; i++) {
		for (int nt = 0; nt < Ntest; nt++) {
			for (int j = 0; j < nt; j++) {
				c_double proj = dot(v_chopped[nt*Nagg+i], v_chopped[j*Nagg+i]);
				v_chopped[nt*Nagg+i] = v_chopped[nt*Nagg+i] - proj * v_chopped[j*Nagg+i];
			}
			normalize(v_chopped[nt*Nagg+i]);
		}
	}

	
	for(int i = 0; i < Ntest; i++){
		for(int j = 0; j < Nagg; j++){
			for(int n = 0; n < LV::Ntot; n++){
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
	int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	using namespace AMGV;
	for (int i = 0; i < Ntest; i++) {
		for (int j = 0; j < LV::Ntot; j++) {
			for (int k = 0; k < 2; k++) {
				test_vectors[i][j][k] = eps * RandomU1();
			}
		}
	}

	for (int i = 0; i < Ntest; i++) {
		//gmres initialization
		//test_vectors[i] = gmres(LV::Ntot,2,GConf.Conf, test_vectors[i], test_vectors[i], m0, 20, 20, 1e-10, false);

		//c_matrix v0(Ntot, c_vector(2, 0)); //Initial guess (zero as right hand side)
		spinor v0 = test_vectors[i]; //Initial guess
		SAP_parallel(GConf.Conf, v0, test_vectors[i], m0, AMGV::SAP_test_vectors_iterations,SAPV::sap_blocks_per_proc); 
		//SAP(GConf.Conf, v0, test_vectors[i], m0, AMGV::SAP_test_vectors_iterations);
	}

	test_vectors_copy = test_vectors; //This step essentially assembles the interpolator
	orthonormalize(); //This modifies test_vectors_copy which is used in the interpolator NOT test_vectors
	if (rank == 0){std::cout << "Improving interpolator" << std::endl;}
	
	for (int n = 0; n < Nit; n++) {
		if (rank == 0){std::cout << "****** Bootstrap iteration " << n << " ******" << std::endl;}
		
		for (int i = 0; i < Ntest; i++) {
			test_vectors[i] = TwoGrid(1,1e-10, test_vectors[i], test_vectors[i], false); 
			//Update the test vectors with the current P
		}
		test_vectors_copy = test_vectors; //"Assemble" interpolator
		orthonormalize();	
	}
	if (rank == 0){std::cout << "Set-up phase finished" << std::endl;}
	

}

//x_i = P_ij v_j. dim(P) = 2 Ntot x Ntest Na, Na = block_x * block_t
//dim(v) = Ntest Na, dim(x) = 2 Ntot
spinor AMG::P_v(const spinor& v) {
	//Prolongation operator times a spinor
	spinor x(LV::Ntot, c_vector(2, 0));
	//Loop over columns
	using namespace AMGV;
	for (int j = 0; j < Ntest * Nagg; j++) {
		int k = j / Nagg; //Number of test vector
		int a = j % Nagg; //Number of aggregate
		for (int i = 0; i < Agg[a].size(); i++) {
			int x_coord = XCoord[Agg[a][i]], t_coord = TCoord[Agg[a][i]], s_coord = SCoord[Agg[a][i]];
			x[Coords[x_coord][t_coord]][s_coord] += test_vectors_copy[k][Coords[x_coord][t_coord]][s_coord] * v[k][a];		
			
		}
	}
	return x;
}

//x_i = P^H_ij v_j. dim(P^H) =  Ntest Na x 2 Ntot, Nagg = block_x * block_t
//dim(v) = 2 NTot, dim(x) = Ntest Nagg
spinor AMG::Pt_v(const spinor& v) {
	//Restriction operator times a spinor
	using namespace AMGV;
	spinor x(Ntest, c_vector(Nagg, 0));
	for (int i = 0; i < Ntest*Nagg; i++) {
		int k = i / Nagg; //number of test vector
		int a = i % Nagg; //number of aggregate
		for (int j = 0; j < Agg[a].size(); j++) {
			int x_coord = XCoord[Agg[a][j]], t_coord = TCoord[Agg[a][j]], s_coord = SCoord[Agg[a][j]];
			x[k][a] += std::conj(test_vectors_copy[k][Coords[x_coord][t_coord]][s_coord]) * v[Coords[x_coord][t_coord]][s_coord];
		}
	}
	return x;
}


//Dc = P^H D P, dim(Dc) = Ntest Nagg x Ntest Nagg, dim(v) = Ntest Nagg, 
spinor AMG::Pt_D_P(const spinor& v){
	return Pt_v(D_phi(GConf.Conf,P_v(v),m0));
}


//x = D^-1 phi
spinor AMG::TwoGrid(const int& max_iter, const double& tol, const spinor& x0, 
	const spinor& phi, const bool& print_message) {
	//nu1 --> pre-smoothing steps
	//nu2 --> post-smoothing steps
	//x0 --> initial guess
	//phi --> right hand side
	//max_iter --> maximum number of iterations
	//tol --> tolerance
	spinor x = x0;
	spinor r;
	double err;
	int k = 0;
	double norm = sqrt(std::real(dot(phi,phi)));

	while(k < max_iter){
		//Pre-smoothing
		if (nu1>0){
			//x = gmres(LV::Ntot,2,GConf.Conf, phi, x, m0, AMGV::gmres_restarts_smoother, nu1, 1e-10, false);
			//SAP(GConf.Conf, phi, x, m0, nu1);
			SAP_parallel(GConf.Conf, phi, x, m0, nu1,SAPV::sap_blocks_per_proc); 
		} 

		double startT, endT;
		startT = MPI_Wtime();
		//*************Coarse grid correction*************//
		//x = x + P*Dc^-1 * P^H * (phi-D*x);  Coarse grid correction
		spinor Pt_r = Pt_v(phi - D_phi(GConf.Conf,x,m0)); //P^H (phi - D x)
		//-- Coarse grid solver --//
		//x = x + P_v(bi_cgstab_Dc(GConf.Conf, Pt_r, Pt_r, m0,AMGV::bi_cgstab_Dc_iterations,AMGV::bi_cgstab_Dc_iterations_tol,false)); 
		x = x + P_v(gmres(AMGV::Ntest,AMGV::Nagg,GConf.Conf, Pt_r, Pt_r, m0,
			AMGV::gmres_restart_length_coarse_level,AMGV::gmres_restarts_coarse_level,AMGV::gmres_tol_coarse_level,false)); 
		//************************************************//
		endT = MPI_Wtime();
		coarse_time += endT - startT; //Add coarse grid time

		//Post-smoothing
		if (nu2>0){
			//x = gmres(LV::Ntot,2,GConf.Conf, phi, x, m0, AMGV::gmres_restarts_smoother, nu2, 1e-10, false);
			//SAP(GConf.Conf, phi, x, m0, nu2);
			double startT, endT;
			startT = MPI_Wtime();
			SAP_parallel(GConf.Conf, phi, x, m0, nu2, SAPV::sap_blocks_per_proc); 
			endT = MPI_Wtime();
			smooth_time += endT - startT; //Add post-smoothing time
		}
		r = phi - D_phi(GConf.Conf, x, m0);
		err = sqrt(std::real(dot(r,r)));

		if (err < tol*norm){
			if (print_message == true){
				std::cout << "Two-grid method converged in " << k+1 << " iterations" << " Error " << err << std::endl;
			}
		return x;
		}
		k++;
	}
	if (print_message == true){
		std::cout << "Two-grid did not converge in " << max_iter << " iterations" << " Error " << err << std::endl;
	}
	return x;
}

//Bi-cgstab for Dc^-1 phi = x
//Coarse grid solver 
spinor AMG::bi_cgstab_Dc(const c_matrix& U, const spinor& phi, const spinor& x0, const double& m0, const int& max_iter, const double& tol, const bool& print_message) {
    //Dc^-1 phi = x
	//Dc = P^T D P 
	int k = 0; //Iteration number
    double err;
	using namespace AMGV;
    spinor r(Ntest, c_vector(Nagg, 0));  //r[coordinate][spin] residual
    spinor r_tilde(Ntest, c_vector(Nagg, 0));  //r[coordinate][spin] residual
    spinor d(Ntest, c_vector(Nagg, 0)); //search direction
    spinor s(Ntest, c_vector(Nagg, 0));
    spinor t(Ntest, c_vector(Nagg, 0));
    spinor Ad(Ntest, c_vector(Nagg, 0)); //DD^dagger*d
    spinor x(Ntest, c_vector(Nagg, 0)); //solution
    c_double alpha, beta, rho_i, omega, rho_i_2;
    x = x0; //initial solution
	
    r = phi - Pt_D_P(x); //r = b - A*x 
    r_tilde = r;
	double norm_phi = sqrt(std::real(dot(phi, phi))); //norm of the right hand side
    while (k<max_iter) {
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

//Solves Dc or D psi = phi using GMRES
spinor AMG::gmres(const int& dim1, const int& dim2,const c_matrix& U, const spinor& phi, const spinor& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message) {
    //GMRES for D^-1 phi
    //phi --> right-hand side
    //x0 --> initial guess  
    //U --> configuration
    //restarts --> number of restarts
    //m --> number of iterations per cycle

    int k = 0; //Iteration number (restart cycle)
    double err;

    spinor r(dim1, c_vector(dim2, 0));  //r[coordinate][spin] residual
   
    //VmT[column vector index][vector arrange in matrix form]
    std::vector<spinor> VmT(m+1, spinor(dim1, c_vector(dim2, 0))); //V matrix transpose-->dimensions exchanged

    spinor Hm(m+1 , c_vector(m, 0)); //H matrix (Hessenberg matrix)
    c_vector gm(m + 1, 0); 

    //Elements of rotation matrix |sn[i]|^2 + |cn[i]|^2 = 1
    c_vector sn(m, 0);
    c_vector cn(m, 0);
    c_vector eta(m, 0);


    spinor w(dim1, c_vector(dim2, 0)); //D*d
    spinor x = x0; //initial solution
    c_double beta;

	r = (dim1 == LV::Ntot) ?  phi - D_phi(U,x,m0) : phi - Pt_D_P(x); //r = b - A*x
	
	
	double norm_phi = sqrt(std::real(dot(phi, phi))); //norm of the right hand side
    while (k < restarts) {
        beta = sqrt(std::real(dot(r, r))) + 0.0 * I_number;
        VmT[0] = 1.0 / beta * r;
        gm[0] = beta; //gm[0] = ||r||
        //-----Arnoldi process to build the Krylov basis and the Hessenberg matrix-----//
        for (int j = 0; j < m; j++) {
			w = (dim1 == LV::Ntot) ? D_phi(U,VmT[j],m0) : Pt_D_P(VmT[j]); //w = D v_j

            for (int i = 0; i <= j; i++) {
                Hm[i][j] = dot(w, VmT[i]); //  (v_i^dagger, w)       
                w = w -  Hm[i][j] * VmT[i];
            }
            
            Hm[j + 1][j] = sqrt(std::real(dot(w, w))); //H[j+1][j] = ||A v_j||
            if (std::real(Hm[j + 1][j]) > 0) {
                VmT[j + 1] = 1.0 / Hm[j + 1][j] * w;
            }
            //----Rotate the matrix----//
            rotation(cn, sn, Hm, j); //Defined in include/gmres.h

            //Rotate gm
            gm[j + 1] = -sn[j] * gm[j];
            gm[j] = std::conj(cn[j]) * gm[j];
        }        
        //Solve the upper triangular system//
	
		eta = solve_upper_triangular(Hm, gm,m);
 
        for (int i = 0; i < dim1 * dim2; i++) {
            int n = i / dim2; int mu = i % dim2;
            for (int j = 0; j < m; j++) {
                x[n][mu] = x[n][mu] + eta[j] * VmT[j][n][mu]; 
            }
        }
        //Compute the residual
		r = (dim1 == LV::Ntot) ?  phi - D_phi(U,x,m0) : phi - Pt_D_P(x);
        err = sqrt(std::real(dot(r, r)));

         if (err < tol* norm_phi) {
             if (print_message == true) {
                 std::cout << "GMRES converged in " << k + 1 << " iterations" << " Error " << err << std::endl;
             }
             return x;
         }

         k++;
    }
	
        //std::cout << "GMRES did not converge in " << restarts << " restarts for tol = " << tol << " Error " << err << std::endl;
	return x;
}


//These functions are for saving the matrices. Useful for testing.
void save_spinor(spinor& phi,char* Name){
    char NameData[500], Data_str[500];
	sprintf(NameData, Name);
	std::ofstream Datfile;
	Datfile.open(NameData);
	int size = phi.size();
	int size2 = phi[0].size();
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size2; j++) {
			sprintf(Data_str, "%-30d%-30d%-30.17g%-30.17g\n", i, j, std::real(phi[i][j]), std::imag(phi[i][j]));
			Datfile << Data_str;
		}
	}
	Datfile.close();
}
