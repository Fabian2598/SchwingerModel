#include "amg.h"

//Aggregates A_j_0 = L_j x {0}, A_j_1 = L_j x {1} (Volume times Spin component)
void Aggregates() {
	using namespace LV; //Lattice parameters namespace
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


void PrintAggregates() {
	for (int i = 0; i < AMGV::Nagg; i++) {
		std::cout << "-------Aggregate-----" << i << std::endl;
		for (int j = 0; j < LV::x_elements * LV::t_elements; j++) {
			std::cout << Agg[i][j] << " ";
		}
		std::cout << std::endl;
	}
}


void normalize(spinor& v){
	c_double norm = sqrt(std::real(dot(v,v))) + 0.0*I_number; 
	scal(1.0/norm, v, v); //v = v / norm
}

spinor canonical_vector(const int& i, const int& N1, const int& N2) {
	spinor e_i(N1, c_vector (N2,0.0));
	int j = i / N2;
	int k = i % N2;
	e_i[j][k] = 1.0;
	return e_i;
}

void AMG::orthonormalize(){
	//Local orthonormalization of the test vectors
	
	//Each test vector is chopped into the Nagg aggregates, which yields Ntest*Nagg columns for the interpolator.
	//Each column is orthonormalized with respect to the others that belong to the same aggregate.
	//This follows the steps from Section 3.1 of A. Frommer et al "Adaptive Aggregation-Based Domain Decomposition 
	//Multigrid for the Lattice Wilson-Dirac Operator", SIAM, 36 (2014).
	

	using namespace AMGV; //AMG parameters
	int x, t, s, n;
	//Orthonormalization by applying Gram-Schmidt
	c_double proj;
	c_double norm;
	for (int a = 0; a < Nagg; a++) {
		for (int nt = 0; nt < Ntest; nt++) {
			for(int ntt=0; ntt < nt; ntt++){
				proj = 0;//dot(v_chopped[nt*Nagg+i], v_chopped[j*Nagg+i]);
				//Addition over the elements of the aggregate a
				for (int j = 0; j < LV::x_elements * LV::t_elements; j++) {
        			x = XCoord[Agg[a][j]], t = TCoord[Agg[a][j]], s = SCoord[Agg[a][j]];
					n = Coords[x][t];
					proj += interpolator_columns[nt][n][s] * std::conj(interpolator_columns[ntt][n][s]);
    			}

				for (int j = 0; j < LV::x_elements * LV::t_elements; j++) {
        			x = XCoord[Agg[a][j]], t = TCoord[Agg[a][j]], s = SCoord[Agg[a][j]];
					n = Coords[x][t];
					interpolator_columns[nt][n][s] -= proj * interpolator_columns[ntt][n][s];
    			}
				
			}

			//normalize
			norm = 0.0;
			for (int j = 0; j < LV::x_elements * LV::t_elements; j++) {
				x = XCoord[Agg[a][j]], t = TCoord[Agg[a][j]], s = SCoord[Agg[a][j]];
				n = Coords[x][t];
				norm += interpolator_columns[nt][n][s] * std::conj(interpolator_columns[nt][n][s]);
			}
			norm = sqrt(std::real(norm)) + 0.0*c_double(0,1); 
			for (int j = 0; j < LV::x_elements * LV::t_elements; j++) {
				x = XCoord[Agg[a][j]], t = TCoord[Agg[a][j]], s = SCoord[Agg[a][j]];
				n = Coords[x][t];
				interpolator_columns[nt][n][s] /= norm;
			}
		}
	}
	
}


void AMG::checkOrthogonality(){
	//Check orthogonality of the test vectors
	//aggregate 
	for(int block = 0; block < LV::Nblocks; block++){
	for (int alf = 0; alf < 2; alf++) {
		//checking orthogonality 
		for (int i = 0; i < AMGV::Ntest; i++) {
		for (int j = 0; j < AMGV::Ntest; j++) {
			c_double dot_product = 0.0;
			for (int n: LatticeBlocks[block]) {
				dot_product += std::conj(interpolator_columns[i][n][alf]) * interpolator_columns[j][n][alf];
			}
			if (std::abs(dot_product) > 1e-8 && i!=j) {
				std::cout << "Block " << block << " spin " << alf << std::endl;
				std::cout << "Test vectors " << i << " and " << j << " are not orthogonal: " << dot_product << std::endl;
				exit(1);
			}
			else if(std::abs(dot_product-1.0) > 1e-8 && i==j){
				std::cout << "Test vector " << i << " not orthonormalized " << dot_product << std::endl;
				exit(1);
			}

		}
		}
	}
	}
	//std::cout << "Test vectors are orthonormalized " << std::endl;

}

void AMG::setUpPhase(const int& Nit) {
	//Call MPI for SAP parallelization
	int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	spinor rhs(LV::Ntot, c_vector(2, 0));
	spinor x0(LV::Ntot, c_vector(2, 0));
	using namespace AMGV; //AMG parameters namespace


  	static std::mt19937 randomInt(50); //Same seed for all the MPI copies
	std::uniform_real_distribution<double> distribution(-1.0, 1.0); //mu, standard deviation

	//Test vectors random initialization


	for (int i = 0; i < Ntest; i++) {
		for (int j = 0; j < LV::Ntot; j++) {
			for (int k = 0; k < 2; k++) {
				interpolator_columns[i][j][k] = distribution(randomInt) + I_number * distribution(randomInt);
			}
		}
	}
	
	for (int i = 0; i < Ntest; i++) {
		// D^{-1} test_vector with SAP
		double startT, endT;
		startT = MPI_Wtime();
		sap.SAP(rhs,interpolator_columns[i],AMGV::SAP_test_vectors_iterations, SAPV::sap_blocks_per_proc,false);
		endT = MPI_Wtime();
		SAP_time += endT - startT; 
			
	}
	
	orthonormalize(); 
	checkOrthogonality();
	initializeCoarseLinks();
	
	spinor D_v(LV::Ntot, c_vector(2, 0));
	int two_grid_iter = 1; 
	double tolerance = 1e-2;
	//Improving the interpolator quality by iterating over the two-grid method defined by the current test vectors
	if (rank == 0){std::cout << "Improving interpolator" << std::endl;}
	for (int n = 0; n < Nit; n++) {
		if (rank == 0){std::cout << "****** Bootstrap iteration " << n << " ******" << std::endl;}
		for (int i = 0; i < Ntest; i++) {
			
			D_phi(GConf.Conf,interpolator_columns[i],D_v,m0); //D times the i-th test vector
			axpy(interpolator_columns[i], D_v, -1.0, rhs); //rhs = v_i - D v_i

			TwoGrid(two_grid_iter,tolerance, x0, interpolator_columns[i],test_vectors[i], false,false); 
			
			for(int n = 0; n < LV::Ntot; n++){
				for(int k = 0; k < 2; k++){
					test_vectors[i][n][k] += interpolator_columns[i][n][k]; 
				}
			}
			
		}
		//"Assemble" the new interpolator
		interpolator_columns = test_vectors; 
		orthonormalize();
		checkOrthogonality();
		initializeCoarseLinks();

	}

	if (rank == 0){std::cout << "Set-up phase finished" << std::endl;}
		
}

/*
	Prolongation operator times a spinor x = P v
	x_i = P_ij v_j. dim(P) = 2 Ntot x Ntest Nagg, 
	dim(v) = Ntest Nagg, dim(x) = 2 Ntot
*/
void AMG::P_v(const spinor& v,spinor& out) {
	//Loop over columns
	for(int n = 0; n < LV::Ntot; n++){
		for(int alf = 0; alf < 2; alf++){
			out[n][alf] = 0.0; //Initialize the output spinor
		}
	}
	using namespace AMGV; //AMG parameters namespace
	int x_coord, t_coord, s_coord; //Coordinates of the lattice point
	int k, a;
	int i, j; //Loop indices
	for (j = 0; j < Ntest * Nagg; j++) {
		k = j % Ntest; //Number of test vector
		a = j / Ntest; //Number of aggregate		
		for (i = 0; i < Agg[a].size(); i++) {
			x_coord = XCoord[Agg[a][i]], t_coord = TCoord[Agg[a][i]], s_coord = SCoord[Agg[a][i]];
			out[Coords[x_coord][t_coord]][s_coord] += interpolator_columns[k][Coords[x_coord][t_coord]][s_coord] * v[k][a];			
		}
	}
}

/*
	Restriction operator times a spinor on the coarse grid, x = P^H v
	x_i = P^H_ij v_j. dim(P^H) =  Ntest Nagg x 2 Ntot, Nagg = block_x * block_t
	dim(v) = 2 Ntot, dim(x) = Ntest Nagg
*/

void AMG::Pt_v(const spinor& v,spinor& out) {
	//Restriction operator times a spinor
	using namespace AMGV;
	for(int n = 0; n < Ntest; n++){
		for(int alf = 0; alf < Nagg; alf++){
			out[n][alf] = 0.0; //Initialize the output spinor
		}
	}
	int k, a;
	int x_coord, t_coord, s_coord; //Coordinates of the lattice point
	int i, j; //Loop indices
	for (i = 0; i < Ntest*Nagg; i++) {
		k = i % Ntest; //Number of test vector
		a = i / Ntest;
		for (j = 0; j < Agg[a].size(); j++) {
			x_coord = XCoord[Agg[a][j]], t_coord = TCoord[Agg[a][j]], s_coord = SCoord[Agg[a][j]];
			out[k][a] += std::conj(interpolator_columns[k][Coords[x_coord][t_coord]][s_coord]) * v[Coords[x_coord][t_coord]][s_coord];
		}
	}

}

/*
	Intialize the coarse gauge links for Dc
*/
void AMG::initializeCoarseLinks(){
	c_double P[2][2][2], M[2][2][2]; 
	
	P[0][0][0] = 1.0; P[0][0][1] = 1.0;
	P[0][1][0] = 1.0; P[0][1][1] = 1.0; 

	P[1][0][0] = 1.0; P[1][0][1] = -I_number;
	P[1][1][0] = I_number; P[1][1][1] = 1.0; 

	M[0][0][0] = 1.0; M[0][0][1] = -1.0;
	M[0][1][0] = -1.0; M[0][1][1] = 1.0; 

	M[1][0][0] = 1.0; M[1][0][1] = I_number;
	M[1][1][0] = -I_number; M[1][1][1] = 1.0; 

	c_matrix &U = GConf.Conf;
	std::vector<spinor> &w = interpolator_columns;
	c_double Lm, Lp, R;
	int block_r, block_l; //Block indices for the right and left periodic boundary
	for(int x=0; x<LV::Nblocks; x++){
	for(int alf=0; alf<2;alf++){
	for(int bet=0; bet<2;bet++){
	for(int p = 0; p<AMGV::Ntest; p++){
	for(int s = 0; s<AMGV::Ntest; s++){

		A_coeff[x][alf][bet][p][s] = 0;
		B_coeff[x][alf][bet][p][s][0] = 0; B_coeff[x][alf][bet][p][s][1] = 0;
		C_coeff[x][alf][bet][p][s][0] = 0; C_coeff[x][alf][bet][p][s][1] = 0;
		for(int n : LatticeBlocks[x]){
		for(int mu : {0,1}){
			getLatticeBlock(RightPB[n][mu], block_r); //Get the block index for the right periodic boundary
			getLatticeBlock(LeftPB[n][mu], block_l); //Get the block index for the right periodic boundary

			Lm = 0.5 * M[mu][alf][bet] * std::conj(w[p][n][alf]) * U[n][mu];
			Lp = 0.5 * P[mu][alf][bet] * std::conj(w[p][n][alf]) * std::conj(U[LeftPB[n][mu]][mu]);	
			//if n+\hat{mu} in Block(x)
			if (block_r == x)
				//[A(x)]^{alf,bet}_{p,s} --> A_coeff[x][alf][bet][p][s] 
				A_coeff[x][alf][bet][p][s] += Lm * w[s][RightPB[n][mu]][bet] * SignR[n][mu];
			//if n+\hat{mu} in Block(x+hat{mu})
			else if (block_r == RightPB_blocks[x][mu])
				//[B_mu(x)]^{alf,bet}_{p,s} --> B_coeff[x][alf][bet][p][s][mu]
				B_coeff[x][alf][bet][p][s][mu] += Lm * w[s][RightPB[n][mu]][bet] * SignR[n][mu];

			//if n-\hat{mu} in Block(x)	
			if (block_l == x)
				A_coeff[x][alf][bet][p][s] += Lp * w[s][LeftPB[n][mu]][bet] * SignL[n][mu];
			//if n-\hat{mu} in Block(x-hat{mu})
			else if (block_l == LeftPB_blocks[x][mu])
				//[C_mu(x)]^{alf,bet}_{p,s} --> C_coeff[x][alf][bet][p][s][mu]
				C_coeff[x][alf][bet][p][s][mu] += Lp * w[s][LeftPB[n][mu]][bet] * SignL[n][mu];
			
			
		}//mu 
		}//n 


	//---------Close loops---------//
	} //s
	} //p
	} //bet
	} //alf
	} //x 

}


/*
	Coarse grid matrix operator Dc = P^H D P times a spinor v
	dim(Dc) = Ntest Nagg x Ntest Nagg, dim(v) = Ntest Nagg,
	Dc_{xy} = A(x) d_{xy} + \sum_{mu} B_mu(x) d_{x+mu,y} + \sum_{mu} C_mu(x)  d_{x-mu,y} 
	With x, y lattice blocks.
*/
void AMG::Pt_D_P(const spinor& v,spinor& out){
	for(int x = 0; x<LV::Nblocks; x++){
	for(int alf = 0; alf<2; alf++){
	for(int p = 0; p<AMGV::Ntest; p++){
		out[p][2*x+alf] = (m0+2)*v[p][2*x+alf]; //Mass term
	for(int bet = 0; bet<2; bet++){
	for(int s = 0; s<AMGV::Ntest; s++){

		//out[test_vec][aggregate]
		out[p][2*x+alf] -= A_coeff[x][alf][bet][p][s] * v[s][2*x+bet];

		for(int mu:{0,1}){
			out[p][2*x+alf] -=  (B_coeff[x][alf][bet][p][s][mu] * v[s][2*RightPB_blocks[x][mu]+bet]
							    +C_coeff[x][alf][bet][p][s][mu] * v[s][2*LeftPB_blocks[x][mu]+bet]);
		}
		

	}
	}
	}
	}
	}
}

/*
	Two-grid method for solving the linear system D x = phi
*/
int AMG::TwoGrid(const int& max_iter, const double& tol, const spinor& x0, 
	const spinor& phi, spinor& x,const bool& save_res, const bool& print_message) {

	x = x0; //Solution spinor
	spinor r(LV::Ntot,c_vector(2,0)); //Residual 
	spinor temp(LV::Ntot,c_vector(2,0));
	spinor Dphi(LV::Ntot,c_vector(2,0));
	double err; //Error = sqrt(dot(r,r))
	int k = 0; //Iteration number
	double norm_phi = sqrt(std::real(dot(phi,phi)));
	std::vector<double> residuals;
	//The convergence criterion is ||r|| < ||phi|| * tol

	while(k < max_iter){
		//Pre-smoothing
		if (nu1>0){
			sap.SAP(phi,x,nu1, SAPV::sap_blocks_per_proc,false); 
		} 

		D_phi(GConf.Conf, x, Dphi,m0); //D x
		axpy(phi,Dphi,-1.0,r); //r = phi - D x
		//----------Coarse grid correction-------------//
		double startT, endT;
		startT = MPI_Wtime();
		//x = x + P*Dc^-1 * P^H * (phi-D*x)  	
		
		//spinor Pt_r = Pt_v(temp); //Pt_r = P^H (phi - D x)
		spinor Pt_r(AMGV::Ntest, c_vector(AMGV::Nagg, 0));
		Pt_v(r,Pt_r);
		
		//Using GMRES for the coarse grid solver 
		spinor gmresResult(AMGV::Ntest, c_vector(AMGV::Nagg, 0));
		gmres_c_level.fgmres(Pt_r,Pt_r,gmresResult,false,false);
		P_v(gmresResult,temp);

		for(int n = 0; n<LV::Ntot; n++){
			for(int alf=0; alf<2; alf++){
				x[n][alf] += temp[n][alf]; 
			}
		}
	
		endT = MPI_Wtime();
		coarse_time += endT - startT; //Measuring time spent for solving the coarse level 
		//------------------------------------------------------//
		
		//Post-smoothing
		if (nu2>0){
			//Measure time spent smoothing
			double startT, endT;
			startT = MPI_Wtime();
			MPI_Barrier(MPI_COMM_WORLD);
			sap.SAP(phi,x,nu2, SAPV::sap_blocks_per_proc,false);
			endT = MPI_Wtime();
			smooth_time += endT - startT; //Add post-smoothing time
			SAP_time += endT - startT; //Add post-smoothing time
		}
		//r = phi - D_phi(GConf.Conf, x, m0);		
		D_phi(GConf.Conf, x, Dphi,m0); //D x
		axpy(phi,Dphi,-1.0,r); //r = phi - D x
		
		
		err = sqrt(std::real(dot(r,r)));

		residuals.push_back(err);

		if (err < tol*norm_phi){
			if (print_message == true){
				std::cout << "Two-grid method converged in " << k+1 << " cycles" << " Error " << err << std::endl;
			}
			if (save_res == true){
                std::ostringstream NameData;
                NameData << "TwoGrid_residual_" << LV::Nx << "x" << LV::Nt << ".txt";
                save_vec(residuals,NameData.str());
            }
		return 1;
		}
		k++;
	}
	if (print_message == true)
		std::cout << "Two-grid did not converge in " << max_iter << " cycles" << " Error " << err << std::endl;
	if (save_res == true){
        std::ostringstream NameData;
        NameData << "TwoGrid_residual_" << LV::Nx << "x" << LV::Nt << ".txt";
        save_vec(residuals,NameData.str());
    }
	return 0;
}


void AMG::testSetUp(){
	//Checking orthogonality

    checkOrthogonality();

	using namespace AMGV;
     // Testing that P^dag D P = D_c 
    spinor in(Ntest,c_vector(Nagg,1)); //in
    spinor temp(LV::Ntot,c_vector(2,0));
    spinor Dphi(LV::Ntot,c_vector(2,0));
    spinor out(Ntest,c_vector(Nagg,0)); //out
    spinor out_v2(Ntest,c_vector(Nagg,0)); //D_c
    //P^H D P
    P_v(in,temp);
    D_phi(GConf.Conf,temp,Dphi,m0);
    Pt_v(Dphi,out);

    Pt_D_P(in,out_v2);
    for(int x = 0; x<Ntest; x++){
        for(int dof = 0; dof<Nagg; dof++){
            if (std::abs(out[x][dof]-out_v2[x][dof]) > 1e-8 ){
            std::cout << "[" << x << "][" << dof << "] " << " different" << std::endl; 
            std::cout << out[x][dof] << "   /=    " << out_v2[x][dof] << std::endl;
            return;
            }
        }
    }
    std::cout << "P^dag D P coincides with Dc" << std::endl;
    std::cout << out[0][0] << "   =    " << out_v2[0][0] << std::endl;
    
     
}