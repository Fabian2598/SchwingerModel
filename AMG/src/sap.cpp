#include "sap.h"

void SchwarzBlocks(){
    using namespace SAPV; //SAP parameters namespace
    schwarz_blocks = true; //Schwarz blocks are initialized after calling this function
    int count, block;
    int x0, t0, x1, t1;
    for (int x = 0; x < sap_block_x; x++) {
        for (int t = 0; t < sap_block_t; t++) {
            x0 = x * sap_x_elements; t0 = t * sap_t_elements;
            x1 = (x + 1) * sap_x_elements; t1 = (t + 1) * sap_t_elements;
            block = x * sap_block_t + t;
            count = 0;  
            //Filling the block with the coordinates of the lattice points
            for(int x = x0; x < x1; x++) {
                for (int t = t0; t < t1; t++) {
                    SAP_Blocks[block][count++] = x * Nt+ t; 
                    //Each block also considers both spin components, 
                    //so we only reference the lattice coordinates here.
                }
            }
            if (count != sap_lattice_sites_per_block) {
                std::cout << "Block " << block << " has " << count << " lattice points" << std::endl;
                std::cout << "Expected " << sap_lattice_sites_per_block << std::endl;
                exit(1);
            }
            //Red-black decomposition for the blocks.
            if (sap_block_t % 2 == 0) {
                if  (x%2 ==0){
                    (block % 2 == 0) ? SAP_RedBlocks[block / 2] = block:SAP_BlackBlocks[block / 2] = block; 
                }
                else{
                    (block % 2 == 0) ? SAP_BlackBlocks[block / 2] = block:SAP_RedBlocks[block / 2] = block; 
                }
            } 
            else {
                (block % 2 == 0) ? SAP_RedBlocks[block / 2] = block:SAP_BlackBlocks[block / 2] = block;                
            }
        }
    }
}

//x = I_B^T v --> Restriction of the vector v to the block B
//dim(v) = 2 Ntot, dim(x) = 2 * sap_lattice_sites_per_block
spinor It_B_v(const spinor& v, const int& block){
    using namespace SAPV;
    //Schwarz blocks have to be initialized first before calling this function
    spinor x(sap_lattice_sites_per_block, c_vector(2, 0)); //Initialize x to zero
    /*
    if (!schwarz_blocks){
        std::cout << "Error: Schwarz blocks not initialized" << std::endl;
        exit(1);
    }
    //Checking dimensions
    if ( checkSize(v, Ntot, 2)== true){
        std::cout << "Error with vector dimensions in It_B_v" << std::endl;
        exit(1);
    } 
    */

    for (int j = 0; j < sap_lattice_sites_per_block; j++){
        //Writing result to x 
        x[j][0] = v[SAP_Blocks[block][j]][0];
        x[j][1] = v[SAP_Blocks[block][j]][1];
    }
    return x;
}

// x = I_B v --> Interpolation of the vector v to the original lattice
//dim(v) = 2 * sap_lattice_sites_per_block, dim(x) = 2 Ntot
spinor I_B_v(const spinor& v,const int& block){
    using namespace SAPV;
    //Schwarz blocks have to be initialized first before calling this function
    /*
    if (!schwarz_blocks){
        std::cout << "Error: Schwarz blocks not initialized" << std::endl;
        exit(1);
    }
    if ( checkSize(v, sap_lattice_sites_per_block, 2) == true){
        std::cout << "Error with dimensions in I_B_v" << std::endl;
        exit(1);
    } 
    */
    spinor x(Ntot, c_vector(2, 0)); //Initialize x to zero

    for (int j = 0; j < sap_lattice_sites_per_block; j++){
        x[SAP_Blocks[block][j]][0] += v[j][0];
        x[SAP_Blocks[block][j]][1] += v[j][1];
    }
    return x;
}

//Restriction of the Dirac operator to the block B
//D_B v = I_B^T D I_B v  
//dim(v) = 2 * sap_lattice_sites_per_block, dim(x) = 2 * sap_lattice_sites_per_block
spinor D_B(const c_matrix& U, const spinor& v, const double& m0,const int& block){
    using namespace SAPV;
    /*
    if (checkSize(v, sap_lattice_sites_per_block, 2)==true){
        std::cout << "Error with vector dimensions in D_B" << std::endl;
        exit(1);
    }
    */
    //return I_B_v(D_phi(U, I_B_v(v, block), m0), block); //D_B v = I_B^T D I_B v

    spinor Dv(sap_lattice_sites_per_block, c_vector(2, 0));
    spinor phi = I_B_v(v, block); //Interpolation of the vector to the original lattice
    for (int m = 0; m < sap_lattice_sites_per_block; m++) {
		//n = x * Nt + t
        int n = SAP_Blocks[block][m]; //n is the index of the lattice point in the block
		Dv[m][0] = (m0 + 2) * phi[n][0] - 0.5 * ( 
			U[n][0] * SignR[n][0] * (phi[RightPB[n][0]][0] - phi[RightPB[n][0]][1]) 
		+   U[n][1] * SignR[n][1] * (phi[RightPB[n][1]][0] + I_number * phi[RightPB[n][1]][1])
		+   std::conj(U[LeftPB[n][0]][0]) * SignL[n][0] * (phi[LeftPB[n][0]][0] + phi[LeftPB[n][0]][1])
		+   std::conj(U[LeftPB[n][1]][1]) * SignL[n][1] * (phi[LeftPB[n][1]][0] - I_number * phi[LeftPB[n][1]][1])
		);

		Dv[m][1] = (m0 + 2) * phi[n][1] - 0.5 * ( 
			U[n][0] * SignR[n][0] * (-phi[RightPB[n][0]][0] + phi[RightPB[n][0]][1]) 
		+   U[n][1] * SignR[n][1] * (-I_number*phi[RightPB[n][1]][0] + phi[RightPB[n][1]][1])
		+   std::conj(U[LeftPB[n][0]][0]) * SignL[n][0] * (phi[LeftPB[n][0]][0] + phi[LeftPB[n][0]][1])
		+   std::conj(U[LeftPB[n][1]][1]) * SignL[n][1] * (I_number*phi[LeftPB[n][1]][0] + phi[LeftPB[n][1]][1])
		);
			
	}

    return Dv;
    
}

//Solves D_B x = phi using GMRES, where D_B is the Dirac matrix restricted to the Schwarz block B 
spinor gmres_D_B(const c_matrix& U, const spinor& phi, const spinor& x0, 
    const double& m0, const int& m, const int& restarts, const double& tol, const int& block,
    const bool& print_message) {
    using namespace SAPV;
    
    if (print_message == true){
        std::cout << "------------------------------------------" << std::endl;
        std::cout << "|GMRES for D_B (block " << block << ") " << std::endl;
        std::cout << "|Restart length " << m << ". Restarts = " << restarts << std::endl;
    }
    //if (checkSize(phi, sap_lattice_sites_per_block, 2) == true){
    //    std::cout << "Error with vector dimensions in gmres_D_B" << std::endl;
    //    exit(1);
    //} 
    if (m> sap_variables_per_block) {
        std::cout << "Error: restart length > sap_variables_per_block" << std::endl;
        exit(1);
    }
    
    int k = 0; //Iteration number (restart cycle)
    double err; //Error = ||r||_2

    spinor r(sap_lattice_sites_per_block, c_vector(2, 0));  //r[coordinate][spin] residual
    
    //VmT[column vector index][vector arrange in matrix form]
    std::vector<spinor> VmT(m+1, spinor(sap_lattice_sites_per_block, c_vector(2, 0))); //V matrix transpose --> dimensions exchanged

    c_matrix Hm(m+1 , c_vector(m, 0)); //H matrix (Hessenberg matrix)
    c_vector gm(m + 1, 0); 

    //Elements of rotation matrix |sn[i]|^2 + |cn[i]|^2 = 1
    c_vector sn(m, 0);
    c_vector cn(m, 0);
    c_vector eta(m, 0);


    spinor w(sap_lattice_sites_per_block, c_vector(2, 0)); //D*d
    spinor x = x0; //initial solution
    c_double beta; //not 1/g^2

    //r = phi - D_B(U, phi, m0, block); //r = b - A*x
    spinor temp =  D_B(U, phi, m0, block);
    for(int i = 0; i < sap_lattice_sites_per_block; i++) {
        r[i][0] = phi[i][0] - temp[i][0];
        r[i][1] = phi[i][1] - temp[i][1];
    }

	
	double norm_phi = sqrt(std::real(dot(phi, phi))); //norm of the right hand side
    err = sqrt(std::real(dot(r, r))); //Initial error
    if (print_message == true){std::cout << "||phi|| * tol = " << norm_phi * tol << std::endl;}
    while (k < restarts) {
        beta = err + 0.0 * I_number;
        VmT[0] = 1.0 / beta * r;
        gm[0] = beta; //gm[0] = ||r||
        //-----Arnoldi process to build the Krylov basis and the Hessenberg matrix-----//
        for (int j = 0; j < m; j++) {
            w = D_B(U, VmT[j], m0, block); //w = D v_j  

            for (int i = 0; i <= j; i++) {
                Hm[i][j] = dot(w, VmT[i]); //  (v_i^dagger, w)

                //w = w -  Hm[i][j] * VmT[i];
				for(int n=0; n<sap_lattice_sites_per_block; n++){
					for(int l=0; l<2; l++){
						w[n][l] -= Hm[i][j] * VmT[i][n][l];
					}
				}
            }
            Hm[j + 1][j] = sqrt(std::real(dot(w, w))); //H[j+1][j] = ||A v_j||
            if (std::real(Hm[j + 1][j]) > 0) {
                VmT[j + 1] = 1.0 / Hm[j + 1][j] * w;
            }
            //----Rotate the matrix----//
            rotation(cn, sn, Hm, j);

            //Rotate gm
            gm[j + 1] = -sn[j] * gm[j];
            gm[j] = std::conj(cn[j]) * gm[j];
        }        
        //Solve the upper triangular system//
	
		eta = solve_upper_triangular(Hm, gm,m);
 
        for (int i = 0; i < sap_variables_per_block; i++) {
            int n = i / 2; int mu = i % 2; //Splitting the spin components
            for (int j = 0; j < m; j++) {
                x[n][mu] = x[n][mu] + eta[j] * VmT[j][n][mu]; 
            }
        }
        //Compute the residual
        //r = phi - D_B(U, x, m0, block); //D_B x; 
        temp =  D_B(U, x, m0, block);
        for(int i = 0; i < sap_lattice_sites_per_block; i++) {
            r[i][0] = phi[i][0] - temp[i][0];
            r[i][1] = phi[i][1] - temp[i][1];
        }

        err = sqrt(std::real(dot(r, r)));

         if (err < tol * norm_phi) {
             if (print_message == true) {
                 std::cout << "|GMRES for D_B (block " << block << ") converged in " << k + 1 << " restarts." << " Error " << err << std::endl;
                 std::cout << "------------------------------------------" << std::endl;
            }
             
             return x;
             
         }
         k++;
    }
    //if (print_message == true) {
      //  std::cout << "|GMRES for D_B (block " << block << ") did not converge in " << restarts << " restarts." << " Error " << err << std::endl;
    //
    //std::cout << "------------------------------------------" << std::endl;
    return x;
}

//A_B = I_B * D_B^-1 * I_B^T v --> Extrapolation of D_B^-1 to the original lattice.
//dim(v) = 2 * Ntot, dim(x) = 2 Ntot
//v: input, x: output
spinor I_D_B_1_It(const c_matrix& U, const spinor& v, const double& m0,const int& block){
    using namespace SAPV;
    bool print_message = false; //good for testing GMRES   
    /*
    if (checkSize(v, Ntot, 2) == true){
        std::cout << "Error with vector dimensions in I_D_B_1_It" << std::endl;
        exit(1);
    } 
    */
    spinor x(Ntot, c_vector(2, 0)); //Initialize x to zero
    spinor temp = It_B_v(v,block); //temp = I_B^T v

    return I_B_v(gmres_D_B(U, temp, temp, m0, 
       sap_gmres_restart_length, sap_gmres_restarts, sap_gmres_tolerance, 
        block, print_message),block); //x = I_B D_B^-1 I_B^T v 
}

/*
phi has the original dimension, the output has block dimension
*/
spinor local_D(const c_matrix& U, const spinor& v, const double& m0,const int& block) {
    using namespace SAPV;

    spinor Dv(sap_lattice_sites_per_block, c_vector(2, 0));
    spinor phi = I_B_v(v, block); //Interpolation of the vector to the original lattice
    for (int m = 0; m < sap_lattice_sites_per_block; m++) {
		//n = x * Nt + t
        int n = SAP_Blocks[block][m]; //n is the index of the lattice point in the block
		Dv[m][0] = (m0 + 2) * phi[n][0] - 0.5 * ( 
			U[n][0] * SignR[n][0] * (phi[RightPB[n][0]][0] - phi[RightPB[n][0]][1]) 
		+   U[n][1] * SignR[n][1] * (phi[RightPB[n][1]][0] + I_number * phi[RightPB[n][1]][1])
		+   std::conj(U[LeftPB[n][0]][0]) * SignL[n][0] * (phi[LeftPB[n][0]][0] + phi[LeftPB[n][0]][1])
		+   std::conj(U[LeftPB[n][1]][1]) * SignL[n][1] * (phi[LeftPB[n][1]][0] - I_number * phi[LeftPB[n][1]][1])
		);

		Dv[m][1] = (m0 + 2) * phi[n][1] - 0.5 * ( 
			U[n][0] * SignR[n][0] * (-phi[RightPB[n][0]][0] + phi[RightPB[n][0]][1]) 
		+   U[n][1] * SignR[n][1] * (-I_number*phi[RightPB[n][1]][0] + phi[RightPB[n][1]][1])
		+   std::conj(U[LeftPB[n][0]][0]) * SignL[n][0] * (phi[LeftPB[n][0]][0] + phi[LeftPB[n][0]][1])
		+   std::conj(U[LeftPB[n][1]][1]) * SignL[n][1] * (I_number*phi[LeftPB[n][1]][0] + phi[LeftPB[n][1]][1])
		);
			
	}

    return Dv; 
}

//Sequential version of the SAP method (used for testing)
spinor SAP(const c_matrix& U, const spinor& v, const double& m0,const int& nu){
    /*
    Solves D x = v using the SAP method
    */
  
    using namespace SAPV;
    spinor x(Ntot, c_vector(2, 0)); //Initialize x to zero
    //if (checkSize(v, Ntot, 2) == true){
    //    std::cout << "Error with vector dimensions in I_KD" << std::endl;
    //    exit(1);
    //}  
    double err;
    double v_norm = sqrt(std::real(dot(v, v))); //norm of the right hand side

    spinor temp(Ntot, c_vector(2, 0)); 
    spinor r(Ntot, c_vector(2, 0)); //residual
    r = v - D_phi(U, x, m0); //r = v - D x
    for (int i = 0; i< nu; i++){
        for (auto block : SAP_RedBlocks){
            temp = I_D_B_1_It(U,r,m0,block);
            //x = x + temp; //x = x + D_B^-1 r
            for(int n = 0; n < Ntot; n++) {
                x[n][0] += temp[n][0];
                x[n][1] += temp[n][1];
            }
        }
        
        r = v - D_phi(U, x, m0); //r = v - D x
        for (auto block : SAP_BlackBlocks){
            temp = I_D_B_1_It(U,r,m0,block);
            //x = x + temp; //x = x + D_B^-1 r
            for(int n = 0; n < Ntot; n++) {
                x[n][0] += temp[n][0];
                x[n][1] += temp[n][1];
            }
           
        }
        r = v - D_phi(U, x, m0); //r = v - D x

        err = sqrt(std::real(dot(r, r))); 
        if (err < sap_tolerance * v_norm) {
            //std::cout << "SAP converged in " << i << " iterations, error: " << err << std::endl;
            return x;
        }
    }
    //std::cout << "SAP did not converge in " << nu << " iterations, error: " << err << std::endl;
    return x; 
}

