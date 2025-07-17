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
void It_B_v(const spinor& v, spinor& x, const int& block){
    using namespace SAPV;
    set_zeros(x,sap_lattice_sites_per_block,2); //Initialize the output vector to zero
    for (int j = 0; j < sap_lattice_sites_per_block; j++){
        //Writing result to x 
        x[j][0] = v[SAP_Blocks[block][j]][0];
        x[j][1] = v[SAP_Blocks[block][j]][1];
    }
}

// x = I_B v --> Interpolation of the vector v to the original lattice
//dim(v) = 2 * sap_lattice_sites_per_block, dim(x) = 2 Ntot
void I_B_v(const spinor& v, spinor& x,const int& block){
    using namespace SAPV;
    set_zeros(x,Ntot,2); //Initialize x to zero
    for (int j = 0; j < sap_lattice_sites_per_block; j++){
        x[SAP_Blocks[block][j]][0] += v[j][0];
        x[SAP_Blocks[block][j]][1] += v[j][1];
    }

}

//Restriction of the Dirac operator to the block B
//D_B v = I_B^T D I_B v  
//dim(v) = 2 * sap_lattice_sites_per_block, dim(x) = 2 * sap_lattice_sites_per_block
void D_B(const c_matrix& U, const spinor& v, spinor& x, const double& m0,const int& block){
    using namespace SAPV;
    int RightPB_0, blockRPB_0; //Right periodic boundary in the 0-direction
    int RightPB_1, blockRPB_1; //Right periodic boundary in the 1-direction
    int LeftPB_0, blockLPB_0; //Left periodic boundary in the 0-direction
    int LeftPB_1, blockLPB_1; //Left periodic boundary in the 1-direction
    for (int m = 0; m < sap_lattice_sites_per_block; m++) {
		//n = x * Nt + t
        int n = SAP_Blocks[block][m]; //n is the index of the lattice point in the block
        
        //Get m and block for the neighbors
        getMandBlock(RightPB[n][0], RightPB_0, blockRPB_0); //I could make an array with this info, but it would have to be for each block
        getMandBlock(RightPB[n][1], RightPB_1, blockRPB_1); 
        getMandBlock(LeftPB[n][0], LeftPB_0, blockLPB_0); 
        getMandBlock(LeftPB[n][1], LeftPB_1, blockLPB_1); 
        
        c_vector phi_RPB_0 = (blockRPB_0 == block) ? v[RightPB_0]: c_vector(2, 0);
        c_vector phi_RPB_1 = (blockRPB_1 == block) ? v[RightPB_1]: c_vector(2, 0);
        c_vector phi_LPB_0 = (blockLPB_0 == block) ? v[LeftPB_0]: c_vector(2, 0);
        c_vector phi_LPB_1 = (blockLPB_1 == block) ? v[LeftPB_1]: c_vector(2, 0);
 

		x[m][0] = (m0 + 2) * v[m][0] - 0.5 * ( 
			U[n][0] * SignR[n][0] * (phi_RPB_0[0] - phi_RPB_0[1]) 
		+   U[n][1] * SignR[n][1] * (phi_RPB_1[0] + I_number * phi_RPB_1[1])
		+   std::conj(U[LeftPB[n][0]][0]) * SignL[n][0] * (phi_LPB_0[0] + phi_LPB_0[1])
		+   std::conj(U[LeftPB[n][1]][1]) * SignL[n][1] * (phi_LPB_1[0] - I_number * phi_LPB_1[1])
		);

		x[m][1] = (m0 + 2) * v[m][1] - 0.5 * ( 
			U[n][0] * SignR[n][0] * (-phi_RPB_0[0] + phi_RPB_0[1]) 
		+   U[n][1] * SignR[n][1] * (-I_number*phi_RPB_1[0] + phi_RPB_1[1])
		+   std::conj(U[LeftPB[n][0]][0]) * SignL[n][0] * (phi_LPB_0[0] + phi_LPB_0[1])
		+   std::conj(U[LeftPB[n][1]][1]) * SignL[n][1] * (I_number*phi_LPB_1[0] + phi_LPB_1[1])
		);
			
	}
    
    

}

//Solves D_B x = phi using GMRES, where D_B is the Dirac matrix restricted to the Schwarz block B 
int gmres_D_B(const c_matrix& U, const spinor& phi, const spinor& x0, spinor& x,
    const double& m0, const int& m, const int& restarts, const double& tol, const int& block,
    const bool& print_message) {
    using namespace SAPV;
   
    if (m> sap_variables_per_block) {
        std::cout << "Error: restart length > sap_variables_per_block" << std::endl;
        exit(1);
    }
    
    int k = 0; //Iteration number (restart cycle)

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
    x = x0; //initial solution
    c_double beta; //not 1/g^2

	spinor temp(sap_lattice_sites_per_block, c_vector(2, 0)); //I_B v
    D_B(U, phi,temp, m0, block);
    r = phi - temp; //r = b - A*x
	
	double norm_phi = sqrt(std::real(dot(phi, phi))); //norm of the right hand side
    double err = sqrt(std::real(dot(r, r))); //Initial error  ||r||_2
    if (print_message == true){std::cout << "||phi|| * tol = " << norm_phi * tol << std::endl;}
    while (k < restarts) {
        beta = err + 0.0 * I_number;
        VmT[0] = 1.0 / beta * r;
        gm[0] = beta; //gm[0] = ||r||
        //-----Arnoldi process to build the Krylov basis and the Hessenberg matrix-----//
        for (int j = 0; j < m; j++) {
            D_B(U, VmT[j], w,m0, block); //w = D v_j  

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
        D_B(U, x,temp, m0, block); //D_B x
        r = phi - temp; 
        err = sqrt(std::real(dot(r, r)));

         if (err < tol * norm_phi) {
             if (print_message == true) {
                 std::cout << "|GMRES for D_B (block " << block << ") converged in " << k + 1 << " restarts." << " Error " << err << std::endl;
                 std::cout << "------------------------------------------" << std::endl;
            }
             return 1;
             
         }
         k++;
    }
    //if (print_message == true) {
      //  std::cout << "|GMRES for D_B (block " << block << ") did not converge in " << restarts << " restarts." << " Error " << err << std::endl;
    //
    //std::cout << "------------------------------------------" << std::endl;
    return 0;
}

//A_B = I_B * D_B^-1 * I_B^T v --> Extrapolation of D_B^-1 to the original lattice.
//dim(v) = 2 * Ntot, dim(x) = 2 Ntot
//v: input, x: output
void I_D_B_1_It(const c_matrix& U, const spinor& v, spinor& x, const double& m0,const int& block){
    using namespace SAPV;
    bool print_message = false; //good for testing GMRES   

    spinor temp(sap_lattice_sites_per_block, c_vector(2, 0)); 
    spinor temp2(sap_lattice_sites_per_block, c_vector(2, 0)); 

    //temp = I_B^T v
    for (int j = 0; j < sap_lattice_sites_per_block; j++){
        //Writing result to x 
        temp[j][0] = v[SAP_Blocks[block][j]][0];
        temp[j][1] = v[SAP_Blocks[block][j]][1];
    }
    
    //temp2 = D_B^-1 I_B^T v 
    gmres_D_B(U, temp, temp,temp2, m0, 
       sap_gmres_restart_length, sap_gmres_restarts, sap_gmres_tolerance, 
        block, print_message);  

    //x = I_B D_B^-1 I_B^T v 
    set_zeros(x,Ntot,2); //Initialize x to zero
    for (int j = 0; j < sap_lattice_sites_per_block; j++){
        x[SAP_Blocks[block][j]][0] += temp2[j][0];
        x[SAP_Blocks[block][j]][1] += temp2[j][1];
    }

    /*
    spinor temp(sap_lattice_sites_per_block, c_vector(2, 0)); 
    spinor temp2(sap_lattice_sites_per_block, c_vector(2, 0)); 
    It_B_v(v,temp,block); //temp = I_B^T v

    gmres_D_B(U, temp, temp,temp2, m0, 
       sap_gmres_restart_length, sap_gmres_restarts, sap_gmres_tolerance, 
        block, print_message);  //temp2 = D_B^-1 I_B^T v 

    I_B_v(temp2,x,block); //x = I_B D_B^-1 I_B^T v 
    */
}

//Sequential version of the SAP method (used for testing)
int SAP(const c_matrix& U, const spinor& v,spinor& x, const double& m0,const int& nu){
    /*
    Solves D x = v using the SAP method
    */
   using namespace SAPV;
    if (checkSize(v, Ntot, 2) == true || checkSize(x, Ntot, 2) == true){
        std::cout << "Error with vector dimensions in I_KD" << std::endl;
        exit(1);
    }  
    double err;
    double v_norm = sqrt(std::real(dot(v, v))); //norm of the right hand side

    spinor temp(Ntot, c_vector(2, 0)); 
    spinor r(Ntot, c_vector(2, 0)); //residual
    r = v - D_phi(U, x, m0); //r = v - D x
    for (int i = 0; i< nu; i++){
        for (auto block : SAP_RedBlocks){
            I_D_B_1_It(U,r,temp,m0,block);
            //x = x + temp; //x = x + D_B^-1 r
            for(int n = 0; n < Ntot; n++) {
                x[n][0] += temp[n][0];
                x[n][1] += temp[n][1];
            }
        }
        
        r = v - D_phi(U, x, m0); //r = v - D x
        for (auto block : SAP_BlackBlocks){
            I_D_B_1_It(U,r,temp,m0,block);
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
            return 1;
        }
    }
    //std::cout << "SAP did not converge in " << nu << " iterations, error: " << err << std::endl;
    return 0; 
}

int SAP_parallel(const c_matrix& U, const spinor& v,spinor &x, const double& m0,const int& nu,const int& blocks_per_proc){
    /*
    Solves D x = v using the SAP method
    */
   using namespace SAPV;
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   /*
    if (checkSize(v, Ntot, 2) == true || checkSize(x, Ntot, 2) == true){
        std::cout << "Error with vector dimensions in SAP" << std::endl;
        exit(1);
    }  
    */

    /*
    if (blocks_per_proc * size != sap_coloring_blocks) {
        std::cout << "blocks_per_proc * no_of_mpi_processes /= sap_coloring_blocks" << std::endl;
        std::cout << "The workload on the processes is not the same. Some processors invert more SAP blocks" << std::endl;
        exit(1);
    }
    */
   /*
    if (size > sap_coloring_blocks) {
        std::cout << "Error: number of processes > number of blocks" << std::endl;
        exit(1);
    }
    */

    double err;
    double v_norm = sqrt(std::real(dot(v, v))); //norm of the right hand side

    //Divide SAP_RedBlocks among processes
    int start = rank * blocks_per_proc;
    int end = std::min(start + blocks_per_proc, sap_coloring_blocks);

    spinor temp(Ntot, c_vector(2, 0)); 
    spinor r(Ntot, c_vector(2, 0)); //residual

    r = v - D_phi(U, x, m0); //r = v - D x

    spinor local_x(Ntot, c_vector(2, 0)); // Local result for this process
    spinor global_x(Ntot, c_vector(2, 0));

    //Prepare buffers for MPI communication
    c_vector local_buffer(Ntot * 2);
    c_vector global_buffer(Ntot * 2);
    spinor Dphi;
    
    for (int i = 0; i< nu; i++){  
        set_zeros(local_x,Ntot,2); //Initialize local_x to zero
        for (int b = start; b < end; b++) {
            int block = SAP_RedBlocks[b];
            I_D_B_1_It(U, r, temp, m0, block);
            //local_x = local_x + temp; // Local computation
            for(int n = 0; n < Ntot; n++) {
                local_x[n][0] += temp[n][0];
                local_x[n][1] += temp[n][1];
            }
        }

        //------MPI communication for red blocks------//
        //1D array for MPI communication
        for (int n = 0; n < Ntot; ++n) {
            local_buffer[2*n]     = local_x[n][0];
            local_buffer[2*n + 1] = local_x[n][1];
        }
        // Perform single allreduce
        MPI_Allreduce(local_buffer.data(), global_buffer.data(), Ntot * 2, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
        // Unpack back into global_x
        for (int n = 0; n < Ntot; ++n) {
            global_x[n][0] = global_buffer[2*n];
            global_x[n][1] = global_buffer[2*n + 1];
        }
        //---------------------------------------------//
        //x = x + global_x;
        for(int n = 0; n < Ntot; n++) {
            x[n][0] += global_x[n][0];
            x[n][1] += global_x[n][1];
        }
        Dphi = D_phi(U, x, m0);
        //r = v - D_phi(U, x, m0); //r = v - D x
        for(int n = 0; n < Ntot; n++) {
            r[n][0] = v[n][0] - Dphi[n][0];
            r[n][1] = v[n][1] - Dphi[n][1];
        }


        set_zeros(local_x,Ntot,2); //Initialize local_x to zero

        for (int b = start; b < end; b++) {
            int block = SAP_BlackBlocks[b];
            I_D_B_1_It(U, r, temp, m0, block);
            //local_x = local_x + temp; // Local computation
            for(int n = 0; n < Ntot; n++) {
                local_x[n][0] += temp[n][0];
                local_x[n][1] += temp[n][1];
            }
        }

        //------MPI communication for black blocks------//
        for (int n = 0; n < Ntot; ++n) {
            local_buffer[2*n]     = local_x[n][0];
            local_buffer[2*n + 1] = local_x[n][1];
        }
        MPI_Allreduce(local_buffer.data(), global_buffer.data(), Ntot * 2, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
        for (int n = 0; n < Ntot; ++n) {
            global_x[n][0] = global_buffer[2*n];
            global_x[n][1] = global_buffer[2*n + 1];
        }
        //---------------------------------------------//
        
         //x = x + global_x;
        for(int n = 0; n < Ntot; n++) {
            x[n][0] += global_x[n][0];
            x[n][1] += global_x[n][1];
        }
        Dphi = D_phi(U, x, m0);
        //r = v - D_phi(U, x, m0); //r = v - D x
        for(int n = 0; n < Ntot; n++) {
            r[n][0] = v[n][0] - Dphi[n][0];
            r[n][1] = v[n][1] - Dphi[n][1];
        }

        err = sqrt(std::real(dot(r, r))); 
        if (err < sap_tolerance * v_norm) {
            //std::cout << "SAP converged in " << i << " iterations, error: " << err << std::endl;
            return 1;
        }
    }
   //std::cout << "SAP did not converge in " << nu << " iterations, error: " << err << std::endl;
    
    return 0; 
}

int SAPV2(const c_matrix& U, const spinor& v,spinor& x, const double& m0,const int& nu){
    /*
    Solves D x = v using the SAP method
    */
   using namespace SAPV;
    if (checkSize(v, Ntot, 2) == true || checkSize(x, Ntot, 2) == true){
        std::cout << "Error with vector dimensions in I_KD" << std::endl;
        exit(1);
    }  
    double err;
    double v_norm = sqrt(std::real(dot(v, v))); //norm of the right hand side

    spinor temp(sap_lattice_sites_per_block, c_vector(2, 0)); 
    spinor r(Ntot, c_vector(2, 0)); //residual
    r = v - D_phi(U, x, m0); //r = v - D x
    for (int i = 0; i< nu; i++){
        for (auto block : SAP_RedBlocks){
            I_D_B_1_ItV2(U,r,temp,m0,block);
            //x = x + temp; //x = x + D_B^-1 r
            for(int n = 0; n < sap_lattice_sites_per_block; n++) {
                x[SAP_Blocks[block][n]][0] += temp[n][0];
                x[SAP_Blocks[block][n]][1] += temp[n][1];
            }
        }
        

        r = v - D_phi(U, x, m0); //r = v - D x
        for (auto block : SAP_BlackBlocks){
            I_D_B_1_ItV2(U,r,temp,m0,block);
            //x = x + temp; //x = x + D_B^-1 r
            for(int n = 0; n < sap_lattice_sites_per_block; n++) {
                x[SAP_Blocks[block][n]][0] += temp[n][0];
                x[SAP_Blocks[block][n]][1] += temp[n][1];
            }
           
        }
        r = v - D_phi(U, x, m0); //r = v - D x

        err = sqrt(std::real(dot(r, r))); 
        if (err < sap_tolerance * v_norm) {
            //std::cout << "SAP converged in " << i << " iterations, error: " << err << std::endl;
            return 1;
        }
    }
    //std::cout << "SAP did not converge in " << nu << " iterations, error: " << err << std::endl;
    return 0; 
}


int SAP_parallelV2(const c_matrix& U, const spinor& v,spinor &x, const double& m0,const int& nu,const int& blocks_per_proc){
    /*
    Solves D x = v using the SAP method
    */
   using namespace SAPV;
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    double err;
    double v_norm = sqrt(std::real(dot(v, v))); //norm of the right hand side

    //Divide SAP_RedBlocks among processes
    int start = rank * blocks_per_proc;
    int end = std::min(start + blocks_per_proc, sap_coloring_blocks);

    spinor temp(sap_lattice_sites_per_block, c_vector(2, 0)); 
    spinor r(Ntot, c_vector(2, 0)); //residual

    r = v - D_phi(U, x, m0); //r = v - D x

    spinor local_x(Ntot, c_vector(2, 0)); // Local result for this process
    spinor global_x(Ntot, c_vector(2, 0));

    //Prepare buffers for MPI communication
    c_vector local_buffer(Ntot * 2, 0);
    c_vector global_buffer(Ntot * 2, 0);
    spinor Dphi;
    
    for (int i = 0; i< nu; i++){  
        //set_zeros(local_x,Ntot,2); //Initialize local_x to zero
        for(int n = 0; n < Ntot * 2; n++) {
            local_buffer[n] = 0.0; //Initialize local_buffer to zero
        }
        for (int b = start; b < end; b++) {
            int block = SAP_RedBlocks[b];
            I_D_B_1_ItV2(U, r, temp, m0, block);
            //local_x = local_x + temp; // Local computation
            for(int n = 0; n < sap_lattice_sites_per_block; n++) {
                //local_x[SAP_Blocks[block][n]][0] += temp[n][0];
                //local_x[SAP_Blocks[block][n]][1] += temp[n][1];
                local_buffer[2*SAP_Blocks[block][n]]     = temp[n][0];
                local_buffer[2*SAP_Blocks[block][n] + 1] = temp[n][1];
            }
        }

        //------MPI communication for red blocks------//
        //1D array for MPI communication
        //for (int n = 0; n < Ntot; ++n) {
        //    local_buffer[2*n]     = local_x[n][0];
        //    local_buffer[2*n + 1] = local_x[n][1];
        //}
        // Perform single allreduce
        MPI_Allreduce(local_buffer.data(), global_buffer.data(), Ntot * 2, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
        // Unpack back into global_x

        /*
        for (int b = start; b < end; b++) {
            int block = SAP_RedBlocks[b];
            for(int n = 0; n < sap_lattice_sites_per_block; n++) {
                global_x[SAP_Blocks[block][n]][0] = global_buffer[2*SAP_Blocks[block][n]];
                global_x[SAP_Blocks[block][n]][1] = global_buffer[2*SAP_Blocks[block][n] + 1];
            }
        }
        */
      

        //for (int n = 0; n < Ntot; ++n) {
        //    global_x[n][0] = global_buffer[2*n];
        //    global_x[n][1] = global_buffer[2*n + 1];
        //}
        //---------------------------------------------//
        //x = x + global_x;

         //Check how to properly do this by looping over the red blocks only ...

        
        for (int b = 0; b < sap_coloring_blocks; b++) {
            int block = SAP_RedBlocks[b];
            for(int m = 0; m < sap_lattice_sites_per_block; m++) {
                int n = SAP_Blocks[block][m];
                x[n][0] += global_buffer[2*n]; 
                x[n][1] += global_buffer[2*n+1];
            }
        }

        /*for(int n = 0; n < Ntot; n++) {
            x[n][0] += global_buffer[2*n]; //global_x[n][0];
            x[n][1] += global_buffer[2*n+1]; //global_x[n][1];
        }*/
        Dphi = D_phi(U, x, m0);
        //r = v - D_phi(U, x, m0); //r = v - D x
        for(int n = 0; n < Ntot; n++) {
            r[n][0] = v[n][0] - Dphi[n][0];
            r[n][1] = v[n][1] - Dphi[n][1];
        }


        //set_zeros(local_x,Ntot,2); //Initialize local_x to zero
        for(int n = 0; n < Ntot * 2; n++) {
            local_buffer[n] = 0.0; //Initialize local_buffer to zero
        }

        for (int b = start; b < end; b++) {
            int block = SAP_BlackBlocks[b];
            I_D_B_1_ItV2(U, r, temp, m0, block);
            //local_x = local_x + temp; // Local computation
            for(int n = 0; n < sap_lattice_sites_per_block; n++) {
                //local_x[SAP_Blocks[block][n]][0] += temp[n][0];
                //local_x[SAP_Blocks[block][n]][1] += temp[n][1];
                local_buffer[2*SAP_Blocks[block][n]]     = temp[n][0];
                local_buffer[2*SAP_Blocks[block][n] + 1] =  temp[n][1];
            }
        }

        //------MPI communication for black blocks------//
        //for (int n = 0; n < Ntot; ++n) {
        //    local_buffer[2*n]     = local_x[n][0];
        //    local_buffer[2*n + 1] = local_x[n][1];
        //}
        MPI_Allreduce(local_buffer.data(), global_buffer.data(), Ntot * 2, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
        /*
        for (int b = start; b < end; b++) {
            int block = SAP_BlackBlocks[b];
            for(int n = 0; n < sap_lattice_sites_per_block; n++) {
                global_x[SAP_Blocks[block][n]][0] = global_buffer[2*SAP_Blocks[block][n]];
                global_x[SAP_Blocks[block][n]][1] = global_buffer[2*SAP_Blocks[block][n] + 1];
            }
        }
        */

        //for (int n = 0; n < Ntot; ++n) {
        //    global_x[n][0] = global_buffer[2*n];
        //    global_x[n][1] = global_buffer[2*n + 1];
        //}
        //---------------------------------------------//
        
         //x = x + global_x;
        //for(int n = 0; n < Ntot; n++) {
        //    x[n][0] += global_x[n][0];
        //    x[n][1] += global_x[n][1];
        //}
        for (int b = 0; b < sap_coloring_blocks; b++) {
            int block = SAP_BlackBlocks[b];
            for(int m = 0; m < sap_lattice_sites_per_block; m++) {
                int n = SAP_Blocks[block][m];
                x[n][0] += global_buffer[2*n]; 
                x[n][1] += global_buffer[2*n+1];
            }
        }
        /*
        for(int n = 0; n < Ntot; n++) {
            x[n][0] += global_buffer[2*n]; //global_x[n][0];
            x[n][1] += global_buffer[2*n+1]; //global_x[n][1];
        }
        */
        Dphi = D_phi(U, x, m0);
        //r = v - D_phi(U, x, m0); //r = v - D x
        for(int n = 0; n < Ntot; n++) {
            r[n][0] = v[n][0] - Dphi[n][0];
            r[n][1] = v[n][1] - Dphi[n][1];
        }

        err = sqrt(std::real(dot(r, r))); 
        if (err < sap_tolerance * v_norm) {
            //std::cout << "SAP converged in " << i << " iterations, error: " << err << std::endl;
            return 1;
        }
    }
   //std::cout << "SAP did not converge in " << nu << " iterations, error: " << err << std::endl;
    
    return 0; 
}

//A_B = I_B * D_B^-1 * I_B^T v --> Extrapolation of D_B^-1 to the original lattice.
//dim(v) = 2 * Ntot, dim(x) = 2 Ntot
//v: input, x: output
void I_D_B_1_ItV2(const c_matrix& U, const spinor& v, spinor& x, const double& m0,const int& block){
    using namespace SAPV;
    bool print_message = false; //good for testing GMRES   

    spinor temp(sap_lattice_sites_per_block, c_vector(2, 0)); 

    //temp = I_B^T v
    for (int j = 0; j < sap_lattice_sites_per_block; j++){
        //Writing result to x 
        temp[j][0] = v[SAP_Blocks[block][j]][0];
        temp[j][1] = v[SAP_Blocks[block][j]][1];
    }
     
    set_zeros(x,sap_lattice_sites_per_block,2); //Initialize x to zero
    gmres_D_B(U, temp, temp,x, m0, 
       sap_gmres_restart_length, sap_gmres_restarts, sap_gmres_tolerance, 
        block, print_message);  


}