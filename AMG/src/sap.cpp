#include "sap.h"

void SchwarzBlocks(){
    schwarz_blocks = true; //Schwarz blocks are initialized
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
                    //but we only reference the lattice coordinates here.
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
            //-----------------------------------//
        }
    }
}

//x = I_B^T v --> Restriction of the vector v to the block B
//dim(v) = 2 Ntot, dim(x) = 2 * sap_lattice_sites_per_block
void It_B_v(const c_matrix& v, c_matrix& x, const int& block){
    //Schwarz blocks have to be initialized first before calling this function
    if (!schwarz_blocks){
        std::cout << "Error: Schwarz blocks not initialized" << std::endl;
        exit(1);
    }
    //Checking dimensions
    if ( checkSize(v, Ntot, 2)== true || checkSize(x, sap_lattice_sites_per_block, 2) == true){
        std::cout << "Error with vector dimensions in It_B_v" << std::endl;
        exit(1);
    } 
    set_zeros(x,sap_lattice_sites_per_block,2); //Initialize x to zero

    for (int j = 0; j < sap_lattice_sites_per_block; j++){
        x[j][0] = v[SAP_Blocks[block][j]][0];
        x[j][1] = v[SAP_Blocks[block][j]][1];
    }
}

// x = I_B v --> Interpolation of the vector v to the original lattice
//dim(v) = 2 * sap_lattice_sites_per_block, dim(x) = 2 Ntot
void I_B_v(const c_matrix& v, c_matrix& x,const int& block){
    //Schwarz blocks have to be initialized first before calling this function
    if (!schwarz_blocks){
        std::cout << "Error: Schwarz blocks not initialized" << std::endl;
        exit(1);
    }
    if ( checkSize(v, sap_lattice_sites_per_block, 2) == true ||  checkSize(x, Ntot, 2) == true ){
        std::cout << "Error with dimensions in I_B_v" << std::endl;
        exit(1);
    } 
    set_zeros(x,Ntot,2); //Initialize x to zero
    for (int j = 0; j < sap_lattice_sites_per_block; j++){
        x[SAP_Blocks[block][j]][0] += v[j][0];
        x[SAP_Blocks[block][j]][1] += v[j][1];
    }

}

//Restriction of the Dirac operator to the block B
//D_B v = I_B^T D I_B v  
//dim(v) = 2 * sap_lattice_sites_per_block, dim(x) = 2 * sap_lattice_sites_per_block
void D_B(const c_matrix& U, const c_matrix& v, c_matrix& x, const double& m0,const int& block){
    if (checkSize(v, sap_lattice_sites_per_block, 2)==true || checkSize(x, sap_lattice_sites_per_block, 2) == true){
        std::cout << "Error with vector dimensions in D_B" << std::endl;
        exit(1);
    }
    c_matrix temp(Ntot, c_vector(2, 0)); //I_B v
    I_B_v(v,temp,block); 
    It_B_v( D_phi(U,temp, m0), x, block); //It_B_v sets x to zero first
}

//Solves D_B x = phi using GMRES, where D_B is the Dirac matrix restricted to the Schwarz block B 
int gmres_D_B(const c_matrix& U, const c_matrix& phi, const c_matrix& x0, c_matrix& x,
    const double& m0, const int& m, const int& restarts, const double& tol, const int& block,
    const bool& print_message) {
    //GMRES for D^-1 phi
    //phi --> right-hand side
    //x0 --> initial guess  
    //x --> solution
    //U --> configuration
    //restarts --> number of restarts
    //m --> number of iterations per cycle
    if (print_message == true){
        std::cout << "------------------------------------------" << std::endl;
        std::cout << "|GMRES for D_B (block " << block << ") " << std::endl;
        std::cout << "|Restart length " << m << ". Restarts = " << restarts << std::endl;
    }
    if (checkSize(phi, sap_lattice_sites_per_block, 2) == true || checkSize(x, sap_lattice_sites_per_block, 2) == true || checkSize(x, sap_lattice_sites_per_block, 2) == true){
        std::cout << "Error with vector dimensions in gmres_D_B" << std::endl;
        exit(1);
    } 
    if (m> sap_variables_per_block) {
        std::cout << "Error: restart length > sap_variables_per_block" << std::endl;
        exit(1);
    }
    
    int k = 0; //Iteration number (restart cycle)
    double err; //Error = ||r||_2

    c_matrix r(sap_lattice_sites_per_block, c_vector(2, 0));  //r[coordinate][spin] residual
    c_matrix r0(sap_lattice_sites_per_block, c_vector(2, 0)); //Initial residual

    //VmT[column vector index][vector arrange in matrix form]
    std::vector<c_matrix> VmT(m+1, c_matrix(sap_lattice_sites_per_block, c_vector(2, 0))); //V matrix transpose --> dimensions exchanged

    c_matrix Hm(m+1 , c_vector(m, 0)); //H matrix (Hessenberg matrix)
    c_vector gm(m + 1, 0); 

    //Elements of rotation matrix |sn[i]|^2 + |cn[i]|^2 = 1
    c_vector sn(m, 0);
    c_vector cn(m, 0);
    c_vector eta(m, 0);


    c_matrix w(sap_lattice_sites_per_block, c_vector(2, 0)); //D*d
    x = x0; //initial solution
    c_double beta; //not 1/g^2

	c_matrix temp(sap_lattice_sites_per_block, c_vector(2, 0)); //I_B v
    D_B(U, phi,temp, m0, block);
    r0 = phi - temp; //r = b - A*x
	
	double norm_phi = sqrt(std::real(dot(phi, phi))); //norm of the right hand side
    if (print_message == true){std::cout << "||phi|| * tol = " << norm_phi * tol << std::endl;}
    while (k < restarts) {
        beta = sqrt(std::real(dot(r0, r0))) + 0.0 * I_number;
        VmT[0] = 1.0 / beta * r0;
        gm[0] = beta; //gm[0] = ||r||
        //-----Arnoldi process to build the Krylov basis and the Hessenberg matrix-----//
        for (int j = 0; j < m; j++) {
            D_B(U, VmT[j], w,m0, block); //w = D v_j  

            for (int i = 0; i <= j; i++) {
                Hm[i][j] = dot(w, VmT[i]); //  (v_i^dagger, w)
                w = w -  Hm[i][j] * VmT[i];
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
			 it_count = k + 1;
             if (print_message == true) {
                 std::cout << "|GMRES for D_B (block " << block << ") converged in " << k + 1 << " restarts." << " Error " << err << std::endl;
                 std::cout << "------------------------------------------" << std::endl;
            }
             
             return 1;
             
         }
         r0 = r;
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
void I_D_B_1_It(const c_matrix& U, const c_matrix& v, c_matrix& x, const double& m0,const int& block){
    bool print_message = false; //good for testing GMRES   
    int rank; 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0){
        print_message = false;
    }
    
    if (checkSize(v, Ntot, 2) == true || checkSize(x, Ntot, 2) == true){
        std::cout << "Error with vector dimensions in I_D_B_1_It" << std::endl;
        exit(1);
    } 
    c_matrix temp(sap_lattice_sites_per_block, c_vector(2, 0)); 
    c_matrix temp2(sap_lattice_sites_per_block, c_vector(2, 0)); 
    It_B_v(v,temp,block); //temp = I_B^T v
    gmres_D_B(U, temp, temp,temp2, m0, 
        sap_gmres_restart_length, sap_gmres_restarts, sap_gmres_tolerance, 
        block, print_message);  //temp2 = D_B^-1 I_B^T v 
    I_B_v(temp2,x,block); //x = I_B D_B^-1 I_B^T v 
}


int SAP(const c_matrix& U, const c_matrix& v,c_matrix &x, const double& m0,const int& nu){
    /*
    Solves D x = v using the SAP method
    */
    if (checkSize(v, Ntot, 2) == true || checkSize(x, Ntot, 2) == true){
        std::cout << "Error with vector dimensions in I_KD" << std::endl;
        exit(1);
    }  
    double err;
    double v_norm = sqrt(std::real(dot(v, v))); //norm of the right hand side

    c_matrix temp(Ntot, c_vector(2, 0)); 
    c_matrix r(Ntot, c_vector(2, 0)); //residual
    //set_zeros(x,Ntot,2); //Initialize x to zero
    r = v - D_phi(U, x, m0); //r = v - D x
    for (int i = 0; i< nu; i++){
        for (auto block : SAP_RedBlocks){
            I_D_B_1_It(U,r,temp,m0,block);
            x = x + temp; //x = x + D_B^-1 r
        }
        
        r = v - D_phi(U, x, m0); //r = v - D x
        for (auto block : SAP_BlackBlocks){
            I_D_B_1_It(U,r,temp,m0,block);
            x = x + temp; //x = x + D_B^-1 r
        }
        r = v - D_phi(U, x, m0); //r = v - D x

        err = sqrt(std::real(dot(r, r))); 
        if (err < sap_tolerance * v_norm) {
            //std::cout << "SAP converged in " << i << " iterations, error: " << err << std::endl;
            return 1;
        }
    }
    //std::cout << "SAP did not converge in " << nu << " iterations, error: " << err << std::endl;
    return 0; //Not converged
}

int SAP_parallel(const c_matrix& U, const c_matrix& v,c_matrix &x, const double& m0,const int& nu,const int& blocks_per_proc){
    /*
    Solves D x = v using the SAP method
    */
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (checkSize(v, Ntot, 2) == true || checkSize(x, Ntot, 2) == true){
        std::cout << "Error with vector dimensions in SAP" << std::endl;
        exit(1);
    }  
    if (blocks_per_proc * size != sap_coloring_blocks) {
        std::cout << "blocks_per_proc * no_of_mpi_processes /= sap_coloring_blocks" << std::endl;
        std::cout << "The workload on the processes is not the same. Some processors invert more SAP blocks" << std::endl;
        exit(1);
    }
    if (size > sap_coloring_blocks) {
        std::cout << "Error: number of processes > number of blocks" << std::endl;
        exit(1);
    }


    double err;
    double v_norm = sqrt(std::real(dot(v, v))); //norm of the right hand side

    //Divide SAP_RedBlocks among processes
    //Red blocks size equal to 
    int start = rank * blocks_per_proc;
    int end = std::min(start + blocks_per_proc, sap_coloring_blocks);

    c_matrix temp(Ntot, c_vector(2, 0)); 
    c_matrix r(Ntot, c_vector(2, 0)); //residual
    //set_zeros(x,Ntot,2); //Initialize x to zero

    r = v - D_phi(U, x, m0); //r = v - D x

    c_matrix local_x(Ntot, c_vector(2, 0)); // Local result for this process
    c_matrix global_x(Ntot, c_vector(2, 0));

    //Prepare buffers for MPI communication
    std::vector<c_double> local_buffer(Ntot * 2);
    std::vector<c_double> global_buffer(Ntot * 2);

    
    for (int i = 0; i< nu; i++){  
        set_zeros(local_x,Ntot,2); //Initialize local_x to zero
        for (int b = start; b < end; b++) {
            int block = SAP_RedBlocks[b];
            I_D_B_1_It(U, r, temp, m0, block);
            local_x = local_x + temp; // Local computation
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

        x = x + global_x;
        r = v - D_phi(U, x, m0); //r = v - D x
        set_zeros(local_x,Ntot,2); //Initialize local_x to zero

        for (int b = start; b < end; b++) {
            int block = SAP_BlackBlocks[b];
            I_D_B_1_It(U, r, temp, m0, block);
            local_x = local_x + temp; // Local computation
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
        
        x =  x + global_x;
        r = v - D_phi(U, x, m0); //r = v - D x

        err = sqrt(std::real(dot(r, r))); 
        if (err < sap_tolerance * v_norm) {
            //std::cout << "SAP converged in " << i << " iterations, error: " << err << std::endl;
            return 1;
        }
    }
   //std::cout << "SAP did not converge in " << nu << " iterations, error: " << err << std::endl;
    
    return 0; //Not converged
}

