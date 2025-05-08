#include "sap.h"

//Verify that input dimensions are correct
bool checkSize(const c_matrix& v, const int& N1, const int& N2) {
    int size1 = v.size(), size2 = v[0].size();
    if (size1!=  N1) {
        std::cout << "Error: v.size() is" << size1 << " and v[0].size is " << size2 << std::endl;
        std::cout << "v.size() has to be " << N1 << ". Instead it is " << size1 << std::endl;
        return true;
    }
    if (size2 != N2) {
        std::cout << "Error: v.size() is" << size1 << " and v[0].size is " << size2 << std::endl;
        std::cout << "v[0].size() has to be " << N2 << ". Instead it is " << size2 << std::endl;
        return true;
    }
    return false;
}

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
c_matrix It_B_v(const c_matrix& v, const int& block){
    //Schwarz blocks have to be initialized first before calling this function
    if (!schwarz_blocks){
        std::cout << "Error: Schwarz blocks not initialized" << std::endl;
        exit(1);
    }
    //Checking dimensions
    if ( checkSize(v, Ntot, 2)== true ){
        std::cout << "Error with vector dimensions in It_B_v" << std::endl;
        exit(1);
    } 
    
    c_matrix x(sap_lattice_sites_per_block, c_vector(2, 0)); 
    for (int j = 0; j < sap_lattice_sites_per_block; j++){
        x[j][0] = v[SAP_Blocks[block][j]][0];
        x[j][1] = v[SAP_Blocks[block][j]][1];
    }
    return x;
}

// x = I_B v --> Interpolation of the vector v to the original lattice
//dim(v) = 2 * sap_lattice_sites_per_block, dim(x) = 2 Ntot
c_matrix I_B_v(const c_matrix& v, const int& block){
    //Schwarz blocks have to be initialized first before calling this function
    if (!schwarz_blocks){
        std::cout << "Error: Schwarz blocks not initialized" << std::endl;
        exit(1);
    }
    if ( checkSize(v, sap_lattice_sites_per_block, 2) == true ){
        std::cout << "Error with dimensions in I_B_v" << std::endl;
        exit(1);
    } 

    c_matrix x(Ntot, c_vector(2, 0)); 
    for (int j = 0; j < sap_lattice_sites_per_block; j++){
        x[SAP_Blocks[block][j]][0] += v[j][0];
        x[SAP_Blocks[block][j]][1] += v[j][1];
    }
    return x;
}

//Restriction of the Dirac operator to the block B
//D_B v = I_B^T D I_B v  
//dim(v) = 2 * sap_lattice_sites_per_block, dim(x) = 2 * sap_lattice_sites_per_block
c_matrix D_B(const c_matrix& U, const c_matrix& v, const double& m0,const int& block){
    if (checkSize(v, sap_lattice_sites_per_block, 2)==true){
        std::cout << "Error with vector dimensions in D_B" << std::endl;
        exit(1);
    }
    return It_B_v( D_phi(U,I_B_v(v,block), m0),  block); 
}

//Solves D_B x = phi using GMRES, where D_B is the Dirac matrix restricted to the Schwarz block B 
c_matrix gmres_D_B(const c_matrix& U, const c_matrix& phi, const c_matrix& x0, const double& m0, const int& m, const int& restarts, const double& tol, const int& block,const bool& print_message) {
    //GMRES for D^-1 phi
    //phi --> right-hand side
    //x0 --> initial guess  
    //U --> configuration
    //restarts --> number of restarts
    //m --> number of iterations per cycle
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "|GMRES for D_B (block " << block << ") " << std::endl;
    std::cout << "|Restart length " << m << ". Restarts = " << restarts << std::endl;
    if (checkSize(phi, sap_lattice_sites_per_block, 2) == true){
        std::cout << "Error with vector dimensions in gmres_D_B" << std::endl;
        exit(1);
    } 
    if (m> sap_variables_per_block) {
        std::cout << "Error: restart length > sap_variables_sites_per_block" << std::endl;
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
    c_matrix x = x0; //initial solution
    c_double beta; //not 1/g^2

	
    r0 = phi - D_B(U, phi, m0, block); //r = b - A*x
	
	double norm_phi = sqrt(std::real(dot(phi, phi))); //norm of the right hand side
    while (k < restarts) {
        beta = sqrt(std::real(dot(r0, r0))) + 0.0 * I_number;
        VmT[0] = 1.0 / beta * r0;
        gm[0] = beta; //gm[0] = ||r||
        //-----Arnoldi process to build the Krylov basis and the Hessenberg matrix-----//
        for (int j = 0; j < m; j++) {
            w = D_B(U, VmT[j], m0, block); //w = D v_j 
            //This for loop is the most time consuming part ...
            //For the values of m that I will use, parallelizing is not really useful due to the overhead
            //#pragma omp parallel for
            for (int i = 0; i <= j; i++) {
                Hm[i][j] = dot(w, VmT[i]); //  (v_i^dagger, w)
                //#pragma omp critical
                w = w -  Hm[i][j] * VmT[i];
            }
            //i.e. the Gramm Schmidt part is highly inefficient. 
            //Could a better implementation of the dot product improve the execution time?
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
        r = phi - D_B(U, x, m0, block); 
        err = sqrt(std::real(dot(r, r)));
        //if (print_message == true) {
        //    std::cout << "GMRES for D_B (block " << block << ") "<< k + 1 << " restart cycle" << " Error " << err << std::endl;
        //}

         if (err < tol * norm_phi) {
			 it_count = k + 1;
             if (print_message == true) {
                 std::cout << "|GMRES for D_B (block " << block << ") converged in " << k + 1 << " restarts." << " Error " << err << std::endl;
             }
             std::cout << "------------------------------------------" << std::endl;
             return x;
             
         }
         r0 = r;
         k++;
    }
    //if (print_message == true) {
        std::cout << "|GMRES for D_B (block " << block << ") did not converge in " << restarts << " restarts." << " Error " << err << std::endl;
    //
    std::cout << "------------------------------------------" << std::endl;
    return x;
}

//I_B * D_B^-1 * I_B^T v --> Extrapolation of D_B^-1 to the original lattice.
//dim(v) = 2 * Ntot, dim(x) = 2 Ntot
c_matrix I_D_B_1_It(const c_matrix& U, const c_matrix& v, const double& m0,const int& block){
    bool print_message = true; //good for testing
    if (checkSize(v, Ntot, 2) == true){
        std::cout << "Error with vector dimensions in I_D_B_1_It" << std::endl;
        exit(1);
    } 
    return I_B_v( 
        gmres_D_B(U, It_B_v(v,block), It_B_v(v,block), m0, sap_gmres_restart_length, sap_gmres_restarts, sap_gmres_tolerance, block, print_message),
        block); 
}


