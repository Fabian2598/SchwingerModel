#include "sap.h"

void SchwarzBlocks(){
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
                    SAP_Blocks[block][++count] = x * Nt+ t; 
                    //Each block also considers both spin components, 
                    //but we only reference the lattice coordinates here.
                }
            }
        
        }
    }
}

//x = I_B^T v --> Restriction of the vector v to the block B
//dim(v) = 2 Ntot, dim(x) = 2 * block_elements
c_matrix It_B_v(const c_matrix& v, const int& block){
    //Schwarz blocks have to be initialized first before calling this function
    int block_elements = sap_x_elements * sap_t_elements; //Number of elements in the block
    c_matrix x(block_elements, c_vector(2, 0)); 
    for (int j = 0; j < block_elements; j++){
        x[j][0] = v[SAP_Blocks[block][j]][0];
        x[j][1] = v[SAP_Blocks[block][j]][1];
    }
    return x;
}

// x = I_B v --> Interpolation of the vector v to the block B
//dim(v) = 2 * block_elements, dim(x) = 2 Ntot
c_matrix I_B_v(const c_matrix& v, const int& block){
    //Schwarz blocks have to be initialized first before calling this function
    int block_elements = sap_x_elements * sap_t_elements; //Number of elements in the block
    c_matrix x(Ntot, c_vector(2, 0)); 
    for (int j = 0; j < block_elements; j++){
        x[SAP_Blocks[block][j]][0] += v[j][0];
        x[SAP_Blocks[block][j]][1] += v[j][1];
    }
    return x;
}

//Restriction of the Dirac operator to the block B
//D_B phi = I_B^T D I_B phi  
c_matrix D_B(const std::vector<std::vector<std::complex<double>>>& U, const std::vector<std::vector<std::complex<double>>>& phi, const double& m0,const int& block){
    return It_B_v( D_phi(U,I_B_v(phi,block), m0),  block); 
}

