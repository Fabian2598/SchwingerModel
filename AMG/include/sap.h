#ifndef SAP_H
#define SAP_H
#include "gmres.h"
#include "mpi.h"

/*
    Build the Schwarz blocks
    Function only has to be called once before using the SAP method.
*/
void SchwarzBlocks();

/*
    Set all elements of a spinor to zero.
    v: input spinor, N1: first index dimension, N2: second index dimension
*/
inline void set_zeros(spinor& v, const int& N1, const int& N2) {
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            v[i][j] = 0;
        }
    }
}

//Schwarz blocks have to be initialized first before calling the I_B operators

/*
    I_B^T v --> Restriction of the spinor v to the block B
    v: input, x: output
*/
void It_B_v(const spinor& v,  spinor& x,const int& block);

/*
    I_B v --> Interpolation of the spinor v to the original lattice
    v: input, x: output
*/
void I_B_v(const spinor& v, spinor& x ,const int& block);

/*
    Restriction of the Dirac operator to the block B. 
    D_B = I_B^T D I_B
    v: input, x: output
*/
void D_B(const c_matrix& U, const spinor& v, spinor& x,const double& m0,const int& block);

/*
    Solves D_B x = phi using GMRES, where D_B is the Dirac matrix restricted to the Schwarz block B     
*/
class GMRES_D_B : public GMRES {
    public:GMRES_D_B(const int& dim1, const int& dim2, const int& m, const int& restarts, const double& tol) :
     GMRES(dim1, dim2, m, restarts, tol) {
        if (m > SAPV::sap_variables_per_block) {
            std::cout << "Error: restart length > sap_variables_per_block" << std::endl;
            exit(1);
        }
    };
    ~GMRES_D_B() { };

    /*
    Set the block index for the GMRES_D_B operator.
    */
    void set_block(const int& block_index) { 
        block = block_index;
    }

    void set_params(const c_matrix& conf, const double& bare_mass){
        U = &conf;
        m0 = &bare_mass;
    }
    
private:
    const c_matrix* U; 
    const double* m0; 
    int block; //block index for the Schwarz block

    /*
    Implementation of the function that computes the matrix-vector product for the fine level
    */
    void func(const spinor& in, spinor& out) override { 
        D_B(*U, in, out, *m0, block);
    }
};

/*
    GMRES solver for D_B
*/
extern GMRES_D_B gmres_DB;

/*
    A_B v = I_B * D_B^-1 * I_B^T v --> Extrapolation of D_B^-1 to the original lattice.
    dim(v) = 2 * Ntot, dim(x) = 2 Ntot
    v: input, x: output 
*/
void I_D_B_1_It(const c_matrix& U, const spinor& v, spinor& x, const double& m0,const int& block);



/*
    Parallel version of the SAP method.
    Solves D x = v using the SAP method.
    U: gauge configuration,
    v: right-hand side,
    x: output,
    m0: mass parameter,
    nu: number of iterations,
    blocks_per_proc: number of blocks per process

    The convergence criterion is ||r|| < ||phi|| * tol
*/
int SAP(const c_matrix& U, const spinor& v,spinor &x, const double& m0,const int& nu,const int& blocks_per_proc);

/*
    Given a lattice point index n, it returns the corresponding 
    SAP block index and the local index m within that block.
*/
inline void getMandBlock(const int& n, int &m, int &block) {
    int x = n / LV::Nt; //x coordinate of the lattice point 
    int t = n % LV::Nt; //t coordinate of the lattice point
    //Reconstructing the block and m index from x and t
    int block_x = x / SAPV::sap_x_elements; //Block index in the x direction
    int block_t = t / SAPV::sap_t_elements; //Block index in the t direction
    block = block_x * SAPV::sap_block_t + block_t; //Block index in the SAP method

    int mx = x % SAPV::sap_x_elements; //x coordinate in the block
    int mt = t % SAPV::sap_t_elements; //t coordinate in the block
    m = mx * SAPV::sap_t_elements + mt; //Index in the block
}


#endif