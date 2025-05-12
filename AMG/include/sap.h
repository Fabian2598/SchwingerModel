//Schwarz alternating method
#ifndef SAP_H
#define SAP_H
#include "gmres.h"
#include "mpi.h"

//Build the blocks for the Schwarz alternating method
void SchwarzBlocks();

//Verify that input dimensions are correct
inline bool checkSize(const c_matrix& v, const int& N1, const int& N2) {
    //DEBUG defined in config.h.in
    #ifdef DEBUG
    int size1 = v.size(), size2 = v[0].size();
    if (size1 != N1 || size2 != N2) {
        std::cerr << "Error: Expected dimensions (" << N1 << ", " << N2 
                  << "), but got (" << size1 << ", " << size2 << ")." << std::endl;
        return true;
    }
    #endif
    return false;
}
inline void set_zeros(c_matrix& v, const int& N1, const int& N2) {
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            v[i][j] = 0;
        }
    }
}

//Schwarz blocks have to be initialized first before calling the I_B operators

// I_B^T v --> Restriction of the vector v to the block B
//v: input, x: output
void It_B_v(const c_matrix& v,  c_matrix& x,const int& block);

//I_B v --> Interpolation of the vector v to the original lattice
//v: input, x: output
void I_B_v(const c_matrix& v, c_matrix& x ,const int& block);

//Restriction of the Dirac operator to the block B. D_B = I_B^T D I_B
//v: input, x: output
void D_B(const c_matrix& U, const c_matrix& v, c_matrix& x,const double& m0,const int& block);

//Solves D_B x = phi using GMRES, where D_B is the Dirac matrix restricted to the Schwarz block B 
//x0: initial guess, x: output
//phi: right-hand side
int gmres_D_B(const c_matrix& U, const c_matrix& phi, const c_matrix& x0, c_matrix& x, const double& m0, 
    const int& m, const int& restarts, const double& tol, const int& block,const bool& print_message);
//returns 1 if converged, 0 if not converged

//A_B v = I_B * D_B^-1 * I_B^T v --> Extrapolation of D_B^-1 to the original lattice.
//dim(v) = 2 * Ntot, dim(x) = 2 Ntot
//v: input, x: output
void I_D_B_1_It(const c_matrix& U, const c_matrix& v, c_matrix& x, const double& m0,const int& block);

//dim(v) = 2 * Ntot, dim(x) = 2 * Ntot
int SAP(const c_matrix& U, const c_matrix& v,c_matrix &x, const double& m0,const int& nu);

int SAP_parallel(const c_matrix& U, const c_matrix& v,c_matrix &x, const double& m0,const int& nu,const int& blocks_per_proc);


//K matrix B_black (I - D B_red) + B_red
//dim(v) = 2 * Ntot, dim(x) = 2 * Ntot
void K_SAP(const c_matrix& U, const c_matrix& v, c_matrix& x, const double& m0);

//I - K D
//dim(v) = 2 * Ntot, dim(x) = 2 * Ntot
void I_KD(const c_matrix& U, const c_matrix& v, c_matrix& x, const double& m0);

//M_SAP^(nu) v = sum_l=0^nu-1 (I - K D)^(l) K v
//dim(v) = 2 * Ntot, dim(x) = 2 * Ntot
int SAP_V2(const c_matrix& U, const c_matrix& v,c_matrix &x, const double& m0,const int& nu);
//This version of SAP implements the M_SAP^(nu) matrix, but 
//according to: Lüscher, M. (2003). Solution of the Dirac equation 
//in lattice QCD using a domain decomposition method, that is not optimal.



#endif