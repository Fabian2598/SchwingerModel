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
    Check if the size of an input spinor is correct.
    Returns true if the size is not correct, false otherwise.
    Only performs the check if DEBUG is defined. Check config.h.
*/
inline bool checkSize(const spinor& v, const int& N1, const int& N2) {
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
    U: gauge configuration,
    phi: right-hand side,
    x0: initial guess, 
    x: output
    phi: right-hand side
    m0: mass parameter,
    m: restart length,
    restarts: number of restarts,
    tol: tolerance for convergence,
    block: Schwarz block index,
    print_message: if true, prints convergence messages
    Returns 1 if converged, 0 if not converged.

    The convergence criterion is ||r|| < ||phi|| * tol
*/
int gmres_D_B(const c_matrix& U, const spinor& phi, const spinor& x0, spinor& x, const double& m0, 
    const int& m, const int& restarts, const double& tol, const int& block,const bool& print_message);

/*
    A_B v = I_B * D_B^-1 * I_B^T v --> Extrapolation of D_B^-1 to the original lattice.
    dim(v) = 2 * Ntot, dim(x) = 2 Ntot
    v: input, x: output 
*/
void I_D_B_1_It(const c_matrix& U, const spinor& v, spinor& x, const double& m0,const int& block);

/*
    Sequential version of the SAP method (used for testing)
    dim(v) = 2 * Ntot, dim(x) = 2 * Ntot
*/
int SAP(const c_matrix& U, const spinor& v,spinor &x, const double& m0,const int& nu);

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
int SAP_parallel(const c_matrix& U, const spinor& v,spinor &x, const double& m0,const int& nu,const int& blocks_per_proc);



#endif