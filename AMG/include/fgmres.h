#ifndef FGMRES_H
#define FGMRES_H

#include "amg.h"

/*
	FGMRES with AMG as preconditioner
    U: gauge configuration
    phi: right hand side
    x0: initial guess
    m0: mass parameter
    m: restart length
    restarts: number of restarts of length m
    tol: tolerance for the solver
    print_message: if true, print the convergence message

    The convergence criterion is ||r|| < ||phi|| * tol
*/
spinor fgmresAMG(const c_matrix& U, const spinor& phi, const spinor& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message);

/*
	FGMRES with SAP as preconditioner
    Parameters same as above.
*/
spinor fgmresSAP(const c_matrix& U, const spinor& phi, const spinor& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message);

/*
    Rotations to transform Hessenberg matrix to upper triangular form
    cn: cosine components of the rotation
    sn: sine components of the rotation
    H: Hessenberg matrix
    j: index of the column being processed
*/
void rotation(c_vector& cn, c_vector& sn, c_matrix& H, const int& j);

/*
    Solves an upper triangular system Ax = b, where A is an upper triangular matrix of dimension n
    A: upper triangular matrix
    b: right-hand side vector
    n: dimension of the matrix
    out: output vector where the solution will be stored
*/
void solve_upper_triangular(const c_matrix& A, const c_vector& b, const int& n, c_vector& out);

#endif