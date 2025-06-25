#ifndef GMRES_H_INCLUDED
#define GMRES_H_INCLUDED

#include "dirac_operator.h"
#include "operator_overloads.h"
#include "variables.h"
#include <cmath>
#include <iostream>

/*
    Implementation of GMRES for inverting D
    func: matrix-vector operation 
    dim1: dimension of the first index of the spinor 
  	dim2: dimension of the second index of the spinor 	
        ->The solution is of the form x[dim1][dim2]
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
spinor gmres(spinor (*func)(const c_matrix&, const c_matrix&, const double&),const int& dim1, const int& dim2,
const c_matrix& U,const spinor& phi, const spinor& x0, const double& m0,const int& m, const int& restarts, const double& tol, const bool& print_message);

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
*/
c_vector solve_upper_triangular(const c_matrix& A, const c_vector& b, const int& n);


#endif 