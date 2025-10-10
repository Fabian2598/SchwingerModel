#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "dirac_operator.h"
#include <cmath>
#include <iostream>

/*
    Conjugate gradient method for computing (DD^dagger)^-1 phi 
    U: gauge configuration
    phi: right-hand side vector
    m0: mass parameter for Dirac matrix 
        
    The convergence criterion is ||r|| < ||phi|| * tol
*/
int conjugate_gradient(const spinor& U, const spinor& phi, spinor &x, const double& m0); 

/*
    Bi-cgstab for inverting the Dirac matrix
    U: gauge configuration
    phi: right-hand side vector
    x0: initial guess
    m0: mass parameter for Dirac matrix 
    max_iter: maximum number of iterations
    tol: tolerance for convergence
    print_message: flag for printing convergence messages
        
    The convergence criterion is ||r|| < ||phi|| * tol
*/
spinor bi_cgstab(const spinor& U, const spinor& phi, const spinor& x0,const double& m0, const int& max_iter, const double& tol); //D^-1 phi

#endif
