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

int conjugate_gradient_v2(const spinor_v2& U, const spinor_v2& phi, spinor_v2 &x); 


#endif
