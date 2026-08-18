#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "dirac_operator.h"

/*
    Conjugate gradient method for computing (DD^dagger)^-1 phi 
    U: gauge configuration
    phi: right-hand side vector
    m0: mass parameter for Dirac matrix 
        
    The convergence criterion is ||r|| < ||phi|| * tol
*/

//Buffers
namespace CG{
    extern spinor r;  //r[coordinate][spin] residual
    extern spinor d; //search direction
    extern spinor Ad; //DD^dagger*d
}
    
int conjugate_gradient(const spinor& U, const spinor& phi, spinor &x); 


#endif
