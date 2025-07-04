#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "dirac_operator.h"
#include <cmath>
#include <iostream>

//Conjugate gradient for computing (DD^dagger)^-1 phi, where phi is a vector represented by a matrix
//phi[Ntot][2]
spinor conjugate_gradient(const c_matrix& U, const spinor& phi, const double& m0); 

//Bi conjugate gradient for inverting D 
spinor bi_cgstab(const c_matrix& U, const spinor& phi, const spinor& x0,const double& m0, const int& max_iter, const double& tol, const bool& print_message); //D^-1 phi

#endif
