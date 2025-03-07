#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "matrix_operations.h"
#include <iostream>

//Conjugate gradient for computing (DD^dagger)^-1 phi, where phi is a vector represented by a matrix
//phi[Ntot][2]
c_matrix conjugate_gradient(const c_matrix& U, const c_matrix& phi, const double& m0); 

#endif
