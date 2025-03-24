#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "matrix_operations.h"
#include "operator_overloads.h"
#include "variables.h"
#include <cmath>
#include <iostream>

c_matrix canonical_vector(const int& i, const int& N1, const int& N2);

//Conjugate gradient for computing (DD^dagger)^-1 phi, where phi is a vector represented by a matrix
//phi[Ntot][2]
c_matrix conjugate_gradient(const c_matrix& U, const c_matrix& phi, const double& m0); 

c_matrix bi_cgstab(const c_matrix& U, const c_matrix& phi, const c_matrix& x0,const double& m0, const int& max_iter, const double& tol, const bool& print_message); //D^-1 phi 

c_matrix kaczmarz(const c_matrix& U, const c_matrix& phi, const c_matrix& x0, const double& m0, const int& max_iter, const double& tol, const bool& print_message); //D^-1 phi
#endif
