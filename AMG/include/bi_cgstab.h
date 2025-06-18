#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "matrix_operations.h"
#include "operator_overloads.h"
#include "variables.h"
#include <cmath>
#include <iostream>

spinor bi_cgstab(const c_matrix& U, const c_matrix& phi, const c_matrix& x0,const double& m0, const int& max_iter, const double& tol, const bool& print_message); //D^-1 phi 

// c_matrix (*func)(const c_matrix&, const c_matrix&, const double&), const std::function<c_matrix(c_matrix, c_matrix, double)>& func
spinor bi_cgstabV2(c_matrix (*func)(const c_matrix&, const c_matrix&, const double&), const int& dim1, const int& dim2,
    const c_matrix& U, const c_matrix& phi, const c_matrix& x0, const double& m0, const int& max_iter, const double& tol, const bool& print_message);

#endif
