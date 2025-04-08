#ifndef GMRES_H_INCLUDED
#define GMRES_H_INCLUDED

#include "matrix_operations.h"
#include "operator_overloads.h"
#include "variables.h"
#include <cmath>
#include <iostream>


//Implementation of GMRES for inverting D

c_matrix gmres(const c_matrix& U, const c_matrix& phi, const c_matrix& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message);

//Rotations to transform Hessenberg matrix to upper triangular form
void rotation(c_vector& cn, c_vector& sn, c_matrix& H, const int& j);

c_vector solve_upper_triangular(const c_matrix& A, const c_vector& b, const int& n);

#endif 