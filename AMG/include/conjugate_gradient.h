#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "matrix_operations.h"
#include <iostream>

typedef std::complex<double> c_double;
typedef std::vector<c_double> c_vector;
typedef std::vector<c_vector> c_matrix;

//Complex dot product (not matrix multiplication)
c_double dot(const c_matrix& x, const c_matrix& y);
//Overload + - and * operators
template <typename T>
c_matrix operator*(const T& lambda, const c_matrix& A);
c_matrix operator+(const c_matrix& A, const c_matrix& B);
c_matrix operator-(const c_matrix& A, const c_matrix& B);


//Conjugate gradient for computing (DD^dagger)^-1 phi, where phi is a vector represented by a matrix
//phi[Ntot][2]
c_matrix conjugate_gradient(const c_matrix& U, const c_matrix& phi, const double& m0); 

c_matrix bi_cgstab(const c_matrix& U, const c_matrix& phi, const c_matrix& x0,const double& m0, const int& max_iter, const double& tol, const bool& print_message); //D^-1 phi 


#endif
