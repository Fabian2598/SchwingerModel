#include "conjugate_gradient.h"

//Conjugate gradient for computing (DD^dagger)^-1 phi, where phi is a vector represented by a matrix
//phi[Ntot][2]
c_matrix conjugate_gradient(const c_matrix& U, const c_matrix& phi, const double& m0) {
    int max_iter = 100;
    double tol = 1e-10; //maybe I lower the tolerance later
    int k = 0; //Iteration number
    double err = 1;

    //D_D_dagger_phi(U, phi, m0); //DD^dagger  
    c_matrix r(Ntot, c_vector(2, 0));  //r[coordinate][spin] residual
    c_matrix d(Ntot, c_vector(2, 0)); //search direction
    c_matrix Ad(Ntot, c_vector(2, 0)); //DD^dagger*d
    c_matrix x(Ntot, c_vector(2, 0)); //solution
    c_double alpha, beta;
	x = phi; //initial solution
    r = phi - D_D_dagger_phi(U, x, m0); //The initial solution can be a vector with zeros 
    d = r; //initial search direction
    c_double r_norm2 = dot(r, r);
    while (k<max_iter && err>tol) {
        Ad = D_D_dagger_phi(U, d, m0); //DD^dagger*d 
        alpha = r_norm2 / dot(d, Ad); //alpha = (r_i,r_i)/(d_i,Ad_i)
        x = x + alpha * d; //x_{i+1} = x_i + alpha*d_i
        r = r - alpha * Ad; //r_{i+1} = r_i - alpha*Ad_i
        err = std::real(dot(r, r)); //err = (r_{i+1},r_{i+1})
        if (err < tol) {
            //std::cout << "Converged in " << k << " iterations" << " Error " << err << std::endl;
            return x;
        }
        beta = err / r_norm2; //beta = (r_{i+1},r_{i+1})/(r_i,r_i)
        d = r + beta * d; //d_{i+1} = r_{i+1} + beta*d_i        
        r_norm2 = err;
        k++;
    }
    std::cout << "Did not converge in " << max_iter << " iterations" << " Error " << err << std::endl;
    return x;
}

/*
// OpenMP parallelization
#include <omp.h>
//Overload + - and * operators
template <typename T>
c_matrix operator*(const T& lambda, const c_matrix& A) {
    c_matrix B(A.size(), c_vector(A[0].size(), 0));
#pragma omp parallel for
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            B[i][j] = lambda * A[i][j];
        }
    }
    return B;
}

c_matrix operator+(const c_matrix& A, const c_matrix& B) {
    c_matrix C(A.size(), c_vector(A[0].size(), 0));
#pragma omp parallel for
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    return C;
}

c_matrix operator-(const c_matrix& A, const c_matrix& B) {
    c_matrix C(A.size(), c_vector(A[0].size(), 0));
#pragma omp parallel for
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

*/