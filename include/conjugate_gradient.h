#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "matrix_operations.h"
#include <iostream>

typedef std::complex<double> c_double; 
typedef std::vector<c_double> c_vector;
typedef std::vector<c_vector> c_matrix;


c_double dot(const c_vector& x,const c_vector& y){
    //Dot product of two vectors
    c_double z = 0;
    for(int i=0;i<x.size();i++){
        z += x[i]*std::conj(y[i]);
    }
    return z;
}

//Overload + - and * operators
template <typename T>
c_vector operator*(const T& v1, const c_vector& v2){
    c_vector y(v2.size(),0);
    for(int i=0;i<v2.size();i++){
        y[i] = v1*v2[i];
    }
    return y;
}

c_vector operator+(const c_vector& v1, const c_vector& v2){
    c_vector y(v1.size(),0);
    for(int i=0;i<v1.size();i++){
        y[i] = v1[i] + v2[i];
    }
    return y;
}

c_vector operator-(const c_vector& v1, const c_vector& v2){
    c_vector y(v1.size(),0);
    for(int i=0;i<v1.size();i++){
        y[i] = v1[i] - v2[i];
    }
    return y;
}

//Conjugate gradient for iverting DD^dagger
c_vector conjugate_gradient(const int& max_iter, const double& tol, const c_matrix& U, const c_matrix& phi, const double& m0){
    //int max_iter = 100;
    //double tol = 1e-10; maybe I lower the tolerance later
    int k = 0; //Iteration number
    double err = 1;
    //c_matrix Dphi(Ntot, c_vector(2, 0)); //Dphi[Ntot][2] 

    //D_D_dagger_phi(U, phi, m0);  
    c_vector r(Ntot,0);
    c_vector d(Ntot,0);
    c_vector Ad(Ntot,0);
    c_vector x(Ntot,0);
    c_double alpha,beta;
    
    //Who is b for this problem?

    r = b - D_D_dagger_phi(U, phi, m0); //residual
    d = r; //initial search direction
    c_double r_norm2;
    while(k<max_iter && err>tol){
        Ad = D_D_dagger_phi(U, d, m0); //c_vector
        r_norm2 = dot(r,r); //c_double

        alpha = r_norm2/dot(d,Ad); //c_double
        x = x + alpha*d; //c_vector
        r = r - alpha*Ad; //c_vector
        beta = dot(r,r)/r_norm2; //c_double
        d = r + beta*d; //c_vector
        err = 0;
        for(int i = 0; i<N; i++){
            err += std::abs(r[i]);
        }
        //std::cout << "Error: " << err << std::endl;
        k++;
    }
    std::cout << "Converged in " << k << "iterations" << " Error " << err << std::endl;
    //std::cout << "Solution: ";
    //for(int i = 0; i<N; i++){
    //    std::cout << x[i] << " ";
   // }
    //std::cout << std::endl;
}

#endif
