#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H
#include <vector>
#include <complex>
#include <iostream>

/*
    dot product between two spinors of the form psi[ntot][2]
    A.B = sum_i A_i conj(B_i) 
*/
inline c_double dot(const spinor& x, const spinor& y) {
    c_double z = 0;
    for (int n = 0; n < LV::Ntot; n++) {
        z += x.mu0[n] * std::conj(y.mu0[n]);
        z += x.mu1[n] * std::conj(y.mu1[n]);
    }
    return z;
}

/*
    "Classic" dot product between two spinors 
    A.B = sum_i sum_j A_ij B_ij (not conjugate)
*/ 
/*
inline c_double dot_v2(const spinor& x, const spinor& y) {
    c_double z = 0;
    int size1 = x.size(), size2 = x[0].size();
    for (int i = 0; i < size1; i++) {
        for (int j = 0; j < size2; j++) {
            z += x[i][j] * y[i][j];
        }
    }
    return z;
}
*/

#endif 
