#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H
#include <vector>
#include <complex>
#include <iostream>
#include "mpi.h"

/*
    dot product between two spinors of the form psi[ntot][2]
    A.B = sum_i A_i conj(B_i) 
*/
inline c_double dot(const spinor& x, const spinor& y) {

    c_double local_z = 0;
    //reduction over all lattice points and spin components
    for (int n = 0; n < mpi::maxSize; n++) {
        local_z += x.mu0[n] * std::conj(y.mu0[n]);
        local_z += x.mu1[n] * std::conj(y.mu1[n]);
    }
    c_double z;
    MPI_Allreduce(&local_z, &z, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi::cart_comm);
    return z;
}

#endif 
