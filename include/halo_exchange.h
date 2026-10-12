#ifndef HALO_EXCHANGE_H
#define HALO_EXCHANGE_H
#include "mpi_setup.h"
#include "utils.h"

//The pointers are field.val, i.e. GPU memory when compiled with CUDA

//Halo exchange of a spinor (rows and columns, the corners are not needed by the Dirac operator)
//Phi has dimension [2*(width_x+2)*(width_t+2)]
void exchange_halo(c_double* phi);

//Halo exchange for a vector (not spinor), corners included
void exchange_halo_vec(c_double* v);

//Halo exchange of the gauge links U_0(n) = U[2n], U_1(n) = U[2n+1], corners included.
//The staples and the clover term need the links at the diagonal neighbors (x+-1, t-+1).
void exchange_halo_gauge(c_double* U);

#endif
