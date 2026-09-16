#ifndef HALO_EXCHANGE_H
#define HALO_EXCHANGE_H
#include "mpi_setup.h"
#include "utils.h"


//Halo exchange
//Phi has dimension [2*(width_x+2)*(width_t+2)]
void exchange_halo(c_double* phi);

//Halo exchange for a vector (not spinor)
void exchange_halo_vec(c_double* v);

#endif