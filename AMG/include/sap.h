//Schwarz alternating method
#ifndef SAP_H
#define SAP_H
#include "variables.h"
#include "operator_overloads.h"
#include "matrix_operations.h"

//Build the blocks for the Schwarz alternating method
void SchwarzBlocks();

// I_B^T v --> Restriction of the vector v to the block B
c_matrix It_B_v(const c_matrix& v);
    //Schwarz blocks have to be initialized first before calling this function


#endif