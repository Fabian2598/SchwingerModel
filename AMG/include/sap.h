//Schwarz alternating method
#ifndef SAP_H
#define SAP_H
#include "gmres.h"

//Build the blocks for the Schwarz alternating method
void SchwarzBlocks();
bool checkSize(const c_matrix& v, const int& N1, const int& N2); //Check vector dimensions

//Schwarz blocks have to be initialized first before calling the I_B operators

// I_B^T v --> Restriction of the vector v to the block B
c_matrix It_B_v(const c_matrix& v, const int& block);
//I_B v --> Interpolation of the vector v to the original lattice
c_matrix I_B_v(const c_matrix& v, const int& block);
//Restriction of the Dirac operator to the block B
c_matrix D_B(const c_matrix& U, const c_matrix& phi, const double& m0,const int& block);

//Solves D_B x = phi using GMRES, where D_B is the Dirac matrix restricted to the Schwarz block B 
c_matrix gmres_D_B(const c_matrix& U, const c_matrix& phi, const c_matrix& x0, const double& m0, 
    const int& m, const int& restarts, const double& tol, const int& block,const bool& print_message);

//I_B * D_B^-1 * I_B^T v --> Extrapolation of D_B^-1 to the original lattice.
c_matrix I_D_B_1_It(const c_matrix& U, const c_matrix& v, const double& m0,const int& block);


#endif