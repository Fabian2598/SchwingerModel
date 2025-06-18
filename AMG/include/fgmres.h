#ifndef FGMRES_H
#define FGMRES_H

#include "gmres.h"
#include "sap.h"
#include "amg.h"

//FGMRES with AMG  as preconditioner
c_matrix fgmresAMG(const c_matrix& U, const spinor& phi, const spinor& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message);

//FGMRES with SAP as preconditioner
c_matrix fgmres(const c_matrix& U, const spinor& phi, const spinor& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message);

c_matrix fgmresParallel(const c_matrix& U, const spinor& phi, const spinor& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message);

#endif