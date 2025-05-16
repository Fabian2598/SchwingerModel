#ifndef FGMRES_H
#define FGMRES_H

#include "gmres.h"
#include "sap.h"

//Implementation of Flexible-GMRES for inverting D

c_matrix fgmres(const c_matrix& U, const c_matrix& phi, const c_matrix& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message);

c_matrix fgmresParallel(const c_matrix& U, const c_matrix& phi, const c_matrix& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message);

#endif