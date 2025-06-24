#ifndef FGMRES_H
#define FGMRES_H

#include "amg.h"

//FGMRES with AMG  as preconditioner
spinor fgmresAMG(const c_matrix& U, const spinor& phi, const spinor& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message);

//FGMRES with SAP as preconditioner
spinor fgmresSAP(const c_matrix& U, const spinor& phi, const spinor& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message);

#endif