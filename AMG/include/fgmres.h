#ifndef FGMRES_H
#define FGMRES_H

#include "amg.h"

/*
	FGMRES with AMG as preconditioner
    U: gauge configuration
    phi: right hand side
    x0: initial guess
    m0: mass parameter
    m: restart length
    restarts: number of restarts of length m
    tol: tolerance for the solver
    print_message: if true, print the convergence message

    The convergence criterion is ||r|| < ||phi|| * tol
*/
spinor fgmresAMG(const c_matrix& U, const spinor& phi, const spinor& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message);

/*
	FGMRES with SAP as preconditioner
    Parameters same as above.
*/
spinor fgmresSAP(const c_matrix& U, const spinor& phi, const spinor& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message);

#endif