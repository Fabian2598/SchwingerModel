#ifndef BI_CGSTAB_H
#define BI_CGSTAB_H

#include "dirac_operator.h"
#include "operator_overloads.h"
#include "variables.h"
#include <cmath>
#include <iostream>

//This function is just for comparison with the AMG method.
spinor bi_cgstab(c_matrix (*func)(const c_matrix&, const c_matrix&, const double&), const int& dim1, const int& dim2,
    const c_matrix& U, const c_matrix& phi, const c_matrix& x0, const double& m0, const int& max_iter, const double& tol, const bool& print_message);

#endif
