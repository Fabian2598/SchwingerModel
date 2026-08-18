#ifndef BI_CGSTAB_H
#define BI_CGSTAB_H

#include "dirac_operator.h"

/*
    Bi-CGstab method for inveting D x = phi.
    U: Gauge conf
    phi: right-hand side
    x0: initial solution
    x: solution buffer
*/
int bi_cgstab(const spinor& U, const spinor& phi, const spinor& x0, spinor& x);

#endif