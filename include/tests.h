#ifndef TESTS_H
#define TESTS_H

#include "conjugate_gradient.h"
#include "utils.h"


//Testing both implementations of the operations
class Tests {

public:
    Tests(double m0) : m0(m0) {
        U     = spinor(mpi::maxSize);
        phi   = spinor(mpi::maxSize);
        left  = spinor(mpi::maxSize);
        right = spinor(mpi::maxSize);
        
        U_v2    = spinor_v2(mpi::maxSizeH);
        phi_v2  = spinor_v2(mpi::maxSizeH);
        left_v2  = spinor_v2(mpi::maxSizeH);
        right_v2 = spinor_v2(mpi::maxSizeH);
    }
    ~Tests() {} 

    
    void initialize_spinors();
    void test_D_operator();
    void test_D_dagger_operator();
    void test_phi_dag_partialD_phi();
    void test_D_D_dagger_phi();
    void test_CG();
    void test_plaquettes();
    void test_staples();
    void test_gauge_action();
    void test_IO();

private:
    double m0;
    spinor U;
    spinor phi;
    spinor left, right;

    spinor_v2 U_v2;
    spinor_v2 phi_v2;
    spinor_v2 left_v2, right_v2;

    GaugeConf GConf;
    GaugeConfV2 GConfV2;
     
};

#endif