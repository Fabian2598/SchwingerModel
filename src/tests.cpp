#include "tests.h"

void Tests::initialize_spinors(){
    for(int x=0; x<mpi::width_x; x++){
        for(int t = 0; t<mpi::width_t; t++){
            int n = x*mpi::width_t+t;
            c_double rand1, rand2;

            rand1 = RandomU1();
            rand2 = RandomU1();
            U.mu0[n] = rand1;
            U.mu1[n] = rand2;
            U_v2.val[idx(x+1,t+1,0)] = rand1;
            U_v2.val[idx(x+1,t+1,1)] = rand2;

            rand1 = RandomU1();
            rand2 = RandomU1();
            phi.mu0[n] = rand1;
            phi.mu1[n] = rand2;
            phi_v2.val[idx(x+1,t+1,0)] = rand1;
            phi_v2.val[idx(x+1,t+1,1)] = rand2;

            rand1 = RandomU1();
            rand2 = RandomU1();
            left.mu0[n] = rand1;
            left.mu1[n] = rand2;
            left_v2.val[idx(x+1,t+1,0)] = rand1;
            left_v2.val[idx(x+1,t+1,1)] = rand2;

            rand1 = RandomU1();
            rand2 = RandomU1();
            right.mu0[n] = rand1;
            right.mu1[n] = rand2;
            right_v2.val[idx(x+1,t+1,0)] = rand1;
            right_v2.val[idx(x+1,t+1,1)] = rand2;

        }
    }
}

void Tests::test_D_operator(){   
    initialize_spinors(); 
    spinor Dphi(mpi::maxSize); 
    spinor_v2 Dphi_v2(mpi::maxSizeH);

    D_phi(U, phi, Dphi, m0);
    D_phi_v2(U_v2, phi_v2, Dphi_v2, m0);
    bool test_passed = true;
    if (mpi::rank2d == 0)
        std::cout << "\n\n\nTesting implementations for D" << std::endl;
    for(int x=0; x<mpi::width_x; x++){
        for(int t = 0; t<mpi::width_t; t++){
            int n = x*mpi::width_t+t;
            if (std::abs(Dphi.mu0[n] - Dphi_v2.val[idx(x+1,t+1,0)])  > 1e-8 
                    || std::abs(Dphi.mu1[n] - Dphi_v2.val[idx(x+1,t+1,1)])  > 1e-8 
                ){
                std::cout << "D operator implementations are different at (x,t)=(" << x << ", " << t << ")\n" << 
                "on rank " << mpi::rank2d << "\n"
                "Dphi.mu0 = " << Dphi.mu0[n] << "  Dphi_v2.mu0 = " << Dphi_v2.val[idx(x+1,t+1,0)] << "\n" <<
                "Dphi.mu1 = " << Dphi.mu1[n] << "  Dphi_v2.mu1 = " << Dphi_v2.val[idx(x+1,t+1,1)] << std::endl; 
                test_passed = false;
            }
        }
    }

    if (test_passed){
        if (mpi::rank2d == 0)
            std::cout << "Both implementations coincide" << std::endl;
    }
    else{
        if (mpi::rank2d == 0)
            std::cout << "Implementations do not coincide" << std::endl;
    }

}


void Tests::test_D_dagger_operator(){
    initialize_spinors();
    spinor Dphi(mpi::maxSize); 
    spinor_v2 Dphi_v2(mpi::maxSizeH);

    D_dagger_phi(U, phi, Dphi, m0);
    D_dagger_phi_v2(U_v2, phi_v2, Dphi_v2, m0);
    bool test_passed = true;
    if (mpi::rank2d == 0)
        std::cout << "\n\n\nTesting implementations for D^+" << std::endl;
    for(int x=0; x<mpi::width_x; x++){
        for(int t = 0; t<mpi::width_t; t++){
            int n = x*mpi::width_t+t;
            if (std::abs(Dphi.mu0[n] - Dphi_v2.val[idx(x+1,t+1,0)])  > 1e-8 
                    || std::abs(Dphi.mu1[n] - Dphi_v2.val[idx(x+1,t+1,1)])  > 1e-8 
                ){
                std::cout << "D^+ operator implementations are different at (x,t)=(" << x << ", " << t << ")\n" << 
                "on rank " << mpi::rank2d << "\n"
                "Dphi.mu0 = " << Dphi.mu0[n] << "  Dphi_v2.mu0 = " << Dphi_v2.val[idx(x+1,t+1,0)] << "\n" <<
                "Dphi.mu1 = " << Dphi.mu1[n] << "  Dphi_v2.mu1 = " << Dphi_v2.val[idx(x+1,t+1,1)] << std::endl; 
                test_passed = false;
            }
        }
    }

    if (test_passed){
        if (mpi::rank2d == 0)
            std::cout << "Both implementations coincide" << std::endl;
    }
    else{
        if (mpi::rank2d == 0)
            std::cout << "Implementations do not coincide" << std::endl;
    }
}


void Tests::test_phi_dag_partialD_phi(){
    initialize_spinors();
    re_field Dphi(mpi::maxSize);
    re_field_v2 Dphi_v2(mpi::maxSizeH);


    Dphi = phi_dag_partialD_phi( U, left, right);
    phi_dag_partialD_phi_v2(U_v2, left_v2,right_v2,Dphi_v2);
    bool test_passed = true;
    if (mpi::rank2d == 0)
        std::cout << "\n\n\nTesting implementations for phi^+ dD phi" << std::endl;
    for(int x=0; x<mpi::width_x; x++){
        for(int t = 0; t<mpi::width_t; t++){
            int n = x*mpi::width_t+t;
            if (std::abs(Dphi.mu0[n] - Dphi_v2.val[idx(x+1,t+1,0)])  > 1e-8 
                    || std::abs(Dphi.mu1[n] - Dphi_v2.val[idx(x+1,t+1,1)])  > 1e-8 
                ){
                std::cout << "phi^+ dD phi operator implementations are different at (x,t)=(" << x << ", " << t << ")\n" << 
                "on rank " << mpi::rank2d << "\n"
                "Dphi.mu0 = " << Dphi.mu0[n] << "  Dphi_v2.mu0 = " << Dphi_v2.val[idx(x+1,t+1,0)] << "\n" <<
                "Dphi.mu1 = " << Dphi.mu1[n] << "  Dphi_v2.mu1 = " << Dphi_v2.val[idx(x+1,t+1,1)] << std::endl; 
                test_passed = false;
            }
        }
    }

    if (test_passed){
        if (mpi::rank2d == 0)
            std::cout << "Both implementations coincide" << std::endl;
    }
    else{
        if (mpi::rank2d == 0)
            std::cout << "Implementations do not coincide" << std::endl;
    }
}

void Tests::test_D_D_dagger_phi(){
    initialize_spinors();
    spinor Dphi(mpi::maxSize); 
    spinor_v2 Dphi_v2(mpi::maxSizeH);

    D_D_dagger_phi(U, phi, Dphi, m0);
    D_D_dagger_phi_v2(U_v2, phi_v2, Dphi_v2, m0);
    bool test_passed = true;
    if (mpi::rank2d == 0)
        std::cout << "\n\n\nTesting implementations for DD^+" << std::endl;
    for(int x=0; x<mpi::width_x; x++){
        for(int t = 0; t<mpi::width_t; t++){
            int n = x*mpi::width_t+t;
            if (std::abs(Dphi.mu0[n] - Dphi_v2.val[idx(x+1,t+1,0)])  > 1e-8 
                    || std::abs(Dphi.mu1[n] - Dphi_v2.val[idx(x+1,t+1,1)])  > 1e-8 
                ){
                std::cout << "DD^+ operator implementations are different at (x,t)=(" << x << ", " << t << ")\n" << 
                "on rank " << mpi::rank2d << "\n"
                "Dphi.mu0 = " << Dphi.mu0[n] << "  Dphi_v2.mu0 = " << Dphi_v2.val[idx(x+1,t+1,0)] << "\n" <<
                "Dphi.mu1 = " << Dphi.mu1[n] << "  Dphi_v2.mu1 = " << Dphi_v2.val[idx(x+1,t+1,1)] << std::endl; 
                test_passed = false;
            }
        }
    }

    if (test_passed){
        if (mpi::rank2d == 0)
            std::cout << "Both implementations coincide" << std::endl;
    }
    else{
        if (mpi::rank2d == 0)
            std::cout << "Implementations do not coincide" << std::endl;
    }
}