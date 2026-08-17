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
            GConf.Conf.mu0[n] = rand1;
            GConf.Conf.mu1[n] = rand2;
            GConfV2.Conf.val[idx(x+1,t+1,0)] = rand1;
            GConfV2.Conf.val[idx(x+1,t+1,1)] = rand2;

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
    D_phi_v2(U_v2, phi_v2, Dphi_v2);
    bool test_passed = true;
    if (mpi::rank2d == 0)
        std::cout << "\n\nTesting implementations for D" << std::endl;
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
    D_dagger_phi_v2(U_v2, phi_v2, Dphi_v2);
    bool test_passed = true;
    if (mpi::rank2d == 0)
        std::cout << "\n\nTesting implementations for D^+" << std::endl;
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
        std::cout << "\n\nTesting implementations for phi^+ dD phi" << std::endl;
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
    D_D_dagger_phi_v2(U_v2, phi_v2, Dphi_v2);
    bool test_passed = true;
    if (mpi::rank2d == 0)
        std::cout << "\n\nTesting implementations for DD^+" << std::endl;
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

void Tests::test_CG(){
    initialize_spinors();
    spinor x(mpi::maxSize); 
    spinor_v2 x_v2(mpi::maxSizeH);

    if (mpi::rank2d == 0)
        std::cout << "Testing conjugate gradient version 1" << std::endl;
    conjugate_gradient(U,phi,x,m0);
    if (mpi::rank2d == 0)
        std::cout << "\n";
    if (mpi::rank2d == 0)
        std::cout << "Testing conjugate gradient version 2" << std::endl;
    conjugate_gradient_v2(U_v2,phi_v2,x_v2);

    
}

void Tests::test_plaquettes(){
    initialize_spinors();

    GConf.Compute_Plaquette01();
    GConfV2.Compute_Plaquette01();

    //Constant pointer to a variable: the pointer can move and change what is pointing to, but variable cannot be modified.
    const c_double* p01 = GConf.Plaquette01;
    const c_double* p01_v2 = GConfV2.Plaquette01;

    bool test_passed = true;
    if (mpi::rank2d == 0)
        std::cout << "\n\nTesting implementations for plaquettes computation" << std::endl;
    for(int x=0; x<mpi::width_x; x++){
        for(int t = 0; t<mpi::width_t; t++){
            int n = x*mpi::width_t+t;
            int n2 = (x+1)*(mpi::width_t+2) + (t+1);
            if (std::abs(p01[n]-p01_v2[n2]) > 1e-8){
                 std::cout << "Plaquette implementations are different at (x,t)=(" << x << ", " << t << ")\n" << 
                "on rank " << mpi::rank2d << "\n"
                "p01 = " << p01[n] << "  p01_v2 = " << p01_v2[n2] << std::endl; 
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

void Tests::test_staples(){
    initialize_spinors();
    GConf.Compute_Staple();
    GConfV2.Compute_Staple();


    const spinor& staples = GConf.Staples;
    const spinor_v2& staples_v2 = GConfV2.Staples;

    bool test_passed = true;
    if (mpi::rank2d == 0)
        std::cout << "\n\nTesting implementations for staples computation" << std::endl;

    for(int x=0; x<mpi::width_x; x++){
        for(int t = 0; t<mpi::width_t; t++){
            int n = x*mpi::width_t+t;
            if (std::abs(staples.mu0[n] - staples_v2.val[idx(x+1,t+1,0)])  > 1e-8 
                    || std::abs(staples.mu1[n] - staples_v2.val[idx(x+1,t+1,1)])  > 1e-8 
                ){
                std::cout << "Staple implementations are different at (x,t)=(" << x << ", " << t << ")\n" << 
                "on rank " << mpi::rank2d << "\n"
                "staples.mu0 = " << staples.mu0[n] << "  staples_v2.mu0 = " << staples_v2.val[idx(x+1,t+1,0)] << "\n" <<
                "staples.mu1 = " << staples.mu1[n] << "  staples_v2.mu1 = " << staples_v2.val[idx(x+1,t+1,1)] << std::endl; 
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


void Tests::test_gauge_action(){
    initialize_spinors();
    GConf.Compute_Plaquette01();
    GConfV2.Compute_Plaquette01();

    double Sp = GConf.MeasureSp_HMC();
    double Sp_v2 = GConfV2.MeasureSp_HMC();
    bool test_passed = true;
    if (mpi::rank2d == 0)
        std::cout << "\n\nTesting implementations for Sp and gauge action" << std::endl;



    if (std::abs(Sp-Sp_v2) > 1e-8){
        test_passed = false;
        if (mpi::rank2d == 0)
            std::cout << "Sp computation differs, Sp " << Sp << ", Sp_v2 " << Sp_v2 << std::endl; 
    }

    double gA, gAv2;
    for(double beta = 1.0; beta<=6.0; beta++){
        gA = GConf.Compute_gaugeAction(beta);
        gAv2 = GConfV2.Compute_gaugeAction(beta);
        if (std::abs(gA-gAv2) > 1e-8){
            test_passed = false;
            if (mpi::rank2d == 0)
                std::cout << "Gauge action computation differs at beta=" << beta <<  ", gA " << gA << ", gAv2 " << gAv2 << std::endl; 
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

void Tests::test_IO(){
    initialize_spinors();
    if (mpi::rank2d == 0)
        std::cout << "\n\nWriting gauge configuration to disk" << std::endl;    
    SaveConf(GConf,"conf_v1.bin");
    GConfV2.SaveConf("conf_v2.bin");

    if (mpi::rank2d == 0)
        std::cout << "Reading gauge configurations and comparing them" << std::endl; 

    GaugeConf Conf2;
    GaugeConfV2 Conf2_v2;
    Conf2.readBinary("conf_v1.bin");
    Conf2_v2.ReadConf("conf_v2.bin");

    const spinor& u1 = Conf2.Conf;
    const spinor_v2&  u2 = Conf2_v2.Conf;

    bool test_passed = true;
    for(int x=0; x<mpi::width_x; x++){
        for(int t = 0; t<mpi::width_t; t++){
            int n = x*mpi::width_t+t;
            if (std::abs(u1.mu0[n] - u2.val[idx(x+1,t+1,0)])  > 1e-8 
                    || std::abs(u1.mu1[n] - u2.val[idx(x+1,t+1,1)])  > 1e-8 
                ){
                std::cout << "I/O implementations are different at (x,t)=(" << x << ", " << t << ")\n" << 
                "on rank " << mpi::rank2d << "\n"
                "u1.mu0 = " << u1.mu0[n] << "  u2.mu0 = " << u2.val[idx(x+1,t+1,0)] << "\n" <<
                "u1.mu1 = " << u1.mu1[n] << "  u2.mu1 = " << u2.val[idx(x+1,t+1,1)] << std::endl; 
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