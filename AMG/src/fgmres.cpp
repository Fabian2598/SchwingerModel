#include "fgmres.h"


//Solves D psi = phi using GMRES
c_matrix fgmres(const c_matrix& U, const c_matrix& phi, const c_matrix& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message) {
    //GMRES for D^-1 phi
    //phi --> right-hand side
    //x0 --> initial guess  
    //U --> configuration
    //restarts --> number of restarts
    //m --> number of iterations per cycle

    int k = 0; //Iteration number (restart cycle)
    double err = 1;


    c_matrix r(Ntot, c_vector(2, 0));  //r[coordinate][spin] residual
    c_matrix r0(Ntot, c_vector(2, 0));

    //VmT[column vector index][vector arrange in matrix form]
    std::vector<c_matrix> VmT(m+1, c_matrix(Ntot, c_vector(2, 0))); //V matrix transpose-->dimensions exchanged
    std::vector<c_matrix> ZmT(m, c_matrix(Ntot, c_vector(2, 0)));  //Z matrix transpose-->dimensions exchanged

    c_matrix Hm(m+1 , c_vector(m, 0)); //H matrix (Hessenberg matrix)
    c_vector gm(m + 1, 0); 

    //Elements of the rotation matrix |sn[i]|^2 + |cn[i]|^2 = 1
    c_vector sn(m, 0);
    c_vector cn(m, 0);
    c_vector eta(m, 0);


    c_matrix w(Ntot, c_vector(2, 0)); //D*d
    c_matrix x = x0; //initial solution
    c_double beta;


    r0 = phi - D_phi(U, x, m0); //r = b - A*x
	double norm_phi = sqrt(std::real(dot(phi, phi))); //norm of the right hand side
    while (k < restarts) {
        beta = sqrt(std::real(dot(r0, r0))) + 0.0 * I_number;
        VmT[0] = 1.0 / beta * r0;
        gm[0] = beta; //gm[0] = ||r||
        //-----Arnoldi process to build the Krylov basis and the Hessenberg matrix-----//
        for (int j = 0; j < m; j++) {
            SAP(U, VmT[j], ZmT[j], m0, 1); //Preconditioner times VmT; // zm = M^-1. vm
            w = D_phi(U, ZmT[j], m0); //w = D v_j
  
            for (int i = 0; i <= j; i++) {
                Hm[i][j] = dot(w, VmT[i]); //  (v_i^dagger, w)
                w = w -  Hm[i][j] * VmT[i];
            }
            //i.e. the Gramm Schmidt part is highly inefficient. 
            //Could a better implementation of the dot product improve the execution time?
            Hm[j + 1][j] = sqrt(std::real(dot(w, w))); //H[j+1][j] = ||A v_j||
            if (std::real(Hm[j + 1][j]) > 0) {
                VmT[j + 1] = 1.0 / Hm[j + 1][j] * w;
            }
            //----Rotate the matrix----//
            rotation(cn, sn, Hm, j);

            //Rotate gm
            gm[j + 1] = -sn[j] * gm[j];
            gm[j] = std::conj(cn[j]) * gm[j];
        }        
        //Solve the upper triangular system//
		eta = solve_upper_triangular(Hm, gm,m);
        
        for (int i = 0; i < 2 * Ntot; i++) {
            int n = i / 2; int mu = i % 2;
            for (int j = 0; j < m; j++) {
                x[n][mu] = x[n][mu] + eta[j] * ZmT[j][n][mu]; 
            }
        }
        //Compute the residual
        r = phi - D_phi(U, x, m0);
        err = sqrt(std::real(dot(r, r)));
        // if (print_message == true) {
        //     std::cout << "GMRES for D " << k + 1 << " restart cycle" << " Error " << err << std::endl;
        // }

         if (err < tol* norm_phi) {
			 it_count = k + 1;
             if (print_message == true) {
                 std::cout << "FGMRES converged in " << k + 1 << " iterations" << " Error " << err << std::endl;
             }
             return x;
         }
         r0 = r;
         k++;
    }
    if (print_message == true) {
        std::cout << "FGMRES for D did not converge in " << restarts << " restarts" << " Error " << err << std::endl;
    }
    return x;
}


//Solves D psi = phi using FGMRES
c_matrix fgmresParallel(const c_matrix& U, const c_matrix& phi, const c_matrix& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message) {
    //GMRES for D^-1 phi
    //phi --> right-hand side
    //x0 --> initial guess  
    //U --> configuration
    //restarts --> number of restarts
    //m --> number of iterations per cycle

    int k = 0; //Iteration number (restart cycle)
    double err = 1;


    c_matrix r(Ntot, c_vector(2, 0));  //r[coordinate][spin] residual
    c_matrix r0(Ntot, c_vector(2, 0));

    //VmT[column vector index][vector arrange in matrix form]
    std::vector<c_matrix> VmT(m+1, c_matrix(Ntot, c_vector(2, 0))); //V matrix transpose-->dimensions exchanged
    std::vector<c_matrix> ZmT(m, c_matrix(Ntot, c_vector(2, 0)));  //Z matrix transpose-->dimensions exchanged

    c_matrix Hm(m+1 , c_vector(m, 0)); //H matrix (Hessenberg matrix)
    c_vector gm(m + 1, 0); 

    //Elements of the rotation matrix |sn[i]|^2 + |cn[i]|^2 = 1
    c_vector sn(m, 0);
    c_vector cn(m, 0);
    c_vector eta(m, 0);


    c_matrix w(Ntot, c_vector(2, 0)); //D*d
    c_matrix x = x0; //initial solution
    c_double beta;


    r0 = phi - D_phi(U, x, m0); //r = b - A*x
	double norm_phi = sqrt(std::real(dot(phi, phi))); //norm of the right hand side
    int blocks_per_proc = 1;
    while (k < restarts) {
        beta = sqrt(std::real(dot(r0, r0))) + 0.0 * I_number;
        VmT[0] = 1.0 / beta * r0;
        gm[0] = beta; //gm[0] = ||r||
        //-----Arnoldi process to build the Krylov basis and the Hessenberg matrix-----//
        for (int j = 0; j < m; j++) {
            //-----Preconditioner-----//
            //zm = M^-1. vm
            //MPI_Barrier(MPI_COMM_WORLD);
            SAP_parallel(U, VmT[j], ZmT[j], m0, 1,blocks_per_proc); 

            //------------------------//
            

            w = D_phi(U, ZmT[j], m0); //w = D v_j
  
            for (int i = 0; i <= j; i++) {
                Hm[i][j] = dot(w, VmT[i]); //  (v_i^dagger, w)
                w = w -  Hm[i][j] * VmT[i];
            }
            //i.e. the Gramm Schmidt part is highly inefficient. 
            //Could a better implementation of the dot product improve the execution time?
            Hm[j + 1][j] = sqrt(std::real(dot(w, w))); //H[j+1][j] = ||A v_j||
            if (std::real(Hm[j + 1][j]) > 0) {
                VmT[j + 1] = 1.0 / Hm[j + 1][j] * w;
            }
            //----Rotate the matrix----//
            rotation(cn, sn, Hm, j);

            //Rotate gm
            gm[j + 1] = -sn[j] * gm[j];
            gm[j] = std::conj(cn[j]) * gm[j];
        }        
        //Solve the upper triangular system//
		eta = solve_upper_triangular(Hm, gm,m);
        
        for (int i = 0; i < 2 * Ntot; i++) {
            int n = i / 2; int mu = i % 2;
            for (int j = 0; j < m; j++) {
                x[n][mu] = x[n][mu] + eta[j] * ZmT[j][n][mu]; 
            }
        }
        //Compute the residual
        r = phi - D_phi(U, x, m0);
        err = sqrt(std::real(dot(r, r)));
        // if (print_message == true) {
        //     std::cout << "GMRES for D " << k + 1 << " restart cycle" << " Error " << err << std::endl;
        // }

         if (err < tol* norm_phi) {
			 it_count = k + 1;
             if (print_message == true) {
                 std::cout << "FGMRES converged in " << k + 1 << " iterations" << " Error " << err << std::endl;
             }
             return x;
         }
         r0 = r;
         k++;
    }
    if (print_message == true) {
        std::cout << "FGMRES for D did not converge in " << restarts << " restarts" << " Error " << err << std::endl;
    }
    return x;
}
