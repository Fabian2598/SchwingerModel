#include "conjugate_gradient.h"

//Conjugate gradient for computing (DD^dagger)^-1 phi, where phi is a vector represented by a matrix
//phi[Ntot][2]
int conjugate_gradient(const c_matrix& U, const spinor& phi, spinor& x,const double& m0) {
    //int max_iter = 10000;
    //double tol = 1e-10; //maybe I lower the tolerance later
    int k = 0; //Iteration number
    double err;
    double err_sqr;

    //D_D_dagger_phi(U, phi, m0); //DD^dagger  
    spinor r(LV::Ntot, c_vector(2, 0));  //r[coordinate][spin] residual
    spinor d(LV::Ntot, c_vector(2, 0)); //search direction
    spinor Ad(LV::Ntot, c_vector(2, 0)); //DD^dagger*d
    //spinor x(LV::Ntot, c_vector(2, 0)); //solution
    c_double alpha, beta;

	x = phi;
    D_D_dagger_phi(U, x, Ad, m0); //DD^dagger*x
    for(int n = 0; n<LV::Ntot; n++){
        r[n][0] = phi[n][0] - Ad[n][0];
        r[n][1] = phi[n][1] - Ad[n][1];
    }

    d = r; //initial search direction
    c_double r_norm2 = dot(r, r);
    double phi_norm2 = sqrt(std::real(dot(phi, phi)));

    while (k<CG::max_iter) {
        D_D_dagger_phi(U, d,Ad, m0); //DD^dagger*d 
        alpha = r_norm2 / dot(d, Ad); //alpha = (r_i,r_i)/(d_i,Ad_i)

        //x = x + alpha * d; //x_{i+1} = x_i + alpha*d_i
        for(int n = 0; n<LV::Ntot; n++){ 
            x[n][0] += alpha*d[n][0];
            x[n][1] += alpha*d[n][1];
        }
        
        //r = r - alpha * Ad; //r_{i+1} = r_i - alpha*Ad_i
        for(int n = 0; n<LV::Ntot; n++){
            r[n][0] -= alpha*Ad[n][0];
            r[n][1] -= alpha*Ad[n][1];
        }
        
        err_sqr = std::real(dot(r, r)); //err_sqr = (r_{i+1},r_{i+1})
		err = sqrt(err_sqr); // err = sqrt(err_sqr)
        if (err < CG::tol*phi_norm2) {
            //std::cout << "Converged in " << k << " iterations" << " Error " << err << std::endl;
            return 1;
        }

        beta = err_sqr / r_norm2; //beta = (r_{i+1},r_{i+1})/(r_i,r_i)

        //d = r + beta * d; //d_{i+1} = r_{i+1} + beta*d_i 
        for(int n = 0; n<LV::Ntot; n++){
            d[n][0] *= beta; 
            d[n][1] *= beta;
            d[n][0] += r[n][0];
            d[n][1] += r[n][1];
        }
               
        r_norm2 = err_sqr;
        k++;
    }
    std::cout << "CG did not converge in " << CG::max_iter << " iterations" << " Error " << err << std::endl;
    return 0;
}

// D x = phi
spinor bi_cgstab(const c_matrix& U, const spinor& phi, const spinor& x0, const double& m0, const int& max_iter, const double& tol, const bool& print_message) {
    //Bi_GCR for D^-1 phi
    //phi --> right-hand side
    //x0 --> initial guess  
    //U --> configuration
    int k = 0; //Iteration number
    double err;

    
    //D_D_dagger_phi(U, phi, m0); //DD^dagger  
    spinor r(LV::Ntot, c_vector(2, 0));  //r[coordinate][spin] residual
    spinor r_tilde(LV::Ntot, c_vector(2, 0));  //r[coordinate][spin] residual
    spinor d(LV::Ntot, c_vector(2, 0)); //search direction
    spinor s(LV::Ntot, c_vector(2, 0));
    spinor t(LV::Ntot, c_vector(2, 0));
    spinor Ad(LV::Ntot, c_vector(2, 0)); //D*d
    spinor x(LV::Ntot, c_vector(2, 0)); //solution
    c_double alpha, beta, rho_i, omega, rho_i_2;
    x = x0; //initial solution
    D_phi(U, x, TEMP, m0);
    r = phi - TEMP; //r = b - A*x
    r = 
    r_tilde = r;
	double norm_phi = sqrt(std::real(dot(phi, phi))); //norm of the right hand side
    while (k<max_iter) {
        rho_i = dot(r, r_tilde); //r . r_dagger
        if (k == 0) {
            d = r; //d_1 = r_0
        }
        else {
            beta = alpha * rho_i / (omega * rho_i_2); //beta_{i-1} = alpha_{i-1} * rho_{i-1} / (omega_{i-1} * rho_{i-2})
            //d = r + beta * (d - omega * Ad); //d_i = r_{i-1} + beta_{i-1} * (d_{i-1} - omega_{i-1} * Ad_{i-1})
            for(int i = 0; i < LV::Ntot; i++) {
                for(int j = 0; j < 2; j++) {
                    d[i][j] = r[i][j] + beta * (d[i][j] - omega * Ad[i][j]); //d_i = r_{i-1} + beta_{i-1} * (d_{i-1} - omega_{i-1} * Ad_{i-1})
                }
            }
        }
        D_phi(U, d, Ad,m0);  //A d_i 
        alpha = rho_i / dot(Ad, r_tilde); //alpha_i = rho_{i-1} / (Ad_i, r_tilde)
        //s = r - alpha * Ad; //s = r_{i-1} - alpha_i * Ad_i
        for(int i = 0; i < LV::Ntot; i++) {
            for (int j = 0; j < 2; j++) {
                s[i][j] = r[i][j] - alpha * Ad[i][j]; //s_i = r_{i-1} - alpha_i * Ad_i
            }
        }
        err = sqrt(std::real(dot(s, s)));
        if (err < tol * norm_phi) {
            x = x + alpha * d;
            if (print_message == true) {
                //std::cout << "Bi-CG-stab for D converged in " << k+1 << " iterations" << " Error " << err << std::endl;
            }
            return x;
        }
        D_phi(U, s,t, m0);   //A s
        omega = dot(s, t) / dot(t, t); //omega_i = t^dagg . s / t^dagg . t
        r = s - omega * t; //r_i = s - omega_i * t
        for(int i = 0; i < LV::Ntot; i++) {
            for (int j = 0; j < 2; j++) {
                r[i][j] = s[i][j] - omega * t[i][j]; //x_i = x_{i-1} + alpha_i * d_i + omega_i * s_i
            }
        }

        //x = x + alpha * d + omega * s; //x_i = x_{i-1} + alpha_i * d_i + omega_i * s
         for(int i = 0; i < LV::Ntot; i++) {
            for (int j = 0; j < 2; j++) {
                x[i][j] = x[i][j] + alpha * d[i][j] + omega * s[i][j]; //x_i = x_{i-1} + alpha_i * d_i + omega_i * s_i
            }
        }
        rho_i_2 = rho_i; //rho_{i-2} = rho_{i-1}
        k++;
    }
    std::cout << "Bi-CG-stab for D did not converge in " << max_iter << " iterations" << " Error " << err << std::endl;
    return x;
}
