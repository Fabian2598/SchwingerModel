#include "conjugate_gradient.h"

//Conjugate gradient for computing (DD^dagger)^-1 phi, where phi is a vector represented by a matrix
//phi[Ntot][2]
c_matrix conjugate_gradient(const c_matrix& U, const c_matrix& phi, const double& m0) {
    int max_iter = 10000;
    double tol = 1e-10; //maybe I lower the tolerance later
    int k = 0; //Iteration number
    double err = 1;

    //D_D_dagger_phi(U, phi, m0); //DD^dagger  
    c_matrix r(Ntot, c_vector(2, 0));  //r[coordinate][spin] residual
    c_matrix d(Ntot, c_vector(2, 0)); //search direction
    c_matrix Ad(Ntot, c_vector(2, 0)); //DD^dagger*d
    c_matrix x(Ntot, c_vector(2, 0)); //solution
    c_double alpha, beta;
	x = phi; //initial solution
    r = phi - D_D_dagger_phi(U, x, m0); //The initial solution can be a vector with zeros 
    d = r; //initial search direction
    c_double r_norm2 = dot(r, r);
    while (k<max_iter && err>tol) {
        Ad = D_D_dagger_phi(U, d, m0); //DD^dagger*d 
        alpha = r_norm2 / dot(d, Ad); //alpha = (r_i,r_i)/(d_i,Ad_i)
        x = x + alpha * d; //x_{i+1} = x_i + alpha*d_i
        r = r - alpha * Ad; //r_{i+1} = r_i - alpha*Ad_i
        err = std::real(dot(r, r)); //err = (r_{i+1},r_{i+1})
        if (err < tol) {
            //std::cout << "Converged in " << k << " iterations" << " Error " << err << std::endl;
            return x;
        }
        beta = err / r_norm2; //beta = (r_{i+1},r_{i+1})/(r_i,r_i)
        d = r + beta * d; //d_{i+1} = r_{i+1} + beta*d_i        
        r_norm2 = err;
        k++;
    }
    std::cout << "Did not converge in " << max_iter << " iterations" << " Error " << err << std::endl;
    return x;
}

//Bi_GCR for D^-1 phi
//phi --> right-hand side
//x0 --> initial guess  
//U --> configuration
c_matrix bi_cgstab(const c_matrix& U, const c_matrix& phi, const c_matrix& x0, const double& m0, const int& max_iter, const double& tol, const bool& print_message) {
    int k = 0; //Iteration number
    double err = 1;

    
    //D_D_dagger_phi(U, phi, m0); //DD^dagger  
    c_matrix r(Ntot, c_vector(2, 0));  //r[coordinate][spin] residual
    c_matrix r_tilde(Ntot, c_vector(2, 0));  //r[coordinate][spin] residual
    c_matrix d(Ntot, c_vector(2, 0)); //search direction
    c_matrix s(Ntot, c_vector(2, 0));
    c_matrix t(Ntot, c_vector(2, 0));
    c_matrix Ad(Ntot, c_vector(2, 0)); //DD^dagger*d
    c_matrix x(Ntot, c_vector(2, 0)); //solution
    c_double alpha, beta, rho_i, omega, rho_i_2;;
    x = x0; //initial solution
    r = phi - D_phi(U, x, m0); //r = b - A*x
    r_tilde = r;
    while (k<max_iter && err>tol) {
        rho_i = dot(r, r_tilde); //r . r_dagger
        if (k == 0) {
            d = r; //d_1 = r_0
        }
        else {
            beta = alpha * rho_i / (omega * rho_i_2); //beta_{i-1} = alpha_{i-1} * rho_{i-1} / (omega_{i-1} * rho_{i-2})
            d = r + beta * (d - omega * Ad); //d_i = r_{i-1} + beta_{i-1} * (d_{i-1} - omega_{i-1} * Ad_{i-1})
        }
        Ad = D_phi(U, d, m0);  //A d_i 
        alpha = rho_i / dot(Ad, r_tilde); //alpha_i = rho_{i-1} / (Ad_i, r_tilde)
        s = r - alpha * Ad; //s = r_{i-1} - alpha_i * Ad_i
        err = std::real(dot(s, s));
        if (err < tol) {
            x = x + alpha * d;
            std::cout << "Converged in " << k << " iterations" << " Error " << err << std::endl;
            return x;
        }
        t = D_phi(U, s, m0);   //A s
        omega = dot(s, t) / dot(t, t); //omega_i = t^dagg . s / t^dagg . t
        r = s - omega * t; //r_i = s - omega_i * t
        x = x + alpha * d + omega * s; //x_i = x_{i-1} + alpha_i * d_i + omega_i * s
        rho_i_2 = rho_i; //rho_{i-2} = rho_{i-1}
        k++;
    }
    if (print_message == true) {
        std::cout << "Did not converge in " << max_iter << " iterations" << " Error " << err << std::endl;
    }
    return x;
}
