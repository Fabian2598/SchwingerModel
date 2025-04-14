#include "conjugate_gradient.h"

//Conjugate gradient for computing (DD^dagger)^-1 phi, where phi is a vector represented by a matrix
//phi[Ntot][2]
c_matrix conjugate_gradient(const c_matrix& U, const c_matrix& phi, const double& m0) {
    int max_iter = 10000;
    double tol = 1e-10; //maybe I lower the tolerance later
    int k = 0; //Iteration number
    double err = 1;
    double err_sqr;

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
	double norm_phi = sqrt(std::real(dot(phi, phi))); //norm of the right hand side
    while (k<max_iter && err>tol * norm_phi) {
        Ad = D_D_dagger_phi(U, d, m0); //DD^dagger*d 
        alpha = r_norm2 / dot(d, Ad); //alpha = (r_i,r_i)/(d_i,Ad_i)
        x = x + alpha * d; //x_{i+1} = x_i + alpha*d_i
        r = r - alpha * Ad; //r_{i+1} = r_i - alpha*Ad_i
        err_sqr = std::real(dot(r, r)); //err_sqrt = (r_{i+1},r_{i+1})
		//err = std::sqrt(err_sqr);
        err = err_sqr;
        if (err < tol * norm_phi) {
            //std::cout << "Converged in " << k << " iterations" << " Error " << err << std::endl;
            return x;
        }
        beta = err_sqr / r_norm2; //beta = (r_{i+1},r_{i+1})/(r_i,r_i)
        d = r + beta * d; //d_{i+1} = r_{i+1} + beta*d_i        
        r_norm2 = err_sqr;
        k++;
    }
    std::cout << "Did not converge in " << max_iter << " iterations" << " Error " << err << std::endl;
    return x;
}

// D x = phi
c_matrix bi_cgstab(const c_matrix& U, const c_matrix& phi, const c_matrix& x0, const double& m0, const int& max_iter, const double& tol, const bool& print_message) {
    //Bi_GCR for D^-1 phi
    //phi --> right-hand side
    //x0 --> initial guess  
    //U --> configuration
    int k = 0; //Iteration number
    double err = 1;

    
    //D_D_dagger_phi(U, phi, m0); //DD^dagger  
    c_matrix r(Ntot, c_vector(2, 0));  //r[coordinate][spin] residual
    c_matrix r_tilde(Ntot, c_vector(2, 0));  //r[coordinate][spin] residual
    c_matrix d(Ntot, c_vector(2, 0)); //search direction
    c_matrix s(Ntot, c_vector(2, 0));
    c_matrix t(Ntot, c_vector(2, 0));
    c_matrix Ad(Ntot, c_vector(2, 0)); //D*d
    c_matrix x(Ntot, c_vector(2, 0)); //solution
    c_double alpha, beta, rho_i, omega, rho_i_2;;
    x = x0; //initial solution
    r = phi - D_phi(U, x, m0); //r = b - A*x
    r_tilde = r;
	double norm_phi = sqrt(std::real(dot(phi, phi))); //norm of the right hand side
    while (k<max_iter && err>tol* norm_phi) {
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
        err = sqrt(std::real(dot(s, s)));
        if (print_message == true) {
            std::cout << "Bi-CG-stab for D " << k + 1 << " iteration" << " Error " << err << std::endl;
        }
        if (err < tol * norm_phi) {
            x = x + alpha * d;
            it_count = k+1;
            //if (print_message == true) {
                std::cout << "Bi-CG-stab for D converged in " << k+1 << " iterations" << " Error " << err << std::endl;
            //}
            return x;
        }
        t = D_phi(U, s, m0);   //A s
        omega = dot(s, t) / dot(t, t); //omega_i = t^dagg . s / t^dagg . t
        r = s - omega * t; //r_i = s - omega_i * t
        x = x + alpha * d + omega * s; //x_i = x_{i-1} + alpha_i * d_i + omega_i * s
        rho_i_2 = rho_i; //rho_{i-2} = rho_{i-1}
        k++;
    }
    it_count = max_iter;
    //if (print_message == true) {
        std::cout << "Bi-CG-stab for D did not converge in " << max_iter << " iterations" << " Error " << err << std::endl;
    //}
    return x;
}

c_matrix canonical_vector(const int& i, const int& N1, const int& N2) {
	c_matrix e_i(N1, c_vector (N2,0.0));
	int j = i / N2;
	int k = i % N2;
	e_i[j][k] = 1.0;
	return e_i;
}


c_matrix kaczmarz(const c_matrix& U, const c_matrix& phi, const c_matrix& x0, const double& m0, const int& max_iter, const double& tol, const bool& print_message){
    //Kaczmarz method for D^-1 phi
    //phi --> right-hand side
    //x0 --> initial guess  
    //U --> configuration
    int k = 0; //Iteration number
    double err = 1;
    double err_sqr;

    c_matrix x(Ntot, c_vector(2, 0)); //solution
    c_matrix r(Ntot, c_vector(2, 0));  //r[coordinate][spin] residual
    c_matrix A(Ntot, c_vector(2, 0)); //D x
    c_matrix en(Ntot, c_vector(2, 0)); //canonical vector
    c_double alpha;
    x = x0; //initial solution
    r = phi - D_phi(U, x, m0); //r = b - A*x
    while (k<max_iter && err>tol) {
        for (int i = 0; i < 2*Ntot; i++) {
            en = canonical_vector(i, Ntot, 2); //e_i
            A = D_phi(U, en, m0); //A_i
            alpha = dot(r, A) / dot(A, A); //alpha = (r,A_i)/(A_i,A_i)
            x = x + alpha * en; //x_{i+1} = x_i + alpha*A_i
            r = phi - D_phi(U, x, m0); //r_{i+1} = b - A*x_{i+1}      
        }
        err_sqr = std::real(dot(r, r)); //err_sqrt = (r_{i+1},r_{i+1})
        err = sqrt(err_sqr);
        if (print_message == true) {
            std::cout << "Kaczmarz for D " << k + 1 << " iteration" << " Error " << err << std::endl;
        }
        if (err < tol) {
            it_count = k+1;
            if (print_message == true) {
                std::cout << "Kaczmarz for D converged in " << k+1 << " iterations" << " Error " << err << std::endl;
            }
            return x;
        }
        k++;
    }
    if (print_message == true){
        std::cout << "Kaczmarz for D did not converge in " << max_iter << " iterations" << " Error " << err << std::endl;
    }
    return x;
}
