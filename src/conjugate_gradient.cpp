#include "conjugate_gradient.h"

//Conjugate gradient for computing (DD^dagger)^-1 phi
int conjugate_gradient(const spinor& U, const spinor& phi, spinor& x,const double& m0) {
    using namespace mpi;
    int k = 0; //Iteration number
    double err;
    double err_sqr;

    spinor r(mpi::maxSize);  //r[coordinate][spin] residual
    spinor d(mpi::maxSize); //search direction
    spinor Ad(mpi::maxSize); //DD^dagger*d
  
    c_double alpha, beta;

	x = phi;
    D_D_dagger_phi(U, x, Ad, m0); //DD^dagger*x
    
    for(int n = 0; n<maxSize; n++){
        r.mu0[n] = phi.mu0[n] - Ad.mu0[n];
        r.mu1[n] = phi.mu1[n] - Ad.mu1[n];
    }

    d = r; //initial search direction
 
    c_double r_norm2 = dot(r, r);
    
    
    double phi_norm2 = sqrt(std::real(dot(phi, phi)));

    while (k<CG::max_iter) {
        D_D_dagger_phi(U, d,Ad, m0); //DD^dagger*d 
        alpha = r_norm2 / dot(d, Ad); //alpha = (r_i,r_i)/(d_i,Ad_i)

        for(int n = 0; n<maxSize; n++){
        //x = x + alpha * d; //x_{i+1} = x_i + alpha*d_i 
            x.mu0[n] += alpha*d.mu0[n];
            x.mu1[n] += alpha*d.mu1[n];
        //r = r - alpha * Ad; //r_{i+1} = r_i - alpha*Ad_i
            r.mu0[n] -= alpha*Ad.mu0[n];
            r.mu1[n] -= alpha*Ad.mu1[n];
        }
        
        err_sqr = std::real(dot(r, r)); //err_sqr = (r_{i+1},r_{i+1})
		err = sqrt(err_sqr); // err = sqrt(err_sqr)
        if (err < CG::tol*phi_norm2) {
            //if (mpi::rank == 0)
            //    std::cout << "Converged in " << k << " iterations" << " Error " << err << std::endl;
            return 1;
        }

        beta = err_sqr / r_norm2; //beta = (r_{i+1},r_{i+1})/(r_i,r_i)

        //d = r + beta * d; //d_{i+1} = r_{i+1} + beta*d_i 
        for(int n = 0; n<maxSize; n++){
            d.mu0[n] *= beta; 
            d.mu1[n] *= beta;
            d.mu0[n] += r.mu0[n];
            d.mu1[n] += r.mu1[n];
        }
               
        r_norm2 = err_sqr;
        k++;
    }
    std::cout << "CG did not converge in " << CG::max_iter << " iterations" << " Error " << err << std::endl;
    return 0;
}

// D x = phi
spinor bi_cgstab(const spinor& U, const spinor& phi, const spinor& x0, const double& m0, const int& max_iter, const double& tol, const bool& print_message) {
    //Bi_GCR for D^-1 phi
    //phi --> right-hand side
    //x0 --> initial guess  
    //U --> configuration
    int k = 0; //Iteration number
    double err;

    spinor r(mpi::maxSize);//r[coordinate][spin] residual
    spinor r_tilde(mpi::maxSize); //r[coordinate][spin] residual
    spinor d(mpi::maxSize);; //search direction
    spinor s(mpi::maxSize);;
    spinor t(mpi::maxSize);;
    spinor Ad(mpi::maxSize);; //D*d
    spinor x(mpi::maxSize);;//solution
    c_double alpha, beta, rho_i, omega, rho_i_2;
    x = x0; //initial solution
    D_phi(U, x, TEMP, m0);
    for(int n = 0; n < mpi::maxSize; n++) {
        r.mu0[n] = phi.mu0[n] - TEMP.mu0[n]; //r = b - A*x
        r.mu1[n] = phi.mu1[n] - TEMP.mu1[n];
    }
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
            for(int n = 0; n < mpi::maxSize; n++) {
                d.mu0[n] = r.mu0[n] + beta * (d.mu0[n] - omega * Ad.mu0[n]); //d_i = r_{i-1} + beta_{i-1} * (d_{i-1} - omega_{i-1} * Ad_{i-1})
                d.mu1[n] = r.mu1[n] + beta * (d.mu1[n] - omega * Ad.mu1[n]);
            }
        }
        D_phi(U, d, Ad,m0);  //A d_i 
        alpha = rho_i / dot(Ad, r_tilde); //alpha_i = rho_{i-1} / (Ad_i, r_tilde)
        //s = r - alpha * Ad; //s = r_{i-1} - alpha_i * Ad_i
        for(int n = 0; n < mpi::maxSize; n++) {
            s.mu0[n] = r.mu0[n] - alpha * Ad.mu0[n]; //s_i = r_{i-1} - alpha_i * Ad_i
            s.mu1[n] = r.mu1[n] - alpha * Ad.mu1[n];
        }
        err = sqrt(std::real(dot(s, s)));
        if (err < tol * norm_phi) {
            for(int n = 0; n <mpi::maxSize; n++) {
                x.mu0[n] = x.mu0[n] + alpha * d.mu0[n];
                x.mu1[n] = x.mu1[n] + alpha * d.mu1[n];
            }
            if (print_message == true) {
                std::cout << "Bi-CG-stab for D converged in " << k+1 << " iterations" << " Error " << err << std::endl;
            }
            return x;
        }
        D_phi(U, s,t, m0);   //A s
        omega = dot(s, t) / dot(t, t); //omega_i = t^dagg . s / t^dagg . t
        for(int n = 0; n < mpi::maxSize; n++) {
            r.mu0[n] = s.mu0[n] - omega * t.mu0[n]; //r_i = s - omega_i * t
            r.mu1[n] = s.mu1[n] - omega * t.mu1[n];
    
        }

         for(int n = 0; n < mpi::maxSize; n++) {
            x.mu0[n] = x.mu0[n] + alpha * d.mu0[n] + omega * s.mu0[n]; //x_i = x_{i-1} + alpha_i * d_i + omega_i * s_i
            x.mu1[n] = x.mu1[n] + alpha * d.mu1[n] + omega * s.mu1[n];
        }
        rho_i_2 = rho_i; //rho_{i-2} = rho_{i-1}
        k++;
    }
    std::cout << "Bi-CG-stab for D did not converge in " << max_iter << " iterations" << " Error " << err << std::endl;
    return x;
}
