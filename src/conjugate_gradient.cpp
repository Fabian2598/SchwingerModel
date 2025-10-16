#include "conjugate_gradient.h"

//Conjugate gradient for computing (DD^dagger)^-1 phi
int conjugate_gradient(const spinor& U, const spinor& phi, spinor& x,const double& m0) {
    using namespace mpi;
    int k = 0; //Iteration number
    double err;
    double err_sqr;

    spinor r(maxSize);  //r[coordinate][spin] residual
    spinor d(maxSize); //search direction
    spinor Ad(maxSize); //DD^dagger*d
  
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
            //if (mpi::rank2d == 0)
            //    std::cout << "Converged in " << k << " iterations" << " Error " << err << std::endl;
            return 1;
        }

        beta = err_sqr / r_norm2; //beta = (r_{i+1},r_{i+1})/(r_i,r_i)

        //d_{i+1} = r_{i+1} + beta*d_i 
        for(int n = 0; n<maxSize; n++){
            d.mu0[n] *= beta; 
            d.mu1[n] *= beta;
            d.mu0[n] += r.mu0[n];
            d.mu1[n] += r.mu1[n];
        }
               
        r_norm2 = err_sqr;
        k++;
    }
    if (rank2d == 0)
        std::cout << "CG for DD^+ did not converge in " << CG::max_iter << " iterations" << " Error " << err << std::endl;
    return 0;
}