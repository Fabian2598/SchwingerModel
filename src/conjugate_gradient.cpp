#include "conjugate_gradient.h"

int conjugate_gradient(const spinor& U, const spinor& phi, spinor &sol){
    using namespace mpi;
    int k = 0; //Iteration number
    double err;
    double err_sqr;

    spinor r(maxSizeH);  //r[coordinate][spin] residual
    spinor d(maxSizeH); //search direction
    spinor Ad(maxSizeH); //DD^dagger*d
  
    c_double alpha, beta;

	sol = phi;
    D_D_dagger_phi(U, sol, Ad); //DD^dagger*x
    
    int n;
    for(int x = 1; x<=width_x; x++){
        for(int t = 1; t<=width_t; t++){
            n = x*(width_t+2)+t;
            r.val[2*n]   = phi.val[2*n]   - Ad.val[2*n]; //mu0
            r.val[2*n+1] = phi.val[2*n+1] - Ad.val[2*n+1]; //mu1
        }
    }


    d = r; //initial search direction
 
    c_double r_norm2 = dot(r, r);
    
    
    double phi_norm2 = sqrt(std::real(dot(phi, phi)));

    while (k<CG::max_iter) {
        D_D_dagger_phi(U, d,Ad); //DD^dagger*d 
        alpha = r_norm2 / dot(d, Ad); //alpha = (r_i,r_i)/(d_i,Ad_i)
        for(int x = 1; x<=width_x; x++){
            for(int t = 1; t<=width_t; t++){
                n =  x*(width_t+2)+t;
                //x = x + alpha * d; //x_{i+1} = x_i + alpha*d_i 
                sol.val[2*n]    += alpha*d.val[2*n];
                sol.val[2*n+1]  += alpha*d.val[2*n+1];
                //r = r - alpha * Ad; //r_{i+1} = r_i - alpha*Ad_i
                r.val[2*n]      -= alpha*Ad.val[2*n];
                r.val[2*n+1]    -= alpha*Ad.val[2*n+1];

            }
        }
        
        err_sqr = std::real(dot(r, r)); //err_sqr = (r_{i+1},r_{i+1})
		err = sqrt(err_sqr); // err = sqrt(err_sqr)
        if (err < CG::tol*phi_norm2) {
            if (mpi::rank2d == 0 && CG::print_convergence_message == true)
                std::cout << "CG converged in " << k+1 << " iterations" << " Error " << err << std::endl;
            return k+1;
        }

        beta = err_sqr / r_norm2; //beta = (r_{i+1},r_{i+1})/(r_i,r_i)

        //d_{i+1} = r_{i+1} + beta*d_i 
       for(int x = 1; x<=width_x; x++){
            for(int t = 1; t<=width_t; t++){
                n =  x*(width_t+2)+t;
                d.val[2*n]      *= beta; 
                d.val[2*n+1]    *= beta;
                d.val[2*n]      += r.val[2*n];
                d.val[2*n+1]    += r.val[2*n+1];
            }
        }
               
        r_norm2 = err_sqr;
        k++;
    }
    if (rank2d == 0)
        std::cout << "CG for DD^+ did not converge in " << CG::max_iter << " iterations" << " Error " << err << std::endl;
    return 0;
}