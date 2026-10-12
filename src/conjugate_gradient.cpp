#include "conjugate_gradient.h"


int conjugate_gradient(const GaugeConf& GConf, const spinor& phi, spinor &sol){
    int k = 0; //Iteration number
    double err;
    double err_sqr;


    spinor r(mpi::maxSizeH);  //r[coordinate][spin] residual
    spinor d(mpi::maxSizeH);  //search direction
    spinor Ad(mpi::maxSizeH); //DD^dagger*d

    c_double alpha, beta;

	sol = phi;
    D_D_dagger_phi(GConf, sol, Ad); //DD^dagger*x

    //Pointers to the data of the spinors, captured by the kernels (check backend.h)
    const c_double* phi_v = phi.val;
    c_double* sol_v = sol.val;
    c_double* r_v   = r.val;
    c_double* d_v   = d.val;
    c_double* Ad_v  = Ad.val;

    for_each_site([=] HOST_DEVICE (int n){
        r_v[2*n]   = phi_v[2*n]   - Ad_v[2*n]; //mu0
        r_v[2*n+1] = phi_v[2*n+1] - Ad_v[2*n+1]; //mu1
    });


    d = r; //initial search direction

    c_double r_norm2 = dot(r, r);


    double phi_norm2 = sqrt(real(dot(phi, phi)));

    while (k<CG::max_iter) {
        D_D_dagger_phi(GConf, d,Ad); //DD^dagger*d
        alpha = r_norm2 / dot(d, Ad); //alpha = (r_i,r_i)/(d_i,Ad_i)
        //x and r are updated and (r_{i+1},r_{i+1}) is computed in the same kernel
        err_sqr = lattice_sum<double>([=] HOST_DEVICE (int n){
            //x = x + alpha * d; //x_{i+1} = x_i + alpha*d_i
            sol_v[2*n]    += alpha*d_v[2*n];
            sol_v[2*n+1]  += alpha*d_v[2*n+1];
            //r = r - alpha * Ad; //r_{i+1} = r_i - alpha*Ad_i
            r_v[2*n]      -= alpha*Ad_v[2*n];
            r_v[2*n+1]    -= alpha*Ad_v[2*n+1];
            return norm(r_v[2*n]) + norm(r_v[2*n+1]); //err_sqr = (r_{i+1},r_{i+1})
        });
		err = sqrt(err_sqr); // err = sqrt(err_sqr)
        if (err < CG::tol*phi_norm2) {
            if (mpi::rank2d == 0 && CG::print_convergence_message == true)
                std::cout << "CG converged in " << k+1 << " iterations" << " Error " << err << std::endl;
            return k+1;
        }

        beta = err_sqr / r_norm2; //beta = (r_{i+1},r_{i+1})/(r_i,r_i)

        //d_{i+1} = r_{i+1} + beta*d_i
        for_each_site([=] HOST_DEVICE (int n){
            d_v[2*n]      *= beta;
            d_v[2*n+1]    *= beta;
            d_v[2*n]      += r_v[2*n];
            d_v[2*n+1]    += r_v[2*n+1];
        });

        r_norm2 = err_sqr;
        k++;
    }
    if (mpi::rank2d == 0){
        std::cout << "CG for DD^+ did not converge in " << CG::max_iter << " iterations with a relative tolerance of "
        << CG::tol*phi_norm2 << ". The final residual norm is " << err << std::endl;
    }
    return 0;
}
