#include "bi_cgstab.h"


//Solves Dx x = phi with the Bi-CGstab method
int bi_cgstab(const GaugeConf& GConf, const spinor& phi, const spinor& x0, spinor& x) {
    spinor r(mpi::maxSizeH);  //r[coordinate][spin] residual
    spinor r_tilde(mpi::maxSizeH);  //r[coordinate][spin] residual
    spinor d(mpi::maxSizeH); //search direction
    spinor s(mpi::maxSizeH);
    spinor t(mpi::maxSizeH);
    spinor Ad(mpi::maxSizeH); //D*d
    spinor Dphi(mpi::maxSizeH); //Temporary spinor for D x

    int k = 0; //Iteration number
    double err; // ||r||
    c_double alpha, beta, rho_i, omega, rho_i_2;

    x = x0; //initial solution

    //Pointers to the data of the spinors, captured by the kernels (check backend.h)
    const c_double* phi_v = phi.val;
    c_double* x_v    = x.val;
    c_double* r_v    = r.val;
    c_double* d_v    = d.val;
    c_double* s_v    = s.val;
    c_double* t_v    = t.val;
    c_double* Ad_v   = Ad.val;
    c_double* Dphi_v = Dphi.val;

    //std::cout << "U[0] from rank " << mpi::rank2d << "   " << U.val[(mpi::width_t+2)+1] << std::endl;
    D_phi(GConf, x, Dphi);
    for_each_site([=] HOST_DEVICE (int n){
        for(int mu=0; mu<2; mu++)
            r_v[2*n+mu] = phi_v[2*n+mu]-Dphi_v[2*n+mu];  //r = b - A*x
    });
    r_tilde = r;
	double norm_phi = sqrt(real(dot(phi, phi))); //norm of the right hand side

    while (k<BiCG::max_iter) {
        rho_i = dot(r, r_tilde); //r . r_dagger

        if (k == 0) {
            d = r; //d_1 = r_0
        }
        else {
            beta = alpha * rho_i / (omega * rho_i_2); //beta_{i-1} = alpha_{i-1} * rho_{i-1} / (omega_{i-1} * rho_{i-2})
            //d = r + beta * (d - omega * Ad);
            for_each_site([=] HOST_DEVICE (int n){
                for(int mu=0; mu<2; mu++)
                    d_v[2*n+mu] = r_v[2*n+mu] + beta * (d_v[2*n+mu] - omega * Ad_v[2*n+mu]); //d_i = r_{i-1} + beta_{i-1} * (d_{i-1} - omega_{i-1} * Ad_{i-1})
            });
        }

        D_phi(GConf, d, Ad);  //A d_i
        alpha = rho_i / dot(Ad, r_tilde); //alpha_i = rho_{i-1} / (Ad_i, r_tilde)
        //s = r - alpha * Ad; //s = r_{i-1} - alpha_i * Ad_i
        //s and its norm are computed in the same kernel
        err = sqrt(lattice_sum<double>([=] HOST_DEVICE (int n){
            for(int mu=0; mu<2; mu++)
                s_v[2*n+mu] = r_v[2*n+mu] - alpha * Ad_v[2*n+mu]; //s_i = r_{i-1} - alpha_i * Ad_i
            return norm(s_v[2*n]) + norm(s_v[2*n+1]);
        }));

        if (err < BiCG::tol * norm_phi) {
            for_each_site([=] HOST_DEVICE (int n){
                for(int mu=0; mu<2; mu++)
                    x_v[2*n+mu] += alpha * d_v[2*n+mu]; //x = x + alpha * d;
            });
            if (BiCG::print_convergence_message == true && mpi::rank2d == 0) {
                std::cout << "Bi-CG-stab for D converged in " << k+1 << " iterations" << " Error " << err << std::endl;
            }
            return k+1;
        }
        D_phi(GConf, s, t);   //A s
        omega = dot(s, t) / dot(t, t); //omega_i = t^dagg . s / t^dagg . t

        for_each_site([=] HOST_DEVICE (int n){
            for(int mu=0; mu<2; mu++){
                r_v[2*n+mu] = s_v[2*n+mu] - omega * t_v[2*n+mu];   //r_i = s - omega_i * t
                x_v[2*n+mu] = x_v[2*n+mu] + alpha * d_v[2*n+mu] + omega * s_v[2*n+mu]; //x_i = x_{i-1} + alpha_i * d_i + omega_i * s_i
            }
        });
        rho_i_2 = rho_i; //rho_{i-2} = rho_{i-1}
        k++;
    }
    if (mpi::rank2d == 0) {
        std::cout << "Bi-CG-stab for D did not converge in " << BiCG::max_iter << " iterations with a relative tolerance of "
        << BiCG::tol*norm_phi << ". The final residual norm is " << err << std::endl;
    }

    return BiCG::max_iter;
}
