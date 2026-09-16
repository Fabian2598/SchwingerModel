#include "bi_cgstab.h"


//Solves Dx x = phi with the Bi-CGstab method
int bi_cgstab(const GaugeConf& GConf, const spinor& phi, const spinor& x0, spinor& x) {
    using namespace BiCG;
    spinor r(mpi::maxSizeH);  //r[coordinate][spin] residual
    spinor r_tilde(mpi::maxSizeH);  //r[coordinate][spin] residual
    spinor d(mpi::maxSizeH); //search direction
    spinor s(mpi::maxSizeH);
    spinor t(mpi::maxSizeH);
    spinor Ad(mpi::maxSizeH); //D*d
    spinor Dphi(mpi::maxSizeH); //Temporary spinor for D x

    int k = 0; //Iteration number
    double err; // ||r||
    int index;
    c_double alpha, beta, rho_i, omega, rho_i_2;

    x = x0; //initial solution
    //std::cout << "U[0] from rank " << mpi::rank2d << "   " << U.val[(mpi::width_t+2)+1] << std::endl; 
    D_phi(GConf, x, Dphi);
    for(int nx = 1; nx<=mpi::width_x; nx++){
    for(int nt = 1; nt<=mpi::width_t; nt++){
    for(int mu=0; mu<2; mu++){
        index = idx(nx,nt,mu);
        r.val[index] = phi.val[index]-Dphi.val[index];  //r = b - A*x
    }
    }
    }
    r_tilde = r;
	double norm_phi = sqrt(std::real(dot(phi, phi))); //norm of the right hand side
    
    while (k<BiCG::max_iter) {
        rho_i = dot(r, r_tilde); //r . r_dagger

        if (k == 0) {
            d = r; //d_1 = r_0
        }
        else {
            beta = alpha * rho_i / (omega * rho_i_2); //beta_{i-1} = alpha_{i-1} * rho_{i-1} / (omega_{i-1} * rho_{i-2})
            //d = r + beta * (d - omega * Ad);
            for(int nx = 1; nx<=mpi::width_x; nx++){
            for(int nt = 1; nt<=mpi::width_t; nt++){
            for(int mu=0; mu<2; mu++){
                index = idx(nx,nt,mu);
                d.val[index] = r.val[index] + beta * (d.val[index] - omega * Ad.val[index]); //d_i = r_{i-1} + beta_{i-1} * (d_{i-1} - omega_{i-1} * Ad_{i-1})
            }
            }
            }
        }

        D_phi(GConf, d, Ad);  //A d_i 
        alpha = rho_i / dot(Ad, r_tilde); //alpha_i = rho_{i-1} / (Ad_i, r_tilde)
        //s = r - alpha * Ad; //s = r_{i-1} - alpha_i * Ad_i
        for(int nx = 1; nx<=mpi::width_x; nx++){
        for(int nt = 1; nt<=mpi::width_t; nt++){
        for(int mu=0; mu<2; mu++){
                index = idx(nx,nt,mu);
                s.val[index] = r.val[index] - alpha * Ad.val[index]; //s_i = r_{i-1} - alpha_i * Ad_i
        }
        }
        }

        err = sqrt(std::real(dot(s, s)));
        
        if (err < BiCG::tol * norm_phi) {
            for(int nx = 1; nx<=mpi::width_x; nx++){
            for(int nt = 1; nt<=mpi::width_t; nt++){
            for(int mu=0; mu<2; mu++){
                    index = idx(nx,nt,mu);
                    x.val[index] += alpha * d.val[index]; //x = x + alpha * d;
            }
            }
            }
            if (BiCG::print_convergence_message == true && mpi::rank2d == 0) {
                std::cout << "Bi-CG-stab for D converged in " << k+1 << " iterations" << " Error " << err << std::endl;
            }
            return k+1;
        }
        D_phi(GConf, s, t);   //A s
        omega = dot(s, t) / dot(t, t); //omega_i = t^dagg . s / t^dagg . t

        for(int nx = 1; nx<=mpi::width_x; nx++){
        for(int nt = 1; nt<=mpi::width_t; nt++){
        for(int mu=0; mu<2; mu++){
            index = idx(nx,nt,mu);
            r.val[index] = s.val[index] - omega * t.val[index];   //r_i = s - omega_i * t
            x.val[index] = x.val[index] + alpha * d.val[index] + omega * s.val[index]; //x_i = x_{i-1} + alpha_i * d_i + omega_i * s_i
        }
        }
        }
        rho_i_2 = rho_i; //rho_{i-2} = rho_{i-1}
        k++;
    }
    if (mpi::rank2d == 0) {
        std::cout << "Bi-CG-stab for D did not converge in " << BiCG::max_iter << " iterations with a relative tolerance of "
        << BiCG::tol*norm_phi << ". The final residual norm is " << err << std::endl;
    }
    
    return BiCG::max_iter;
}