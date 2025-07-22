#include "bi_cgstab.h"


//Solves Dx x = phi with the Bi-CGstab method
spinor bi_cgstab(void (*func)(const c_matrix&, const spinor&, spinor&,const double&), const int& dim1, const int& dim2,
const c_matrix& U, const spinor& phi, const spinor& x0, const double& m0, const int& max_iter, const double& tol, const bool& print_message) {

    int k = 0; //Iteration number
    double err; // ||r||

    spinor r(dim1, c_vector(dim2, 0));  //r[coordinate][spin] residual
    spinor r_tilde(dim1, c_vector(dim2, 0));  //r[coordinate][spin] residual
    spinor d(dim1, c_vector(dim2, 0)); //search direction
    spinor s(dim1, c_vector(dim2, 0));
    spinor t(dim1, c_vector(dim2, 0));
    spinor Ad(dim1, c_vector(dim2, 0)); //D*d
    spinor x(dim1, c_vector(dim2, 0)); //solution
    c_double alpha, beta, rho_i, omega, rho_i_2;;
    x = x0; //initial solution
    spinor Dphi(dim1, c_vector(dim2, 0)); //Temporary spinor for D x
    func(U, x, Dphi, m0);
    r = phi - Dphi; //r = b - A*x
    r_tilde = r;
	double norm_phi = sqrt(std::real(dot(phi, phi))); //norm of the right hand side

    while (k<max_iter) {
        rho_i = dot(r, r_tilde); //r . r_dagger
        if (k == 0) {
            d = r; //d_1 = r_0
        }
        else {
            beta = alpha * rho_i / (omega * rho_i_2); //beta_{i-1} = alpha_{i-1} * rho_{i-1} / (omega_{i-1} * rho_{i-2})
            d = r + beta * (d - omega * Ad); //d_i = r_{i-1} + beta_{i-1} * (d_{i-1} - omega_{i-1} * Ad_{i-1})
        }
        func(U, d, Ad, m0);  //A d_i 
        alpha = rho_i / dot(Ad, r_tilde); //alpha_i = rho_{i-1} / (Ad_i, r_tilde)
        s = r - alpha * Ad; //s = r_{i-1} - alpha_i * Ad_i
        err = sqrt(std::real(dot(s, s)));
        
        if (err < tol * norm_phi) {
            x = x + alpha * d;
            if (print_message == true) {
                std::cout << "Bi-CG-stab for D converged in " << k+1 << " iterations" << " Error " << err << std::endl;
            }
            return x;
        }
        func(U, s, t,m0);   //A s
        omega = dot(s, t) / dot(t, t); //omega_i = t^dagg . s / t^dagg . t
        r = s - omega * t; //r_i = s - omega_i * t
        x = x + alpha * d + omega * s; //x_i = x_{i-1} + alpha_i * d_i + omega_i * s
        rho_i_2 = rho_i; //rho_{i-2} = rho_{i-1}
        k++;
    }
    if (print_message == true) {
        std::cout << "Bi-CG-stab for D did not converge in " << max_iter << " iterations" << " Error " << err << std::endl;
    }
    return x;
}




