#include "fgmres.h"
#include "iomanip"

//Solves D psi = phi with FGMRES using a parallelized version of the SAP preconditioner
spinor fgmresSAP(const c_matrix& U, const spinor& phi, const spinor& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message) {
    using namespace LV; //Use the lattice variables namespace
    int k = 0; //Iteration number (restart cycle)
    double err; // ||r|| residual norm

    spinor r(Ntot, c_vector(2, 0));  //r[coordinate][spin] residual
    //VmT[column vector index][vector arrange in matrix form]
    std::vector<spinor> VmT(m+1, spinor(Ntot, c_vector(2, 0))); //V matrix transpose
    std::vector<spinor> ZmT(m, spinor(Ntot, c_vector(2, 0)));  //Z matrix transpose

    c_matrix Hm(m+1 , c_vector(m, 0)); //H matrix (Hessenberg matrix)
    c_vector gm(m + 1, 0); 

    //Elements of the Givens rotation matrix |sn[i]|^2 + |cn[i]|^2 = 1
    c_vector sn(m, 0);
    c_vector cn(m, 0);
    c_vector eta(m, 0);


    spinor w(Ntot, c_vector(2, 0)); 
    spinor x = x0; //initial solution
    c_double beta; //not 1/g^2 from lattice simulations 
 
    spinor Dphi(Ntot, c_vector(2, 0)); //Temporary spinor for D x
    D_phi(U, x,Dphi, m0);
    axpy(phi,Dphi, -1.0, r); //r = b - A*x
	double norm_phi = sqrt(std::real(dot(phi, phi))); //norm of the right hand side


    err = sqrt(std::real(dot(r, r))); //Initial error
    while (k < restarts) {
        beta = err + 0.0 * I_number;
        scal(1.0/beta, r,VmT[0]); //VmT[0] = r / ||r||
        gm[0] = beta; //gm[0] = ||r||
        //-----Arnoldi process to build the Krylov basis and the Hessenberg matrix-----//
        for (int j = 0; j < m; j++) {
            //-----Preconditioner-----//
            //zm = M^-1. vm
            set_zeros(ZmT[j], Ntot, 2); //Initialize ZmT[j] to zero
            SAP(U, VmT[j], ZmT[j], m0, 1,SAPV::sap_blocks_per_proc); //One SAP iteration

        
            D_phi(U, ZmT[j],w, m0); //w = D v_j

            //----Gram-Schmidt process----//
            for (int i = 0; i <= j; i++) {
                Hm[i][j] = dot(w, VmT[i]); //  (v_i^dagger, w)
                //w = w -  Hm[i][j] * VmT[i];
				for(int n=0; n<Ntot; n++){
					for(int l=0; l<2; l++){
						w[n][l] -= Hm[i][j] * VmT[i][n][l];
					}
				}
            }

            Hm[j + 1][j] = sqrt(std::real(dot(w, w))); //H[j+1][j] = ||A v_j||
            if (std::real(Hm[j + 1][j]) > 0) {
                scal(1.0 / Hm[j + 1][j], w, VmT[j + 1]); //VmT[j + 1] = w / ||A v_j||
            }
            //----Rotate the matrix----//
            rotation(cn, sn, Hm, j);

            //Rotate gm
            gm[j + 1] = -sn[j] * gm[j];
            gm[j] = std::conj(cn[j]) * gm[j];
        }        
        //Solve the upper triangular system//
		solve_upper_triangular(Hm, gm,m,eta);
        
        for (int i = 0; i < 2 * Ntot; i++) {
            int n = i / 2; int mu = i % 2;
            for (int j = 0; j < m; j++) {
                x[n][mu] = x[n][mu] + eta[j] * ZmT[j][n][mu]; 
            }
        }
        //Compute the residual
        D_phi(U, x,Dphi, m0);
        axpy(phi,Dphi, -1.0, r); //r = b - A*x
        err = sqrt(std::real(dot(r, r)));
         if (err < tol* norm_phi) {
             if (print_message == true) {
                 std::cout << "FGMRES with SAP preconditioner converged in " << k + 1 << " iterations" << " Error " << err << std::endl;
             }
             return x;
         }
         k++;
    }
    if (print_message == true) {
        std::cout << "FGMRES with SAP preconditioner did not converge in " << restarts << " restarts" << " Error " << err << std::endl;
    }
    return x;
}



//Solves D psi = phi with FGMRES using a two-grid preconditioner
spinor fgmresAMG(const c_matrix& U, const spinor& phi, const spinor& x0, const double& m0, const int& m, const int& restarts, const double& tol, const bool& print_message) {
    using namespace LV; //Use the lattice variables namespace
    int k = 0; //Iteration number (restart cycle)
    double err;


    spinor r(Ntot, c_vector(2, 0));  //r[coordinate][spin] residual

    //VmT[column vector index][vector arrange in matrix form]
    std::vector<spinor> VmT(m+1, spinor(Ntot, c_vector(2, 0))); //V matrix transpose
    std::vector<spinor> ZmT(m, spinor(Ntot, c_vector(2, 0)));  //Z matrix transpose

    c_matrix Hm(m+1 , c_vector(m, 0)); //H matrix (Hessenberg matrix)
    c_vector gm(m + 1, 0); 

    //Elements of the rotation matrix |sn[i]|^2 + |cn[i]|^2 = 1
    c_vector sn(m, 0);
    c_vector cn(m, 0);
    c_vector eta(m, 0);


    spinor w(Ntot, c_vector(2, 0)); 
    spinor x = x0; //initial solution
    c_double beta;

    spinor Dphi(Ntot, c_vector(2, 0)); //Temporary spinor for D x
    D_phi(U, x, Dphi, m0);
    axpy(phi,Dphi, -1.0, r); //r = b - A*x
	double norm_phi = sqrt(std::real(dot(phi, phi))); //norm of the right hand side

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    //       Set up phase for the two-grid method         //
    GaugeConf Gconf = GaugeConf(Nx, Nt); //Gauge configuration
    Gconf.setGconf(U); //Set the gauge configuration
    AMG amg = AMG(Gconf, m0,AMGV::nu1,AMGV::nu2);   //nu1 pre-smoothing it, nu2 post-smoothing it

    double elapsed_time;
    double startT, endT;     
    startT = MPI_Wtime();
    amg.setUpPhase(1, AMGV::Nit); //test vectors intialization
    endT = MPI_Wtime();
    elapsed_time = endT - startT;
    std::cout << "[MPI Process " << rank << "] Elapsed time for Set-up phase = " << elapsed_time << " seconds" << std::endl;   
    //--------------------------------------------------//

    err = sqrt(std::real(dot(r, r))); //Initial error
    while (k < restarts) {
        beta = err + 0.0 * I_number;
        scal(1.0/beta, r,VmT[0]); //VmT[0] = r / ||r||
        gm[0] = beta; //gm[0] = ||r||
        //-----Arnoldi process to build the Krylov basis and the Hessenberg matrix-----//
        for (int j = 0; j < m; j++) {
            set_zeros(ZmT[j], Ntot, 2); //Initialize ZmT[j] to zero
            ZmT[j] = amg.TwoGrid(1, 1e-10, ZmT[j], VmT[j], false); //One two-grid iteration
           
            D_phi(U, ZmT[j],w, m0); //w = D v_j
            //Gram-Schmidt process to orthogonalize the vectors
            for (int i = 0; i <= j; i++) {
                Hm[i][j] = dot(w, VmT[i]); //  (v_i^dagger, w)
                //w = w -  Hm[i][j] * VmT[i];
                for(int n=0; n<Ntot; n++){
					for(int l=0; l<2; l++){
						w[n][l] -= Hm[i][j] * VmT[i][n][l];
					}
				}
            }

            Hm[j + 1][j] = sqrt(std::real(dot(w, w))); //H[j+1][j] = ||A v_j||
            if (std::real(Hm[j + 1][j]) > 0) {
                scal(1.0 / Hm[j + 1][j], w, VmT[j + 1]); //VmT[j + 1] = w / ||A v_j||
            }
            //----Rotate the matrix----//
            rotation(cn, sn, Hm, j);

            //Rotate gm
            gm[j + 1] = -sn[j] * gm[j];
            gm[j] = std::conj(cn[j]) * gm[j];
        }        
        //Solve the upper triangular system//
		solve_upper_triangular(Hm, gm,m,eta);
        
        for (int i = 0; i < 2 * Ntot; i++) {
            int n = i / 2; int mu = i % 2;
            for (int j = 0; j < m; j++) {
                x[n][mu] = x[n][mu] + eta[j] * ZmT[j][n][mu]; 
            }
        }
        //Compute the residual
        D_phi(U, x, Dphi, m0);
        axpy(phi,Dphi, -1.0, r); //r = b - A*x
        err = sqrt(std::real(dot(r, r)));
        
        //std::cout << "FGMRES iteration " << k + 1 << " Error " << std::setprecision(17) << err << "  from process " << rank << std::endl;
        
         if (err < tol* norm_phi) {
             if (print_message == true) {
                 std::cout << "FGMRES with AMG preconditioning converged in " << k + 1 << " iterations" << " Error " << err << std::endl;
             }
             return x;
         }
         k++;
    }
    if (print_message == true) {
        std::cout << "FGMRES with AMG preconditioning for D did not converge in " << restarts << " restarts" << " Error " << err << std::endl;
    }
    return x;
}


//FGMRES still uses these functions
void rotation(c_vector& cn, c_vector& sn, c_matrix& H, const int& j) {
    //Rotation of the column elements that are <j
    for (int i = 0; i < j; i++) {
		c_double temp = std::conj(cn[i]) * H[i][j] + std::conj(sn[i]) * H[i + 1][j];
		H[i + 1][j] = -sn[i] * H[i][j] + cn[i] * H[i + 1][j];
		H[i][j] = temp;
    }
    //Rotation of the diagonal and element right below the diagonal
    c_double den = sqrt(std::conj(H[j][j] ) * H[j][j] + std::conj(H[j + 1][j]) * H[j + 1][j]);
	sn[j] = H[j + 1][j] / den; cn[j] = H[j][j] / den;
	H[j][j] = std::conj(cn[j]) * H[j][j] + std::conj(sn[j]) * H[j + 1][j];
    H[j + 1][j] = 0.0;

}

//x = A^-1 b, A an upper triangular matrix of dimension n
void solve_upper_triangular(const c_matrix& A, const c_vector& b, const int& n, c_vector& out) {
	for (int i = n - 1; i >= 0; i--) {
		out[i] = b[i];
		for (int j = i + 1; j < n; j++) {
			out[i] -= A[i][j] * out[j];
		}
		out[i] /= A[i][i];
	}
}