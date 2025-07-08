#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "bi_cgstab.h"
#include "fgmres.h"
#include "mpi.h"
#include <thread>
#include <chrono>


int main(int argc, char **argv) {

    using namespace SAPV;
    //srand(18);
    srand(time(0));

    initialize_matrices(); //Initialize gamma matrices, identity and unit vectors
    Coordinates(); //Builds array with coordinates of the lattice points x * Nt + t 
    periodic_boundary(); //Builds LeftPB and RightPB (periodic boundary for U_mu(n))

    double m0 = -0.2; 
    
    GaugeConf GConf = GaugeConf(Nx, Nt);
    GConf.initialize(); //Initialize a random gauge configuration

    spinor rhs(Ntot, c_vector(2, 0)); //random right hand side
    spinor X(Ntot, c_vector(2, 0)); //solution vector 
    spinor Xold (Ntot, c_vector(2, 0));

    clock_t start, end;
    double elapsed_time;
    //Random right hand side
    for(int i = 0; i < Ntot; i++) {
        rhs[i][0] = RandomU1();
        rhs[i][1] = RandomU1();
    }

    int tries = 100000;

 

    start = clock();
    for(int i = 0; i<tries; i++) {
     X = D_phi(GConf.Conf, rhs, m0); //Apply the operator D to the right hand side
    }
    end = clock();
    elapsed_time = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed time for new = " << elapsed_time << " seconds" << std::endl;   
    
    
     start = clock();
    for(int i = 0; i<tries; i++) {
     Xold = D_phi_old(GConf.Conf, rhs, m0); //Apply the old operator D to the right hand side
    }
    end = clock();
    elapsed_time = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed time for old = " << elapsed_time << " seconds" << std::endl;   

    /*
    for(int x = 0; x < Nx; x++) {
        for(int t = 0; t < Nt; t++) {
            int n = x * Nt + t;
            std::cout << "n=" << n << ", x[n][0]=" << X[n][0] << ", x[n][1]=" << X[n][1] << std::endl;
            std::cout << "n=" << n << ", xold[n][0]=" << Xold[n][0] << ", xold[n][1]=" << Xold[n][1] << std::endl;
            std::cout << std::endl;
        }
    }
    for(int i = 0; i < Ntot; i++) {
        if(std::abs(X[i][0] - Xold[i][0]) > 1e-10 || std::abs(X[i][1] - Xold[i][1]) > 1e-10) {
            std::cout << "Error in the operator D_phi " << std::endl;
            std::cout << "n " << i << std::endl;
            return 1;
        }
    }
    */

    //Bi-cgstab inversion for comparison
    
    /*
    std::cout << "--------------Bi-CGstab inversion--------------" << std::endl;
    start = clock();
    spinor x0(Ntot, c_vector(2, 0)); //Initial guess
    spinor x_bi = bi_cgstab(&D_phi,Ntot,2,GConf.Conf, rhs, x0, m0, 100000, 1e-15, true);
    end = clock();
    elapsed_time = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed time for Bi-CGstab = " << elapsed_time << " seconds" << std::endl;    

    std::cout << "--------------Bi-CGstab old inversion--------------" << std::endl;
    start = clock();
    spinor x_bi_old = bi_cgstab(&D_phi_old,Ntot,2,GConf.Conf, rhs, x0, m0, 100000, 1e-15, true);
    end = clock();
    elapsed_time = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed time for Bi-CGstab = " << elapsed_time << " seconds" << std::endl;    

    for(int i = 0; i < Ntot; i++) {
        if(std::abs(x_bi[i][0] - x_bi_old[i][0]) > 1e-6 || std::abs(x_bi[i][1] - x_bi_old[i][1]) > 1e-6) {
            std::cout << "Error in the operator D_phi" << std::endl;
            return 1;
        }
    }
    std::cout << "D_phi operator is working correctly" << std::endl;
    */

    return 0;
}