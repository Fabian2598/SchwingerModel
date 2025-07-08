#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "bi_cgstab.h"
#include "fgmres.h"
#include "mpi.h"


int main(int argc, char **argv) {

    using namespace SAPV;
    srand(19);
    //srand(time(0));
    initialize_matrices(); //Initialize gamma matrices, identity and unit vectors
    Coordinates(); //Builds array with coordinates of the lattice points x * Nt + t 
    periodic_boundary(); //Builds LeftPB and RightPB (periodic boundary for U_mu(n))
    double m0 = -0.7; 
    
    GaugeConf GConf = GaugeConf(Nx, Nt);
    GConf.initialize(); //Initialize a random gauge configuration

    spinor rhs(Ntot, c_vector(2, 0)); //random right hand side 
    spinor x(Ntot, c_vector(2, 0)); //solution vector 
    //Random right hand side
    for(int i = 0; i < Ntot; i++) {
        rhs[i][0] = RandomU1();
        rhs[i][1] = RandomU1();
    }

    clock_t start, end;
    double elapsed_time;
    int apps = 1000;

   
    //Bi-cgstab inversion for comparison
    std::cout << "Applying the Dirac operator " << apps << " times" << std::endl;
    start = clock();
    for (int i = 0; i < apps; i++) {
        x = D_phi(GConf.Conf, rhs, m0); //Apply the Dirac operator
    }
    end = clock();
    elapsed_time = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed time = " << elapsed_time << " seconds" << std::endl;    
    

    return 0;
}