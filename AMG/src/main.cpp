#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include "amg.h"

int main() {
    srand(time(0));//srand(0); //srand(time(0));
    initialize_matrices(); //Initialize gamma matrices, identity and unit vectors
    Coordinates(); //Vectorized coordinates
    periodic_boundary(); //Builds LeftPB and RightPB (periodic boundary for U_mu(n))
    std::cout << "Ns = " << Ns << " Nt = " << Nt << std::endl;
    std::cout << "block_x = " << block_x << " block_t = " << block_t << std::endl;
    std::cout << "x_elements = " << x_elements << " t_elements = " << t_elements << std::endl;
    Aggregates(); //Aggregates
    std::cout << "Number of test vectors = " << Ntest << std::endl;
    std::cout << "Dirac matrix dimension = " << (2 * Ntot) << std::endl;
    std::cout << "Dc dimension = " << Ntest * Nagg << std::endl;

    GaugeConf GConf = GaugeConf(Ns, Nt);
    GConf.initialization(); //Initialize the gauge configuration
    //PrintAggregates();
    double m0 = -0.5; //for negative masses close to zero this gets very ill-conditioned for V = 8^2

    c_matrix PHI(Ntot, c_vector(2, 0));
    for (int i = 0; i < Ntot; i++) {
        for (int j = 0; j < 2; j++) {
            PHI[i][j] = 1.0 * RandomU1();
        }
    }
  
    std::cout << "################Bi-CGstab inversion##############" << std::endl;
    c_matrix x = bi_cgstab(GConf.Conf, PHI, PHI, m0, 1000, 1e-10, true);//bi_cgstab(GConf.Conf, phi,m0); //D^-1 phi
    std::cout << "------D^-1 phi------" << std::endl;
    std::cout << x[0][0] << " " << x[0][1] << std::endl;

    c_matrix phi(Ntest, c_vector(Nagg, 0));
    for (int i = 0; i < Ntest; i++) {
        for (int j = 0; j < Nagg; j++) {
            phi[i][j] = 1.0 * RandomU1();
        }
    }

    std::cout << "################Set-up phase##############" << std::endl;
    AMG amg = AMG(GConf, Ns, Nt, Ntest, m0);
    amg.tv_init(1, 6); //test vectors intialization
    std::cout << "################Two-grid inversion of D##############" << std::endl;
    int nu1 = 0, nu2 = 2;
    std::cout << "Pre-smoothing steps " << nu1 << " Post-smoothing steps " << nu2 << std::endl;
    c_matrix x0;
    x0 = amg.TwoGrid(nu1, nu2, 100, 1e-10,  PHI, PHI, true);
    std::cout << "------PHI------" << std::endl;
    std::cout << PHI[0][0] << " " << PHI[0][1] << std::endl;
    std::cout << "------D^-1 phi------" << std::endl;
    std::cout << x0[0][0] << " " << x0[0][1] << std::endl;
    std::cout << "------D D^-1 phi------" << std::endl;
    x0 = D_phi(GConf.Conf, x0, m0);
    std::cout << x0[0][0] << " " << x0[0][1] << std::endl;
    

    return 0;
}