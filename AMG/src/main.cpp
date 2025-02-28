#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include "amg.h"

int main() {
    srand(time(0));
    initialize_matrices(); //Initialize gamma matrices, identity and unit vectors
    Coordinates(); //Vectorized coordinates
    periodic_boundary(); //Builds LeftPB and RightPB (periodic boundary for U_mu(n))
    Aggregates(); //Aggregates
    std::cout << "Ns = " << Ns << " Nt = " << Nt << std::endl;
    std::cout << "block_x = " << block_x << " block_t = " << block_t << std::endl;
    std::cout << "x_elements = " << x_elements << " t_elements = " << t_elements << std::endl;
    std::cout << "Number of test vectors = " << Ntest << std::endl;
    GaugeConf GConf = GaugeConf(Ns, Nt);
    GConf.initialization(); //Initialize the gauge configuration
    PrintAggregates();
    double m0 = 1;

    c_matrix phi(Ntot, c_vector(2, 0));
    for (int i = 0; i < Ntot; i++) {
        for (int j = 0; j < 2; j++) {
            phi[i][j] = 1.0 * RandomU1();
        }
    }
    std::cout << "------phi------" << std::endl;
    //PrintVector(phi);
    std::cout << phi[0][0] << " " << phi[0][1] << std::endl;
    std::cout << "################Bi-CGstab inversion##############" << std::endl;
    c_matrix x = bi_cgstab(GConf.Conf, phi, phi, m0, 1000, 1e-10, true);//bi_cgstab(GConf.Conf, phi,m0); //D^-1 phi
    std::cout << "------D^-1 phi------" << std::endl;
    //PrintVector(x);
    std::cout << x[0][0] << " " << x[0][1] << std::endl;
    std::cout << "------D D^-1 phi------" << std::endl;
    x = D_phi(GConf.Conf, x, m0); //D D^-1 phi
    //PrintVector(x);
    std::cout << x[0][0] << " " << x[0][1] << std::endl;
    std::cout << std::endl;
    std::cout << "################Multigrid inversion##############" << std::endl;
    AMG amg = AMG(GConf, Ns, Nt, Ntest, m0);

    amg.tv_init(1, 3); //test vectors intialization
    std::vector<std::vector<std::complex<double>>> x0(Ntot, std::vector<std::complex<double>>(2, 0));
    x0 = amg.TwoGrid(1, 1, x0, phi, true);
    std::cout << "------D^-1 phi------" << std::endl;
    //PrintVector(x0);
    std::cout << x0[0][0] << " " << x0[0][1] << std::endl;
    std::cout << "------D D^-1 phi------" << std::endl;
    x0 = D_phi(GConf.Conf, x0, m0); //D D^-1 phi
    std::cout << x0[0][0] << " " << x0[0][1] << std::endl;
    //PrintVector(x0);

    return 0;
}



/*
    AMG amg = AMG(GConf,Ns,Nt,Ntest,m0);
    amg.tv_init(0.5); //test vectors intialization
	//Print P e_i where e_i are the basis vectors. The result should return the test vectors chopped over the aggregates.
    std::vector<std::vector<std::vector<std::complex<double>>>> tv = amg.test_vectors;
    for(int i =0; i<Ntest; i++){
        std::cout << "------tv" << i << "------" << std::endl;
        PrintVector(tv[i]);
    }
    for(int i = 0; i < Ntest*block_x*block_t; i++){
        std::cout << "------ P e_" << i << "------" << std::endl;
        std::vector<std::vector<std::complex<double>>> e_i = canonical_vector(i, Ntest,block_x*block_t);
        PrintVector(e_i);
        std::vector<std::vector<std::complex<double>>> Pe = amg.P_v(e_i); 
        PrintVector(Pe);    
    }
    for(int i = 0; i < 2*Ntot;i++){
        std::cout << "------P^T e_" << i << "------" << std::endl;
        std::vector<std::vector<std::complex<double>>> e_i = canonical_vector(i, Ntot, 2);
        //PrintVector(e_i);
        std::vector<std::vector<std::complex<double>>> Pe = amg.Pt_v(e_i); 
        PrintVector(Pe);
    }
    
*/