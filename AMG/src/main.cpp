#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include "amg.h"

int main() {
    srand(0);//srand(time(0));
    initialize_matrices(); //Initialize gamma matrices, identity and unit vectors
    Coordinates(); //Vectorized coordinates
    periodic_boundary(); //Builds LeftPB and RightPB (periodic boundary for U_mu(n))
    Aggregates(); //Aggregates
    std::cout << "Ns = " << Ns << " Nt = " << Nt << std::endl;
    std::cout << "block_x = " << block_x << " block_t = " << block_t << std::endl;
    std::cout << "x_elements = " << x_elements << " t_elements = " << t_elements << std::endl;
    std::cout << "Number of test vectors = " << Ntest << std::endl;
    std::cout << "Dirac matrix dimension = " << (2 * Ntot) << std::endl;
    std::cout << "Dc dimension = " << Ntest * Nagg << std::endl;

    GaugeConf GConf = GaugeConf(Ns, Nt);
    GConf.initialization(); //Initialize the gauge configuration
    //PrintAggregates();
    double m0 = 1;

    c_matrix PHI(Ntot, c_vector(2, 0));
    for (int i = 0; i < Ntot; i++) {
        for (int j = 0; j < 2; j++) {
            PHI[i][j] = 1.0 * RandomU1();
        }
    }

    std::cout << "################Bi-CGstab inversion##############" << std::endl;
    clock_t begin = clock();
    c_matrix x = bi_cgstab(GConf.Conf, PHI, PHI, m0, 1000, 1e-10, true);//bi_cgstab(GConf.Conf, phi,m0); //D^-1 phi
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "Time = " << elapsed_secs << " s" << std::endl;

    c_matrix phi(Ntest, c_vector(Nagg, 0));
    for (int i = 0; i < Ntest; i++) {
        for (int j = 0; j < Nagg; j++) {
            phi[i][j] = 1.0 * RandomU1();
        }
    }




    std::cout << "################Bi-CGstab_DC inversion##############" << std::endl;
    AMG amg = AMG(GConf, Ns, Nt, Ntest, m0);
    amg.tv_init(1, 3); //test vectors intialization
    //Testing execution time
  

    c_matrix Dc(Ntest * Nagg, c_vector(Ntest * Nagg, 0));

    for (int col = 0; col < Ntest * Nagg; col++) {
        c_matrix v = amg.Pt_D_P(canonical_vector(col, Ntest, Nagg)); //column
        int count = 0;
        for (int i = 0; i < v.size(); i++) {
            for (int j = 0; j < v[i].size(); j++) {
                Dc[count][col] = v[i][j];
                count += 1;
            }
        }
    }

    //Flatten phi
    c_vector phi_flat(Ntest * Nagg, 0);
    int count = 0;
    for (int i = 0; i < phi.size(); i++) {
        for (int j = 0; j < phi[i].size(); j++) {
            phi_flat[count] = phi[i][j];
            count += 1;
        }
    }
    //PrintVector(Dc);
    save_matrix(Dc, "Dc.dat");
    save_vector(phi_flat, "phi.dat");
    c_vector Dc_phi_flat = Dc * phi_flat;
    std::cout << std::endl;
    std::cout << "------phi_flat------" << std::endl;
    std::cout << phi_flat[0] << "   " << phi_flat[1] << std::endl;
    std::cout << "------Dc phi_flat------" << std::endl;
    std::cout << Dc_phi_flat[0] << "   " << Dc_phi_flat[1] << std::endl;

    c_matrix Dc_phi = amg.Pt_D_P(phi);
    std::cout << "------Dc phi------" << std::endl;
    std::cout << Dc_phi[0][0] << " " << Dc_phi[0][1] << std::endl;

    std::cout << "------Dc^-1 phi------" << std::endl;
    begin = clock();
    c_matrix x1 = amg.bi_cgstab_Dc(GConf.Conf, phi, phi, m0, 100000, 1e-10, true);
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "Time = " << elapsed_secs << " s" << std::endl;

    std::cout << x1[0][0] << " " << x1[0][1] << std::endl;
    std::cout << "------Dc Dc^-1 phi------" << std::endl;
    c_matrix Dc_x = amg.Pt_D_P(x1);
    std::cout << Dc_x[0][0] << " " << Dc_x[0][1] << std::endl;

    return 0;
}

