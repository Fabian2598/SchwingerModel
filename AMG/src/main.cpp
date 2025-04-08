#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "amg.h"
#include "gmres.h"
#include "statistics.h"

//Formats decimal numbers
static std::string format(const double& number) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << number;
    std::string str = oss.str();
    str.erase(str.find('.'), 1); //Removes decimal dot
    return str;
}

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
    //GConf.initialization(); //Initialize the gauge configuration
    int nu1 = 0, nu2 = 2;
    std::cout << "Pre-smoothing steps " << nu1 << " Post-smoothing steps " << nu2 << std::endl;

    double m0 = -0.65;
    double beta = 1;
    int n_conf = 3;
    std::cout << "m0 " << m0 << " beta " << beta << std::endl;
    std::vector<double> bi_cg_it(n_conf);
    std::vector<double> multigrid_it(n_conf);
    for (int n = 0; n < n_conf; n++) {
        char NameData[500];
        sprintf(NameData, "../confs/2D_U1_Ns%d_Nt%d_b%s_m%s_%d.txt", Ns, Nt, format(beta).c_str(), format(m0).c_str(), n);
        std::ifstream infile(NameData);
        if (!infile) {
            std::cerr << "File not found" << std::endl;
            return 1;
        }
        int x, t, mu;
        double re, im;

        while (infile >> x >> t >> mu >> re >> im) {
            GConf.Conf[Coords[x][t]][mu] = c_double(re, im);
        }
        infile.close();

        c_matrix PHI(Ntot, c_vector(2, 0));
        for (int i = 0; i < Ntot; i++) {
            for (int j = 0; j < 2; j++) {
                PHI[i][j] = 1.0 * RandomU1();
            }
        }

        //Store the matrix//
        c_matrix D(2 * Ntot, c_vector(2 * Ntot, 0));
        for (int col = 0; col < 2 * Ntot; col++) {
            c_matrix v = canonical_vector(col, Ntot, 2);  //column
            c_matrix Dv = D_phi(GConf.Conf, v, m0);
            int count = 0;
            for (int i = 0; i < Dv.size(); i++) {
                for (int j = 0; j < Dv[i].size(); j++) {
                    D[count][col] = Dv[i][j];
                    count += 1;
                }
            }
        }
        char NameD[500], NamePhi[500];
        sprintf(NameD, "D%d.dat", n);  sprintf(NamePhi, "Phi%d.dat", n);
        save_matrix(D, NameD);
		save_matrix(PHI, NamePhi);

        std::cout << "##### Conf " << n << "#####" << std::endl;
        /*
        std::cout << "--Bi-CGstab inversion--" << std::endl;
        c_matrix X = bi_cgstab(GConf.Conf, PHI, PHI, m0, 100000, 1e-10, true);
        bi_cg_it[n] = it_count;
        std::cout << "-------Bi inversion done-------" << std::endl;
        */
        std::cout << "-------GMRES inversion-------" << std::endl;
        
        c_matrix X_GMRES = gmres(GConf.Conf, PHI, PHI, m0, 80, 50, 1e-10, true);
        std::cout << "GMRES" << X_GMRES[0][0] << std::endl;

        //std::cout << "Comparison : bi-cg" << X[0][0] << "    GMRES" << X_GMRES[0][0] << std::endl;
        //std::cout << "-------Kaczmarz inversion done-------" << std::endl;
        //std::cout << "Comparison : bi-cg" << X[0][0] << "    Kaczmarz" << X_K[0][0] << std::endl;

        //std::cout << "-------Two-grid inversion with Kaczmarz as a smoother-------" << std::endl;
        //AMG amg = AMG(GConf, Ns, Nt, Ntest, m0,nu1,nu2);
        //amg.tv_init(1, 6); //test vectors intialization
        
        //c_matrix x0;
        //x0 = amg.TwoGrid(200, 1e-10, PHI, PHI, true);
        //multigrid_it[n] = it_count;
        
    }

    std::cout << "Mean number of iterations for bi-cg " << mean(bi_cg_it) << " +- " << Jackknife_error(bi_cg_it, 5) << std::endl;
    std::cout << "Mean number of iterations for two-grid " << mean(multigrid_it) << " +- " << Jackknife_error(multigrid_it, 5) << std::endl;


    return 0;
}