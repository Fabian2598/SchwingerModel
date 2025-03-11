#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "amg.h"
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
    
    double m0 = -0.6; 
    double beta = 1;
    int n_conf = 50;
    std::cout << "m0 " << m0 << " beta " << beta << std::endl;
    std::vector<double> bi_cg_it(n_conf);
    std::vector<double> multigrid_it(n_conf);
    int non_convergent_conf = 0;
    for (int n = 29; n < 31; n++) {
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

        srand(time(0));
        c_matrix PHI(Ntot, c_vector(2, 0));
        for (int i = 0; i < Ntot; i++) {
            for (int j = 0; j < 2; j++) {
                PHI[i][j] = 1.0 * RandomU1();
            }
        }
        
        std::cout << "##### Conf " << n << "#####" << std::endl;
        std::cout << "Bi-CGstab inversion" << std::endl;
        c_matrix X = bi_cgstab(GConf.Conf, PHI, PHI, m0, 100000, 1e-10, true);
        std::cout << "----------------" << std::endl;   
        bi_cg_it[n] = it_count;
      

        AMG amg = AMG(GConf, Ns, Nt, Ntest, m0);
        amg.tv_init(1, 2); //test vectors intialization
        std::cout << "Two-grid inversion" << std::endl;
        c_matrix x0;
        c_matrix zero = c_matrix(Ntot, c_vector(2, 0));
        x0 = amg.TwoGrid(nu1, nu2, 300, 1e-10, PHI, PHI, true);
        std::cout << "--------------------" << std::endl;
        multigrid_it[n] = it_count;
        if (it_count == 300) { 
            non_convergent_conf += 1; 
            std::cout << "***** Conf " << n << "  ****" << std::endl;
            std::cout << "***** did not converge  ****" << std::endl;
            std::cout << "*****----------****" << std::endl;
            std::cout << "*******---------- ******" << std::endl;
            std::cout << "*******---------- ******" << std::endl;
        }
    }
    std::cout << "Number of confs that did not converge " << non_convergent_conf << std::endl;
    std::cout << "Mean number of iterations for bi-cg " << mean(bi_cg_it) << " +- " << Jackknife(bi_cg_it, { 2,25,10,5 }) << std::endl;
    std::cout << "Mean number of iterations for two-grid " << mean(multigrid_it) << " +- " << Jackknife(multigrid_it, { 2,25,10,5 }) << std::endl;
    

    return 0;
}