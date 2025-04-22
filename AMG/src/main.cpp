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
    srand(time(0));
    initialize_matrices(); //Initialize gamma matrices, identity and unit vectors
    Coordinates(); //Vectorized coordinates
    periodic_boundary(); //Builds LeftPB and RightPB (periodic boundary for U_mu(n))
    std::cout << "******************* Two-grid method for the Dirac matrix in the Schwinger model *******************" << std::endl;
    std::cout << "Ns = " << Ns << " Nt = " << Nt << std::endl;
    std::cout << "block_x = " << block_x << " block_t = " << block_t << std::endl;
    std::cout << "x_elements = " << x_elements << " t_elements = " << t_elements << std::endl;
    Aggregates(); //Aggregates
    std::cout << "Number of test vectors = " << Ntest << std::endl;
    std::cout << "Dirac matrix dimension = " << (2 * Ntot) << std::endl;
    std::cout << "Dc dimension = " << Ntest * Nagg << std::endl;
    std::cout << "*****************************************************************************************************" << std::endl;


    GaugeConf GConf = GaugeConf(Ns, Nt);
    //GConf.initialization(); //Initialize a random gauge configuration
    int nu1 = 0, nu2 = 2;
    std::cout << "Pre-smoothing steps " << nu1 << " Post-smoothing steps " << nu2 << std::endl;

    double m0 = -0.6;
    double beta = 1;
    int n_conf = 10;
    std::cout << "m0 = " << m0 << " beta = " << beta << std::endl;
    std::vector<double> bi_cg_it(n_conf), bi_extime(n_conf);
    std::vector<double> multigrid_it(n_conf), amg_extime(n_conf);
    std::vector<double> gmres_it(n_conf), gmres_extime(n_conf);
	c_matrix x_bi, x_gmres, x_multigrid;
    for (int n = 0; n < n_conf; n++) {
        //***Open Conf File***//
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

		//***Generate random right hand side***//
        c_matrix PHI(Ntot, c_vector(2, 0));
        for (int i = 0; i < Ntot; i++) {
            for (int j = 0; j < 2; j++) {
                PHI[i][j] = 1.0 * RandomU1(); //Each element with norm=1
            }
        }

        //****Assemble the matrix****// This is only necessary for comparing with python ...
        /*
        c_matrix D(2 * Ntot, c_vector(2 * Ntot, 0));
        for (int col = 0; col < 2 * Ntot; col++) {
            c_matrix v = canonical_vector(col, Ntot, 2); 
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
        */

        std::cout << "########################################## Conf " << n << " ##########################################" << std::endl;
        
        std::cout << "--------------Bi-CGstab inversion--------------" << std::endl;
        clock_t begin = clock();
        x_bi = bi_cgstab(GConf.Conf, PHI, PHI, m0, 10000, 1e-10, false);//bi_cgstab(GConf.Conf, PHI, PHI, m0, 100000, 1e-10, false);
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        std::cout << "Bi-CGstab total computing time = " << elapsed_secs << " s" << std::endl;
        bi_cg_it[n] = it_count;
		bi_extime[n] = elapsed_secs;

        std::cout << "--------------GMRES inversion--------------" << std::endl;
        int m = 200;
        int restarts = 100;
        std::cout << "iterations per restart = " << m << " restarts = " << restarts << std::endl;

        begin = clock();
        x_gmres = gmres(GConf.Conf, PHI, PHI, m0, m, restarts, 1e-10, true);
        end = clock();
        
        elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        std::cout << "GMRES total computing time = " << elapsed_secs << " s" << std::endl;
		gmres_it[n] = it_count;
		gmres_extime[n] = elapsed_secs;

        std::cout << "--------------Two-grid inversion with GMRES as a smoother--------------" << std::endl;
        begin = clock();
        int rpc = 20; //iterations per GMRES cycle
        std::cout << "iterations per GMRES cycle for AMG = " << rpc << std::endl;
        AMG amg = AMG(GConf, Ns, Nt, Ntest, m0,nu1,nu2,rpc);
        amg.tv_init(1, 3); //test vectors intialization

        c_matrix x_ini(Ntot, c_vector(2, 0));
        for (int j = 0; j < Ntot; j++) {
			for (int k = 0; k < 2; k++) {
				x_ini[j][k] = RandomU1();
            }
		}
        
        c_matrix x0;
        x_multigrid = amg.TwoGrid(100,rpc, 1e-10, x_ini, PHI, true);
        end = clock();
        elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        std::cout << "Two-grid total computing time = " << elapsed_secs << " s" << std::endl;
        std::cout << "----------------------------------------------" << std::endl;
        multigrid_it[n] = it_count;
		amg_extime[n] = elapsed_secs;
        std::cout << "Two-grid " << std::setprecision(10) << x_multigrid[0][0] << std::endl;
        std::cout << "BiCGstab " << std::setprecision(10) << x_bi[0][0] << std::endl;
        std::cout << "GMRES " << std::setprecision(10) << x_gmres[0][0] << std::endl;
    }

    std::cout << "Mean number of iterations for bi-cg " << mean(bi_cg_it) << " +- " << Jackknife_error(bi_cg_it, 2) << std::endl;
    std::cout << "Mean execution time for bi-cg " << mean(bi_extime) << " +- " << Jackknife_error(bi_extime, 2) << std::endl;
    std::cout << " " << std::endl;

    std::cout << "Mean number of iterations for gmres " << mean(gmres_it) << " +- " << Jackknife_error(gmres_it, 2) << std::endl;
    std::cout << "Mean execution time for gmres " << mean(gmres_extime) << " +- " << Jackknife_error(gmres_extime, 2) << std::endl;
    std::cout << " " << std::endl;

    std::cout << "Mean number of iterations for two-grid " << mean(multigrid_it) << " +- " << Jackknife_error(multigrid_it, 2) << std::endl;
    std::cout << "Mean execution time for two-grid " << mean(amg_extime) << " +- " << Jackknife_error(amg_extime, 2) << std::endl;
    std::cout << " " << std::endl;

    return 0;
}