#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "sap.h"
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
    std::cout << "Lattice dimension = " << (Ns * Nt) << std::endl;
    std::cout << "Number of entries of the Dirac matrix = (" << (2 * Ns * Nt) << ")^2 = " << (2 * Ns * Nt) * (2 * Ns * Nt) << std::endl;
    std::cout << "-----------------------------------" << std::endl;
    std::cout << "| Lattice blocking for the aggregates" << std::endl;
    std::cout << "| block_x = " << block_x << " block_t = " << block_t << std::endl;
    std::cout << "| x_elements = " << x_elements << " t_elements = " << t_elements << std::endl;
    std::cout << "| Each aggregate has x_elements * t_elements = " <<  x_elements * t_elements << " elements" << std::endl;
    std::cout << "| Number of aggregates = " << Nagg << std::endl;
    Aggregates(); //Aggregates
    std::cout << "-----------------------------------" << std::endl;
    std::cout << "| Variable blocking for SAP" << std::endl;
    std::cout << "| sap_block_x = " << sap_block_x << " sap_block_t = " << sap_block_t << std::endl;
    std::cout << "| Lattice sites in the x direction = " << sap_x_elements<< " and in the t direction = " << sap_t_elements<< std::endl;
    std::cout << "| Each Schwarz block has " <<  sap_lattice_sites_per_block << " lattice points and " << sap_variables_per_block << " variables" << std::endl;
    std::cout << "| D restricted to each block has (" << 2 * sap_lattice_sites_per_block << ")^2 = " << sap_variables_per_block*sap_variables_per_block << " entries" << std::endl;
    std::cout << "| Number of Schwarz blocks = " << N_sap_blocks << std::endl;
    CheckBlocks(); //Check blocks dimensions
    SchwarzBlocks(); //Builds the blocks for the Schwarz alternating method
     
    std::cout << "Number of test vectors = " << Ntest << std::endl;
    std::cout << "Dc dimension = " << Ntest * Nagg << std::endl;
    std::cout << "*****************************************************************************************************" << std::endl;


    GaugeConf GConf = GaugeConf(Ns, Nt);
    GConf.initialization(); //Initialize a random gauge configuration
    int nu1 = 0, nu2 = 2;
    std::cout << "Pre-smoothing steps " << nu1 << " Post-smoothing steps " << nu2 << std::endl;

    double m0 = -0.1;
    double beta = 1;
    int n_conf = 10;
    std::cout << "m0 = " << m0 << " beta = " << beta << std::endl;
    
   

    //Assemble matrix for checking Schwarz decomposition. Only try with small latices.
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

    std::cout << "#### Dirac Matrix ####" << std::endl;
    PrintVector(D);
    */



    
 
    c_matrix rhs(Ntot, c_vector(2, 0)); //random right hand side 
    c_matrix x_B(Ntot, c_vector(2, 0)); //solution vector on the original lattice
    for(int i = 0; i < sap_variables_per_block; i++) {
        rhs[i][0] = RandomU1();
        rhs[i][1] = RandomU1();
    }
    //Default values in variables.cpp
    sap_gmres_restart_length = 5; //GMRES restart length for the Schwarz blocks. Set to 20 by default
    sap_gmres_restarts = 10; //GMRES iterations for the Schwarz blocks. Set to 10 by default.
    sap_gmres_tolerance = 1e-10; //GMRES tolerance for the Schwarz blocks
    for(int block = 0; block<N_sap_blocks; block++){
        x_B = I_D_B_1_It(GConf.Conf, rhs, m0,block);
        std::cout << "right-hand side rhs[0][0] = " << rhs[0][0] << "  rhs[0][1] = " << rhs[0][1]  << std::endl;
        std::cout << "solution" << std::endl;
        for(int i = 0; i < Ntot; i++) {
            std::cout << "x_B[" << i << "][0] = " << x_B[i][0] << "  ";
            std::cout << "x_B[" << i << "][1] = " << x_B[i][1];
            std::cout << std::endl;
        }
        std::cout << "**********************" << std::endl;
    }
    

    return 0;
}