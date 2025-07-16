#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "bi_cgstab.h"
#include "fgmres.h"


//Formats decimal numbers
//For opening file with confs 
static std::string format(const double& number) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << number;
    std::string str = oss.str();
    str.erase(str.find('.'), 1); //Removes decimal dot 
    return str;
}

int main(int argc, char **argv) {

    using namespace SAPV;

    
    srand(19);

    //srand(time(0));
    
    Coordinates(); //Builds array with coordinates of the lattice points x * Nt + t 
    periodic_boundary(); //Builds LeftPB and RightPB (periodic boundary for U_mu(n))
    double m0 = -0.18840579710144945; 

    //Default values in variables.cpp
    sap_gmres_restart_length = 50; //GMRES restart length for the Schwarz blocks. Set to 20 by default
    sap_gmres_restarts = 1; //GMRES iterations for the Schwarz blocks. Set to 10 by default.
    sap_gmres_tolerance = 1e-3; //GMRES tolerance for the Schwarz blocks
    sap_tolerance = 1e-10; //Tolerance for the SAP method
    sap_blocks_per_proc = 1; //Number of blocks per process for the parallel SAP method

    AMGV::SAP_test_vectors_iterations = 1;
    AMGV::gmres_restarts_coarse_level = 10; 
    AMGV::gmres_restart_length_coarse_level = 100; //GMRES restart length for the coarse level
    AMGV::gmres_tol_coarse_level = 0.1; //GMRES tolerance for the coarse level
    AMGV::nu1 = 0; //Pre-smoothing iterations
    AMGV::nu2 = 2; //Post-smoothing iterations
    AMGV::Nit = 3; //Number of iterations for improving the interpolator

    FGMRESV::fgmres_tolerance = 1e-10; //Tolerance for FGMRES
    FGMRESV::fgmres_restart_length = 20; //Restart length for FGMRES
    FGMRESV::fgmres_restarts = 50; //Number of restarts for FGMRES

        std::cout << "********** Critical masses ********** " << std::endl;
        std::cout << "* beta = 2, m0_crit = -0.1968(9)"  << std::endl;
        std::cout << "* beta = 3, m0_crit = 0.1351(2)" << std::endl;
        std::cout << "* beta = 4, m0_crit = 0.1033(1)" << std::endl;
        std::cout << "* beta = 5, m0_crit = 0.0840(1)" << std::endl;
        std::cout << "************************************* " << std::endl;

        std::cout << "******************* Two-grid method for the Dirac matrix in the Schwinger model *******************" << std::endl;
        std::cout << " Nx = " << Nx << " Nt = " << Nt << std::endl;
        std::cout << " Lattice dimension = " << (Nx * Nt) << std::endl;
        std::cout << " Number of entries of the Dirac matrix = (" << (2 * Nx * Nt) << ")^2 = " << (2 * Nx * Nt) * (2 * Nx * Nt) << std::endl;
        std::cout << " Bare mass parameter m0 = " << m0 << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        std::cout << " Lattice blocking for the aggregates" << std::endl;
        std::cout << "| block_x = " << block_x << " block_t = " << block_t << std::endl;
        std::cout << "| x_elements = " << x_elements << " t_elements = " << t_elements << std::endl;
        std::cout << "| Each aggregate has x_elements * t_elements = " <<  x_elements * t_elements << " elements" << std::endl;
        std::cout << "| Number of aggregates = " << AMGV::Nagg << std::endl;
    
    Aggregates(); //build aggregates
    CheckAggregates();

        std::cout << "----------------------------------" << std::endl;
        std::cout << " Variable blocking for SAP" << std::endl;
        std::cout << "| sap_block_x = " << sap_block_x << " sap_block_t = " << sap_block_t << std::endl;
        std::cout << "| Lattice sites in the x direction = " << sap_x_elements<< " and in the t direction = " << sap_t_elements<< std::endl;
        std::cout << "| Each Schwarz block has " <<  sap_lattice_sites_per_block << " lattice points and " << sap_variables_per_block << " variables" << std::endl;
        std::cout << "| D restricted to each block has (" << 2 * sap_lattice_sites_per_block << ")^2 = " << sap_variables_per_block*sap_variables_per_block << " entries" << std::endl;
        std::cout << "| Number of Schwarz blocks = " << N_sap_blocks << std::endl;
        std::cout << "| Red/Black blocks = " << sap_coloring_blocks << std::endl;
        std::cout << "| Number of processes = " << 1 << std::endl;
        std::cout << "| Number of blocks per process = " << sap_blocks_per_proc << std::endl;
        std::cout << "| GMRES restart length for SAP blocks = " << sap_gmres_restart_length << std::endl;
        std::cout << "| GMRES iterations for SAP blocks = " << sap_gmres_restarts << std::endl;
        std::cout << "| GMRES tolerance for SAP blocks = " << sap_gmres_tolerance << std::endl;
    
    SchwarzBlocks(); //Builds the blocks for the Schwarz alternating method
    CheckBlocks(); //Check blocks dimensions

        std::cout << "----------------------------------" << std::endl;
        std::cout << " Two-grid parameters" << std::endl;
        std::cout << "| Number of test vectors = " << AMGV::Ntest << std::endl;
        std::cout << "| nu1 (pre-smoothing) = " << AMGV::nu1 << " nu2 (post-smoothing) = " << AMGV::nu2 << std::endl;
        std::cout << "| Number of iterations for improving the interpolator = " << AMGV::Nit << std::endl;
        std::cout << "| Coarse grid matrix (Dc) dimension = " << AMGV::Ntest * AMGV::Nagg << std::endl;
        std::cout << "| Number of SAP iterations to smooth test vectors = " << AMGV::SAP_test_vectors_iterations << std::endl; 
        std::cout << "| Restart length of GMRES at the coarse level = " << AMGV::gmres_restart_length_coarse_level << std::endl;
        std::cout << "| Restarts of GMRES at the coarse level = " << AMGV::gmres_restarts_coarse_level << std::endl;
        std::cout << "| GMRES tolerance for the coarse level solution = " << AMGV::gmres_tol_coarse_level << std::endl;
        std::cout << "----------------------------------" << std::endl;
        std::cout << " FGMRES with AMG preconditioning parameters" << std::endl;
        std::cout << "| FGMRES restart length = " << FGMRESV::fgmres_restart_length << std::endl;
        std::cout << "| FGMRES restarts = " << FGMRESV::fgmres_restarts << std::endl;
        std::cout << "| FGMRES tolerance = " << FGMRESV::fgmres_tolerance << std::endl;
        std::cout << "*****************************************************************************************************" << std::endl;
    
    
    GaugeConf GConf = GaugeConf(Nx, Nt);
    GConf.initialize(); //Initialize a random gauge configuration

    //Open conf from file//
    
    {
        double beta = 2;
        int nconf = 1;
        std::ostringstream NameData;
        NameData << "../confs/b" << beta << "_" << LV::Nx << "x" << LV::Nt << "/m-018/2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt << "_b" << 
        format(beta).c_str() << "_m" << format(m0).c_str() << "_" << nconf << ".ctxt";
        //std::cout << "Reading conf from file: " << NameData.str() << std::endl;
        std::ifstream infile(NameData.str());
        if (!infile) {
           std::cerr << "File not found" << std::endl;
       }
        int x, t, mu;
        double re, im;
        
        c_matrix CONF(Ntot,c_vector(2,0)); 
        while (infile >> x >> t >> mu >> re >> im) {
            CONF[Coords[x][t]][mu] = c_double(re, im); 
        }
        GConf.setGconf(CONF);
        infile.close();
    }


    spinor rhs(Ntot, c_vector(2, 0)); //random right hand side 
    spinor x(Ntot, c_vector(2, 0)); //solution vector 
    //Random right hand side
    for(int i = 0; i <Ntot; i++) {
        rhs[i][0] = RandomU1();
        rhs[i][1] = RandomU1();
    }

    clock_t start, end;
    double elapsed_time;
    double startT, endT;

            
    //Bi-cgstab inversion for comparison
    std::cout << "--------------Bi-CGstab inversion--------------" << std::endl;
    start = clock();
    spinor x0(Ntot, c_vector(2, 0)); //Initial guess
    spinor x_bi = bi_cgstab(&D_phi,Ntot,2,GConf.Conf, rhs, x0, m0, 100000, 1e-10, true);
    end = clock();
    elapsed_time = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed time for Bi-CGstab = " << elapsed_time << " seconds" << std::endl;    
    
    std::cout << "\n\n";
    

    /*
    std::cout << "--------------Flexible GMRES with SAP preconditioning --------------" << std::endl; 
    start = clock();
    spinor xfgmres = fgmresSAP(GConf.Conf, rhs, x, m0, FGMRESV::fgmres_restart_length,FGMRESV::fgmres_restarts, FGMRESV::fgmres_tolerance , true);
    end = clock();
    elapsed_time = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed time for FGMRES with SAP = " << elapsed_time << " seconds" << std::endl;

    std::cout << "\n\n";
    */

    
    std::cout << "--------------Flexible GMRES with AMG preconditioning--------------" << std::endl;
    start = clock();
    spinor xAMG = fgmresAMG(GConf.Conf, rhs, x, m0, FGMRESV::fgmres_restart_length,FGMRESV::fgmres_restarts, FGMRESV::fgmres_tolerance , true);
    end = clock();
    elapsed_time = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed time for FGMRES with AMG = " << elapsed_time << " seconds" << std::endl;
    std::cout << "coarse time = " << coarse_time << " seconds" << std::endl;
    std::cout << "smooth time = " << smooth_time << " seconds" << std::endl;
    std::cout << "SAP time = " << SAP_time << " seconds" << std::endl;
    



    return 0;
}





    /*
    spinor RHS(sap_lattice_sites_per_block, c_vector(2, 0)); //random right hand side 

    for(int i = 0; i <sap_lattice_sites_per_block; i++) {
        RHS[i][0] = RandomU1();
        RHS[i][1] = RandomU1();
    }

    int block = 0;
    spinor xv1 = D_B(GConf.Conf, RHS, m0,block);
    for(int n = 0; n < 2; n++) {
        std::cout << "xv1[" << n << "] = (" << xv1[n][0] << ", " << xv1[n][1] << ")" << std::endl;
    }

    spinor xv2 = local_D(GConf.Conf, RHS, m0,block);
    for(int n = 0; n < 2; n++) {
        std::cout << "xv2[" << n << "] = (" << xv2[n][0] << ", " << xv2[n][1] << ")" << std::endl;
   }

   //check if both are the same
   for(int n = 0; n < sap_lattice_sites_per_block; n++) {
        if (std::abs(xv1[n][0] - xv2[n][0]) > 1e-10 || std::abs(xv1[n][1] - xv2[n][1]) > 1e-10) {
            std::cout << "Error: xv1 and xv2 are not the same!" << std::endl;
            return 0;
        } 
    }
     std::cout << "xv1 and xv2 are the same" << std::endl;
     */