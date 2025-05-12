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
#include "mpi.h"

//Formats decimal numbers
static std::string format(const double& number) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << number;
    std::string str = oss.str();
    str.erase(str.find('.'), 1); //Removes decimal dot 
    return str;
}

int main(int argc, char **argv) {


    MPI_Init(&argc, &argv);
    int rank, size; 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //srand(time(0));
    srand(0);
    initialize_matrices(); //Initialize gamma matrices, identity and unit vectors
    Coordinates(); //Vectorized coordinates
    periodic_boundary(); //Builds LeftPB and RightPB (periodic boundary for U_mu(n))
    double m0 = -0.2;
    double beta = 1;
    if (rank == 0){
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
    } 
    Aggregates(); //Aggregates
    if (rank == 0){
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "| Variable blocking for SAP" << std::endl;
        std::cout << "| sap_block_x = " << sap_block_x << " sap_block_t = " << sap_block_t << std::endl;
        std::cout << "| Lattice sites in the x direction = " << sap_x_elements<< " and in the t direction = " << sap_t_elements<< std::endl;
        std::cout << "| Each Schwarz block has " <<  sap_lattice_sites_per_block << " lattice points and " << sap_variables_per_block << " variables" << std::endl;
        std::cout << "| D restricted to each block has (" << 2 * sap_lattice_sites_per_block << ")^2 = " << sap_variables_per_block*sap_variables_per_block << " entries" << std::endl;
        std::cout << "| Number of Schwarz blocks = " << N_sap_blocks << std::endl;
        std::cout << "| Red/Black blocks = " << sap_coloring_blocks << std::endl;
        CheckBlocks(); //Check blocks dimensions
    }
    
    SchwarzBlocks(); //Builds the blocks for the Schwarz alternating method
     
    if (rank == 0){
        std::cout << "Number of test vectors = " << Ntest << std::endl;
        std::cout << "Dc dimension = " << Ntest * Nagg << std::endl;
        std::cout << "*****************************************************************************************************" << std::endl;
        std::cout << "m0 = " << m0 << " beta = " << beta << std::endl;
        std::cout << "-------------------------" << std::endl;    
    }
    
    GaugeConf GConf = GaugeConf(Ns, Nt);
    GConf.initialization(); //Initialize a random gauge configuration


     //Default values in variables.cpp
    sap_gmres_restart_length = 20; //GMRES restart length for the Schwarz blocks. Set to 20 by default
    sap_gmres_restarts = 10; //GMRES iterations for the Schwarz blocks. Set to 10 by default.
    sap_gmres_tolerance = 1e-6; //GMRES tolerance for the Schwarz blocks
    sap_tolerance = 1e-10; //Tolerance for the SAP method
    c_matrix rhs(Ntot, c_vector(2, 0)); //random right hand side 
    c_matrix x(Ntot, c_vector(2, 0)); //solution vector 
    for(int i = 0; i < Ntot; i++) {
        rhs[i][0] = RandomU1();
        rhs[i][1] = RandomU1();
    }
    c_matrix D_x(Ntot, c_vector(2, 0)); //This should be the rhs again
    int nu_sap = 20;

    clock_t start, end;
    double elapsed_time;
    if (rank == 0){
        std::cout << "SAP Serial version" << std::endl;
        start = clock();
        SAP(GConf.Conf, rhs, x, m0, nu_sap);
        end = clock();
        elapsed_time = double(end - start) / CLOCKS_PER_SEC;
        std::cout << "Elapsed time for sequential SAP = " << elapsed_time << " seconds" << std::endl;
    }
    
    /*
    D_x = D_phi(GConf.Conf, xv1, m0);
    std::cout << "-------------------------" << std::endl;
    for(int i = 0; i < 10; i++) {
        std::cout << "rhs[" << i << "][0] = " << rhs[i][0] << " | ";
        std::cout  << D_x[i][0]  << " = D x_B[" << i << "][0]" << std::endl;
        std::cout << "rhs[" << i << "][1] = " << rhs[i][1] << " | ";
        std::cout  << D_x[i][1]  << " = D x_B[" << i << "][0]" << std::endl;
        std::cout << std::endl;
    }
    */

   double startT, endT;

   MPI_Barrier(MPI_COMM_WORLD);
   startT = MPI_Wtime();
   SAP_parallel(GConf.Conf, rhs, x, m0, nu_sap,1);
   MPI_Barrier(MPI_COMM_WORLD);
   endT = MPI_Wtime();

   printf("[MPI process %d] time elapsed during the job: %.4fs.\n", rank, endT - startT);
   MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0){
        std::cout << "**********************" << std::endl;
        std::cout << "--------------Bi-CGstab inversion--------------" << std::endl;
        start = clock();
        c_matrix x_bi = bi_cgstab(GConf.Conf, rhs, rhs, m0, 10000, 1e-10, false);
        end = clock();
        elapsed_time = double(end - start) / CLOCKS_PER_SEC;
        std::cout << "Elapsed time for Bi-CGstab = " << elapsed_time << " seconds" << std::endl;
        for(int i = 0; i < 10; i++) {
            std::cout << "x[" << i << "][0] = " << x[i][0] << " | ";
            std::cout  << x_bi[i][0]  << " = x_bi[" << i << "][0]" << std::endl;
            std::cout << "x[" << i << "][1] = " << x[i][1] << " | ";
            std::cout  << x_bi[i][1]  << " = x_bi[" << i << "][1]" << std::endl;
            std::cout << std::endl;
        }
    }
    

    MPI_Finalize();
    return 0;
}