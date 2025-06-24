#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "bi_cgstab.h"
#include "statistics.h"
#include "fgmres.h"
#include "mpi.h"

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
    MPI_Init(&argc, &argv);
    int rank, size; 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //srand(19);

    srand(time(0));
    initialize_matrices(); //Initialize gamma matrices, identity and unit vectors
    Coordinates(); //Builds array with coordinates of the lattice points x * Nt + t 
    periodic_boundary(); //Builds LeftPB and RightPB (periodic boundary for U_mu(n))
    double m0 = -0.65;
    //double m0 = -0.3;
    if (rank == 0){
        std::cout << "******************* Two-grid method for the Dirac matrix in the Schwinger model *******************" << std::endl;
        std::cout << "Nx = " << Nx << " Nt = " << Nt << std::endl;
        std::cout << "Lattice dimension = " << (Nx * Nt) << std::endl;
        std::cout << "Number of entries of the Dirac matrix = (" << (2 * Nx * Nt) << ")^2 = " << (2 * Nx * Nt) * (2 * Nx * Nt) << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "| Lattice blocking for the aggregates" << std::endl;
        std::cout << "| block_x = " << block_x << " block_t = " << block_t << std::endl;
        std::cout << "| x_elements = " << x_elements << " t_elements = " << t_elements << std::endl;
        std::cout << "| Each aggregate has x_elements * t_elements = " <<  x_elements * t_elements << " elements" << std::endl;
        std::cout << "| Number of aggregates = " << AMGV::Nagg << std::endl;
    } 
    Aggregates(); //build aggregates
    CheckAggregates();
    if (rank == 0){
        std::cout << "----------------------------------" << std::endl;
        std::cout << "| Variable blocking for SAP" << std::endl;
        std::cout << "| sap_block_x = " << sap_block_x << " sap_block_t = " << sap_block_t << std::endl;
        std::cout << "| Lattice sites in the x direction = " << sap_x_elements<< " and in the t direction = " << sap_t_elements<< std::endl;
        std::cout << "| Each Schwarz block has " <<  sap_lattice_sites_per_block << " lattice points and " << sap_variables_per_block << " variables" << std::endl;
        std::cout << "| D restricted to each block has (" << 2 * sap_lattice_sites_per_block << ")^2 = " << sap_variables_per_block*sap_variables_per_block << " entries" << std::endl;
        std::cout << "| Number of Schwarz blocks = " << N_sap_blocks << std::endl;
        std::cout << "| Red/Black blocks = " << sap_coloring_blocks << std::endl;
        
    }
    SchwarzBlocks(); //Builds the blocks for the Schwarz alternating method
    CheckBlocks(); //Check blocks dimensions
    AMGV::gmres_restarts_coarse_level = 20; 
    AMGV::gmres_restart_length_coarse_level = 140; //GMRES restart length for the coarse level
    AMGV::gmres_tol_coarse_level = 1e-3; //GMRES tolerance for the coarse level
    AMGV::nu1 = 0; //Pre-smoothing iterations
    AMGV::nu2 = 2; //Post-smoothing iterations
    AMGV::Nit = 3; //Number of iterations for improving the interpolator
    if (rank == 0){
        std::cout << "----------------------------------" << std::endl;
        std::cout << "|Other parameters" << std::endl;
        std::cout << "|Number of test vectors = " << AMGV::Ntest << std::endl;
        std::cout << "|nu1 (pre-smoothing) = " << AMGV::nu1 << " nu2 (post-smoothing) = " << AMGV::nu2 << std::endl;
        std::cout << "|Number of iterations for improving the interpolator = " << AMGV::Nit << std::endl;
        std::cout << "|Dc dimension = " << AMGV::Ntest * AMGV::Nagg << std::endl;
        std::cout << "|m0 = " << m0 << std::endl;
        std::cout << "|Number of SAP iterations to smooth test vectors = " << AMGV::SAP_test_vectors_iterations << std::endl; 
        std::cout << "|Restart length of GMRES at the coarse level = " << AMGV::gmres_restart_length_coarse_level << std::endl;
        std::cout << "|Restarts of GMRES at the coarse level = " << AMGV::gmres_restarts_coarse_level << std::endl;
        std::cout << "|GMRES tolerance for the coarse level solution = " << AMGV::gmres_tol_coarse_level << std::endl;
        std::cout << "*****************************************************************************************************" << std::endl;
    }
    
    GaugeConf GConf = GaugeConf(Nx, Nt);
    GConf.initialization(); //Initialize a random gauge configuration

    //Default values in variables.cpp
    sap_gmres_restart_length = 20; //GMRES restart length for the Schwarz blocks. Set to 20 by default
    sap_gmres_restarts = 10; //GMRES iterations for the Schwarz blocks. Set to 10 by default.
    sap_gmres_tolerance = 1e-3; //GMRES tolerance for the Schwarz blocks
    sap_tolerance = 1e-10; //Tolerance for the SAP method
    spinor rhs(Ntot, c_vector(2, 0)); //random right hand side 
    spinor x(Ntot, c_vector(2, 0)); //solution vector 
    for(int i = 0; i < Ntot; i++) {
        rhs[i][0] = RandomU1();
        rhs[i][1] = RandomU1();
    }

    clock_t start, end;
    double elapsed_time;
    double startT, endT;

   
   int gmres_restarts = 50, gmres_restart_length = 20; //for fgmres and gmres
    if (rank == 0){
        std::cout << "--------------Bi-CGstab inversion--------------" << std::endl;
        start = clock();
        spinor x0(Ntot, c_vector(2, 0)); //Initial guess
        spinor x_bi = bi_cgstab(&D_phi,Ntot,2,GConf.Conf, rhs, x0, m0, 100000, 1e-10, true);
        end = clock();
        elapsed_time = double(end - start) / CLOCKS_PER_SEC;
        std::cout << "Elapsed time for Bi-CGstab = " << elapsed_time << " seconds" << std::endl;    
    }

    MPI_Barrier(MPI_COMM_WORLD);
        std::cout << "--------------Flexible GMRES with SAP preconditioning --------------" << std::endl;   
    startT = MPI_Wtime();
    spinor xfgmres = fgmresSAP(GConf.Conf, rhs, x, m0, gmres_restart_length,gmres_restarts, 1e-10, true);
    endT = MPI_Wtime();
    printf("[MPI process %d] time elapsed during the job: %.4fs.\n", rank, endT - startT);

    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "--------------Flexible GMRES with AMG preconditioning--------------" << std::endl;
    startT = MPI_Wtime();
    spinor xAMG = fgmresAMG(GConf.Conf, rhs, x, m0, gmres_restart_length,gmres_restarts, 1e-10, true);
    endT = MPI_Wtime();
    printf("[MPI process %d] time elapsed during the job: %.4fs.\n", rank, endT - startT);
    
    MPI_Finalize();
    return 0;
}