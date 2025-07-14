#include "variables.h"

typedef std::complex<double> c_double;
double coarse_time = 0.0; //Time spent in the coarse grid solver
double smooth_time = 0.0; //Time spent in the smoother
int nonzero = 0; //Count the number of non-zero elements in the coarse grid operator

std::vector<std::vector<int>>Coords = std::vector<std::vector<int>>(LV::Nx, std::vector<int>(LV::Nt, 0));
void Coordinates() {
	for (int x = 0; x < LV::Nx; x++) {
		for (int t = 0; t < LV::Nt; t++) {
			Coords[x][t] = x * LV::Nt+ t;
		}
	}
}


//Aggregates A_j_0 = L_j x {0}, A_j_1 = L_j x {1}
std::vector<std::vector<int>>Agg = std::vector<std::vector<int>>(2*LV::block_x*LV::block_t, std::vector<int>(LV::x_elements*LV::t_elements, 0));
std::vector<int> XCoord = std::vector<int>(2*LV::Ntot, 0);
std::vector<int> TCoord = std::vector<int>(2*LV::Ntot, 0);
std::vector<int> SCoord = std::vector<int>(2*LV::Ntot, 0);

//--Coordinates of the neighbors to avoid recomputing them each time the operator D is called--//
std::vector<std::vector<int>>LeftPB = std::vector<std::vector<int>>(LV::Ntot, std::vector<int>(2,0)); //LeftPB[x][t][mu]
std::vector<std::vector<int>>RightPB = std::vector<std::vector<int>>(LV::Ntot, std::vector<int>(2,0)); //RightPB[x][t][mu]



std::vector<std::vector<c_double>>SignL =std::vector<std::vector<c_double>>(LV::Ntot, std::vector<c_double>(2,0)); //SignL[x][t][mu]
std::vector<std::vector<c_double>>SignR = std::vector<std::vector<c_double>>(LV::Ntot, std::vector<c_double>(2,0)); ////SignR[x][t][mu]

//--SAP blocks--//
std::vector<std::vector<int>>SAP_Blocks = std::vector<std::vector<int>>(SAPV::sap_block_x*SAPV::sap_block_t, 
    std::vector<int>(SAPV::sap_x_elements*SAPV::sap_t_elements, 0));
std::vector<int> SAP_RedBlocks = std::vector<int>(SAPV::sap_coloring_blocks, 0); //Red blocks
std::vector<int> SAP_BlackBlocks = std::vector<int>(SAPV::sap_coloring_blocks, 0); //Black blocks

//--Coarse grid operator--//
std::vector<std::vector<std::complex<double>>> DcMatrix = std::vector<std::vector<std::complex<double>>> (AMGV::Ntest*AMGV::Nagg, 
   std::vector<std::complex<double> >(AMGV::Ntest*AMGV::Nagg,0));


namespace SAPV {
    bool schwarz_blocks = false; //Schwarz blocks are not initialized by default
    int sap_gmres_restart_length = 20; //GMRES restart length for the Schwarz blocks. Set to 20 by default
    int sap_gmres_restarts = 10; //GMRES iterations for the Schwarz blocks. Set to 10 by default.
    double sap_gmres_tolerance = 1e-10; //GMRES tolerance for the Schwarz blocks
    double sap_tolerance = 1e-10; //Tolerance for the SAP method
    int sap_blocks_per_proc = 1; //Number of blocks per process for the parallel SAP method
}

namespace AMGV {
    int SAP_test_vectors_iterations = 2; //Number of SAP iterations to smooth test vectors
    bool aggregates_initialized = false;  //Aggregates are not initialized by default
    //Parameters for the coarse level solver. They can be changed in the main function
    int gmres_restarts_coarse_level = 10; 
    int gmres_restart_length_coarse_level = 250; //GMRES restart length for the coarse level
    double gmres_tol_coarse_level = 1e-10; //GMRES tolerance for the coarse level

    int gmres_restarts_smoother = 20; //Iterations for GMRES as a smoother (SAP is the default)

    int bi_cgstab_Dc_iterations= 1000; //Number of iterations for the bi-cgstab method
    double bi_cgstab_Dc_iterations_tol = 1e-10; //Tolerance for the bi-cgstab method
    int nu1 = 0; //Pre-smoothing iterations
    int nu2 = 2; //Post-smoothing iterations
    int Nit = 3; //Number of iterations for improving the interpolator

    bool SetUpDone = false; //Set to true when the setup is done
}

//--------------Parameters for FGMRES--------------//
namespace FGMRESV {
    double fgmres_tolerance = 1e-10; //Tolerance for FGMRES
    int fgmres_restart_length = 20; //Restart length for FGMRES
    int fgmres_restarts = 50; //Number of restarts for FGMRES
}


void CheckBlocks(){
    bool check = true;
    using namespace SAPV;
    if (Nx % block_x != 0) {
        std::cout << "Error: Ns/block_x is not an integer" << std::endl;
        check = false;
    }
    if (Nt % block_t != 0) {
        std::cout << "Error: Nt/block_t is not an integer" << std::endl;
        check = false;
    }
    if (Nx % sap_block_x != 0) {
        std::cout << "Error: Ns/sap_block_x is not an integer" << std::endl;
        check = false;
    }
    if (Nt % sap_block_t != 0) {
        std::cout << "Error: Nt/sap_block_t is not an integer" << std::endl;
        check = false;
    }
    if (sap_block_t*sap_block_x % 2 != 0 ){
        std::cout << "Error: sap_block_t*sap_block_x is not even" << std::endl;
        std::cout << "Expected an even number of SAP blocks" << std::endl;
        check = false;
    }
    if (check == false){
        exit(1);
    }

}

void CheckAggregates(){
    bool check = true;
    using namespace AMGV;
    if (aggregates_initialized == false) {
        std::cout << "Error: Aggregates are not initialized" << std::endl;
        check = false;
    }
    if (Nagg*Ntest > 2*LV::Ntot) {
        std::cout << "Error: Nagg*Ntest > 2*Ntot" << std::endl;
        check = false;
    }
    if (check == false) {
        exit(1);
    }
}
