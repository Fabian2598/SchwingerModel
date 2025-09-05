#include "variables.h"

typedef std::complex<double> c_double;
double coarse_time = 0.0; //Time spent in the coarse grid solver
double smooth_time = 0.0; //Time spent in the smoother
double SAP_time = 0.0; //Time spent in the SAP method
double m0 = 0.0;


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


namespace SAPV {
    bool schwarz_blocks = false; //Schwarz blocks are not initialized by default
    int sap_gmres_restart_length = 5; //GMRES restart length for the Schwarz blocks. Set to 20 by default
    int sap_gmres_restarts = 5; //GMRES iterations for the Schwarz blocks. Set to 10 by default.
    double sap_gmres_tolerance = 1e-3; //GMRES tolerance for the Schwarz blocks
    double sap_tolerance = 1e-10; //Tolerance for the SAP method
    int sap_blocks_per_proc = 1; //Number of blocks per process for the parallel SAP method
}

namespace AMGV {
    int SAP_test_vectors_iterations = 1; //Number of SAP iterations to smooth test vectors
    bool aggregates_initialized = false;  //Aggregates are not initialized by default
    //Parameters for the coarse level solver. They can be changed in the main function
    int gmres_restarts_coarse_level = 10; 
    int gmres_restart_length_coarse_level = 20; //GMRES restart length for the coarse level
    double gmres_tol_coarse_level = 0.1; //GMRES tolerance for the coarse level

    int nu1 = 0; //Pre-smoothing iterations
    int nu2 = 2; //Post-smoothing iterations
    int Nit = 1; //Number of iterations for improving the interpolator

    bool SetUpDone = false; //Set to true when the setup is done
}

//--------------Parameters for FGMRES--------------//
namespace FGMRESV {
    double fgmres_tolerance = 1e-10; //Tolerance for FGMRES
    int fgmres_restart_length = 20; //Restart length for FGMRES
    int fgmres_restarts = 50; //Number of restarts for FGMRES
}

namespace CG{
    int max_iter = 100000;
    double tol = 1e-10;
}

std::vector<std::vector<c_double>>D_TEMP = std::vector<std::vector<c_double>>(LV::Ntot, std::vector<c_double>(2,0)); 


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

std::vector<std::vector<int>> LatticeBlocks = std::vector<std::vector<int>> (LV::Nblocks, std::vector<int>(LV::lattice_sites_per_block));
int RightPB_blocks[LV::Nblocks][2];
int LeftPB_blocks[LV::Nblocks][2];

void MakeBlocks(){
	using namespace LV; //Lattice parameters namespace
	int count; 
	for (int x = 0; x < block_x; x++) {
		for (int t = 0; t < block_t; t++) {
				int block = x * block_t + t; //block index
				int x0 = x * x_elements, t0 = t * t_elements;
				int x1 = (x + 1) * x_elements, t1 = (t + 1) * t_elements;
            	count = 0;  
            	for(int x = x0; x < x1; x++) {
                	for (int t = t0; t < t1; t++) {
                    	LatticeBlocks[block][count++] = x * Nt+ t; 
                    	
					}
            	}
		}
	}
}

void save_vec(const std::vector<double>& vec,const std::string& Name){
    std::ofstream Datfile(Name);
    if (!Datfile.is_open()) {
        std::cerr << "Error opening file: " << Name << std::endl;
        return;
    }
    int size = vec.size();
    for (int n = 0; n < size; n++) {
        Datfile << n
                << std::setw(30) << std::setprecision(17) << std::scientific << vec[n]
                << "\n";
    }
        
    
}

void read_rhs(std::vector<std::vector<c_double>>& vec,const std::string& name){
    std::ifstream infile(name);
    if (!infile) {
        std::cerr << "File " << name << " not found" << std::endl;
    }
    int x, t, mu;
    double re, im;
    //x, t, mu, real part, imaginary part
    while (infile >> x >> t >> mu >> re >> im) {
        vec[Coords[x][t]][mu] = c_double(re, im); 
    }
    infile.close();
  
}

void save_rhs(std::vector<std::vector<c_double>>& rhs,const std::string& name){
    std::ofstream rhsfile(name);
    if (!rhsfile.is_open()) {
        std::cerr << "Error opening rhs.txt for writing." << std::endl;
    } 
    else {
        int x,t;
        //x, t, mu, real part, imaginary part
        for (int n = 0; n < LV::Ntot; ++n) {
            x = n/LV::Nt;
            t = n%LV::Nt;
            rhsfile << x << std::setw(30) << t << std::setw(30) << 0 << std::setw(30)
                    << std::setprecision(17) << std::scientific << std::real(rhs[n][0]) << std::setw(30)
                    << std::setprecision(17) << std::scientific << std::imag(rhs[n][0]) << "\n";

            rhsfile << x << std::setw(30) << t << std::setw(30) << 1 << std::setw(30)
                    << std::setprecision(17) << std::scientific << std::real(rhs[n][1]) << std::setw(30)
                    << std::setprecision(17) << std::scientific << std::imag(rhs[n][1]) << "\n";
          
        }
        rhsfile.close();
    }

}


void random_rhs(std::vector<std::vector<c_double>>& vec,const int seed){
    c_double I_number(0,1);
    static std::mt19937 randomInt(seed);
	std::uniform_real_distribution<double> distribution(-1.0, 1.0); //mu, standard deviation
    for(int i = 0; i < LV::Ntot; i++) {
        vec[i][0] = distribution(randomInt) + I_number * distribution(randomInt); //RandomU1();
        vec[i][1] = distribution(randomInt) + I_number * distribution(randomInt);
    }

}


void print_parameters(){
    using namespace SAPV;
    using namespace AMGV;
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
     
    
    std::cout << "----------------------------------" << std::endl;
    std::cout << " Variable blocking for SAP" << std::endl;
    std::cout << "| sap_block_x = " << sap_block_x << " sap_block_t = " << sap_block_t << std::endl;
    std::cout << "| Lattice sites in the x direction = " << sap_x_elements<< " and in the t direction = " << sap_t_elements<< std::endl;
    std::cout << "| Each Schwarz block has " <<  sap_lattice_sites_per_block << " lattice points and " << sap_variables_per_block << " variables" << std::endl;
    std::cout << "| D restricted to each block has (" << 2 * sap_lattice_sites_per_block << ")^2 = " << sap_variables_per_block*sap_variables_per_block << " entries" << std::endl;
    std::cout << "| Number of Schwarz blocks = " << N_sap_blocks << std::endl;
    std::cout << "| Red/Black blocks = " << sap_coloring_blocks << std::endl;
    std::cout << "| Number of blocks per process = " << sap_blocks_per_proc << std::endl;
    std::cout << "| GMRES restart length for SAP blocks = " << sap_gmres_restart_length << std::endl;
    std::cout << "| GMRES iterations for SAP blocks = " << sap_gmres_restarts << std::endl;
    std::cout << "| GMRES tolerance for SAP blocks = " << sap_gmres_tolerance << std::endl;
    
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
    

}