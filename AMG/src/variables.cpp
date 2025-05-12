#include "variables.h"

//We define the variables that we declared in variables.h
double pi=3.14159265359;
double it_count = 0;

std::vector<std::vector<int>>Coords = std::vector<std::vector<int>>(Ns, std::vector<int>(Nt, 0));

std::vector<int> XCoord = std::vector<int>(2*Ntot, 0);
std::vector<int> TCoord = std::vector<int>(2*Ntot, 0);
std::vector<int> SCoord = std::vector<int>(2*Ntot, 0);

std::vector<std::vector<int>>x_1_t1 = std::vector<std::vector<int>>(Ns, std::vector<int>(Nt, 0));
std::vector<std::vector<int>>x1_t_1 = std::vector<std::vector<int>>(Ns, std::vector<int>(Nt, 0));
std::vector<std::vector<std::vector<int>>>LeftPB = std::vector<std::vector<std::vector<int>>>(Ns, std::vector<std::vector<int>>(Nt,std::vector<int>(2, 0)));
std::vector<std::vector<std::vector<int>>>RightPB = std::vector<std::vector<std::vector<int>>>(Ns, std::vector<std::vector<int>>(Nt, std::vector<int>(2, 0)));
std::vector<std::vector<int>>Agg = std::vector<std::vector<int>>(2*block_x*block_t, std::vector<int>(x_elements*t_elements, 0));

std::vector<std::vector<int>>SAP_Blocks = std::vector<std::vector<int>>(sap_block_x*sap_block_t, std::vector<int>(sap_x_elements*sap_t_elements, 0));
std::vector<int> SAP_RedBlocks = std::vector<int>(sap_coloring_blocks, 0); //Red blocks
std::vector<int> SAP_BlackBlocks = std::vector<int>(sap_coloring_blocks, 0); //Red blocks

bool schwarz_blocks = false; //Schwarz blocks are not initialized by default
int sap_gmres_restart_length = 20; //GMRES restart length for the Schwarz blocks. Set to 20 by default
int sap_gmres_restarts = 10; //GMRES iterations for the Schwarz blocks. Set to 10 by default.
double sap_gmres_tolerance = 1e-10; //GMRES tolerance for the Schwarz blocks
double sap_tolerance = 1e-10; //Tolerance for the SAP method


void CheckBlocks(){
    bool check = true;
    if (Ns % block_x != 0) {
        std::cout << "Error: Ns/block_x is not an integer" << std::endl;
        check = false;
    }
    if (Nt % block_t != 0) {
        std::cout << "Error: Nt/block_t is not an integer" << std::endl;
        check = false;
    }
    if (Ns % sap_block_x != 0) {
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
