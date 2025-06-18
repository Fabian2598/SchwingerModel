#ifndef VARIABLES_H_INCLUDED
#define VARIABLES_H_INCLUDED
#include "config.h"
#include <iostream>
#include <complex>
#include <vector>

extern double pi;
extern double it_count;
//Lattice variables
namespace LV {
    constexpr int Nx=NS; //We extract this value from config.h
    constexpr int Nt = NT; //We extract this value from config.h
    constexpr int block_x = BLOCK_X; //We extract this value from config.h
    constexpr int block_t = BLOCK_T; //We extract this value from config.h
    constexpr int x_elements = Nx/block_x; //Number of elements in the x direction
    constexpr int t_elements = Nt/block_t; //Number of elements in the t direction
    constexpr int Ntot = Nx*Nt; //Total number of lattice points
}
//------------Schwarz alternating procedure parameters--------------//
namespace SAPV {
    using namespace LV; 
    constexpr int sap_block_x = SAP_BLOCK_X; //We extract this value from config.h
    constexpr int sap_block_t = SAP_BLOCK_T; //We extract this value from config.h
    constexpr int sap_x_elements = Nx/sap_block_x; //Number of lattice points in the x direction (without the spin index)
    constexpr int sap_t_elements = Nt/sap_block_t; //Number of lattice points in the t direction (without the spin index)
    constexpr int N_sap_blocks = sap_block_x * sap_block_t; //Number of Schwarz blocks
    constexpr int sap_lattice_sites_per_block = sap_x_elements * sap_t_elements; //Number of lattice points in the block
    constexpr int sap_variables_per_block = 2 * sap_lattice_sites_per_block; //Number of variables in the block
    constexpr int sap_coloring_blocks = N_sap_blocks/2; //Number of red or black blocks 
    extern bool schwarz_blocks; //True if the Schwarz blocks are initialized
    extern int sap_gmres_restart_length; //GMRES restart length for the Schwarz blocks.
    extern int sap_gmres_restarts; //GMRES iterations for the Schwarz blocks
    extern double sap_gmres_tolerance; //GMRES tolerance for the Schwarz blocks 
    extern double sap_tolerance; //tolerance for the SAP method
}

constexpr int Ntest = NTEST; //Number of test vectors
constexpr int Nagg = 2*LV::block_t*LV::block_x; //Number of aggregates (one spin per aggregate)
extern bool aggregates_initialized; //True if the aggregates are initialized

//Vectorize coordinates for the lattice points
extern std::vector<std::vector<int>>Coords; 
void Coordinates();

 //Vectorized Coordinates for the aggregates. 
 //Intialized in Aggregates()
extern std::vector<int> XCoord; 
extern std::vector<int> TCoord;
extern std::vector<int> SCoord;
extern std::vector<std::vector<int>> Agg; //Aggregates[number_of_aggregate][vectorized_coordinate of the lattice point]
//The vectorization also takes into account the spin index


//--Coordinates of the neighbors to avoid recomputing them each time the operator D is called--//
//Check matrix_operations.h for the definition of RightPB and LeftPB
extern std::vector<std::vector<int>>x_1_t1;
extern std::vector<std::vector<int>>x1_t_1;
extern std::vector<std::vector<std::vector<int>>>RightPB; //Right periodic boundary
extern std::vector<std::vector<std::vector<int>>>LeftPB; //Left periodic boundary
//---------------------------------------------------------------------------------------------//


extern std::vector<std::vector<int>> SAP_Blocks; //SAP_Blocks[number_of_block][vectorized_coordinate of the lattice point]
//The vectorization does not take into account the spin index, since both spin indices are in the same block.
extern std::vector<int> SAP_RedBlocks; //Block index for the red blocks
extern std::vector<int> SAP_BlackBlocks; //Block index for the black blocks

void CheckBlocks(); //Check that Nx/block_x and Nt/block_t are integers, the same for Schwarz blocks

#endif 