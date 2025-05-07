#ifndef VARIABLES_H_INCLUDED
#define VARIABLES_H_INCLUDED
#include "config.h"
#include <iostream>
#include <complex>
#include <vector>

extern double pi;
extern double it_count;

constexpr int Ns=NS; //We extract this value from config.h
constexpr int Nt = NT; //We extract this value from config.h
constexpr int block_x = BLOCK_X; //We extract this value from config.h
constexpr int block_t = BLOCK_T; //We extract this value from config.h

constexpr int sap_block_x = SAP_BLOCK_X; //We extract this value from config.h
constexpr int sap_block_t = SAP_BLOCK_T; //We extract this value from config.h
constexpr int sap_x_elements = Ns/sap_block_x; //Number of elements in the x direction
constexpr int sap_t_elements = Nt/sap_block_t; //Number of elements in the t direction
constexpr int N_sap_blocks = sap_block_x * sap_block_t; 

constexpr int Ntest = NTEST; //Number of test vectors
constexpr int x_elements = Ns/block_x; //Number of elements in the x direction
constexpr int t_elements = Nt/block_t; //Number of elements in the t direction
constexpr int Nagg = NAGG;//2*block_x * block_t; //Number of aggregates (one spin per aggregate)
constexpr int Ntot = Ns*Nt; //Total number of lattice points

extern std::vector<std::vector<int>>Coords; //Vectorize coordinates for the lattice points
 //Vectorized Coordinates for the aggregates
extern std::vector<int> XCoord; 
extern std::vector<int> TCoord;
extern std::vector<int> SCoord;


extern std::vector<std::vector<int>>x_1_t1;
extern std::vector<std::vector<int>>x1_t_1;
extern std::vector<std::vector<std::vector<int>>>RightPB; //Right periodic boundary
extern std::vector<std::vector<std::vector<int>>>LeftPB; //Left periodic boundary
extern std::vector<std::vector<int>> Agg; //Aggregates[number_of_aggregate][vectorized_coordinate of the lattice point]
//The vectorization also takes into account the spin index

extern std::vector<std::vector<int>> SAP_Blocks; //SAP_Blocks[number_of_block][vectorized_coordinate of the lattice point]
//The vectorization does not take into account the spin index, since both spin indices are in the same block.

void CheckBlocks(); //Check that Ns/block_x and Nt/block_t are integers, the same for Schwarz blocks

#endif 