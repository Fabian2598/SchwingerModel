#ifndef VARIABLES_H_INCLUDED
#define VARIABLES_H_INCLUDED
#include "config.h"
#include <complex>
#include <vector>

extern double pi;
extern double it_count;

constexpr int Ns=NS; //We extract this value from config.h
constexpr int Nt = NT; //We extract this value from config.h
constexpr int block_x = BLOCK_X; //We extract this value from config.h
constexpr int block_t = BLOCK_T; //We extract this value from config.h
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


#endif 