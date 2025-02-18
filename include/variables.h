#ifndef VARIABLES_H_INCLUDED
#define VARIABLES_H_INCLUDED
#include "config.h"
#include <vector>

extern double pi;

constexpr int Ns=NS; //We extract this value from config.h
constexpr int Nt = NT; //We extract this value from config.h

extern std::vector<std::vector<int>>Coords;
extern std::vector<std::vector<int>>x_1_t1;
extern std::vector<std::vector<int>>x1_t_1;
extern std::vector<std::vector<std::vector<int>>>RightPB; //Right periodic boundary
extern std::vector<std::vector<std::vector<int>>>LeftPB; //Left periodic boundary



#endif 