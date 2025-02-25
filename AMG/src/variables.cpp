#include "variables.h"
//We define the variables that we declared in variables.h
double pi=3.14159265359;

std::vector<std::vector<int>>Coords = std::vector<std::vector<int>>(Ns, std::vector<int>(Nt, 0));
std::vector<std::vector<int>>x_1_t1 = std::vector<std::vector<int>>(Ns, std::vector<int>(Nt, 0));
std::vector<std::vector<int>>x1_t_1 = std::vector<std::vector<int>>(Ns, std::vector<int>(Nt, 0));
std::vector<std::vector<std::vector<int>>>LeftPB = std::vector<std::vector<std::vector<int>>>(Ns, std::vector<std::vector<int>>(Nt,std::vector<int>(2, 0)));
std::vector<std::vector<std::vector<int>>>RightPB = std::vector<std::vector<std::vector<int>>>(Ns, std::vector<std::vector<int>>(Nt, std::vector<int>(2, 0)));
std::vector<std::vector<int>>Agg = std::vector<std::vector<int>>(block_x*block_t, std::vector<int>(x_elements*t_elements, 0));
std::vector<std::unordered_set<int>> Agg_sets(Nagg); //aggregate sets