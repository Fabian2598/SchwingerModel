#include "variables.h"

double pi=3.14159265359;
/*
	Vectorized coordinates
*/
void Coordinates() {
	for (int x = 0; x < LV::Nx; x++) {
		for (int t = 0; t < LV::Nt; t++) {
			Coords[x][t] = x * LV::Nx + t;
		}
	}
}

std::vector<std::vector<int>>Coords = std::vector<std::vector<int>>(LV::Nx, std::vector<int>(LV::Nt, 0));
std::vector<std::vector<int>>x_1_t1 = std::vector<std::vector<int>>(LV::Nx, std::vector<int>(LV::Nt, 0));
std::vector<std::vector<int>>x1_t_1 = std::vector<std::vector<int>>(LV::Nx, std::vector<int>(LV::Nt, 0));
std::vector<std::vector<std::vector<int>>>LeftPB = std::vector<std::vector<std::vector<int>>>(LV::Nx, std::vector<std::vector<int>>(LV::Nt,std::vector<int>(2, 0)));
std::vector<std::vector<std::vector<int>>>RightPB = std::vector<std::vector<std::vector<int>>>(LV::Nx, std::vector<std::vector<int>>(LV::Nt, std::vector<int>(2, 0)));



std::vector<std::vector<int>>LeftPBT = std::vector<std::vector<int>>(LV::Ntot, 
   std::vector<int>(2,0)); //LeftPB[x][t][mu]
std::vector<std::vector<int>>RightPBT = std::vector<std::vector<int>>(LV::Ntot, 
   std::vector<int>(2,0)); //RightPB[x][t][mu]


std::vector<std::vector<c_double>>SignLT =std::vector<std::vector<c_double>>(LV::Ntot, 
    std::vector<c_double>(2,0)); //SignL[x][t][mu]
std::vector<std::vector<c_double>>SignRT = std::vector<std::vector<c_double>>(LV::Ntot, 
    std::vector<c_double>(2,0)); ////SignR[x][t][mu]
