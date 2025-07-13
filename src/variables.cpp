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

namespace CG{
	int max_iter = 1000; //Maximum number of iterations for the conjugate gradient method
	double tol = 1e-10; //Tolerance for convergence
}

std::vector<std::vector<int>>Coords = std::vector<std::vector<int>>(LV::Nx, std::vector<int>(LV::Nt, 0));
std::vector<std::vector<int>>LeftPB = std::vector<std::vector<int>>(LV::Ntot, std::vector<int>(2,0)); 
std::vector<std::vector<int>>RightPB = std::vector<std::vector<int>>(LV::Ntot, std::vector<int>(2,0)); 
std::vector<std::vector<c_double>>SignL =std::vector<std::vector<c_double>>(LV::Ntot,std::vector<c_double>(2,0)); 
std::vector<std::vector<c_double>>SignR = std::vector<std::vector<c_double>>(LV::Ntot,std::vector<c_double>(2,0)); 
std::vector<int>x_1_t1 = std::vector<int>(LV::Ntot, 0);
std::vector<int>x1_t_1 =  std::vector<int>(LV::Ntot, 0);
