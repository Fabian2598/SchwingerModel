#include "variables.h"

double pi=3.14159265359;

/*
	Vectorized lattice coords.*/
int Coords(const int& x, const int& t){
	return x*LV::Nt + t;
}

namespace CG{
	int max_iter = 1000; //Maximum number of iterations for the conjugate gradient method
	double tol = 1e-10; //Tolerance for convergence
}

int LeftPB[LV::Ntot*2]; 
int RightPB[LV::Ntot*2]; 
c_double SignL[LV::Ntot*2]; 
c_double SignR[LV::Ntot*2]; 
int x_1_t1[LV::Ntot]; 
int x1_t_1[LV::Ntot]; 


//Memory preallocation
std::vector<std::vector<c_double>> TEMP = std::vector<std::vector<c_double>>(LV::Ntot, std::vector<c_double>(2, 0));
std::vector<std::vector<c_double>> DTEMP = std::vector<std::vector<c_double>>(LV::Ntot, std::vector<c_double>(2, 0)); 