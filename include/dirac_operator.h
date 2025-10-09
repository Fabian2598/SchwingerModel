#ifndef DIRAC_OPERATOR_INCLUDED
#define DIRAC_OPERATOR_INCLUDED
#include <complex>
#include "gauge_conf.h"
#include "matrix_operations.h"
#include "mpi.h"


extern c_double I_number; //imaginary number


/*
	Periodic boundary conditions used for the link variables U_mu(n).
	This function builds th#include "omp.h"e arrays x_1_t1, x1_t_1, RightPB and LeftPB, which
	store the neighbor coordinates for the periodic boundary conditions.
	This prevents recalculation every time we call the operator D.
	The function is only called once at the beginning of the program.

	right periodic boundary x+hat{mu}
	left periodic boundary x-hat{mu}
	hat_mu[0] = { 1, 0 } --> hat_t
	hat_mu[1] = { 0, 1 } --> hat_x
*/
inline void periodic_boundary() {
	using namespace LV; //Lattice parameters namespace
	std::vector<std::vector<int>>hat_mu(2, std::vector<int>(2, 0));
	hat_mu[0] = { 1, 0 }; //hat_t
	hat_mu[1] = { 0, 1 }; //hat_x
	int x, t, n;
	for(x = 0; x<mpi::width; x++){
	for(t = 0; t<mpi::width; t++){
		n = x * mpi::width + t;
		x_1_t1[n] = Coords(mod(x - 1, mpi::width),mod(t + 1, mpi::width));
		x1_t_1[n] = Coords(mod(x + 1, mpi::width),mod(t - 1, mpi::width)); 
		for (int mu = 0; mu < 2; mu++) {
			RightPB[2*n+mu] = Coords(mod(x + hat_mu[mu][1], mpi::width),mod(t + hat_mu[mu][0], mpi::width)); 
			LeftPB[2*n+mu] = Coords(mod(x - hat_mu[mu][1], mpi::width),mod(t - hat_mu[mu][0], mpi::width));

			SignR[2*n+mu] = 1; //sign for the right boundary in time
			SignL[2*n+mu] = 1;
			if ((mpi::rank+1) % (mpi::size/2) == 0){
				SignR[2*n+mu] = (mu == 0 && t == mpi::width - 1) ? -1 : 1; //sign for the right boundary in time
			} 
			if (mpi::rank % (mpi::size/2) == 0){
				SignL[2*n+mu] = (mu == 0 && t == 0) ? -1 : 1; //sign for the left boundary in time	
			}
		}
	}
	}	
}


/*
	Dirac operator application D phi
	U: gauge configuration
	phi: spinor to apply the operator to
	m0: mass parameter
*/
void D_phi(const spinor& U, const spinor& phi, spinor &Dphi, const double& m0);


/*
	Dirac dagger operator application D^+ phi
	U: gauge configuration
	phi: spinor to apply the operator to
	m0: mass parameter
*/
void D_dagger_phi(const spinor& U, const spinor& phi, spinor &Dphi, const double& m0);


/*
	Application of D D^+
	It just calls the previous functions
*/
void D_D_dagger_phi(const spinor& U, const spinor& phi, spinor &Dphi,const double& m0);

/*
	2* Re ( left^+ d D / d omega(z) right )
	This derivative is needed for the fermion force
*/
re_field phi_dag_partialD_phi(const spinor& U, const spinor& left, const spinor& right);

#endif