#ifndef DIRAC_OPERATOR_INCLUDED
#define DIRAC_OPERATOR_INCLUDED
#include <complex>
#include "mpi.h"
#include "halo_exchange.h"
#include "boundary.h"


/*
	Dirac operator application D phi
	U: gauge configuration
	phi: spinor to apply the operator to
	m0: mass parameter
*/
void D_phi(const spinor& U, const spinor&  phi, spinor&  Dphi);

/*
	Dirac dagger operator application D^+ phi
	U: gauge configuration
	phi: spinor to apply the operator to
	m0: mass parameter
*/
void D_dagger_phi(const spinor& U, const spinor&  phi, spinor&  Dphi);

/*
	Application of D D^+
	It just calls the previous functions
*/
void D_D_dagger_phi(const spinor& U, const spinor& phi, spinor &Dphi);

/*
	2* Re ( left^+ d D / d omega(z) right )
	This derivative is needed for the fermion force
*/
void phi_dag_partialD_phi(const spinor& U, const spinor& left_term,const spinor& right_term, re_field& Dphi);


#endif