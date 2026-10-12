#include "dirac_operator.h"


/*
 *                t                    2D parallelization
 *   0  +-------------------+  Nt   +---------------------+
 *      |                   |       |  rank 0  |  rank 1  |
 *      |                   |       |---------------------|
 *      |                   |       |  rank 2  |  rank 3  |
 *   x  |                   |       |---------------------|
 *      |                   |       |  rank 4  |  rank 5  |
 *      |                   |       |---------------------|
 *      |                   |       |  rank 6  |  rank 7  |
 *   Nx +-------------------+ Nt    +---------------------+
 *                Nx
 * x+1, t down
 * x-1, t up
 * x, t+1 right
 * x, t-1 left
 * n = x * Nt + t = (x,t) coordinates
 *
 * The loops over the lattice are executed by for_each_site (GPU kernels when compiled with CUDA). The kernels
 * only capture raw pointers and plain values (check backend.h), which is why the spinors are called phi_in,
 * Dphi_out, ... and the pointers to their data phi, Dphi, ...
 */
//Eqs (34) of the documentation

 void D_phi(const GaugeConf& GConf, const spinor&  phi_in, spinor&  Dphi_out){
	//Communicate halos
	exchange_halo(phi_in.val);
	//exchange_halo(U.val);
	const LocalLattice lat = local_lattice();
	const c_double* U = GConf.Conf.val;
	const c_double* phi = phi_in.val;
	c_double* Dphi = Dphi_out.val;
	const double m0 = sim_params::m0;
	#ifdef CLOVER
	const c_double* Q01 = GConf.Q01.val;
	const double csw = sim_params::csw;
	#endif
	for_each_site([=] HOST_DEVICE (int n){
		int right, down, left, up;
		double rsign, lsign;
		//get coordinates of the neighbors and boundary sign
		get_neighbors(lat, n, right, down, left, up, rsign, lsign); //check boundary.h for conventions

		#ifdef CLOVER
		const c_double DeltaQ = Q01[n]-conj(Q01[n]);
		const c_double diag0 = m0 + 2 - I_number*l_16*csw*DeltaQ; //mu = 0
		const c_double diag1 = m0 + 2 + I_number*l_16*csw*DeltaQ; //mu = 1
		#else
		const double diag0 = m0 + 2;
		const double diag1 = m0 + 2;
		#endif

		//mu = 0
		Dphi[2*n] = diag0  * phi[2*n] - 0.5 * (
				U[2*n] 	 			* rsign  	* (phi[2*right] - phi[2*right+1])
			+	U[2*n+1] 			*  			  (phi[2*down] + I_number * phi[2*down+1])
			+ conj(U[2*left])  		* lsign		* (phi[2*left] + phi[2*left+1])
			+ conj(U[2*up+1]) 		*  			  (phi[2*up] - I_number*phi[2*up+1])
		);
		//mu = 1
		Dphi[2*n+1] = diag1 * phi[2*n+1] - 0.5 * (
				U[2*n] 	 			* rsign 	* (-phi[2*right] + phi[2*right+1])
			+	U[2*n+1] 			* 			  (-I_number*phi[2*down] + phi[2*down+1])
			+ conj(U[2*left])  		* lsign 	* (phi[2*left] + phi[2*left+1])
			+ conj(U[2*up+1]) 		* 			  (I_number*phi[2*up] + phi[2*up+1])
		);
	});
}

void D_dagger_phi(const GaugeConf& GConf, const spinor&  phi_in, spinor&  Dphi_out){
	//Communicate halos
	exchange_halo(phi_in.val);
	//exchange_halo(U.val);
	const LocalLattice lat = local_lattice();
	const c_double* U = GConf.Conf.val;
	const c_double* phi = phi_in.val;
	c_double* Dphi = Dphi_out.val;
	const double m0 = sim_params::m0;
	#ifdef CLOVER
	const c_double* Q01 = GConf.Q01.val;
	const double csw = sim_params::csw;
	#endif
	for_each_site([=] HOST_DEVICE (int n){
		int right, down, left, up;
		double rsign, lsign;
		//get coordinates of the neighbors and boundary sign
		get_neighbors(lat, n, right, down, left, up, rsign, lsign); //check boundary.h for conventions

		#ifdef CLOVER
		//1/16=0.0625
		const c_double DeltaQ = Q01[n]-conj(Q01[n]);
		const c_double diag0 = m0 + 2.0 - I_number*l_16*csw*DeltaQ; //mu = 0
		const c_double diag1 = m0 + 2.0 + I_number*l_16*csw*DeltaQ; //mu = 1
		#else
		const double diag0 = m0 + 2.0;
		const double diag1 = m0 + 2.0;
		#endif

		//mu = 0
		Dphi[2*n] = diag0 * phi[2*n] -0.5 * (
			conj(U[2*left]) 		* lsign 	* (phi[2*left] - phi[2*left+1])
		+   conj(U[2*up+1]) 	 	* (phi[2*up] + I_number * phi[2*up+1])
		+   U[2*n] 					* rsign 	* (phi[2*right] + phi[2*right+1])
		+	U[2*n+1] 				* (phi[2*down] - I_number * phi[2*down+1])
		);
		//mu = 1
		Dphi[2*n+1] = diag1 * phi[2*n+1] -0.5 * (
			conj(U[2*left]) 		* lsign 	* (-phi[2*left] + phi[2*left+1])
		+   conj(U[2*up+1]) 	 	* (-I_number*phi[2*up] + phi[2*up+1])
		+   U[2*n] 					* rsign 	* (phi[2*right] + phi[2*right+1])
		+	U[2*n+1] 	 			* (I_number * phi[2*down] + phi[2*down+1])
		);
	});
}

//Buffer for D^+ phi. It is allocated at the first call, once initializeMPI() has selected the GPU.
static spinor& ddagg_buffer(){
	static spinor buffer(mpi::maxSizeH);
	return buffer;
}

void D_D_dagger_phi(const GaugeConf& GConf, const spinor& phi_in, spinor &Dphi_out){
	D_dagger_phi(GConf, phi_in, ddagg_buffer());
	D_phi(GConf,  ddagg_buffer(), Dphi_out);

	#ifdef TWISTED_MASS
	//We add the twisted mass term
	//Due to the way things are implemented in the code, we don't need to modify anything else.
	const double tm2 = sim_params::tm*sim_params::tm;
	const c_double* phi = phi_in.val;
	c_double* Dphi = Dphi_out.val;
	for_each_site([=] HOST_DEVICE (int n){
		Dphi[2*n]   += tm2*phi[2*n];
		Dphi[2*n+1] += tm2*phi[2*n+1];
	});
	#endif

}


//Fermion force of the hopping terms
//2* Re ( left^dag \partial D / \partial omega(z) right )
//Eqs (37) and (38) of the documentation
void phi_dag_partialD_phi(const GaugeConf& GConf, const spinor& left_in,const spinor& right_in,re_field& Dphi_out){
	//Communicate halos
	exchange_halo(left_in.val);
	exchange_halo(right_in.val);
	//exchange_halo(U.val);
	const LocalLattice lat = local_lattice();
	const c_double* U = GConf.Conf.val;
	const c_double* left_term = left_in.val;
	const c_double* right_term = right_in.val;
	double* Dphi = Dphi_out.val;
	for_each_site([=] HOST_DEVICE (int n){
		int right, down, left, up;
		double rsign, lsign;
		//get coordinates of the neighbors and boundary sign
		get_neighbors(lat, n, right, down, left, up, rsign, lsign); //check boundary.h for conventions
		//n = x * Nt + t

		//mu = 0
		Dphi[2*n] = imag(
		U[2*n] * rsign * ( conj(left_term[2*n] - left_term[2*n+1]) ) * (right_term[2*right] - right_term[2*right+1])
		- conj(U[2*n]) * rsign * ( conj(left_term[2*right] + left_term[2*right+1]) ) * (right_term[2*n] + right_term[2*n+1])
		);
		//mu = 1
		Dphi[2*n+1] = imag(
		U[2*n+1] * ( conj(left_term[2*n]) - I_number*conj(left_term[2*n+1]) ) * (right_term[2*down] + I_number * right_term[2*down+1])
		+ conj(U[2*n+1]) * ( conj(left_term[2*down]) + I_number*conj(left_term[2*down+1]) )
		* (-right_term[2*n] + I_number * right_term[2*n+1])
		);
	});
}
