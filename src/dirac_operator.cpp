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
 */
//Eqs (34) of the documentation

 void D_phi(const GaugeConf& GConf, const spinor&  phi, spinor&  Dphi){
	using namespace mpi;
	using namespace sim_params;

	int n, right, down, left, up;
	double rsign, lsign;
	c_double DeltaQ;
	const spinor& U = GConf.Conf;
	//Communicate halos 
	exchange_halo(phi.val);
	//exchange_halo(U.val);
	for(int x = 1; x<=width_x; x++){
		for(int t = 1; t<=width_t; t++){
			n = x*(width_t+2)+t;
			//get coordinates of the neighbors and boundary sign 
			get_neighbors(x, t,right, down, left, up, rsign, lsign); //check boundary.h for conventions
			
			#ifdef CLOVER
			//mu = 0
			DeltaQ = GConf.Q01[n]-std::conj(GConf.Q01[n]);
			Dphi.val[2*n] = (m0 + 2 - I_number*l_16*csw*DeltaQ)  * phi.val[2*n] - 0.5 * ( 
					U.val[2*n] 	 			* rsign  	* (phi.val[2*right] - phi.val[2*right+1])
				+	U.val[2*n+1] 			*  			  (phi.val[2*down] + I_number * phi.val[2*down+1])
				+ std::conj(U.val[2*left])  * lsign		* (phi.val[2*left] + phi.val[2*left+1])
				+ std::conj(U.val[2*up+1]) 	*  			  (phi.val[2*up] - I_number*phi.val[2*up+1])
			);
			//mu = 1
			Dphi.val[2*n+1] = (m0 + 2 + I_number*l_16*csw*DeltaQ) * phi.val[2*n+1] - 0.5 * ( 
					U.val[2*n] 	 			* rsign 	* (-phi.val[2*right] + phi.val[2*right+1])
				+	U.val[2*n+1] 			* 			  (-I_number*phi.val[2*down] + phi.val[2*down+1])
				+ std::conj(U.val[2*left])  * lsign 	* (phi.val[2*left] + phi.val[2*left+1])
				+ std::conj(U.val[2*up+1]) 	* 			  (I_number*phi.val[2*up] + phi.val[2*up+1])
			);
			#endif


			#ifndef CLOVER
			//mu = 0
			Dphi.val[2*n] = (m0 + 2) * phi.val[2*n] - 0.5 * ( 
					U.val[2*n] 	 			* rsign  	* (phi.val[2*right] - phi.val[2*right+1])
				+	U.val[2*n+1] 			*  			  (phi.val[2*down] + I_number * phi.val[2*down+1])
				+ std::conj(U.val[2*left])  * lsign		* (phi.val[2*left] + phi.val[2*left+1])
				+ std::conj(U.val[2*up+1]) 	*  			  (phi.val[2*up] - I_number*phi.val[2*up+1])
			);
			//mu = 1
			Dphi.val[2*n+1] = (m0 + 2) * phi.val[2*n+1] - 0.5 * ( 
					U.val[2*n] 	 			* rsign 	* (-phi.val[2*right] + phi.val[2*right+1])
				+	U.val[2*n+1] 			* 			  (-I_number*phi.val[2*down] + phi.val[2*down+1])
				+ std::conj(U.val[2*left])  * lsign 	* (phi.val[2*left] + phi.val[2*left+1])
				+ std::conj(U.val[2*up+1]) 	* 			  (I_number*phi.val[2*up] + phi.val[2*up+1])
			);
			#endif
		}
	}		
}

void D_dagger_phi(const GaugeConf& GConf, const spinor&  phi, spinor&  Dphi){
	using namespace mpi;
	using namespace sim_params;

	int n, right, down, left, up;
	double rsign, lsign;
	c_double DeltaQ;
	const spinor& U = GConf.Conf;
	//Communicate halos 
	exchange_halo(phi.val);
	//exchange_halo(U.val);
	for(int x = 1; x<=width_x; x++){
		for(int t = 1; t<=width_t; t++){
			n = x*(width_t+2)+t;
			//get coordinates of the neighbors and boundary sign 
			get_neighbors(x, t,right, down, left, up, rsign, lsign); //check boundary.h for conventions

			#ifdef CLOVER
			//mu = 0 
			//1/16=0.0625
			DeltaQ = GConf.Q01[n]-std::conj(GConf.Q01[n]);
			Dphi.val[2*n] = (m0 + 2 + I_number*l_16*csw*DeltaQ ) * phi.val[2*n] -0.5 * ( 
				std::conj(U.val[2*left]) 		* lsign 	* (phi.val[2*left] - phi.val[2*left+1])
			+   std::conj(U.val[2*up+1]) 	 	* (phi.val[2*up] + I_number * phi.val[2*up+1])
			+   U.val[2*n] 						* rsign 		* (phi.val[2*right] + phi.val[2*right+1])
			+	U.val[2*n+1] 					* (phi.val[2*down] - I_number * phi.val[2*down+1])
			);
			//mu = 1
			Dphi.val[2*n+1] = (m0 + 2 - I_number*l_16*csw*DeltaQ) * phi.val[2*n+1] -0.5 * ( 
				std::conj(U.val[2*left]) 		* lsign 	* (-phi.val[2*left] + phi.val[2*left+1])
			+   std::conj(U.val[2*up+1]) 	 	* (-I_number*phi.val[2*up] + phi.val[2*up+1])
			+   U.val[2*n] 						* rsign 		* (phi.val[2*right] + phi.val[2*right+1])
			+	U.val[2*n+1] 	 				* (I_number * phi.val[2*down] + phi.val[2*down+1])
			);
			#endif

			#ifndef CLOVER
			//mu = 0
			Dphi.val[2*n] = (m0 + 2) * phi.val[2*n] -0.5 * ( 
				std::conj(U.val[2*left]) 		* lsign 	* (phi.val[2*left] - phi.val[2*left+1])
			+   std::conj(U.val[2*up+1]) 	 	* (phi.val[2*up] + I_number * phi.val[2*up+1])
			+   U.val[2*n] 						* rsign 		* (phi.val[2*right] + phi.val[2*right+1])
			+	U.val[2*n+1] 					* (phi.val[2*down] - I_number * phi.val[2*down+1])
			);
			//mu = 1
			Dphi.val[2*n+1] = (m0 + 2) * phi.val[2*n+1] -0.5 * ( 
				std::conj(U.val[2*left]) 		* lsign 	* (-phi.val[2*left] + phi.val[2*left+1])
			+   std::conj(U.val[2*up+1]) 	 	* (-I_number*phi.val[2*up] + phi.val[2*up+1])
			+   U.val[2*n] 						* rsign 		* (phi.val[2*right] + phi.val[2*right+1])
			+	U.val[2*n+1] 	 				* (I_number * phi.val[2*down] + phi.val[2*down+1])
			);
			#endif
		}
		
	}	
}

spinor ddagg_buffer(mpi::maxSizeH);//Note that maxSizeH = 2*(Nx+2)*(Nt+2) by default, so this buffer is quite large ... 
void D_D_dagger_phi(const GaugeConf& GConf, const spinor& phi, spinor &Dphi){
	D_dagger_phi(GConf, phi, ddagg_buffer);
	D_phi(GConf,  ddagg_buffer, Dphi);
	
	#ifdef TWISTED_MASS
	//We add the twisted mass term
	//Due to the way things are implemented in the code, we don't need to modify anything else. 
	int n;
	for(int x = 1; x<=mpi::width_x; x++){
		for(int t = 1; t<=mpi::width_t; t++){
			n = x*(mpi::width_t+2)+t;
			Dphi.val[2*n]   += sim_params::tm*sim_params::tm*phi.val[2*n];
			Dphi.val[2*n+1] += sim_params::tm*sim_params::tm*phi.val[2*n+1];
		}
	}
	#endif

}


//Fermion force of the hopping terms
//2* Re ( left^dag \partial D / \partial omega(z) right )
//Eqs (37) and (38) of the documentation
void phi_dag_partialD_phi(const GaugeConf& GConf, const spinor& left_term,const spinor& right_term,re_field& Dphi){
	using namespace mpi;

	int n, right, down, left, up;
	double rsign, lsign;
	const spinor& U = GConf.Conf;
	//Communicate halos 
	exchange_halo(left_term.val);
	exchange_halo(right_term.val);
	//exchange_halo(U.val);
	for(int x = 1; x<=width_x; x++){
		for(int t = 1; t<=width_t; t++){
			n = x*(width_t+2)+t;
			//get coordinates of the neighbors and boundary sign 
			get_neighbors(x, t,right, down, left, up, rsign, lsign); //check boundary.h for conventions
			//n = x * Nt + t

			//mu = 0
			Dphi.val[2*n] = std::imag(
			U.val[2*n] * rsign * ( std::conj(left_term.val[2*n] - left_term.val[2*n+1]) ) * (right_term.val[2*right] - right_term.val[2*right+1])
			- std::conj(U.val[2*n]) * rsign * ( std::conj(left_term.val[2*right] + left_term.val[2*right+1]) ) * (right_term.val[2*n] + right_term.val[2*n+1])
			);
			//mu = 1
			Dphi.val[2*n+1] = std::imag(
			U.val[2*n+1] * ( std::conj(left_term.val[2*n]) - I_number*std::conj(left_term.val[2*n+1]) ) * (right_term.val[2*down] + I_number * right_term.val[2*down+1])
			+ std::conj(U.val[2*n+1]) * ( std::conj(left_term.val[2*down]) + I_number*std::conj(left_term.val[2*down+1]) ) 
			* (-right_term.val[2*n] + I_number * right_term.val[2*n+1])
			);
		}
	}
}