#include "dirac_operator.h"

c_double I_number(0, 1); //imaginary number


void D_phi(const spinor& U, const spinor& phi, spinor &Dphi, const double& m0) {
	using namespace LV;
	using namespace mpi;
	MPI_Status status;

	//n = x * Nt + t
	for(int n = maxSize-Nt; n < maxSize; n++){
		BottomRow.mu0[n - (maxSize-Nt)] = (phi.mu0[n] + I_number * phi.mu1[n]);
		BottomRow.mu1[n - (maxSize-Nt)] = (-I_number*phi.mu0[n] + phi.mu1[n]);
	}
	for(int n = 0; n < Nt; n++){
		TopRow.mu0[n] = (phi.mu0[n] - I_number * phi.mu1[n]);
		TopRow.mu1[n] = (I_number*phi.mu0[n] + phi.mu1[n]);
	}

	if (size > 1){
	//Send last row to rank+1 
	//std::cout << "Rank " << rank << " send last row to rank " << mod(rank+1,size) << BottomRow.mu0[0] << std::endl;
	MPI_Send(BottomRow.mu0, Nt, MPI_DOUBLE_COMPLEX, mod(rank+1,size), tagTop, MPI_COMM_WORLD);
	MPI_Send(BottomRow.mu1, Nt, MPI_DOUBLE_COMPLEX, mod(rank+1,size), tagTop, MPI_COMM_WORLD);
	
	//Send first row to rank-1
	MPI_Send(TopRow.mu0, Nt, MPI_DOUBLE_COMPLEX, mod(rank-1,size), tagBottom, MPI_COMM_WORLD);
	MPI_Send(TopRow.mu1, Nt, MPI_DOUBLE_COMPLEX, mod(rank-1,size), tagBottom, MPI_COMM_WORLD);

	//Recieve last row from rank-1
	MPI_Recv(BottomRow.mu0, Nt, MPI_DOUBLE_COMPLEX, mod(rank-1,size), tagTop, MPI_COMM_WORLD, &status);
	MPI_Recv(BottomRow.mu1, Nt, MPI_DOUBLE_COMPLEX, mod(rank-1,size), tagTop, MPI_COMM_WORLD, &status);
	//std::cout << "Rank " << rank << " received last row from rank " << mod(rank-1,size) << BottomRow.mu0[0] << std::endl;
	//Recieve first row from rank+1 
	MPI_Recv(TopRow.mu0, Nt, MPI_DOUBLE_COMPLEX, mod(rank+1,size), tagBottom, MPI_COMM_WORLD, &status);
	MPI_Recv(TopRow.mu1, Nt, MPI_DOUBLE_COMPLEX, mod(rank+1,size), tagBottom, MPI_COMM_WORLD, &status);
	}
	
	//Update interior points (no communication needed here)
	for (int n = Nt; n < maxSize-Nt; n++) {
		Dphi.mu0[n] = (m0 + 2) * phi.mu0[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * (phi.mu0[RightPB[2*n]] - phi.mu1[RightPB[2*n]])
		+   U.mu1[n] * SignR[2*n+1] * (phi.mu0[RightPB[2*n+1]] + I_number * phi.mu1[RightPB[2*n+1]])
		+   std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
		+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (phi.mu0[LeftPB[2*n+1]] - I_number * phi.mu1[LeftPB[2*n+1]])
		);

		Dphi.mu1[n] = (m0 + 2) * phi.mu1[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * (-phi.mu0[RightPB[2*n]] + phi.mu1[RightPB[2*n]])
		+   U.mu1[n] * SignR[2*n+1] * (-I_number*phi.mu0[RightPB[2*n+1]] + phi.mu1[RightPB[2*n+1]])
		+   std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
		+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (I_number*phi.mu0[LeftPB[2*n+1]] + phi.mu1[LeftPB[2*n+1]])
		);	
	}
		
	//First row: Receives last row from rank-1
	for(int n = 0; n<Nt; n++){
		Dphi.mu0[n] = (m0 + 2) * phi.mu0[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * (phi.mu0[RightPB[2*n]] - phi.mu1[RightPB[2*n]])
			//Receives (phi.mu0[RightPB[2*n+1]] + I_number * phi.mu1[RightPB[2*n+1]]) from rank-1
			+   U.mu1[n] * BottomRow.mu0[n]	
			+   std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
			+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (phi.mu0[LeftPB[2*n+1]] - I_number * phi.mu1[LeftPB[2*n+1]])
			);

		Dphi.mu1[n] = (m0 + 2) * phi.mu1[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * (-phi.mu0[RightPB[2*n]] + phi.mu1[RightPB[2*n]])
			//Receives  (-I_number*phi.mu0[RightPB[2*n+1]] + phi.mu1[RightPB[2*n+1]]) from rank-1
			+   U.mu1[n] * BottomRow.mu1[n]	
			+   std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
			+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (I_number*phi.mu0[LeftPB[2*n+1]] + phi.mu1[LeftPB[2*n+1]])
		);
	}

	//Last row: Receives first row from rank+1
	for(int n = maxSize-Nt; n<maxSize; n++){
		Dphi.mu0[n] = (m0 + 2) * phi.mu0[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * (phi.mu0[RightPB[2*n]] - phi.mu1[RightPB[2*n]])
		+   U.mu1[n] * SignR[2*n+1] * (phi.mu0[RightPB[2*n+1]] + I_number * phi.mu1[RightPB[2*n+1]])
		+   std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
		+   std::conj(U.mu1[LeftPB[2*n+1]]) * TopRow.mu0[n]
		);

		Dphi.mu1[n] = (m0 + 2) * phi.mu1[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * (-phi.mu0[RightPB[2*n]] + phi.mu1[RightPB[2*n]])
		+   U.mu1[n] * SignR[2*n+1] * (-I_number*phi.mu0[RightPB[2*n+1]] + phi.mu1[RightPB[2*n+1]])
		+   std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
		+   std::conj(U.mu1[LeftPB[2*n+1]]) * TopRow.mu1[n]
		);	
	}
	
}

void D_dagger_phi(const spinor& U, const spinor& phi, spinor &Dphi,const double& m0) {
	using namespace LV;

	#pragma omp parallel for
	for (int n = 0; n < Ntot; n++) {
		//n = x * Nt + t
		Dphi.mu0[n] = (m0 + 2) * phi.mu0[n] -0.5 * ( 
			std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] - phi.mu1[LeftPB[2*n]])
		+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (phi.mu0[LeftPB[2*n+1]] + I_number * phi.mu1[LeftPB[2*n+1]])
		+   U.mu0[n] * SignR[2*n] * (phi.mu0[RightPB[2*n]] + phi.mu1[RightPB[2*n]])
		+   U.mu1[n] * SignR[2*n+1] * (phi.mu0[RightPB[2*n+1]] - I_number * phi.mu1[RightPB[2*n+1]])
		);

		Dphi.mu1[n] = (m0 + 2) * phi.mu1[n] -0.5 * ( 
			std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (-phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
		+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (-I_number*phi.mu0[LeftPB[2*n+1]] + phi.mu1[LeftPB[2*n+1]])
		+   U.mu0[n] * SignR[2*n] * (phi.mu0[RightPB[2*n]] + phi.mu1[RightPB[2*n]])
		+   U.mu1[n] * SignR[2*n+1] * (I_number*phi.mu0[RightPB[2*n+1]] + phi.mu1[RightPB[2*n+1]])
		);
		
	}

}


//D D^dagger phi 
void D_D_dagger_phi(const spinor& U, const spinor& phi, spinor &Dphi,const double& m0) {
	D_dagger_phi(U, phi, DTEMP, m0);
	D_phi(U,  DTEMP, Dphi, m0);
}



//2* Re ( left^dag \partial D / \partial omega(z) right )
re_field phi_dag_partialD_phi(const spinor& U, const spinor& left,const spinor& right){
	using namespace LV;
	re_field Dphi; 

	//#pragma omp parallel for
	for (int n = 0; n < Ntot; n++) {
		//mu = 0
		Dphi.mu0[n] = std::imag(
		U.mu0[n] * SignR[2*n] * ( std::conj(left.mu0[n] - left.mu1[n]) ) * (right.mu0[RightPB[2*n]] - right.mu1[RightPB[2*n]])
		- std::conj(U.mu0[n]) * SignR[2*n] * ( std::conj(left.mu0[RightPB[2*n]] + left.mu1[RightPB[2*n]]) ) * (right.mu0[n] + right.mu1[n])
		);

		//mu = 1
		Dphi.mu1[n] = std::imag(
		U.mu1[n] * SignR[2*n+1] * ( std::conj(left.mu0[n]) - I_number*std::conj(left.mu1[n]) ) * (right.mu0[RightPB[2*n+1]] + I_number * right.mu1[RightPB[2*n+1]])
		+ std::conj(U.mu1[n]) * SignR[2*n+1] * ( std::conj(left.mu0[RightPB[2*n+1]]) + I_number*std::conj(left.mu1[RightPB[2*n+1]]) ) 
				* (-right.mu0[n] + I_number * right.mu1[n])
		);

	}

	return Dphi;
 }

