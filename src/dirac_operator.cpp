#include "dirac_operator.h"

c_double I_number(0, 1); //imaginary number
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
 * RightPB[2*n+1] = x+1, t (towards down)
 * LeftPB[2*n+1]  = x-1, t (towards up)
 * RightPB[2*n]   = x, t+1 (towards right)
 * LeftPB[2*n]    = x, t-1 (towards left)
 * n = x * Nt + t = (x,t) coordinates
 * 
 */

void D_phi(const spinor& U, const spinor& phi, spinor &Dphi, const double& m0) {
	using namespace LV;
	using namespace mpi;
	MPI_Status status;
	if (size == 1){	
		for (int n = 0; n < maxSize; n++) {
			//n = x * Nt + t
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
	}

	else{

	for(int n = maxSize-width_t; n < maxSize; n++){
		BottomRow.mu0[n - (maxSize-width_t)] = std::conj(U.mu1[n]) * (phi.mu0[n] - I_number * phi.mu1[n]);
		BottomRow.mu1[n - (maxSize-width_t)] = std::conj(U.mu1[n]) * (I_number*phi.mu0[n] + phi.mu1[n]);
	}
	for(int n = 0; n < width_t; n++){
		TopRow.mu0[n] = (phi.mu0[n] + I_number * phi.mu1[n]);
		TopRow.mu1[n] = (-I_number*phi.mu0[n] + phi.mu1[n]);
	}
	for(int n = width_t - 1; n<maxSize; n+=width_t){
		RightCol.mu0[n/width_t] = std::conj(U.mu0[n]) * (phi.mu0[n] + phi.mu1[n]);
		RightCol.mu1[n/width_t] = std::conj(U.mu0[n]) * (phi.mu0[n] + phi.mu1[n]);
	}
	for(int n = 0; n<maxSize; n+=width_t ){
		LeftCol.mu0[n/width_t] = phi.mu0[n] - phi.mu1[n];
		LeftCol.mu1[n/width_t] = -phi.mu0[n] + phi.mu1[n];
	}

	MPI_Send(TopRow.mu0, width_t, MPI_DOUBLE_COMPLEX, top, 0, MPI_COMM_WORLD);
	MPI_Send(TopRow.mu1, width_t, MPI_DOUBLE_COMPLEX, top, 1, MPI_COMM_WORLD);

	MPI_Recv(TopRow.mu0, width_t, MPI_DOUBLE_COMPLEX, bot, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(TopRow.mu1, width_t, MPI_DOUBLE_COMPLEX, bot, 1, MPI_COMM_WORLD, &status);
	
	MPI_Send(BottomRow.mu0, width_t, MPI_DOUBLE_COMPLEX, bot, 2, MPI_COMM_WORLD);
	MPI_Send(BottomRow.mu1, width_t, MPI_DOUBLE_COMPLEX, bot, 3, MPI_COMM_WORLD);
		
	MPI_Recv(BottomRow.mu0, width_t, MPI_DOUBLE_COMPLEX, top, 2, MPI_COMM_WORLD, &status);
	MPI_Recv(BottomRow.mu1, width_t, MPI_DOUBLE_COMPLEX, top, 3, MPI_COMM_WORLD, &status);

	MPI_Send(RightCol.mu0, width_x, MPI_DOUBLE_COMPLEX, right, 5, MPI_COMM_WORLD);
	MPI_Send(RightCol.mu1, width_x, MPI_DOUBLE_COMPLEX, right, 6, MPI_COMM_WORLD);

	MPI_Recv(RightCol.mu0, width_x, MPI_DOUBLE_COMPLEX, left, 5, MPI_COMM_WORLD, &status);
	MPI_Recv(RightCol.mu1, width_x, MPI_DOUBLE_COMPLEX, left, 6, MPI_COMM_WORLD, &status);

	MPI_Send(LeftCol.mu0, width_x, MPI_DOUBLE_COMPLEX, left, 7, MPI_COMM_WORLD);
	MPI_Send(LeftCol.mu1, width_x, MPI_DOUBLE_COMPLEX, left, 8, MPI_COMM_WORLD);

	MPI_Recv(LeftCol.mu0, width_x, MPI_DOUBLE_COMPLEX, right, 7, MPI_COMM_WORLD, &status);
	MPI_Recv(LeftCol.mu1, width_x, MPI_DOUBLE_COMPLEX, right, 8, MPI_COMM_WORLD, &status);

	
	//Update interior points (no communication needed here)
	for (int x = 1; x<width_x-1; x++) {
	for(int t = 1; t<width_t-1; t++){
		int n = x * width_t + t;
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
	}
		
	
	//First row: Receives last row from top
	for(int n = 1; n<width_t-1; n++){
		Dphi.mu0[n] = (m0 + 2) * phi.mu0[n] - 0.5 * ( 
				U.mu0[n] * SignR[2*n] * (phi.mu0[RightPB[2*n]] - phi.mu1[RightPB[2*n]])
			+   U.mu1[n] * SignR[2*n+1] * (phi.mu0[RightPB[2*n+1]] + I_number * phi.mu1[RightPB[2*n+1]])
			+   std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
			+   BottomRow.mu0[n]
			);

		Dphi.mu1[n] = (m0 + 2) * phi.mu1[n] - 0.5 * ( 
				U.mu0[n] * SignR[2*n] * (-phi.mu0[RightPB[2*n]] + phi.mu1[RightPB[2*n]])
			+   U.mu1[n] * SignR[2*n+1] * (-I_number*phi.mu0[RightPB[2*n+1]] + phi.mu1[RightPB[2*n+1]])
			+   std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
			+   BottomRow.mu1[n]
		);		
	}

	//Last row: Receives first row from bot
	for(int n = maxSize-width_t+1; n<maxSize-1; n++){
		Dphi.mu0[n] = (m0 + 2) * phi.mu0[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * (phi.mu0[RightPB[2*n]] - phi.mu1[RightPB[2*n]])
		+   U.mu1[n] * SignR[2*n+1] *TopRow.mu0[n-(maxSize-width_t)]
		+   std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
		+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (phi.mu0[LeftPB[2*n+1]] - I_number * phi.mu1[LeftPB[2*n+1]])
		);

		Dphi.mu1[n] = (m0 + 2) * phi.mu1[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * (-phi.mu0[RightPB[2*n]] + phi.mu1[RightPB[2*n]])
		+   U.mu1[n] * SignR[2*n+1] * TopRow.mu1[n-(maxSize-width_t)]
		+   std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
		+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (I_number*phi.mu0[LeftPB[2*n+1]] + phi.mu1[LeftPB[2*n+1]])
		);	
	}

	//Right column: receives left column from right
	for(int n = 2*width_t - 1 ; n<maxSize-width_t; n+=width_t){
		Dphi.mu0[n] = (m0 + 2) * phi.mu0[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * LeftCol.mu0[n/width_t]
		+   U.mu1[n] * SignR[2*n+1] * (phi.mu0[RightPB[2*n+1]] + I_number * phi.mu1[RightPB[2*n+1]])
		+   std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
		+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (phi.mu0[LeftPB[2*n+1]] - I_number * phi.mu1[LeftPB[2*n+1]])
		);

		Dphi.mu1[n] = (m0 + 2) * phi.mu1[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * LeftCol.mu1[n/width_t]
		+   U.mu1[n] * SignR[2*n+1] * (-I_number*phi.mu0[RightPB[2*n+1]] + phi.mu1[RightPB[2*n+1]])
		+   std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
		+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (I_number*phi.mu0[LeftPB[2*n+1]] + phi.mu1[LeftPB[2*n+1]])
		);	
	}
	//Left column: receives right column from left
	for(int n = width_t; n<maxSize-width_t; n+=width_t ){
		Dphi.mu0[n] = (m0 + 2) * phi.mu0[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * (phi.mu0[RightPB[2*n]] - phi.mu1[RightPB[2*n]])
		+   U.mu1[n] * SignR[2*n+1] * (phi.mu0[RightPB[2*n+1]] + I_number * phi.mu1[RightPB[2*n+1]])
		+   SignL[2*n] * RightCol.mu0[n/width_t]
		+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (phi.mu0[LeftPB[2*n+1]] - I_number * phi.mu1[LeftPB[2*n+1]])
		);

		Dphi.mu1[n] = (m0 + 2) * phi.mu1[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * (-phi.mu0[RightPB[2*n]] + phi.mu1[RightPB[2*n]])
		+   U.mu1[n] * SignR[2*n+1] * (-I_number*phi.mu0[RightPB[2*n+1]] + phi.mu1[RightPB[2*n+1]])
		+   SignL[2*n] * RightCol.mu1[n/width_t]
		+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (I_number*phi.mu0[LeftPB[2*n+1]] + phi.mu1[LeftPB[2*n+1]])
		);	
	}

	//Top-left corner
	int n = 0;
	Dphi.mu0[n] = (m0 + 2) * phi.mu0[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * (phi.mu0[RightPB[2*n]] - phi.mu1[RightPB[2*n]])
		+   U.mu1[n] * SignR[2*n+1] * (phi.mu0[RightPB[2*n+1]] + I_number * phi.mu1[RightPB[2*n+1]])
		+   SignL[2*n] * RightCol.mu0[n/width_t]
		+   BottomRow.mu0[n]
		);

	Dphi.mu1[n] = (m0 + 2) * phi.mu1[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * (-phi.mu0[RightPB[2*n]] + phi.mu1[RightPB[2*n]])
		+   U.mu1[n] * SignR[2*n+1] * (-I_number*phi.mu0[RightPB[2*n+1]] + phi.mu1[RightPB[2*n+1]])
		+   SignL[2*n] * RightCol.mu1[n/width_t]
		+   BottomRow.mu1[n]
	);	

	//Top-right corner
	n = width_t - 1;
	Dphi.mu0[n] = (m0 + 2) * phi.mu0[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * LeftCol.mu0[n/width_t]
		+   U.mu1[n] * SignR[2*n+1] * (phi.mu0[RightPB[2*n+1]] + I_number * phi.mu1[RightPB[2*n+1]])
		+   std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
		+   BottomRow.mu0[n]
	);
	Dphi.mu1[n] = (m0 + 2) * phi.mu1[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * LeftCol.mu1[n/width_t]
		+   U.mu1[n] * SignR[2*n+1] * (-I_number*phi.mu0[RightPB[2*n+1]] + phi.mu1[RightPB[2*n+1]])
		+   std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
		+   BottomRow.mu1[n]
	);	

	//Bottom-left corner
	n = maxSize-width_t;
	Dphi.mu0[n] = (m0 + 2) * phi.mu0[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * (phi.mu0[RightPB[2*n]] - phi.mu1[RightPB[2*n]])
		+   U.mu1[n] * SignR[2*n+1] * TopRow.mu0[n-(maxSize-width_t)]
		+   SignL[2*n] * RightCol.mu0[n/width_t]
		+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (phi.mu0[LeftPB[2*n+1]] - I_number * phi.mu1[LeftPB[2*n+1]])
	);

	Dphi.mu1[n] = (m0 + 2) * phi.mu1[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * (-phi.mu0[RightPB[2*n]] + phi.mu1[RightPB[2*n]])
		+   U.mu1[n] * SignR[2*n+1] * TopRow.mu1[n-(maxSize-width_t)]
		+   SignL[2*n] * RightCol.mu1[n/width_t]
		+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (I_number*phi.mu0[LeftPB[2*n+1]] + phi.mu1[LeftPB[2*n+1]])
	);
	n = maxSize-1;
	//Bottom-right corner
	Dphi.mu0[n] = (m0 + 2) * phi.mu0[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * LeftCol.mu0[n/width_t]
		+   U.mu1[n] * SignR[2*n+1] * TopRow.mu0[n-(maxSize-width_t)]
		+   std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
		+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (phi.mu0[LeftPB[2*n+1]] - I_number * phi.mu1[LeftPB[2*n+1]])
	);

	Dphi.mu1[n] = (m0 + 2) * phi.mu1[n] - 0.5 * ( 
			U.mu0[n] * SignR[2*n] * LeftCol.mu0[n/width_t]
		+   U.mu1[n] * SignR[2*n+1] * TopRow.mu1[n-(maxSize-width_t)]
		+   std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
		+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (I_number*phi.mu0[LeftPB[2*n+1]] + phi.mu1[LeftPB[2*n+1]])
	);	


	} //end else (size>1)


}

void D_dagger_phi(const spinor& U, const spinor& phi, spinor &Dphi,const double& m0) {
	using namespace LV;
	using namespace mpi;
	MPI_Status status;
	
	if (size == 1){
		for (int n = 0; n < Ntot; n++) {
			//n = x * Nt + t
			Dphi.mu0[n] = (m0 + 2) * phi.mu0[n] -0.5 * ( 
				std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] - phi.mu1[LeftPB[2*n]])
			+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (phi.mu0[LeftPB[2*n+1]] + I_number * phi.mu1[LeftPB[2*n+1]])
			+   U.mu0[n] * SignR[2*n] * (phi.mu0[RightPB[2*n]] + phi.mu1[RightPB[2*n]])
			+	U.mu1[n] * SignR[2*n+1] * (phi.mu0[RightPB[2*n+1]] - I_number * phi.mu1[RightPB[2*n+1]])
			);

			Dphi.mu1[n] = (m0 + 2) * phi.mu1[n] -0.5 * ( 
				std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (-phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
			+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (-I_number*phi.mu0[LeftPB[2*n+1]] + phi.mu1[LeftPB[2*n+1]])
			+   U.mu0[n] * SignR[2*n] * (phi.mu0[RightPB[2*n]] + phi.mu1[RightPB[2*n]])
			+   U.mu1[n] * SignR[2*n+1] * (I_number*phi.mu0[RightPB[2*n+1]] + phi.mu1[RightPB[2*n+1]])
			);
		}
	}

	else{

	for(int n = maxSize-width_t; n < maxSize; n++){
		BottomRow.mu0[n - (maxSize-Nt)] = std::conj(U.mu1[n]) * (phi.mu0[n] + I_number * phi.mu1[n]);
		BottomRow.mu1[n - (maxSize-Nt)] = std::conj(U.mu1[n]) * (-I_number*phi.mu0[n] + phi.mu1[n]);
	}
	for(int n = 0; n < width_t; n++){
		TopRow.mu0[n] = (phi.mu0[n] - I_number * phi.mu1[n]);
		TopRow.mu1[n] = (I_number*phi.mu0[n] + phi.mu1[n]);
	}
	for(int n = width_t - 1; n<maxSize; n+=width_t){
		RightCol.mu0[n/width_t] = std::conj(U.mu0[n]) * (phi.mu0[n] + phi.mu1[n]);
		RightCol.mu1[n/width_t] = std::conj(U.mu0[n]) * (phi.mu0[n] + phi.mu1[n]);
	}
	for(int n = 0; n<maxSize; n+=width_t ){
		LeftCol.mu0[n/width_t] = phi.mu0[n] - phi.mu1[n];
		LeftCol.mu1[n/width_t] = -phi.mu0[n] + phi.mu1[n];
	}

	MPI_Send(TopRow.mu0, width_t, MPI_DOUBLE_COMPLEX, top, 0, MPI_COMM_WORLD);
	MPI_Send(TopRow.mu1, width_t, MPI_DOUBLE_COMPLEX, top, 1, MPI_COMM_WORLD);

	MPI_Recv(TopRow.mu0, width_t, MPI_DOUBLE_COMPLEX, bot, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(TopRow.mu1, width_t, MPI_DOUBLE_COMPLEX, bot, 1, MPI_COMM_WORLD, &status);
	
	MPI_Send(BottomRow.mu0, width_t, MPI_DOUBLE_COMPLEX, bot, 2, MPI_COMM_WORLD);
	MPI_Send(BottomRow.mu1, width_t, MPI_DOUBLE_COMPLEX, bot, 3, MPI_COMM_WORLD);
		
	MPI_Recv(BottomRow.mu0, width_t, MPI_DOUBLE_COMPLEX, top, 2, MPI_COMM_WORLD, &status);
	MPI_Recv(BottomRow.mu1, width_t, MPI_DOUBLE_COMPLEX, top, 3, MPI_COMM_WORLD, &status);

	MPI_Send(RightCol.mu0, width_x, MPI_DOUBLE_COMPLEX, right, 5, MPI_COMM_WORLD);
	MPI_Send(RightCol.mu1, width_x, MPI_DOUBLE_COMPLEX, right, 6, MPI_COMM_WORLD);

	MPI_Recv(RightCol.mu0, width_x, MPI_DOUBLE_COMPLEX, left, 5, MPI_COMM_WORLD, &status);
	MPI_Recv(RightCol.mu1, width_x, MPI_DOUBLE_COMPLEX, left, 6, MPI_COMM_WORLD, &status);

	MPI_Send(LeftCol.mu0, width_x, MPI_DOUBLE_COMPLEX, left, 7, MPI_COMM_WORLD);
	MPI_Send(LeftCol.mu1, width_x, MPI_DOUBLE_COMPLEX, left, 8, MPI_COMM_WORLD);

	MPI_Recv(LeftCol.mu0, width_x, MPI_DOUBLE_COMPLEX, right, 7, MPI_COMM_WORLD, &status);
	MPI_Recv(LeftCol.mu1, width_x, MPI_DOUBLE_COMPLEX, right, 8, MPI_COMM_WORLD, &status);

	//Update interior points (no communication needed here)
	for (int n = Nt; n < maxSize-Nt; n++) {
		Dphi.mu0[n] = (m0 + 2) * phi.mu0[n] -0.5 * ( 
				std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] - phi.mu1[LeftPB[2*n]])
			+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (phi.mu0[LeftPB[2*n+1]] + I_number * phi.mu1[LeftPB[2*n+1]])
			+   U.mu0[n] * SignR[2*n] * (phi.mu0[RightPB[2*n]] + phi.mu1[RightPB[2*n]])
			+	U.mu1[n] * SignR[2*n+1] * (phi.mu0[RightPB[2*n+1]] - I_number * phi.mu1[RightPB[2*n+1]])
		);

			Dphi.mu1[n] = (m0 + 2) * phi.mu1[n] -0.5 * ( 
				std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (-phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
			+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (-I_number*phi.mu0[LeftPB[2*n+1]] + phi.mu1[LeftPB[2*n+1]])
			+   U.mu0[n] * SignR[2*n] * (phi.mu0[RightPB[2*n]] + phi.mu1[RightPB[2*n]])
			+   U.mu1[n] * SignR[2*n+1] * (I_number*phi.mu0[RightPB[2*n+1]] + phi.mu1[RightPB[2*n+1]])
		);
	}
		
	
	//First row: Receives last row from rank-1
	for(int n = 0; n<Nt; n++){
		Dphi.mu0[n] = (m0 + 2) * phi.mu0[n] -0.5 * ( 
				std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] - phi.mu1[LeftPB[2*n]])
			+   BottomRow.mu0[n]
			+   U.mu0[n] * SignR[2*n] * (phi.mu0[RightPB[2*n]] + phi.mu1[RightPB[2*n]])
			+	U.mu1[n] * SignR[2*n+1] * (phi.mu0[RightPB[2*n+1]] - I_number * phi.mu1[RightPB[2*n+1]])
		);

		Dphi.mu1[n] = (m0 + 2) * phi.mu1[n] -0.5 * ( 
				std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (-phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
			+   BottomRow.mu1[n]
			+   U.mu0[n] * SignR[2*n] * (phi.mu0[RightPB[2*n]] + phi.mu1[RightPB[2*n]])
			+   U.mu1[n] * SignR[2*n+1] * (I_number*phi.mu0[RightPB[2*n+1]] + phi.mu1[RightPB[2*n+1]])
		);
	}

	//Last row: Receives first row from rank+1
	for(int n = maxSize-Nt; n<maxSize; n++){
		Dphi.mu0[n] = (m0 + 2) * phi.mu0[n] -0.5 * ( 
				std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (phi.mu0[LeftPB[2*n]] - phi.mu1[LeftPB[2*n]])
			+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (phi.mu0[LeftPB[2*n+1]] + I_number * phi.mu1[LeftPB[2*n+1]])
			+   U.mu0[n] * SignR[2*n] * (phi.mu0[RightPB[2*n]] + phi.mu1[RightPB[2*n]])
			+	U.mu1[n] * TopRow.mu0[n - (maxSize-Nt)]
		);

		Dphi.mu1[n] = (m0 + 2) * phi.mu1[n] -0.5 * ( 
				std::conj(U.mu0[LeftPB[2*n]]) * SignL[2*n] * (-phi.mu0[LeftPB[2*n]] + phi.mu1[LeftPB[2*n]])
			+   std::conj(U.mu1[LeftPB[2*n+1]]) * SignL[2*n+1] * (-I_number*phi.mu0[LeftPB[2*n+1]] + phi.mu1[LeftPB[2*n+1]])
			+   U.mu0[n] * SignR[2*n] * (phi.mu0[RightPB[2*n]] + phi.mu1[RightPB[2*n]])
			+   U.mu1[n] * TopRow.mu1[n - (maxSize-Nt)]
		);
		
	}


	} //end else




}


//D D^dagger phi 
void D_D_dagger_phi(const spinor& U, const spinor& phi, spinor &Dphi,const double& m0) {
	D_dagger_phi(U, phi, DTEMP, m0);
	D_phi(U,  DTEMP, Dphi, m0);
}



//2* Re ( left^dag \partial D / \partial omega(z) right )
re_field phi_dag_partialD_phi(const spinor& U, const spinor& left,const spinor& right){
	using namespace LV;
	using namespace mpi;
	re_field Dphi(mpi::maxSize); 
	MPI_Status status;

	if (size == 1){
		#pragma omp parallel for
		for (int n = 0; n < maxSize; n++) {
			//n = x * Nt + t
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
	}

	else{
		//Update mu0 component (does not need communication)
		for (int n = 0; n < maxSize; n++) {
			Dphi.mu0[n] = std::imag(
				U.mu0[n] * SignR[2*n] * ( std::conj(left.mu0[n] - left.mu1[n]) ) * (right.mu0[RightPB[2*n]] - right.mu1[RightPB[2*n]])
				- std::conj(U.mu0[n]) * SignR[2*n] * ( std::conj(left.mu0[RightPB[2*n]] + left.mu1[RightPB[2*n]]) ) * (right.mu0[n] + right.mu1[n])
			);
		}	
	

	for(int n = 0; n < Nt; n++){
		TopRow.mu0[n] = right.mu0[n] + I_number * right.mu1[n];
		TopRow.mu1[n] = std::conj(left.mu0[n]) + I_number*std::conj(left.mu1[n]); 
	}

	MPI_Send(TopRow.mu0, Nt, MPI_DOUBLE_COMPLEX, mod(rank-1,size), 0, MPI_COMM_WORLD);
	MPI_Send(TopRow.mu1, Nt, MPI_DOUBLE_COMPLEX, mod(rank-1,size), 1, MPI_COMM_WORLD);

	MPI_Recv(TopRow.mu0, Nt, MPI_DOUBLE_COMPLEX, mod(rank+1,size), 0, MPI_COMM_WORLD, &status);
	MPI_Recv(TopRow.mu1, Nt, MPI_DOUBLE_COMPLEX, mod(rank+1,size), 1, MPI_COMM_WORLD, &status);
	


	//Interior points and first row: no communication needed here
	for (int n = 0; n < maxSize-Nt; n++) {
		Dphi.mu1[n] = std::imag(
			U.mu1[n] * SignR[2*n+1] * ( std::conj(left.mu0[n]) - I_number*std::conj(left.mu1[n]) ) * (right.mu0[RightPB[2*n+1]] + I_number * right.mu1[RightPB[2*n+1]])
			+ std::conj(U.mu1[n]) * SignR[2*n+1] * ( std::conj(left.mu0[RightPB[2*n+1]]) + I_number*std::conj(left.mu1[RightPB[2*n+1]]) ) 
			* (-right.mu0[n] + I_number * right.mu1[n])
		);
	}	

	//Last row: Receives first row from rank+1
	for(int n = maxSize-Nt; n<maxSize; n++){
		//mu = 1
		Dphi.mu1[n] = std::imag(
			U.mu1[n] * SignR[2*n+1] * ( std::conj(left.mu0[n]) - I_number*std::conj(left.mu1[n]) ) * TopRow.mu0[n-(maxSize-Nt)]
			+ std::conj(U.mu1[n]) * SignR[2*n+1] * TopRow.mu1[n-(maxSize-Nt)] 
			* (-right.mu0[n] + I_number * right.mu1[n])
		);
	}

	} //end else (size>1)

	return Dphi;
 }