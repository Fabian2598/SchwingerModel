#include "gauge_conf.h"

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

c_double RandomU1() {
	//Random angle in (0,2*pi) with uniform distribution 
	double cociente = ((double) rand() / (RAND_MAX));
    double theta = 2.0*pi * cociente;
	c_double z(cos(theta), sin(theta));
	return z;
}

void GaugeConf::initialization() {
	for (int n = 0; n < mpi::maxSize; n++) {
		Conf.mu0[n] = RandomU1(); 
		Conf.mu1[n] = RandomU1(); 
	}
}


//U_mv(x) = U_m(x) U_v(x+m) U*_m(x+v) U*_v(x)
//mu = 0 time direction, mu = 1 space direction
void GaugeConf::Compute_Plaquette01() {
    MPI_Status status;
	//U_01(n) = U_0(n) U_1(n+0) U*_0(n+1) U*_1(n)
    if (mpi::size == 1) {
        for (int n = 0; n < mpi::maxSize; n++){
            Plaquette01[n] = Conf.mu0[n] * Conf.mu1[RightPB[2*n]] * std::conj(Conf.mu0[RightPB[2*n+1]]) * std::conj(Conf.mu1[n]);
        }
    }
    else{

        for(int n = 0; n<mpi::width_t; n++)
            TopRow.mu0[n] = Conf.mu0[n];

        for(int n = 0; n < mpi::maxSize; n+=mpi::width_t)
            LeftCol.mu0[n/mpi::width_t] = Conf.mu1[n];
        
        MPI_Send(TopRow.mu0, mpi::width_t, MPI_DOUBLE_COMPLEX, mpi::top, 0, mpi::cart_comm);
	    MPI_Recv(TopRow.mu0, mpi::width_t, MPI_DOUBLE_COMPLEX, mpi::bot, 0, mpi::cart_comm, &status);

        MPI_Send(LeftCol.mu0, mpi::width_x, MPI_DOUBLE_COMPLEX, mpi::left, 1, mpi::cart_comm);
	    MPI_Recv(LeftCol.mu0, mpi::width_x, MPI_DOUBLE_COMPLEX, mpi::right, 1, mpi::cart_comm, &status);
        
        //Interior points (without last row and right column)
        for(int x = 0; x < mpi::width_x-1; x++){
        for(int t = 0; t < mpi::width_t-1;t++){
            int n = x * mpi::width_t + t;
            Plaquette01[n] = Conf.mu0[n] * Conf.mu1[RightPB[2*n]] * std::conj(Conf.mu0[RightPB[2*n+1]]) * std::conj(Conf.mu1[n]);
        }
        }
        //U_01(n) = U_0(n) U_1(n+0) conj(TopRow) U*_1(n)
        //last row 
        for(int n = mpi::maxSize-mpi::width_t; n<mpi::maxSize; n++){
            Plaquette01[n] = Conf.mu0[n] * Conf.mu1[RightPB[2*n]] * std::conj(TopRow.mu0[n-(mpi::maxSize-mpi::width_t)]) * std::conj(Conf.mu1[n]);
        }
        //U_01(n) = U_0(n) LeftCol U*_0(n+1)  U*_1(n)
        //right column
        for(int n = mpi::width_t-1; n<mpi::maxSize; n+=mpi::width_t){
            Plaquette01[n] = Conf.mu0[n] * LeftCol.mu0[n/mpi::width_t] * std::conj(Conf.mu0[RightPB[2*n+1]]) * std::conj(Conf.mu1[n]);
        }
        //Bottom-right corner
        int n = mpi::maxSize-1;
        Plaquette01[n] = Conf.mu0[n] * LeftCol.mu0[n/mpi::width_t]  * std::conj(TopRow.mu0[n-(mpi::maxSize-mpi::width_t)]) * std::conj(Conf.mu1[n]);

    }
}


//Compute staple at coordinate (x,t) in the mu-direction
void GaugeConf::Compute_Staple() {
    using namespace LV;
    using namespace mpi;
    //U_v(x) U_m(x+v) U*_v(x+m) + U*_v(x-v) U_m(x-v) U_v(x+m-v)
    //mu = 0 time direction, mu = 1 space direction
    int x1, x_1, t1, t_1;
    if (size == 1){
        for (int n = 0; n < LV::Ntot; n++) {
            //These coordinates could change depending on the conventions 
		    x1 = RightPB[2*n+1];  //Coords[modulo(x + 1, Ns) ,t]
		    x_1 = LeftPB[2*n+1];  //Coords[modulo(x - 1, Ns) ,t]
		    t1 = RightPB[2*n];  //Coords[x, modulo(t + 1, Nt)]
		    t_1 = LeftPB[2*n];  //Coords[x, modulo(t - 1, Nt)]
            for (int mu = 0; mu < 2; mu++) {
                if (mu == 0) {
                    //U_1(n) U_0(n+1) U*_1(n+0) + U*_1(n-1) U_0(n-1) U_1(n-1+0)
                    const c_double& conf1 = Conf.mu1[n]; 
                    const c_double& conf2 = Conf.mu0[x1]; 
                    const c_double& conf3 = Conf.mu1[t1];
                    const c_double& conf4 = Conf.mu1[x_1];
                    const c_double& conf5 = Conf.mu0[x_1];
                    const c_double& conf6 = Conf.mu1[x_1_t1[n]]; //Coords[mod(x - 1, Ns)][mod(t + 1, Nt)];
                    Staples.mu0[n] = conf1 * conf2 * std::conj(conf3) +
                        std::conj(conf4) * conf5 * conf6;
                }
                else {
                    //U_0(n) U_1(n+0) U*_0(n+1) + U*_0(n-0) U_1(n-0) U_0(n+1-0)
                    const c_double& conf1 = Conf.mu0[n];
                    const c_double& conf2 = Conf.mu1[t1];
                    const c_double& conf3 = Conf.mu0[x1];
                    const c_double& conf4 = Conf.mu0[t_1];
                    const c_double& conf5 = Conf.mu1[t_1];
                    const c_double& conf6 = Conf.mu0[x1_t_1[n]]; //Coords[mod(x + 1, Ns)][mod(t - 1, Nt)];
                    Staples.mu1[n] = conf1 * conf2 * std::conj(conf3) +
                        std::conj(conf4) * conf5 * conf6;
                }
            }
        
        }
    }
    else{
        MPI_Status status;

        //-----------------mu = 0 direction--------------------//
        //Interior points
        //Includes left column
        for (int x = 1; x<width_x-1; x++) {
        for(int t = 0; t<width_t-1; t++){
            int n = x * width_t + t;
		    x1 = RightPB[2*n+1];  
		    x_1 = LeftPB[2*n+1]; 
		    t1 = RightPB[2*n];  
		    t_1 = LeftPB[2*n];  
            //U_1(n) U_0(n+1) U*_1(n+0) + U*_1(n-1) U_0(n-1) U_1(n-1+0)
            const c_double& conf1 = Conf.mu1[n]; 
            const c_double& conf2 = Conf.mu0[x1]; 
            const c_double& conf3 = Conf.mu1[t1];
            const c_double& conf4 = Conf.mu1[x_1];
            const c_double& conf5 = Conf.mu0[x_1];
            const c_double& conf6 = Conf.mu1[x_1_t1[n]];
            Staples.mu0[n] = conf1 * conf2 * std::conj(conf3) +
                std::conj(conf4) * conf5 * conf6;
        }
        }

        //Sending top row to previous rank
        for( int n = 0; n < width_t; n++){
            TopRow.mu0[n] = Conf.mu0[n]; //U_0(n+1)
        }
        //Sending bottom row.
        for( int n = maxSize-width_t; n < maxSize; n++){
            BottomRow.mu0[n-(maxSize-width_t)] = std::conj(Conf.mu1[n]) * Conf.mu0[n]; // U*_1(n-1) U_0(n-1) 
            BottomRow.mu1[n-(maxSize-width_t)] = Conf.mu1[RightPB[2*n]]; //U_1(n-1+0)
        }
        //Sending left column
        for( int n = 0; n < maxSize; n+=width_t){
            LeftCol.mu0[n/width_t] = std::conj(Conf.mu1[n]); //U*_1(n+0)
            LeftCol.mu1[n/width_t] = Conf.mu1[LeftPB[2*n+1]]; //U_1(n-1+0)
        }

        MPI_Send(TopRow.mu0, width_t, MPI_DOUBLE_COMPLEX, top, 0, cart_comm);
        MPI_Recv(TopRow.mu0, width_t, MPI_DOUBLE_COMPLEX, bot, 0, cart_comm, &status);

        MPI_Send(BottomRow.mu0, width_t, MPI_DOUBLE_COMPLEX, bot, 1, cart_comm);
        MPI_Recv(BottomRow.mu0, width_t, MPI_DOUBLE_COMPLEX, top, 1, cart_comm, &status);

        MPI_Send(BottomRow.mu1, width_t, MPI_DOUBLE_COMPLEX, bot, 2, cart_comm);
        MPI_Recv(BottomRow.mu1, width_t, MPI_DOUBLE_COMPLEX, top, 2, cart_comm, &status);

        MPI_Send(LeftCol.mu0, width_x, MPI_DOUBLE_COMPLEX, left, 3, cart_comm);
        MPI_Recv(LeftCol.mu0, width_x, MPI_DOUBLE_COMPLEX, right, 3, cart_comm, &status);

        MPI_Send(LeftCol.mu1, width_x, MPI_DOUBLE_COMPLEX, left, 4, cart_comm);
        MPI_Recv(LeftCol.mu1, width_x, MPI_DOUBLE_COMPLEX, right, 4, cart_comm, &status);
        
        //Update first row (except for last point)
        for(int n = 0; n < width_t-1; n++){
            //U_1(n) U_0(n+1) U*_1(n+0) + U*_1(n-1) U_0(n-1) U_1(n-1+0)
            x1 = RightPB[2*n+1];
		    t1 = RightPB[2*n];  //n+0
            const c_double& conf1 = Conf.mu1[n]; 
            const c_double& conf2 = Conf.mu0[x1];
            const c_double& conf3 = Conf.mu1[t1];
            Staples.mu0[n] = conf1 * conf2 * std::conj(conf3) +
                BottomRow.mu0[n] * BottomRow.mu1[n];
        }
        
        //Update last row (except for last point)
        for(int n = maxSize-width_t; n < maxSize-1; n++){
            //U_1(n) U_0(n+1) U*_1(n+0) + U*_1(n-1) U_0(n-1) U_1(n-1+0)
		    x_1 = LeftPB[2*n+1]; 
		    t1 = RightPB[2*n];              
            const c_double& conf1 = Conf.mu1[n]; 
            const c_double& conf3 = Conf.mu1[t1];
            const c_double& conf4 = Conf.mu1[x_1];
            const c_double& conf5 = Conf.mu0[x_1];
            const c_double& conf6 = Conf.mu1[x_1_t1[n]]; 
            Staples.mu0[n] = conf1 * TopRow.mu0[n-(maxSize-width_t)] * std::conj(conf3) +
                std::conj(conf4) * conf5 * conf6;
        }

        //Update right column (except for first point)
        for (int n = 2*width_t-1; n<maxSize-width_t; n+=width_t) {
            //U_1(n) U_0(n+1) U*_1(n+0) + U*_1(n-1) U_0(n-1) U_1(n-1+0)
		    x_1 = LeftPB[2*n+1]; 
            x1 = RightPB[2*n+1];   
            const c_double& conf1 = Conf.mu1[n]; 
            const c_double& conf2 = Conf.mu0[x1]; 
            const c_double& conf4 = Conf.mu1[x_1];
            const c_double& conf5 = Conf.mu0[x_1];
            Staples.mu0[n] = conf1 * conf2 * LeftCol.mu0[n/width_t] +
                 std::conj(conf4) * conf5 * LeftCol.mu1[n/width_t];
        }
        //------------------------------------------------------------//

        //Update top-right corner
        c_double bottom_left = Conf.mu1[maxSize-width_t]; //U_1(n-1+0)
        MPI_Send(&bottom_left, 1, MPI_DOUBLE_COMPLEX, bot_left, 5, cart_comm);
        MPI_Recv(&bottom_left, 1, MPI_DOUBLE_COMPLEX, top_right, 5, cart_comm, &status);
        int n = width_t - 1;
        x1 = RightPB[2*n+1];  
		x_1 = LeftPB[2*n+1]; 
		t1 = RightPB[2*n];  
		t_1 = LeftPB[2*n];  
        //U_1(n) U_0(n+1) U*_1(n+0) + U*_1(n-1) U_0(n-1) U_1(n-1+0)
        {
        const c_double& conf1 = Conf.mu1[n]; 
        const c_double& conf2 = Conf.mu0[x1];
        Staples.mu0[n] = conf1 * conf2 * LeftCol.mu0[n/width_t] +
            BottomRow.mu0[n] * bottom_left;
        }

        //Update bottom-right corner
        n = maxSize - 1;
        x1 = RightPB[2*n+1];  
		x_1 = LeftPB[2*n+1]; 
		t1 = RightPB[2*n];  
		t_1 = LeftPB[2*n];  
        //U_1(n) U_0(n+1) U*_1(n+0) + U*_1(n-1) U_0(n-1) U_1(n-1+0)
        {
        const c_double& conf1 = Conf.mu1[n]; 
        const c_double& conf4 = Conf.mu1[x_1];
        const c_double& conf5 = Conf.mu0[x_1];
        Staples.mu0[n] = conf1 * TopRow.mu0[n-(maxSize-width_t)] * LeftCol.mu0[n/width_t] +
            std::conj(conf4) * conf5 * LeftCol.mu1[n/width_t];
        }


//----------------------------mu = 1 direction-----------------------------//
        //Interior points (includes first row)
        for(int x = 0; x < width_x-1; x++) {
        for(int t = 1; t < width_t-1; t++) {  
            int n = x * width_t + t;
            //U_0(n) U_1(n+0) U*_0(n+1) + U*_0(n-0) U_1(n-0) U_0(n+1-0)
            x1 = RightPB[2*n+1];  
		    t1 = RightPB[2*n];  
		    t_1 = LeftPB[2*n];  
            const c_double& conf1 = Conf.mu0[n];
            const c_double& conf2 = Conf.mu1[t1];
            const c_double& conf3 = Conf.mu0[x1];
            const c_double& conf4 = Conf.mu0[t_1];
            const c_double& conf5 = Conf.mu1[t_1];
            const c_double& conf6 = Conf.mu0[x1_t_1[n]]; 
            Staples.mu1[n] = conf1 * conf2 * std::conj(conf3) +
                std::conj(conf4) * conf5 * conf6;
        }
        }

        for(int n = 0; n < width_t; n++){
            TopRow.mu0[n] = std::conj(Conf.mu0[n]); //U*_0(n+1)
            TopRow.mu1[n] = Conf.mu0[LeftPB[2*n]]; //U_0(n+1-0)
        }
        for(int n = width_t-1; n<maxSize; n+=width_t){
            RightCol.mu0[n/width_t] = std::conj(Conf.mu0[n]) * Conf.mu1[n]; //U*_0(n-0) U_1(n-0)
            RightCol.mu1[n/width_t] = Conf.mu0[RightPB[2*n+1]]; //U_0(n+1-0)
        }

        for(int n = 0; n<maxSize; n+=width_t){
            LeftCol.mu0[n/width_t] = Conf.mu1[n]; //U_1(n+0)
        }


        MPI_Send(TopRow.mu0, width_t, MPI_DOUBLE_COMPLEX, top, 0, cart_comm);
        MPI_Recv(TopRow.mu0, width_t, MPI_DOUBLE_COMPLEX, bot, 0, cart_comm, &status);

        MPI_Send(TopRow.mu1, width_t, MPI_DOUBLE_COMPLEX, top, 1, cart_comm);
        MPI_Recv(TopRow.mu1, width_t, MPI_DOUBLE_COMPLEX, bot, 1, cart_comm, &status);

        MPI_Send(RightCol.mu0, width_x, MPI_DOUBLE_COMPLEX, right, 2, cart_comm);
        MPI_Recv(RightCol.mu0, width_x, MPI_DOUBLE_COMPLEX, left, 2, cart_comm, &status);

        MPI_Send(RightCol.mu1, width_x, MPI_DOUBLE_COMPLEX, right, 3, cart_comm);
        MPI_Recv(RightCol.mu1, width_x, MPI_DOUBLE_COMPLEX, left, 3, cart_comm, &status);

        MPI_Send(LeftCol.mu0, width_x, MPI_DOUBLE_COMPLEX, left, 4, cart_comm);
        MPI_Recv(LeftCol.mu0, width_x, MPI_DOUBLE_COMPLEX, right, 4, cart_comm, &status);

        //Last row (except first and last points)
        for(int n = maxSize-width_t + 1; n<maxSize - 1; n++) {
            //U_0(n) U_1(n+0) U*_0(n+1) + U*_0(n-0) U_1(n-0) U_0(n+1-0)
            t1 = RightPB[2*n];  
		    t_1 = LeftPB[2*n];  
            const c_double& conf1 = Conf.mu0[n];
            const c_double& conf2 = Conf.mu1[t1];
            const c_double& conf4 = Conf.mu0[t_1];
            const c_double& conf5 = Conf.mu1[t_1];
            Staples.mu1[n] = conf1 * conf2 * TopRow.mu0[n-(maxSize-width_t)] +
                std::conj(conf4) * conf5 * TopRow.mu1[n-(maxSize-width_t)];
        }
        //Left column (except last point)
        for(int n = 0; n<maxSize-width_t; n+=width_t) {
            //U_0(n) U_1(n+0) U*_0(n+1) + U*_0(n-0) U_1(n-0) U_0(n+1-0)
		    x1 = RightPB[2*n+1];
            t1 = RightPB[2*n];  
            const c_double& conf1 = Conf.mu0[n];
            const c_double& conf2 = Conf.mu1[t1];
            const c_double& conf3 = Conf.mu0[x1];
            Staples.mu1[n] = conf1 * conf2 * std::conj(conf3) +
                RightCol.mu0[n/width_t] * RightCol.mu1[n/width_t];
        }
        //Right column (except for last point)
        for(int n = width_t-1; n<maxSize-width_t; n+=width_t) {
            //U_0(n) U_1(n+0) U*_0(n+1) + U*_0(n-0) U_1(n-0) U_0(n+1-0)
		    x1 = RightPB[2*n+1]; 
		    t_1 = LeftPB[2*n];  
            const c_double& conf1 = Conf.mu0[n];
            const c_double& conf3 = Conf.mu0[x1];
            const c_double& conf4 = Conf.mu0[t_1];
            const c_double& conf5 = Conf.mu1[t_1];
            const c_double& conf6 = Conf.mu0[x1_t_1[n]]; 
            Staples.mu1[n] = conf1 * LeftCol.mu0[n/width_t] * std::conj(conf3) +
                std::conj(conf4) * conf5 * conf6;
        }

        //Bottom right corner
        {
        n = maxSize - 1;
        //U_0(n) U_1(n+0) U*_0(n+1) + U*_0(n-0) U_1(n-0) U_0(n+1-0)
		t1 = RightPB[2*n];  
		t_1 = LeftPB[2*n];  
        const c_double& conf1 = Conf.mu0[n];
        const c_double& conf4 = Conf.mu0[t_1];
        const c_double& conf5 = Conf.mu1[t_1];
        Staples.mu1[n] = conf1 * LeftCol.mu0[n/width_t] * TopRow.mu0[n-(maxSize-width_t)] +
            std::conj(conf4) * conf5 * TopRow.mu1[n-(maxSize-width_t)];
        }

        //Bottom left corner
        {
        c_double topright = Conf.mu0[width_t-1]; //U_0(n+1-0)
        MPI_Send(&topright, 1, MPI_DOUBLE_COMPLEX, top_right, 6, cart_comm);
        MPI_Recv(&topright, 1, MPI_DOUBLE_COMPLEX, bot_left, 6, cart_comm, &status);
        n = maxSize - width_t;
        //U_0(n) U_1(n+0) U*_0(n+1) + U*_0(n-0) U_1(n-0) U_0(n+1-0)
		t1 = RightPB[2*n];  
        const c_double& conf1 = Conf.mu0[n];
        const c_double& conf2 = Conf.mu1[t1];
        Staples.mu1[n] = conf1 * conf2 * TopRow.mu0[n-(maxSize-width_t)] +
           RightCol.mu0[n/width_t] * topright;
        }


    } //end else

}

/*
Save Gauge configuration
*/   
void SaveConf(const GaugeConf& GConf, const std::string& Name) {
    using namespace LV;
    int counts[mpi::size];
    int displs[mpi::size];
    int offset;
    for(int i = 0; i < mpi::size; i++){
        counts[i] = 1;
        offset = (i/mpi::ranks_t);
        displs[i] = (i < mpi::ranks_t) ? i : (i-mpi::ranks_t*offset) + offset * mpi::ranks_t * mpi::width_x; 
    }

    spinor GlobalConf(LV::Ntot); //Temporary variable to store the full configuration
    MPI_Gatherv(GConf.Conf.mu0, mpi::maxSize, MPI_DOUBLE_COMPLEX,
            GlobalConf.mu0, counts, displs, sub_block_resized,
            0, mpi::cart_comm);
    MPI_Gatherv(GConf.Conf.mu1, mpi::maxSize, MPI_DOUBLE_COMPLEX,
            GlobalConf.mu1, counts, displs, sub_block_resized,
            0, mpi::cart_comm);

    if (mpi::rank == 0){
        std::ofstream Datfile(Name,std::ios::binary);
        //std::ofstream Datfile(Name);
        if (!Datfile.is_open()) {
            std::cerr << "Error opening file: " << Name << std::endl;
            return;
        }
        for (int x = 0; x < Nx; x++) {
        for (int t = 0; t < Nt; t++) {
        int n = x * Nx + t;
            for (int mu = 0; mu < 2; mu++) {
                const double& re = std::real(mu == 0 ? GlobalConf.mu0[n] : GlobalConf.mu1[n]);
                const double& im = std::imag(mu == 0 ? GlobalConf.mu0[n] : GlobalConf.mu1[n]);
                                
                Datfile.write(reinterpret_cast<const char*>(&x), sizeof(int));
                Datfile.write(reinterpret_cast<const char*>(&t), sizeof(int));
                Datfile.write(reinterpret_cast<const char*>(&mu), sizeof(int));
                Datfile.write(reinterpret_cast<const char*>(&re), sizeof(double));
                Datfile.write(reinterpret_cast<const char*>(&im), sizeof(double));
                
            }
        }
        }
        Datfile.close();
    } 
    
}



double GaugeConf::MeasureSp_HMC() {
	//Plaquettes have to be computed during the HMC update
    double local_Sp = 0.0;
    //reduction over all lattice points and spin components
    for (int n = 0; n < mpi::maxSize; n++) {
        local_Sp += std::real(Plaquette01[n]);
    }
    double Sp;
    MPI_Allreduce(&local_Sp, &Sp, 1, MPI_DOUBLE, MPI_SUM, mpi::cart_comm);
	return Sp;
}



double GaugeConf::Compute_gaugeAction(const double& beta) {
	double local_action = 0.0;
	for (int n = 0; n < mpi::maxSize; n++) {
        local_action += beta * std::real(1.0-Plaquette01[n]);
	}
    double action;
    MPI_Allreduce(&local_action, &action, 1, MPI_DOUBLE, MPI_SUM, mpi::cart_comm);
	return action;
}



void GaugeConf::read_conf(const std::string& name){
    std::ifstream infile(name);
    if (!infile) {
        std::cerr << "File " << name << " not found " << std::endl;
        exit(1);
    }
   
    spinor GlobalConf(LV::Ntot); //Temporary variable to store the full configuration
    int counts[mpi::size];
    int displs[mpi::size];
    int offset;
    for(int i = 0; i < mpi::size; i++){
        counts[i] = 1;
        offset = (i/mpi::ranks_t);
        displs[i] = (i < mpi::ranks_t) ? i : (i-mpi::ranks_t*offset) + offset * mpi::ranks_t * mpi::width_x; 
    }

    if (mpi::rank2d == 0){
        int x, t, mu;
        double re, im; 
        while (infile >> x >> t >> mu >> re >> im) {
            if (mu == 0)
                GlobalConf.mu0[x*LV::Nt+t] = c_double(re, im); //x*LV::Nt+t;
            else
                GlobalConf.mu1[x*LV::Nt+t] = c_double(re, im); //x*LV::Nt+t;
        }
        infile.close();
        std::cout << "Conf read from " << name << std::endl;    
        
        
    }

    MPI_Scatterv(GlobalConf.mu0, counts, displs, sub_block_resized,
                 Conf.mu0, mpi::maxSize, MPI_DOUBLE_COMPLEX, 0, mpi::cart_comm);

    MPI_Scatterv(GlobalConf.mu1, counts, displs, sub_block_resized,
                 Conf.mu1, mpi::maxSize, MPI_DOUBLE_COMPLEX, 0, mpi::cart_comm);


}


void GaugeConf::readBinary(const std::string& name){
    using namespace LV;
    std::ifstream infile(name, std::ios::binary);
    if (!infile) {
       std::cerr << "File " << name << " not found " << std::endl;
        exit(1);
    }
    spinor GlobalConf(LV::Ntot); //Temporary variable to store the full configuration

    int counts[mpi::size];
    int displs[mpi::size];
    int offset;
    for(int i = 0; i < mpi::size; i++){
        counts[i] = 1;
        offset = (i/mpi::ranks_t);
        displs[i] = (i < mpi::ranks_t) ? i : (i-mpi::ranks_t*offset) + offset * mpi::ranks_t * mpi::width_x; 
    }

    //In case I forget, check this presentation: https://engineering.purdue.edu/~smidkiff/KKU/files/MPI2.pdf
    if (mpi::rank2d == 0){
        for (int x = 0; x < Nx; x++) {
        for (int t = 0; t < Nt; t++) {
            int n = x * Nt + t;
            for (int mu = 0; mu < 2; mu++) {
                int x_read, t_read, mu_read;
                double re, im;
                infile.read(reinterpret_cast<char*>(&x_read), sizeof(int));
                infile.read(reinterpret_cast<char*>(&t_read), sizeof(int));
                infile.read(reinterpret_cast<char*>(&mu_read), sizeof(int));
                infile.read(reinterpret_cast<char*>(&re), sizeof(double));
                infile.read(reinterpret_cast<char*>(&im), sizeof(double));
                if (mu == 0)
                    GlobalConf.mu0[n] = c_double(re, im); //n --> for testing
                else
                    GlobalConf.mu1[n] = c_double(re, im);
            }
        }
        }
        infile.close();
        //std::cout << "Binary conf read from " << name << std::endl;     
    }

    MPI_Scatterv(GlobalConf.mu0, counts, displs, sub_block_resized,
                 Conf.mu0, mpi::maxSize, MPI_DOUBLE_COMPLEX, 0, mpi::cart_comm);

    MPI_Scatterv(GlobalConf.mu1, counts, displs, sub_block_resized,
                 Conf.mu1, mpi::maxSize, MPI_DOUBLE_COMPLEX, 0, mpi::cart_comm);




}
