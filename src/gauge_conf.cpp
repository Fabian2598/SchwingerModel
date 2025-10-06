#include "gauge_conf.h"

/*
 *                t                  Stripes parallelization
 *   0  +-------------------+  Nt   +-------------------+
 *      |                   |       |       rank 0      |
 *      |                   |       |-------------------|
 *      |                   |       |       rank 1      |
 *   x  |                   |       |-------------------|
 *      |                   |       |       rank 2      |
 *      |                   |       |-------------------|
 *      |                   |       |       rank 3      |
 *   Nx +-------------------+ Nt    +-------------------+
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

        for(int n = 0; n<LV::Nt; n++){
            TopRow.mu0[n] = Conf.mu0[n];
        }
        MPI_Send(TopRow.mu0, LV::Nt, MPI_DOUBLE_COMPLEX, mod(mpi::rank-1,mpi::size), 0, MPI_COMM_WORLD);
	    MPI_Recv(TopRow.mu0, LV::Nt, MPI_DOUBLE_COMPLEX, mod(mpi::rank+1,mpi::size), 0, MPI_COMM_WORLD, &status);
        
        for (int n = 0; n<mpi::maxSize-LV::Nt; n++){
            Plaquette01[n] = Conf.mu0[n] * Conf.mu1[RightPB[2*n]] * std::conj(Conf.mu0[RightPB[2*n+1]]) * std::conj(Conf.mu1[n]);
        }
        for(int n = mpi::maxSize-LV::Nt; n<mpi::maxSize; n++){
            Plaquette01[n] = Conf.mu0[n] * Conf.mu1[RightPB[2*n]] * std::conj(TopRow.mu0[n-(mpi::maxSize-LV::Nt)]) * std::conj(Conf.mu1[n]);
        }		
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
        for (int n = Nt; n<maxSize-Nt; n++) {
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

        //Sending top row to previous rank
        for( int n = 0; n < Nt; n++){
            TopRow.mu0[n] = Conf.mu0[n];
        }
        MPI_Send(TopRow.mu0, Nt, MPI_DOUBLE_COMPLEX, mod(rank-1,size), 0, MPI_COMM_WORLD);
        MPI_Recv(TopRow.mu0, Nt, MPI_DOUBLE_COMPLEX, mod(rank+1,size), 0, MPI_COMM_WORLD, &status);

        //Last row
        for(int n = maxSize-Nt; n < maxSize; n++){
            //U_1(n) U_0(n+1) U*_1(n+0) + U*_1(n-1) U_0(n-1) U_1(n-1+0)
		    x_1 = LeftPB[2*n+1]; 
		    t1 = RightPB[2*n];  
		    t_1 = LeftPB[2*n];  
            
            const c_double& conf1 = Conf.mu1[n]; 
            const c_double& conf3 = Conf.mu1[t1];
            const c_double& conf4 = Conf.mu1[x_1];
            const c_double& conf5 = Conf.mu0[x_1];
            const c_double& conf6 = Conf.mu1[x_1_t1[n]]; 
            Staples.mu0[n] = conf1 * TopRow.mu0[n-(maxSize-Nt)] * std::conj(conf3) +
                std::conj(conf4) * conf5 * conf6;
        }

        //Sending bottom row to next rank
        for(int n = maxSize-Nt; n < maxSize; n++){
            BottomRow.mu0[n - (maxSize-Nt)] = std::conj(Conf.mu1[n]) * Conf.mu0[n]; //U*_1(n-1) U_0(n-1)
            BottomRow.mu1[n - (maxSize-Nt)] = Conf.mu1[RightPB[2*n]]; //U_1(n-1+0)
        }
        MPI_Send(BottomRow.mu0, Nt, MPI_DOUBLE_COMPLEX, mod(rank+1,size), 2, MPI_COMM_WORLD);
        MPI_Send(BottomRow.mu1, Nt, MPI_DOUBLE_COMPLEX, mod(rank+1,size), 3, MPI_COMM_WORLD);

        MPI_Recv(BottomRow.mu0, Nt, MPI_DOUBLE_COMPLEX, mod(rank-1,size), 2, MPI_COMM_WORLD, &status);
        MPI_Recv(BottomRow.mu1, Nt, MPI_DOUBLE_COMPLEX, mod(rank-1,size), 3, MPI_COMM_WORLD, &status);
        for (int n = 0; n<Nt; n++) {
            //U_1(n) U_0(n+1) U*_1(n+0) + U*_1(n-1) U_0(n-1) U_1(n-1+0)
		    x1 = RightPB[2*n+1];   
		    t1 = RightPB[2*n];  
            const c_double& conf1 = Conf.mu1[n]; 
            const c_double& conf2 = Conf.mu0[x1]; 
            const c_double& conf3 = Conf.mu1[t1];
            Staples.mu0[n] = conf1 * conf2 * std::conj(conf3) +
                BottomRow.mu0[n] * BottomRow.mu1[n];
        }
        //------------------------------------------------------------//

        //-----------------mu = 1 direction--------------------//
        //Interior points
        for(int n = 0; n<maxSize-Nt; n++) {
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

        //Last row
        for(int n = 0; n < Nt; n++){
            TopRow.mu0[n] = std::conj(Conf.mu0[n]); //U*_0(n+1)
            TopRow.mu1[n] = Conf.mu0[LeftPB[2*n]]; //U_0(n+1-0)
        }
        MPI_Send(TopRow.mu0, Nt, MPI_DOUBLE_COMPLEX, mod(rank-1,size), 0, MPI_COMM_WORLD);
        MPI_Recv(TopRow.mu0, Nt, MPI_DOUBLE_COMPLEX, mod(rank+1,size), 0, MPI_COMM_WORLD, &status);

        MPI_Send(TopRow.mu1, Nt, MPI_DOUBLE_COMPLEX, mod(rank-1,size), 1, MPI_COMM_WORLD);
        MPI_Recv(TopRow.mu1, Nt, MPI_DOUBLE_COMPLEX, mod(rank+1,size), 1, MPI_COMM_WORLD, &status);
        
        for(int n = maxSize-Nt; n<maxSize; n++) {
            //U_0(n) U_1(n+0) U*_0(n+1) + U*_0(n-0) U_1(n-0) U_0(n+1-0)
		    t1 = RightPB[2*n];  
		    t_1 = LeftPB[2*n];  
            const c_double& conf1 = Conf.mu0[n];
            const c_double& conf2 = Conf.mu1[t1];
            const c_double& conf4 = Conf.mu0[t_1];
            const c_double& conf5 = Conf.mu1[t_1];
            Staples.mu1[n] = conf1 * conf2 * TopRow.mu0[n-(maxSize-Nt)] +
                std::conj(conf4) * conf5 * TopRow.mu1[n-(maxSize-Nt)];
        }


    } //end else

}

/*
Save Gauge configuration
*/   

void SaveConf(const GaugeConf& GConf, const std::string& Name) {
    using namespace LV;

    spinor GlobalConf(LV::Ntot); //Temporary variable to store the full configuration
    int counts[mpi::size], displs[mpi::size];
    for(int i = 0; i < mpi::size; i++) {
        counts[i] = (i != mpi::size-1) ?  (LV::Nx/mpi::size) * LV::Nt :  (LV::Nx/mpi::size) * LV::Nt + (LV::Nx%mpi::size)*LV::Nt;
        displs[i] = i * (LV::Nx/mpi::size) * LV::Nt;
    }
    MPI_Gatherv(GConf.Conf.mu0, mpi::maxSize, MPI_DOUBLE_COMPLEX,
            GlobalConf.mu0, counts, displs, MPI_DOUBLE_COMPLEX,
            0, MPI_COMM_WORLD);
    MPI_Gatherv(GConf.Conf.mu1, mpi::maxSize, MPI_DOUBLE_COMPLEX,
            GlobalConf.mu1, counts, displs, MPI_DOUBLE_COMPLEX,
            0, MPI_COMM_WORLD);
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
                
                /*
                Datfile << x
                        << std::setw(10) << t
                        << std::setw(10) << mu
                        << std::setw(30) << std::setprecision(17) << std::scientific << re
                        << std::setw(30) << std::setprecision(17) << std::scientific << im
                        << "\n";
                */
                
                
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
    #pragma omp parallel for reduction(+:local_Sp)
    for (int n = 0; n < mpi::maxSize; n++) {
        local_Sp += std::real(Plaquette01[n]);
    }
    double Sp;
    MPI_Allreduce(&local_Sp, &Sp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return Sp;
}



double GaugeConf::Compute_gaugeAction(const double& beta) {
	double local_action = 0.0;
    #pragma omp parallel for reduction(+:local_action)
	for (int n = 0; n < mpi::maxSize; n++) {
        local_action += beta * std::real(1.0-Plaquette01[n]);
	}
    double action;
    MPI_Allreduce(&local_action, &action, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return action;
}



void GaugeConf::read_conf(const std::string& name){
    std::ifstream infile(name);
    if (!infile) {
        std::cerr << "File " << name << " not found " << std::endl;
        exit(1);
    }
    spinor GlobalConf(LV::Ntot); //Temporary variable to store the full configuration
    int counts[mpi::size], displs[mpi::size];
    for(int i = 0; i < mpi::size; i++) {
        counts[i] = (i != mpi::size-1) ?  (LV::Nx/mpi::size) * LV::Nt :  (LV::Nx/mpi::size) * LV::Nt + (LV::Nx%mpi::size)*LV::Nt;
        displs[i] = i * (LV::Nx/mpi::size) * LV::Nt;
    }
    if (mpi::rank == 0){
        int x, t, mu;
        double re, im; 
        while (infile >> x >> t >> mu >> re >> im) {
            if (mu == 0)
                GlobalConf.mu0[x*LV::Nt+t] = c_double(re, im); 
            else
                GlobalConf.mu1[x*LV::Nt+t] = c_double(re, im); 
        }
        infile.close();
        std::cout << "Conf read from " << name << std::endl;     
    }

    MPI_Scatterv(GlobalConf.mu0, counts, displs, MPI_DOUBLE_COMPLEX,
                 Conf.mu0, mpi::maxSize, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Scatterv(GlobalConf.mu1, counts, displs, MPI_DOUBLE_COMPLEX,
                 Conf.mu1, mpi::maxSize, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

}


void GaugeConf::readBinary(const std::string& name){
    using namespace LV;
    std::ifstream infile(name, std::ios::binary);
    if (!infile) {
        std::cerr << "File " << name << " not found " << std::endl;
        exit(1);
    }
    spinor GlobalConf(LV::Ntot); //Temporary variable to store the full configuration
    int counts[mpi::size], displs[mpi::size];
    for(int i = 0; i < mpi::size; i++) {
        counts[i] = (i != mpi::size-1) ?  (LV::Nx/mpi::size) * LV::Nt :  (LV::Nx/mpi::size) * LV::Nt + (LV::Nx%mpi::size)*LV::Nt;
        displs[i] = i * (LV::Nx/mpi::size) * LV::Nt;
    }

    if (mpi::rank == 0){
        for (int x = 0; x < Nx; x++) {
        for (int t = 0; t < Nt; t++) {
            int n = x * Nx + t;
            for (int mu = 0; mu < 2; mu++) {
                int x_read, t_read, mu_read;
                double re, im;
                infile.read(reinterpret_cast<char*>(&x_read), sizeof(int));
                infile.read(reinterpret_cast<char*>(&t_read), sizeof(int));
                infile.read(reinterpret_cast<char*>(&mu_read), sizeof(int));
                infile.read(reinterpret_cast<char*>(&re), sizeof(double));
                infile.read(reinterpret_cast<char*>(&im), sizeof(double));
                if (mu_read == 0)
                    GlobalConf.mu0[n] = c_double(re, im);
                else
                    GlobalConf.mu1[n] = c_double(re, im);
            }
        }
        }
        infile.close();
        //std::cout << "Binary conf read from " << name << std::endl;     
    }

    MPI_Scatterv(GlobalConf.mu0, counts, displs, MPI_DOUBLE_COMPLEX,
                 Conf.mu0, mpi::maxSize, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Scatterv(GlobalConf.mu1, counts, displs, MPI_DOUBLE_COMPLEX,
                 Conf.mu1, mpi::maxSize, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
}
