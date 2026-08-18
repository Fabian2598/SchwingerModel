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
 * x+1, t down
 * x-1, t up
 * x, t+1 right
 * x, t-1 left
 * n = x * Nt + t = (x,t) coordinates
 * 
 */

void GaugeConf::initialization() {
    int n;
	for (int x = 1; x <= mpi::width_x; x++) {
        for (int t = 1; t <= mpi::width_t; t++) {
            n = x*(mpi::width_t+2)+t;
		    Conf.val[2*n]   = RandomU1(); 
            Conf.val[2*n+1] = RandomU1(); 
	    }
    }
}


//U_mv(x) = U_m(x) U_v(x+m) U*_m(x+v) U*_v(x)
//mu = 0 time direction, mu = 1 space direction
void GaugeConf::Compute_Plaquette01() {
	//U_01(n) = U_0(n) U_1(n+0) U*_0(n+1) U*_1(n)
    //Halo must be communicated externally
    //exchange_halo(Conf.val);
    int n, right, down, left, up;
    double lsign, rsign;
    for(int x = 1; x<=mpi::width_x; x++){
		for(int t = 1; t<=mpi::width_t; t++){
			n = x*(mpi::width_t+2)+t;
			//get coordinates of the neighbors and boundary sign 
			get_neighbors(x, t,right, down, left, up, rsign, lsign); 
            Plaquette01[n] = Conf.val[2*n] * Conf.val[2*right+1] * std::conj(Conf.val[2*down]) * std::conj(Conf.val[2*n+1]);
        }
    }
   
}

//Compute staple at coordinate (x,t) in the mu-direction
void GaugeConf::Compute_Staple() {
    MPI_Status status;
    //U_v(x) U_m(x+v) U*_v(x+m) + U*_v(x-v) U_m(x-v) U_v(x+m-v)
    //mu = 0 time direction, mu = 1 space direction
    int n;
    int x1, x_1, t1, t_1; //Nearest neighbors
    int x_1_t1, x1_t_1;   //Diagonal neighbors
    double rsign, lsign;
    //Halo must be communicated externally
    //exchange_halo(Conf.val);
    //Corners we have to communicate manually 
    //Update top-right corner (needs bottom-left corner from diagonal rank)
    {
        int x0 = mpi::width_x, t0 = 1;
        int n0 = x0*(mpi::width_t+2)+t0;
        c_double bottom_left = Conf.val[2*n0+1]; //U_1(n-1+0)
        MPI_Send(&bottom_left, 1, MPI_DOUBLE_COMPLEX, mpi::bot_left, 0, mpi::cart_comm);
        MPI_Recv(&bottom_left, 1, MPI_DOUBLE_COMPLEX, mpi::top_right, 0, mpi::cart_comm, &status);
        n0 = mpi::width_t+1; 
        Conf.val[2*n0+1] = bottom_left;    
    }

    //Update bottom-left corner (needs top-right corner from diagonal rank)
    {
        int x0 = 1, t0 = mpi::width_t;
        int n0 = x0*(mpi::width_t+2)+t0;
        c_double top_right = Conf.val[2*n0]; //U_0(n+1-0)
        MPI_Send(&top_right, 1, MPI_DOUBLE_COMPLEX, mpi::top_right, 1, mpi::cart_comm);
        MPI_Recv(&top_right, 1, MPI_DOUBLE_COMPLEX, mpi::bot_left, 1, mpi::cart_comm, &status);

        x0 = mpi::width_x+1; t0 = 0;
        n0 = x0*(mpi::width_t+2)+t0;
        Conf.val[2*n0] = top_right;
    }

    for(int x = 1; x<=mpi::width_x; x++){
	    for(int t = 1; t<=mpi::width_t; t++){
            n = x*(mpi::width_t+2)+t;
            //These coordinates could change depending on the conventions 
            get_neighbors(x,t,t1,x1,t_1,x_1,rsign,lsign);
            x_1_t1 = (x-1)*(mpi::width_t+2)+(t+1);
            x1_t_1 = (x+1)*(mpi::width_t+2)+(t-1);
            for (int mu = 0; mu < 2; mu++) {
                if (mu == 0) {
                    //U_1(n) U_0(n+1) U*_1(n+0) + U*_1(n-1) U_0(n-1) U_1(n-1+0)
                    const c_double& conf1 = Conf.val[2*n+1]; 
                    const c_double& conf2 = Conf.val[2*x1]; 
                    const c_double& conf3 = Conf.val[2*t1+1];
                    const c_double& conf4 = Conf.val[2*x_1+1];
                    const c_double& conf5 = Conf.val[2*x_1];
                    const c_double& conf6 = Conf.val[2*x_1_t1+1]; 
                    Staples.val[2*n] = conf1 * conf2 * std::conj(conf3) +
                        std::conj(conf4) * conf5 * conf6;
                }
                else {
                    //U_0(n) U_1(n+0) U*_0(n+1) + U*_0(n-0) U_1(n-0) U_0(n+1-0)
                    const c_double& conf1 = Conf.val[2*n];
                    const c_double& conf2 = Conf.val[2*t1+1];
                    const c_double& conf3 = Conf.val[2*x1];
                    const c_double& conf4 = Conf.val[2*t_1];
                    const c_double& conf5 = Conf.val[2*t_1+1];
                    const c_double& conf6 = Conf.val[2*x1_t_1]; 
                    Staples.val[2*n+1] = conf1 * conf2 * std::conj(conf3) +
                        std::conj(conf4) * conf5 * conf6;
                }
            }   
        }
    }
    
}

double GaugeConf::MeasureSp_HMC() {
	//Plaquettes have to be computed during the HMC update
    double local_Sp = 0.0;
    int n;
    //reduction over all lattice points and spin components
    for(int x = 1; x<=mpi::width_x; x++){
	    for(int t = 1; t<=mpi::width_t; t++){
            n = x*(mpi::width_t+2)+t;
            local_Sp += std::real(Plaquette01[n]);
        }
    }
    double Sp;
    MPI_Allreduce(&local_Sp, &Sp, 1, MPI_DOUBLE, MPI_SUM, mpi::cart_comm);
	return Sp;
}


double GaugeConf::Compute_gaugeAction(const double& beta) {
	double local_action = 0.0;
    int n;
	for(int x = 1; x<=mpi::width_x; x++){
	    for(int t = 1; t<=mpi::width_t; t++){
            n = x*(mpi::width_t+2)+t;
            local_action += beta * std::real(1.0-Plaquette01[n]);
        }
	}
    double action;
    MPI_Allreduce(&local_action, &action, 1, MPI_DOUBLE, MPI_SUM, mpi::cart_comm);
	return action;
}

void GaugeConf::SaveConf(const std::string& Name){
    using namespace LV;

    int counts[mpi::size];
    int displs[mpi::size];
    for (int r = 0; r < mpi::size; ++r) {
        counts[r] = 1;
        int rx = r / mpi::ranks_t;
        int rt = r % mpi::ranks_t;
        int global_x_start = rx * mpi::width_x + 1;
        int global_t_start = rt * mpi::width_t + 1;
        displs[r] = (global_x_start * (LV::Nt + 2) + global_t_start)*2; // in complex-element units
    }

    // index of first inner element (skip halo) in local Conf.val
    int input_ini_local = 2 * (mpi::width_t + 2 + 1); // 2*(1*(width_t+2) + 1) -> 2*(width_t+3)
    spinor GlobalConf((Nt+2)*(Nx+2)*2);
    // gather local inner blocks into GlobalConf at root
    MPI_Gatherv(&Conf.val[input_ini_local], 1, mpi::local_conf_resized,
                &GlobalConf.val[0], counts, displs, mpi::global_conf_resized,
                0, mpi::cart_comm);

    if (mpi::rank2d == 0){

        std::ofstream Datfile(Name, std::ios::binary);
        if (!Datfile.is_open()) {
            std::cerr << "Error opening file: " << Name << std::endl;
            return;
        }

        for (int x = 0; x < LV::Nx; x++) {
        for (int t = 0; t < LV::Nt; t++) {
        int n = (x+1) * (LV::Nt+2) + (t+1);
        for (int mu = 0; mu < 2; mu++) {
            const double& re = std::real(GlobalConf.val[2*n+mu]);
            const double& im = std::imag(GlobalConf.val[2*n+mu]);
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



void GaugeConf::ReadConf(const std::string& Name){
    using namespace LV;
    std::ifstream infile(Name, std::ios::binary);
    if (!infile) {
       std::cerr << "File " << Name << " not found " << std::endl;
        exit(1);
    }
    spinor GlobalConf((Nt+2)*(Nx+2)*2); //Temporary variable to store the full configuration

    int counts[mpi::size];
    int displs[mpi::size];
    for(int r = 0; r < mpi::size; r++){
        counts[r] = 1;
        int rx = r / mpi::ranks_t; 
        int rt = r % mpi::ranks_t; 

        // Global starting position inside the buffer including halo (halo at index 0)
        int global_x_start = rx * mpi::width_x + 1; // +1 to skip halo
        int global_t_start = rt * mpi::width_t + 1; // +1 to skip halo
        // Displacement in complex-element units into buffer.val (including halo padding)
        displs[r] = (global_x_start * (LV::Nt + 2) + global_t_start) * 2;

        //if (mpi::rank2d == 0)
        //    std::cout << "displ " << displs[r] << std::endl;
    }

    if (mpi::rank2d == 0){
        for (int x = 1; x <= LV::Nx; x++) {
        for (int t = 1; t <= LV::Nt; t++) {
            int n = x * (LV::Nt+2) + t;
            for (int mu = 0; mu < 2; mu++) {
                int x_read, t_read, mu_read;
                double re, im;
                infile.read(reinterpret_cast<char*>(&x_read), sizeof(int));
                infile.read(reinterpret_cast<char*>(&t_read), sizeof(int));
                infile.read(reinterpret_cast<char*>(&mu_read), sizeof(int));
                infile.read(reinterpret_cast<char*>(&re), sizeof(double));
                infile.read(reinterpret_cast<char*>(&im), sizeof(double));
                GlobalConf.val[2*n+mu] = c_double(re, im); 
            }
        }
        }
        infile.close();
        //std::cout << "Binary conf read from " << Name << std::endl;     
    }

    int input_ini_local = 2 * (mpi::width_t + 2 + 1);
    MPI_Scatterv(&GlobalConf.val[0], counts, displs, mpi::global_conf_resized, &Conf.val[input_ini_local],1,
        	mpi::local_conf_resized, 0, mpi::cart_comm);

}