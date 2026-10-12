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
 * The loops over the lattice are executed by for_each_site / lattice_sum (GPU kernels when compiled with CUDA).
 * The lambdas only capture raw pointers (field.val) and plain values, check backend.h.
 */

void GaugeConf::initialization() {
    //Random numbers are generated on the CPU and then copied to the GPU
    std::vector<c_double> U = Conf.toHost();
    int n;
	for (int x = 1; x <= mpi::width_x; x++) {
        for (int t = 1; t <= mpi::width_t; t++) {
            n = x*(mpi::width_t+2)+t;
		    U[2*n]   = RandomU1();
            U[2*n+1] = RandomU1();
	    }
    }
    Conf.fromHost(U);
}


//U_mv(x) = U_m(x) U_v(x+m) U*_m(x+v) U*_v(x)
//mu = 0 time direction, mu = 1 space direction
void GaugeConf::Compute_Plaquette01() {
	//U_01(n) = U_0(n) U_1(n+0) U*_0(n+1) U*_1(n)
    //Halo must be communicated externally
    const LocalLattice lat = local_lattice();
    const c_double* U = Conf.val;
    c_double* plaquette = Plaquette01.val;
    for_each_site([=] HOST_DEVICE (int n){
        int right, down, left, up;
        double rsign, lsign;
        //get coordinates of the neighbors and boundary sign
        get_neighbors(lat, n, right, down, left, up, rsign, lsign);
        plaquette[n] = U[2*n] * U[2*right+1] * conj(U[2*down]) * conj(U[2*n+1]);
    });
}

//mu = 0 time direction, mu = 1 space direction
void GaugeConf::Compute_Q() {
    //Halo, corners included, must be communicated externally (exchange_halo_gauge)
    const LocalLattice lat = local_lattice();
    const c_double* U = Conf.val;
    c_double* q01 = Q01.val;
    c_double* q10 = Q10.val;
    c_double* p1 = P1.val;
    c_double* p2 = P2.val;
    c_double* p3 = P3.val;
    c_double* p4 = P4.val;
    for_each_site([=] HOST_DEVICE (int n){
        int right, down, left, up;
        double lsign, rsign;
        int x1_t_1, x_1_t_1, x_1_t1, x1_t1; //n-0+1, n-0-1, n+0-1, n+0+1
        c_double Umv, Uv_m, U_m_v, U_vm;  //U_{m,v}(n) + U_{v,-m}(n) + U_{-m,-v}(n) + U_{-v,m}(n)
        get_neighbors(lat, n, right, down, left, up, rsign, lsign);
        get_corners(lat, n, x1_t_1, x_1_t_1, x_1_t1, x1_t1);

        //Q01

        //U_01(n) = U_0(n) U_1(n+0) U*_0(n+1) U*_1(n)
        Umv = U[2*n] * U[2*right+1] * conj(U[2*down]) * conj(U[2*n+1]);
        p1[n] = Umv;

        //U_{1-0}(n) = U_1(n) U*_0(n-0+1) U*_1(n-0) U_0(n-0)
        Uv_m = U[2*n+1] * conj(U[2*x1_t_1]) * conj(U[2*left+1]) * U[2*left];
        p2[n] = Uv_m;

        //U_{-0,-1}(n) = U*_0(n-0) U*_1(n-0-1) U_0(n-0-1) U_1(n-1)
        U_m_v = conj(U[2*left]) * conj(U[2*x_1_t_1+1]) * U[2*x_1_t_1] * U[2*up+1];
        p3[n] = U_m_v;

        //U_{-10}(n) = U*_1(n-1) U_0(n-1) U_1(n+0-1) U*_0(n)
        U_vm = conj(U[2*up+1]) * U[2*up] * U[2*x_1_t1+1] * conj(U[2*n]);
        p4[n] = U_vm;

        q01[n] = Umv+Uv_m+U_m_v+U_vm;

        //This is only for testing, I should get rid of it eventually
        //Q10
        //U_10(n) = U_1(n) U_0(n+1) U*_1(n+0) U*_0(n)
        Umv = U[2*n+1] * U[2*down] * conj(U[2*right+1]) * conj(U[2*n]);

        //U_{0-1}(n) = U_0(n) U*_1(n-1+0) U*_0(n-1) U_1(n-1)
        Uv_m = U[2*n] * conj(U[2*x_1_t1+1]) * conj(U[2*up]) * U[2*up+1];

        //U_{-1,-0}(n) = U*_1(n-1) U*_0(n-0-1) U_1(n-0-1) U_0(n-0)
        U_m_v = conj(U[2*up+1]) * conj(U[2*x_1_t_1]) * U[2*x_1_t_1+1] * U[2*left];

        //U_{-01}(n) = U*_0(n-0) U_1(n-0) U_0(n+1-0) U*_1(n)
        U_vm = conj(U[2*left]) * U[2*left+1] * U[2*x1_t_1] * conj(U[2*n+1]);

        q10[n] = Umv+Uv_m+U_m_v+U_vm;
    });
}


//Compute staple at coordinate (x,t) in the mu-direction
void GaugeConf::Compute_Staple() {
    //U_v(x) U_m(x+v) U*_v(x+m) + U*_v(x-v) U_m(x-v) U_v(x+m-v)
    //mu = 0 time direction, mu = 1 space direction
    //Halo, corners included, must be communicated externally (exchange_halo_gauge)
    const LocalLattice lat = local_lattice();
    const c_double* U = Conf.val;
    c_double* staple = Staples.val;
    for_each_site([=] HOST_DEVICE (int n){
        int x1, x_1, t1, t_1; //Nearest neighbors
        int x1_t_1, x_1_t_1, x_1_t1, x1_t1; //Diagonal neighbors
        double rsign, lsign;
        //These coordinates could change depending on the conventions
        get_neighbors(lat, n, t1, x1, t_1, x_1, rsign, lsign);
        get_corners(lat, n, x1_t_1, x_1_t_1, x_1_t1, x1_t1);
        {
            //mu = 0
            //U_1(n) U_0(n+1) U*_1(n+0) + U*_1(n-1) U_0(n-1) U_1(n-1+0)
            const c_double conf1 = U[2*n+1];
            const c_double conf2 = U[2*x1];
            const c_double conf3 = U[2*t1+1];
            const c_double conf4 = U[2*x_1+1];
            const c_double conf5 = U[2*x_1];
            const c_double conf6 = U[2*x_1_t1+1];
            staple[2*n] = conf1 * conf2 * conj(conf3) +
                conj(conf4) * conf5 * conf6;
        }
        {
            //mu = 1
            //U_0(n) U_1(n+0) U*_0(n+1) + U*_0(n-0) U_1(n-0) U_0(n+1-0)
            const c_double conf1 = U[2*n];
            const c_double conf2 = U[2*t1+1];
            const c_double conf3 = U[2*x1];
            const c_double conf4 = U[2*t_1];
            const c_double conf5 = U[2*t_1+1];
            const c_double conf6 = U[2*x1_t_1];
            staple[2*n+1] = conf1 * conf2 * conj(conf3) +
                conj(conf4) * conf5 * conf6;
        }
    });
}

double GaugeConf::MeasureSp_HMC() {
	//Plaquettes have to be computed during the HMC update
    const c_double* plaquette = Plaquette01.val;
    //reduction over all lattice points
    return lattice_sum<double>([=] HOST_DEVICE (int n){
        return real(plaquette[n]);
    });
}


double GaugeConf::Compute_gaugeAction(const double& beta) {
    const c_double* plaquette = Plaquette01.val;
    const double b = beta;
    return lattice_sum<double>([=] HOST_DEVICE (int n){
        return b * real(1.0-plaquette[n]);
    });
}

void GaugeConf::SaveConf(const std::string& Name){
    using namespace LV;

    std::vector<int> counts(mpi::size);
    std::vector<int> displs(mpi::size);
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
    std::vector<c_double> LocalConf = Conf.toHost(); //Copy of the links on the CPU
    std::vector<c_double> GlobalConf((Nt+2)*(Nx+2)*2);
    // gather local inner blocks into GlobalConf at root
    MPI_Gatherv(&LocalConf[input_ini_local], 1, mpi::local_conf_resized,
                GlobalConf.data(), counts.data(), displs.data(), mpi::global_conf_resized,
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
            const double re = real(GlobalConf[2*n+mu]);
            const double im = imag(GlobalConf[2*n+mu]);
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
    std::vector<c_double> GlobalConf((Nt+2)*(Nx+2)*2); //Temporary variable to store the full configuration

    std::vector<int> counts(mpi::size);
    std::vector<int> displs(mpi::size);
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
                GlobalConf[2*n+mu] = c_double(re, im);
            }
        }
        }
        infile.close();
        //std::cout << "Binary conf read from " << Name << std::endl;
    }

    int input_ini_local = 2 * (mpi::width_t + 2 + 1);
    std::vector<c_double> LocalConf = Conf.toHost();
    MPI_Scatterv(GlobalConf.data(), counts.data(), displs.data(), mpi::global_conf_resized, &LocalConf[input_ini_local],1,
        	mpi::local_conf_resized, 0, mpi::cart_comm);
    Conf.fromHost(LocalConf); //Copy the links to the GPU

}
