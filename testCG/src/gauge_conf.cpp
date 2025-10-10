#include "gauge_conf.h"

/*
 *                t                  Strips parallelization
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
        std::cout << "Binary conf read from " << name << std::endl;     
    }

    MPI_Scatterv(GlobalConf.mu0, counts, displs, MPI_DOUBLE_COMPLEX,
                 Conf.mu0, mpi::maxSize, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Scatterv(GlobalConf.mu1, counts, displs, MPI_DOUBLE_COMPLEX,
                 Conf.mu1, mpi::maxSize, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
}
