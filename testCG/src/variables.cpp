#include "variables.h"

double pi=3.14159265359;

/*
	Vectorized lattice coords.*/
int Coords(const int& x, const int& t){
	return x*LV::Nt + t;
}

namespace mpi{
    int rank = 0;
    int size = 1; 
    int maxSize = LV::Ntot; //Default value, will be updated in main
}

namespace CG{
	int max_iter = 1000; //Maximum number of iterations for the conjugate gradient method
	double tol = 1e-10; //Tolerance for convergence
}

int* LeftPB = nullptr;
int* RightPB = nullptr;
c_double* SignL = nullptr;
c_double* SignR = nullptr;
int* x_1_t1 = nullptr;
int* x1_t_1 = nullptr;

void allocate_lattice_arrays() {
    using namespace mpi;
    LeftPB  = new int[maxSize * 2];
    RightPB = new int[maxSize * 2];
    SignL   = new c_double[maxSize * 2];
    SignR   = new c_double[maxSize * 2];
    x_1_t1  = new int[maxSize];
    x1_t_1  = new int[maxSize];
}

void free_lattice_arrays() {
    delete[] LeftPB;
    delete[] RightPB;
    delete[] SignL;
    delete[] SignR;
    delete[] x_1_t1;
    delete[] x1_t_1;
}

//Memory preallocation I need to give them the right dimension mpi::maxSize, but I only know it in main after MPI_Init ...
spinor DTEMP;
spinor TEMP; 


//Buffers for MPI communication
spinor TopRow(LV::Nt);
spinor BottomRow(LV::Nt);
spinor TopRowSend(LV::Nt);
spinor BottomRowSend(LV::Nt);

void save_rhs(spinor& rhs,const std::string& name){
    using namespace LV;
    spinor GlobalRhs(LV::Ntot); //Temporary variable to store the full configuration
    int counts[mpi::size], displs[mpi::size];
    for(int i = 0; i < mpi::size; i++) {
        counts[i] = (i != mpi::size-1) ?  (LV::Nx/mpi::size) * LV::Nt :  (LV::Nx/mpi::size) * LV::Nt + (LV::Nx%mpi::size)*LV::Nt;
        displs[i] = i * (LV::Nx/mpi::size) * LV::Nt;
    }

    MPI_Gatherv(rhs.mu0, mpi::maxSize, MPI_DOUBLE_COMPLEX,
            GlobalRhs.mu0, counts, displs, MPI_DOUBLE_COMPLEX,
            0, MPI_COMM_WORLD);
    MPI_Gatherv(rhs.mu1, mpi::maxSize, MPI_DOUBLE_COMPLEX,
            GlobalRhs.mu1, counts, displs, MPI_DOUBLE_COMPLEX,
            0, MPI_COMM_WORLD);

    if (mpi::rank == 0){
        std::ofstream Datfile(name,std::ios::binary);
        //std::ofstream Datfile(Name);
        if (!Datfile.is_open()) {
            std::cerr << "Error opening file: " << name << std::endl;
            return;
        }
        for (int x = 0; x < Nx; x++) {
        for (int t = 0; t < Nt; t++) {
        int n = x * Nx + t;
        for (int mu = 0; mu < 2; mu++) {
            const double& re = std::real(mu == 0 ? GlobalRhs.mu0[n] : GlobalRhs.mu1[n]);
            const double& im = std::imag(mu == 0 ? GlobalRhs.mu0[n] : GlobalRhs.mu1[n]);
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

void read_rhs(spinor& rhs, const std::string& name){
    using namespace LV;
    std::ifstream infile(name, std::ios::binary);
    if (!infile) {
        std::cerr << "File " << name << " not found " << std::endl;
        exit(1);
    }
    spinor Globalrhs(LV::Ntot); //Temporary variable to store the full configuration
    int counts[mpi::size], displs[mpi::size];
    for(int i = 0; i < mpi::size; i++) {
        counts[i] = (i != mpi::size-1) ?  (LV::Nx/mpi::size) * LV::Nt :  (LV::Nx/mpi::size) * LV::Nt + (LV::Nx%mpi::size)*LV::Nt;
        displs[i] = i * (LV::Nx/mpi::size) * LV::Nt;
    }

    if (mpi::rank == 0){
        int n;
        for (int x = 0; x < Nx; x++) {
        for (int t = 0; t < Nt; t++) { 
            for (int mu = 0; mu < 2; mu++) {
                int x_read, t_read, mu_read;
                double re, im;
                infile.read(reinterpret_cast<char*>(&x_read), sizeof(int));
                infile.read(reinterpret_cast<char*>(&t_read), sizeof(int));
                infile.read(reinterpret_cast<char*>(&mu_read), sizeof(int));
                infile.read(reinterpret_cast<char*>(&re), sizeof(double));
                infile.read(reinterpret_cast<char*>(&im), sizeof(double));
                n = x_read * Nt + t_read;
                if (mu_read == 0)
                    Globalrhs.mu0[n] = c_double(re, im);
                else
                    Globalrhs.mu1[n] = c_double(re, im);
            }
        }
        }
        infile.close();
        std::cout << "Rhs read from " << name << std::endl;     
    }

    MPI_Scatterv(Globalrhs.mu0, counts, displs, MPI_DOUBLE_COMPLEX,
                 rhs.mu0, mpi::maxSize, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Scatterv(Globalrhs.mu1, counts, displs, MPI_DOUBLE_COMPLEX,
                 rhs.mu1, mpi::maxSize, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
}