#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include "conjugate_gradient.h"
#include <format>

//Formats decimal numbers
//For opening file with confs 
static std::string format(const double& number) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << number;
    std::string str = oss.str();
    str.erase(str.find('.'), 1); //Removes decimal dot 
    return str;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi::size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi::rank);

    mpi::maxSize = (mpi::rank != mpi::size-1) ? LV::Nx/mpi::size * LV::Nt :  LV::Nx/mpi::size * LV::Nt + (LV::Nx%mpi::size)*LV::Nt;
    if (mpi::size == 1) mpi::maxSize = LV::Ntot;
    std::cout << "MPI size: " << mpi::size << " Rank: " << mpi::rank << " maxSize: " << mpi::maxSize << std::endl;
    allocate_lattice_arrays();

    

    srand(mpi::rank*19);

    double beta = 2; 
    double m0 = -0.18840579710144945;
    int nconf=0;
    if (LV::Nx == 64)
        nconf = 0;
    else if (LV::Nx == 128)
        nconf = 3;
    else if (LV::Nx == 256)
        nconf = 20;  

    CG::max_iter = 10000; //Maximum number of iterations for the conjugate gradient method
    CG::tol = 1e-10; //Tolerance for convergence
    
    periodic_boundary(); //Compute right and left periodic boundary
    
	GaugeConf GConf = GaugeConf();  //Gauge configuration
    GConf.initialization(); //Random initialization of the gauge configuration 

    std::ostringstream NameData;
    NameData << "../confs/b" << beta << "_" << LV::Nx << "x" << LV::Nt << "/m-018/2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt << "_b" << 
    format(beta).c_str() << "_m" << format(m0).c_str() << "_" << nconf << ".ctxt";
    GConf.read_conf(NameData.str());
 
    spinor sol(mpi::maxSize), rhs(mpi::maxSize);

    for(int n = 0; n <mpi::maxSize; n++) {
        rhs.mu0[n] = 1;//RandomU1(); //spin up
        rhs.mu1[n] = 1;//RandomU1(); //spin down
    }
  
   D_phi(GConf.Conf, rhs, sol, m0); //Applies Dirac operator to rhs and stores the result in sol

    
    
    for(int i = 0; i < mpi::size; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (i == mpi::rank) {
            //printf("Rank %d\n", mpi::rank);
            for(int n = 0; n < mpi::maxSize; n++) {
                //std::cout << "Rank " << mpi::rank << " sol.mu0[" << n << "] = " << sol.mu0[n] << std::endl;
                //std::cout << "Rank " << mpi::rank << " sol.mu1[" << n << "] = " << sol.mu1[n] << std::endl;
                std::cout << sol.mu0[n] << std::endl;
                std::cout << sol.mu1[n] << std::endl;
            }
        }
    }
    

    /*
    for(int i = 0; i < mpi::size; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (i == mpi::rank) {
            //printf("Rank %d\n", mpi::rank);
            for(int n = 0; n < mpi::maxSize; n++) {
                //std::cout << "Rank " << mpi::rank << " U.mu0[" << n << "] = " << GConf.Conf.mu0[n] << std::endl;
                //std::cout << "Rank " << mpi::rank << " U.mu1[" << n << "] = " << GConf.Conf.mu1[n] << std::endl;
                std::cout << GConf.Conf.mu0[n] << std::endl;
                std::cout << GConf.Conf.mu1[n] << std::endl;
            }
        }
    }
    */
    
        
    
        
        
        
    


    //Free coordinate arrays
    free_lattice_arrays();

     MPI_Finalize();

	return 0;
}

