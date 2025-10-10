#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include <format>
#include <sstream>
#include "conjugate_gradient.h"


//Formats decimal numbers
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
    
    allocate_lattice_arrays();
    //srand(mpi::rank*time(0));
    srand(mpi::rank);

    double beta, m0; 
    int nconf;

    CG::max_iter = 10000; //Maximum number of iterations for the conjugate gradient method
    CG::tol = 1e-10; //Tolerance for convergence

    if (mpi::rank == 0){
         //---Input data---//
        std::cout << "  -----------------------------" << std::endl;
        std::cout << "  |          Testing CG       |" << std::endl;
        std::cout << "  -----------------------------" << std::endl;
        std::cout << "Nx " << LV::Nx << " Nt " << LV::Nt << std::endl;
        std::cout << "beta : ";
        std::cin >> beta;
        std::cout << "m0: ";
        std::cin >> m0;
        std::cout << "Configuration id";
        std::cin >> nconf;
        std::cout << " " << std::endl;
        
    }
   
    MPI_Bcast(&beta, 1, MPI_DOUBLE,  0, MPI_COMM_WORLD);
    MPI_Bcast(&m0, 1, MPI_DOUBLE,  0, MPI_COMM_WORLD);
    MPI_Bcast(&nconf, 1, MPI_INT,  0, MPI_COMM_WORLD);
    
    periodic_boundary(); //Compute right and left periodic boundary
    
	GaugeConf GConf = GaugeConf();  //Gauge configuration
    
    std::ostringstream file;
    file << "../../confs/b" << beta << "_" <<  LV::Nx << "x" << LV::Nt << "/m-018/2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt
        << "_b" << format(beta)
        << "_m" << format(m0)
        << "_" << nconf << ".ctxt";

    //GConf.read_conf(file.str());
    GConf.readBinary(file.str());

    if (mpi::rank == 0){
        std::cout << "**********************************************************************" << std::endl;
        std::cout << "*                              PARAMETERS" << std::endl;
        std::cout << "* Nx = " << LV::Nx << ", Nt = " << LV::Nt << std::endl;
        std::cout << "* m0 = " << m0 << ", kappa = " << 1/(2*(m0+2)) << std::endl;
        std::cout << "* beta = " << beta << std::endl;
        std::cout << "* Conf ID = " << nconf << std::endl;
        std::cout << "* CG max iterations = " << CG::max_iter << ", CG tolerance = " << CG::tol << std::endl;
        std::cout << "* Number of MPI ranks = " << mpi::size << std::endl;
        std::cout << "**********************************************************************" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    spinor rhs(mpi::maxSize), solCG(mpi::maxSize), solBi(mpi::maxSize);
    for(int n = 0; n<mpi::maxSize; n++){
        rhs.mu0[n] = RandomU1();
        rhs.mu1[n] = RandomU1();
    }
    

    double begin = MPI_Wtime();
    conjugate_gradient(GConf.Conf, rhs, solCG, m0);
    double end = MPI_Wtime();
    double elapsed_secs = end - begin;
    if (mpi::rank == 0)
    std::cout << "Cg time = " << elapsed_secs << " s" << std::endl;
        
    MPI_Barrier(MPI_COMM_WORLD);

    begin = MPI_Wtime();
    bi_cgstab(GConf.Conf, rhs, solBi, m0, 10*CG::max_iter, CG::tol);
    end = MPI_Wtime();
    elapsed_secs = end - begin;
    if (mpi::rank == 0)
    std::cout << "BiCG time = " << elapsed_secs << " s" << std::endl;
    
    //Free coordinate arrays
    free_lattice_arrays();
    
    MPI_Finalize();

	return 0;
}

