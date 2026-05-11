#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include <format>
#include <cstdint>
#include <cstring>
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
    srand(mpi::rank*time(0));
    //srand(mpi::rank);

    double beta, m0; 
    int nconf;

    CG::max_iter = 15000; //Maximum number of iterations for the conjugate gradient method
    CG::tol = 1e-10; //Tolerance for convergence
    std::string fileName;
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
        std::cout << "Configuration id: ";
        std::cin >> nconf;
        std::cout << "file path: ";
        std::cin >> fileName;
        std::cout << " " << std::endl;
    }
   
    MPI_Bcast(&beta, 1, MPI_DOUBLE,  0, MPI_COMM_WORLD);
    MPI_Bcast(&m0, 1, MPI_DOUBLE,  0, MPI_COMM_WORLD);
    MPI_Bcast(&nconf, 1, MPI_INT,  0, MPI_COMM_WORLD);
    //MPI_Bcast(&fileName, 1, MPI_INT,  0, MPI_COMM_WORLD);

    // Broadcast std::string fileName:
    int filename_len = 0;
    if (mpi::rank == 0) {
        filename_len = static_cast<int>(fileName.size()) + 1; // include null terminator
    }
    MPI_Bcast(&filename_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<char> filename_buf(filename_len);
    if (mpi::rank == 0) {
        std::memcpy(filename_buf.data(), fileName.c_str(), filename_len);
    }
    MPI_Bcast(filename_buf.data(), filename_len, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (mpi::rank != 0) {
        fileName.assign(filename_buf.data());
    }

    periodic_boundary(); //Compute right and left periodic boundary
    
	GaugeConf GConf = GaugeConf();  //Gauge configuration
    
    std::ostringstream file;


    //file << "../../confs/b" << beta << "_" <<  LV::Nx << "x" << LV::Nt << "/m-018/2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt
    //file << "../../confs/b" << beta << "_" <<  LV::Nx << "x" << LV::Nt << "/m-0709/2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt
    //file << "../../confs/b" << beta << "_" <<  LV::Nx << "x" << LV::Nt << "/m-01023/2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt
    //    << "_b" << format(beta)
    //    << "_m" << format(m0)
    //    << "_" << nconf << ".ctxt";

    //GConf.read_conf(file.str());
    //GConf.readBinary(file.str());
    if (mpi::rank == 0){
        std::cout << "**********************************************************************" << std::endl;
        std::cout << "*                              PARAMETERS" << std::endl;
        std::cout << "* Nx = " << LV::Nx << ", Nt = " << LV::Nt << std::endl;
        std::cout << "* m0 = " << m0 << ", kappa = " << 1/(2*(m0+2)) << std::endl;
        std::cout << "* beta = " << beta << std::endl;
        std::cout << "* Conf ID = " << nconf << std::endl;
        std::cout << "* CG max iterations = " << CG::max_iter << ", CG tolerance = " << CG::tol << std::endl;
        std::cout << "* Number of MPI ranks = " << mpi::size << std::endl;
        std::cout << "* File = " << fileName << std::endl;
        std::cout << "**********************************************************************" << std::endl;
    }

    GConf.readBinary(fileName);
    MPI_Barrier(MPI_COMM_WORLD);

    spinor rhs(mpi::maxSize), solCG(mpi::maxSize), solBi(mpi::maxSize);
    spinor x0(mpi::maxSize);
    //rhs.mu0[0] = 1;
    for(int n = 0; n<mpi::maxSize; n++){
        rhs.mu0[n] = RandomU1();
        rhs.mu1[n] = RandomU1();
    }
    
    std::ostringstream rhsName;
    rhsName << "../../confs/rhs/rhs_conf" << nconf << "_" << LV::Nx << "x" << LV::Nt
            << "_b" << format(beta)
            << "_m" << format(m0)
            << ".rhs";
    //save_rhs(rhs,rhsName.str());
    
    read_rhs(rhs,rhsName.str());


    double begin = MPI_Wtime();
    conjugate_gradient(GConf.Conf, rhs, solCG, m0);
    double end = MPI_Wtime();
    double elapsed_secs = end - begin;
    if (mpi::rank == 0)
    std::cout << "Cg time = " << elapsed_secs << " s" << std::endl;
        
    MPI_Barrier(MPI_COMM_WORLD);

    
    begin = MPI_Wtime();
    solBi = bi_cgstab(GConf.Conf, rhs, x0, m0, 10000, CG::tol);
    end = MPI_Wtime();
    elapsed_secs = end - begin;
    if (mpi::rank == 0)
    std::cout << "BiCG time = " << elapsed_secs << " s" << std::endl;
    
    //Free coordinate arrays
    free_lattice_arrays();
    
    MPI_Finalize();

	return 0;
}

