#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include "mpi_setup.h"
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
    srand(mpi::rank*1);
    if (mpi::rank == 0){
        std::cout << "ranks_x: ";
        std::cin >> mpi::ranks_x;
        std::cout << "ranks_t: ";
        std::cin >> mpi::ranks_t;
    }
    MPI_Bcast(&mpi::ranks_x, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(&mpi::ranks_t, 1, MPI_INT,  0, MPI_COMM_WORLD);

    initializeMPI(); //2D rank topology
    allocate_lattice_arrays(); //Allocates memory for coordinates array
    periodic_boundary(); //Computes arrays with neighbors
   

    double m0 = -0.18840579710144945;
    int nconf = 12;
    double beta = 2;
	GaugeConf GConf = GaugeConf();  //Gauge configuration
    std::ostringstream NameData;
    NameData << "../confs/b" << beta << "_" << LV::Nx << "x" << LV::Nt << "/m-018/2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt << "_b" << 
    format(beta).c_str() << "_m" << format(m0).c_str() << "_" << nconf << ".ctxt";
    GConf.read_conf(NameData.str());
    GConf.Compute_Staple();

    
    c_double local_staple = 0.0;
	for (int n = 0; n < mpi::maxSize; n++) {
        local_staple += GConf.Staples.mu1[n];
	}
    c_double staple;
    MPI_Allreduce(&local_staple, &staple, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi::cart_comm);

    if (mpi::rank2d == 0) 
        std::cout << "Staple sum " << staple << std::endl; 
     

    /*
    spinor Staple(LV::Ntot);    
    int counts[mpi::size];
    int displs[mpi::size];
    int offset;
    for(int i = 0; i < mpi::size; i++){
        counts[i] = 1;
        offset = (i/mpi::ranks_t);
        displs[i] = (i < mpi::ranks_t) ? i : (i-mpi::ranks_t*offset) + offset * mpi::ranks_t * mpi::width_x; 
    }
    
    MPI_Gatherv(GConf.Staples.mu1, mpi::maxSize, MPI_DOUBLE_COMPLEX,
            Staple.mu1, counts, displs, sub_block_resized,
            0, mpi::cart_comm);

    if (mpi::rank2d == 0) {
        std::cout << "Reconstructed Staple on rank 0" << std::endl;
        for(int n = 0; n<LV::Ntot; n++){
            std::cout << Staple.mu1[n] << "\n";
        }
        std::cout << std::endl;
    }
   */ 
    

    /*
    for(int i = 0; i < mpi::size; i++) {
        MPI_Barrier(mpi::cart_comm);
        if (i == mpi::rank2d) {
            std::cout << "rank " << mpi::rank2d << std::endl;
            for(int n = 0; n<mpi::maxSize; n++){
                std::cout << GConf.Conf.mu0[n] << " ";
                if ((n+1) % mpi::width_t == 0 )
                    std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
    
    if (mpi::rank2d == 0) {
        std::cout << "Reconstructed GlobalConf on rank 0" << std::endl;
        for(int n = 0; n<LV::Ntot; n++){
            std::cout << GlobalConf.mu0[n] << " ";
            if ((n+1) % LV::Nt == 0 )
                std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    */
    
    

 
    spinor sol(mpi::maxSize), rhs(mpi::maxSize);

    for(int n = 0; n <mpi::maxSize; n++) {
        rhs.mu0[n] = 1;//RandomU1(); //spin up
        rhs.mu1[n] = 1;//RandomU1(); //spin down
    }



    //SaveConf(GConf, "binaryConf");
    
    //GaugeConf GConfBinary = GaugeConf();
    //GConf.readBinary(NameData.str());

    
    double startT, endT;
    startT = MPI_Wtime();
    conjugate_gradient(GConf.Conf, rhs, sol,m0);
    endT = MPI_Wtime();
    printf("[rank %d] time elapsed during CG implementation: %.4fs.\n", mpi::rank, endT - startT);

    
    /*
    spinor x0(mpi::maxSize);
    startT = MPI_Wtime();
    bi_cgstab(GConf.Conf, rhs, x0, m0, CG::max_iter, CG::tol, true);
    endT = MPI_Wtime();
    printf("[rank %d] time elapsed during BiCGstab implementation: %.4fs.\n", mpi::rank, endT - startT);
    */

    //Free coordinate arrays
    free_lattice_arrays();

    MPI_Finalize();

	return 0;
}

