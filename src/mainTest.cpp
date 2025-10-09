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

void RankErrorMessage(){
    if (mpi::ranks_t * mpi::ranks_x != mpi::size){
        std::cout << "ranks_t * ranks_x != total number of ranks" << std::endl;
        exit(1);
    }
    if (LV::Nx % mpi::ranks_x!= 0 ||LV::Nt % mpi::ranks_t != 0){
        std::cout << "Nx (Nt) is not exactly divisible by rank_x (rank_t)" << std::endl;
        exit(1);
    }
    mpi::width_x = LV::Nx/mpi::ranks_x;
    mpi::width_t = LV::Nt/mpi::ranks_t;
    mpi::maxSize = mpi::width_t * mpi::width_x;
    if (mpi::rank == 0){
        std::cout << "width_x " << mpi::width_x << std::endl;
        std::cout << "width_t " << mpi::width_t << std::endl;
    }
}
 
void buildCartesianTopology(){
    //if (mpi::rank == 0){
    //    std::cout << "ranks_x " << mpi::ranks_x << "  ranks_t  " << mpi::ranks_t << std::endl;
    //}
    
    int periods[2] = {true, true}; // Make both dimensions periodic
    int reorder = true; // Let MPI assign arbitrary ranks if it deems it necessary
    int dims[2] = {mpi::ranks_x,mpi::ranks_t};
    // Create a communicator given the 2D torus topology.
    MPI_Comm new_communicator;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &new_communicator);
 
    // My rank in the new communicator
    MPI_Comm_rank(new_communicator, &mpi::rank2d);
 
    // Get my coordinates in the new communicatormpi::
    MPI_Cart_coords(new_communicator, mpi::rank2d, 2, mpi::coords);

    //Does this coincide with the ranks from the new communicator?
    mpi::top = mod(coords[0]-1,mpi::ranks_x) * mpi::ranks_t + coords[1]; //rank above
	mpi::bot = mod(coords[0]+1,mpi::ranks_x) * mpi::ranks_t + coords[1]; //rank below
	mpi::right = coords[0] * mpi::ranks_t + mod(coords[1]+1,mpi::ranks_t); //rank to the right
	mpi::left = coords[0] * mpi::ranks_t + mod(coords[1]-1,mpi::ranks_t); //rank to the left
 
    // Print my location in the 2D torus.
    //printf("[Original Comm %d][MPI process %d] I am located at (%d, %d).\n",mpi::rank ,mpi::rank2d, mpi::coords[0],mpi::coords[1]);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi::size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi::rank);

    mpi::ranks_x = 2;
    mpi::ranks_t = 4;
    RankErrorMessage();

    allocate_lattice_arrays();
    srand(mpi::rank*1);

    periodic_boundary(); //Compute right and left periodic boundary
    
    using namespace mpi;
    /*
    if (rank == 0){
    std::cout << "BottomRow" << std::endl;
    for(int n = maxSize-width_t; n < maxSize; n++){
		std::cout << n << " ";
	}
    std::cout << std::endl;
    for(int n = maxSize-width_t; n < maxSize; n++){
		std::cout << n - (maxSize-width_t) << " ";
	}
    std::cout << "----------------------------------" << std::endl;
    std::cout << std::endl;

    std::cout << "TopRow" << std::endl;
	for(int n = 0; n < width_t; n++){
		std::cout << n << " ";
	}
    std::cout << std::endl;
    for(int n = 0; n < width_t; n++){
		std::cout << n << " ";
	}
    std::cout << "----------------------------------" << std::endl;
    std::cout << std::endl;

    std::cout << "RightCol" << std::endl;
	for(int n = width_t - 1; n<maxSize; n+=width_t){
		std::cout << n << " ";
	}
    std::cout << std::endl;
    for(int n = width_t - 1; n<maxSize; n+=width_t){
		std::cout << n/width_t << " ";
	}
    std::cout << "----------------------------------" << std::endl;
    std::cout << std::endl;


    std::cout << "LeftCol" << std::endl;
	for(int n = 0; n<maxSize; n+=width_t ){
		std::cout << n << " ";
	}
    std::cout << std::endl;
    for(int n = 0; n<maxSize; n+=width_t ){
		std::cout << n/width_t << " ";
	}
    std::cout << std::endl;


    }
    */
    buildCartesianTopology();
   



    /*
	GaugeConf GConf = GaugeConf();  //Gauge configuration
    //GConf.initialization(); //Random initialization of the gauge configuration 

    std::ostringstream NameData;
    NameData << "../confs/b" << beta << "_" << LV::Nx << "x" << LV::Nt << "/m-018/2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt << "_b" << 
    //NameData << "/wsgjsc/home/nietocastellanos1/Downloads/" << "2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt << "_b" << 
    format(beta).c_str() << "_m" << format(m0).c_str() << "_" << nconf << ".ctxt";
    GConf.read_conf(NameData.str());
 
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

