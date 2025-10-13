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
    //Check that number of ranks on the x and t direction match the number of total ranks called.
    if (mpi::ranks_t * mpi::ranks_x != mpi::size){
        std::cout << "ranks_t * ranks_x != total number of ranks" << std::endl;
        std::cout << mpi::ranks_t * mpi::ranks_x << " != " << mpi::size << std::endl;
        exit(1);
    }
    if (LV::Nx % mpi::ranks_x!= 0 ||LV::Nt % mpi::ranks_t != 0){
        std::cout << "Nx (Nt) is not exactly divisible by rank_x (rank_t)" << std::endl;
        exit(1);
    }
    mpi::width_x = LV::Nx/mpi::ranks_x;
    mpi::width_t = LV::Nt/mpi::ranks_t;
    mpi::maxSize = mpi::width_t * mpi::width_x;
    //if (mpi::rank == 0){
    //    std::cout << "width_x " << mpi::width_x << std::endl;
    //    std::cout << "width_t " << mpi::width_t << std::endl;
    //}
}
 
/*
 *                t                    2D parallelization
 *   0  +-------------------+  Nt   +-----------------------+
 *      |                   |       |                       |
 *      |                   |       |          top          |
 *      |                   |       |           |           |
 *   x  |                   |       |--left--rank2d--right--|
 *      |                   |       |           |           |
 *      |                   |       |          bot          |
 *      |                   |       |                       |
 *   Nx +-------------------+ Nt    +-----------------------+
 *                Nx
*/
void buildCartesianTopology(){
    //if (mpi::rank == 0){
    //    std::cout << "ranks_x " << mpi::ranks_x << "  ranks_t  " << mpi::ranks_t << std::endl;
    //}
    int dims[2] = {mpi::ranks_x, mpi::ranks_t};
    int periods[2] = {1, 1}; // periodic in both dims
    int reorder = 1;         // allow rank reordering
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &mpi::cart_comm);

    //rank in the Cartesian communicator and its coordinates
    MPI_Comm_rank(mpi::cart_comm, &mpi::rank2d);
    MPI_Cart_coords(mpi::cart_comm, mpi::rank2d, 2, mpi::coords); // mpi::coords[0]=x coord, [1]=t coord
    
    //MPI_Cart_shift(cart_comm, Direction, Displacement, - direction,  +direction);
    //Along t direction
    MPI_Cart_shift(mpi::cart_comm, 1, 1, &mpi::left, &mpi::right);
    //Along x direction
    MPI_Cart_shift(mpi::cart_comm, 0, 1, &mpi::top , &mpi::bot);

    //In case I want to know the rank2d of a particular set of coordinates.
    //int coords_query[2] = {0, 3};
    //int rank_q;
    //MPI_Cart_rank(mpi::cart_comm, coords_query, &rank_q);
    //printf("[Rank %d] coords (%d, %d)",rank_q,coords_query[0],coords_query[1]);
    
    //printf("[MPI process %d] I am located at (%d, %d). Top %d bot %d right %d left %d \n",
    //       mpi::rank2d, mpi::coords[0], mpi::coords[1], mpi::top, mpi::bot, mpi::right, mpi::left);
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

    RankErrorMessage();
    buildCartesianTopology();
    allocate_lattice_arrays();
    periodic_boundary(); //Compute right and left periodic boundary
   

    double m0 = -0.18840579710144945;
    int nconf = 12;
    double beta = 2;
	GaugeConf GConf = GaugeConf();  //Gauge configuration
    std::ostringstream NameData;
    NameData << "../confs/b" << beta << "_" << LV::Nx << "x" << LV::Nt << "/m-018/2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt << "_b" << 
    format(beta).c_str() << "_m" << format(m0).c_str() << "_" << nconf << ".ctxt";
    GConf.read_conf(NameData.str());

    spinor left_phi(mpi::maxSize), Dphi(mpi::maxSize), temp(mpi::maxSize); 
    spinor right_phi(mpi::maxSize);
    spinor Dphi_check(LV::Ntot); //For checking the result on rank 0
    re_field real_field(mpi::maxSize);
    re_field real_check(LV::Ntot);
    for(int n = 0; n <mpi::maxSize; n++) {
        left_phi.mu0[n] = 1; //RandomU1(); //spin up
        left_phi.mu1[n] = 2; //RandomU1(); //spin down
        right_phi.mu0[n] = 3; //RandomU1(); //spin up
        right_phi.mu1[n] = 4; //RandomU1(); //spin down
        
    }

    real_field = phi_dag_partialD_phi(GConf.Conf,left_phi,right_phi);
    double local_z = 0;
    //reduction over all lattice points and spin components
    for (int n = 0; n < mpi::maxSize; n++) {
        local_z += real_field.mu0[n];
        local_z += real_field.mu1[n];
    }
    double z;
    MPI_Allreduce(&local_z, &z, 1, MPI_DOUBLE, MPI_SUM, mpi::cart_comm);
    if (mpi::rank2d == 0){
       std::cout << "< , > = " << z << std::endl;
    }



  
    //Create a new data type for the blocks corresponding to each rank  
    /*
    MPI_Datatype sub_block_type;
    //int MPI_Type_vector(int block_count, int block_length, int stride, MPI_Datatype old_datatype, MPI_Datatype* new_datatype);
    MPI_Type_vector(mpi::width_x, mpi::width_t, LV::Nt, MPI_DOUBLE, &sub_block_type);
    MPI_Type_commit(&sub_block_type);

    //Resize the data type to use scatterV properly
    int extent = mpi::width_t;
    MPI_Datatype sub_block_resized;
    MPI_Type_create_resized(sub_block_type, 0, extent * sizeof(double), &sub_block_resized);
    MPI_Type_commit(&sub_block_resized);

    int counts[mpi::size];
    int displs[mpi::size];
    int offset;
    for(int i = 0; i < mpi::size; i++){
        counts[i] = 1;
        offset = (i/mpi::ranks_t);
        displs[i] = (i < mpi::ranks_t) ? i : (i-mpi::ranks_t*offset) + offset * mpi::ranks_t * mpi::width_x; 
    }
    
    MPI_Gatherv(real_field.mu0, mpi::maxSize, MPI_DOUBLE,
            real_check.mu0, counts, displs, sub_block_resized,
            0, mpi::cart_comm);

    MPI_Gatherv(real_field.mu1, mpi::maxSize, MPI_DOUBLE,
            real_check.mu1, counts, displs, sub_block_resized,
            0, mpi::cart_comm);
    

    if (mpi::rank2d == 0) {
        std::cout << "Reconstructed Dphi on rank 0" << std::endl;
        for(int n = 0; n<LV::Ntot; n++){
            std::cout << real_check.mu0[n] << "\n";
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
    
    

/* 
    spinor sol(mpi::maxSize), rhs(mpi::maxSize);

    for(int n = 0; n <mpi::maxSize; n++) {
        rhs.mu0[n] = 1;//RandomU1(); //spin up
        rhs.mu1[n] = 1;//RandomU1(); //spin down
    }
*/

    //SaveConf(GConf, "binaryConf");
    
    //GaugeConf GConfBinary = GaugeConf();
    //GConf.readBinary(NameData.str());

    /*
    double startT, endT;
    startT = MPI_Wtime();
    conjugate_gradient(GConf.Conf, rhs, sol,m0);
    endT = MPI_Wtime();
    printf("[rank %d] time elapsed during CG implementation: %.4fs.\n", mpi::rank, endT - startT);
    */
    
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

