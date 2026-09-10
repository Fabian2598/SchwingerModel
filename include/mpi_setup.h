#ifndef MPI_SETUP_H
#define MPI_SETUP_H
#include "variables.h"

//Check that number of ranks on the x and t direction match the number of total ranks called.
inline void assignWidth(){
    if (mpi::ranks_t * mpi::ranks_x != mpi::size){
        if (mpi::rank == 0){
            std::cerr << "ranks_t * ranks_x != total number of ranks" << std::endl;
            std::cerr << mpi::ranks_t * mpi::ranks_x << " != " << mpi::size << std::endl;
        }
        exit(1);
    }
    //We do this to enforce an equal workload on each rank
    if (LV::Nx % mpi::ranks_x!= 0 ||LV::Nt % mpi::ranks_t != 0){
        if (mpi::rank == 0)
            std::cerr << "Nx (Nt) is not exactly divisible by rank_x (rank_t)" << std::endl;
        exit(1);
    }

    if ((mpi::ranks_x == 1 && mpi::ranks_t != 1) || (mpi::ranks_x != 1 && mpi::ranks_t == 1)){
        if (mpi::rank == 0){
            std::cerr << "Unsupported MPI topology: this code requires both Cartesian dimensions to have more than one rank. "
                      << "Got ranks_x = " << mpi::ranks_x << " and ranks_t = " << mpi::ranks_t << ". "
                      << "Please run with a 2D grid such as (1,1), (2,2), (2,3), etc., not (2,1) or (1,2)."
                      << std::endl;
        }
        exit(1);
    }

    mpi::width_x = LV::Nx/mpi::ranks_x;
    mpi::width_t = LV::Nt/mpi::ranks_t;
    mpi::maxSizeH = 2*(mpi::width_x+2)*(mpi::width_t+2); //With halos included
    mpi::sitesH = (mpi::width_x+2)*(mpi::width_t+2); 
}
 
/*
 * Two-dimensional cartesian topology
 *                t                    2D parallelization
 *   0  +-------------------+  Nt   +-----------------------+
 *      |                   |       |                       |
 *      |                   |       | top-left top top-right|
 *      |                   |       |           |           |
 *   x  |                   |       |--left--rank2d--right--|
 *      |                   |       |           |           |
 *      |                   |       | bot-left bot bot-right|
 *      |                   |       |                       |
 *   Nx +-------------------+ Nt    +-----------------------+
 *                Nx
*/
inline void buildCartesianTopology(){
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

    //Diagonal ranks (needed for the staples)
    int coords_bot_left[2] = {mod(mpi::coords[0]+1,mpi::ranks_x), mod(mpi::coords[1]-1,mpi::ranks_t)}; //bot-left
    MPI_Cart_rank(mpi::cart_comm, coords_bot_left, &mpi::bot_left);

    int coords_bot_right[2] = {mod(mpi::coords[0]+1,mpi::ranks_x), mod(mpi::coords[1]+1,mpi::ranks_t)}; //bot-right
    MPI_Cart_rank(mpi::cart_comm, coords_bot_right, &mpi::bot_right);

    int coords_top_left[2] = {mod(mpi::coords[0]-1,mpi::ranks_x), mod(mpi::coords[1]-1,mpi::ranks_t)}; //top-left
    MPI_Cart_rank(mpi::cart_comm, coords_top_left, &mpi::top_left);

    int coords_top_right[2] = {mod(mpi::coords[0]-1,mpi::ranks_x), mod(mpi::coords[1]+1,mpi::ranks_t)}; //top-right
    MPI_Cart_rank(mpi::cart_comm, coords_top_right, &mpi::top_right);
    
    //printf("[MPI process %d] I am located at (%d, %d). Top %d bot %d right %d left %d bot-left %d bot-right %d top-left %d top-right %d  \n",
    //       mpi::rank2d, mpi::coords[0], mpi::coords[1], mpi::top, mpi::bot, mpi::right, mpi::left,
    //        mpi::bot_left,mpi::bot_right,mpi::top_left,mpi::top_right);
}

inline void defineDataTypes(){
    //Create a new data type for the blocks corresponding to each rank
    /*  
    *              width_t
    *          ---------------     
    *          |             |
    *          |             |
    * width_x  |             |
    *          |             |
    *          |             |
    *          ---------------
    */
    //int MPI_Type_vector(int block_count, int block_length, int stride, MPI_Datatype old_datatype, MPI_Datatype* new_datatype);
   
    MPI_Type_vector(mpi::width_x,mpi::width_t*2,2*(mpi::width_t+2),MPI_DOUBLE_COMPLEX, &mpi::local_conf_type);
    MPI_Type_commit(&mpi::local_conf_type);

    //The displacement of local_domain_resized is in units of std::complex<double>
    MPI_Type_create_resized(mpi::local_conf_type, 0, sizeof(std::complex<double>), &mpi::local_conf_resized);
    MPI_Type_commit(&mpi::local_conf_resized);

    // Gather inner domains from all ranks in the coarse communicator
    // Buffer has size (Nx_coarse_rank+2)*(Nt_coarse_rank+2)*DOF
    // Create a type that matches the global buffer layout (strided by full global row including halo)
    MPI_Type_vector(mpi::width_x,                 // number of rows to place per rank
        2 * mpi::width_t,               		// elements per row (complex numbers)
        2 * (LV::Nt + 2),  	// stride between rows in global buffer (complex elements) including halo
            MPI_DOUBLE_COMPLEX,
            &mpi::global_conf_type);
    MPI_Type_commit(&mpi::global_conf_type);

    // Resize type so displacements are specified in units of one complex element
    MPI_Type_create_resized(mpi::global_conf_type, 0, sizeof(std::complex<double>), &mpi::global_conf_resized);
    MPI_Type_commit(&mpi::global_conf_resized);


    //Datatype for the halo exchange (spinor: 2 components per site)
    MPI_Type_vector(mpi::width_x, 2, 2*(mpi::width_t+2), MPI_DOUBLE_COMPLEX, &mpi::column_type);
    MPI_Type_commit(&mpi::column_type);

    //Datatype for the halo exchange of vectors (1 component per site, no spinor factor)
    MPI_Type_vector(mpi::width_x, 1, (mpi::width_t+2), MPI_DOUBLE_COMPLEX, &mpi::column_type_vec);
    MPI_Type_commit(&mpi::column_type_vec);
}

inline void initializeMPI(){
    assignWidth();
    buildCartesianTopology();
    defineDataTypes();
}


#endif