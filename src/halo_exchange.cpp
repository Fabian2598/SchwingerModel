#include "halo_exchange.h"


void exchange_halo(c_double* phi){
	using namespace mpi;
    int row_size = 2 * width_t;
    //Send top row to top rank. Receive top row from bot rank.
    MPI_Sendrecv(&phi[idx(1,1,0)], row_size, MPI_DOUBLE_COMPLEX, top, 0,
        &phi[idx(width_x+1,1,0)], row_size, MPI_DOUBLE_COMPLEX, bot, 0,
        cart_comm, MPI_STATUS_IGNORE);

    //Send bot row to bot rank. Receive bot row from top rank.
    MPI_Sendrecv(&phi[idx(width_x,1,0)], row_size, MPI_DOUBLE_COMPLEX, bot, 1,
        &phi[idx(0,1,0)], row_size, MPI_DOUBLE_COMPLEX, top, 1,
        cart_comm, MPI_STATUS_IGNORE);

    //Send left column to left rank. Receive left column from right rank. 
    MPI_Sendrecv(&phi[idx(1,1,0)], 1, column_type, left, 2,
    &phi[idx(1,width_t+1,0)], 1, column_type, right, 2,
    cart_comm, MPI_STATUS_IGNORE);

    //Send right column to right rank. Receive right column from left rank. 
    MPI_Sendrecv(&phi[idx(1,width_t,0)], 1, column_type, right, 3,
    &phi[idx(1,0,0)], 1, column_type, left, 3,
    cart_comm, MPI_STATUS_IGNORE);
}

void exchange_halo_vec(c_double* v){
    using namespace mpi;
    int row_size = width_t;
    //Send top row to top rank. Receive top row from bot rank.
    MPI_Sendrecv(&v[idx_vec(1,1)], row_size, MPI_DOUBLE_COMPLEX, top, 0,
        &v[idx_vec(width_x+1,1)], row_size, MPI_DOUBLE_COMPLEX, bot, 0,
        cart_comm, MPI_STATUS_IGNORE);

    //Send bot row to bot rank. Receive bot row from top rank.
    MPI_Sendrecv(&v[idx_vec(width_x,1)], row_size, MPI_DOUBLE_COMPLEX, bot, 1,
        &v[idx_vec(0,1)], row_size, MPI_DOUBLE_COMPLEX, top, 1,
        cart_comm, MPI_STATUS_IGNORE);

    //Send left column to left rank. Receive left column from right rank. 
    MPI_Sendrecv(&v[idx_vec(1,1)], 1, column_type_vec, left, 2,
    &v[idx_vec(1,width_t+1)], 1, column_type_vec, right, 2,
    cart_comm, MPI_STATUS_IGNORE);

    //Send right column to right rank. Receive right column from left rank. 
    MPI_Sendrecv(&v[idx_vec(1,width_t)], 1, column_type_vec, right, 3,
    &v[idx_vec(1,0)], 1, column_type_vec, left, 3,
    cart_comm, MPI_STATUS_IGNORE);

    //Update top-right corner (needs bottom-left corner from diagonal rank)
    {
        c_double bottom_left = v[idx_vec(width_x,1)];   
        MPI_Send(&bottom_left, 1, MPI_DOUBLE_COMPLEX, mpi::bot_left, 0, mpi::cart_comm);
        MPI_Recv(&bottom_left, 1, MPI_DOUBLE_COMPLEX, mpi::top_right, 0, mpi::cart_comm, &status);
        v[idx_vec(0,width_t+1)] = bottom_left;  
    }

    //Update bottom-left corner (needs top-right corner from diagonal rank)
    {
        c_double top_right = v[idx_vec(1,width_t)];     //U_0(n+1-0)
        MPI_Send(&top_right, 1, MPI_DOUBLE_COMPLEX, mpi::top_right, 1, mpi::cart_comm);
        MPI_Recv(&top_right, 1, MPI_DOUBLE_COMPLEX, mpi::bot_left, 1 mpi::cart_comm, &status);
        v[idx_vec(width_x+1,0)] = top_right;
    }
    //Update top-left corner (needs bot-right corner from diagonal rank)
    {
        c_double bot_right = v[idx_vec(width_x,width_t)];   //U_0(n-1-0)
        MPI_Send(&bot_right, 1, MPI_DOUBLE_COMPLEX, mpi::bot_right, 2, mpi::cart_comm);
        MPI_Recv(&bot_right, 1, MPI_DOUBLE_COMPLEX, mpi::top_left, 2, mpi::cart_comm, &status);
        v[0]   = bot_right;

    }
    //Update bottom-right corner (needs top-left corner from diagonal rank)
    {
        c_double top_left = v[idx_vec(1,1)];   //U_(n+1+0)
        MPI_Send(&top_left, 1, MPI_DOUBLE_COMPLEX, mpi::top_left, 2, mpi::cart_comm);
        MPI_Recv(&top_left, 1, MPI_DOUBLE_COMPLEX, mpi::bot_right, 2, mpi::cart_comm, &status);
        v[idx_vec(width_x+1,width_t+1)]  = top_left;
    }

}