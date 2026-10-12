#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "mpi_setup.h"

/*
    Geometry of the local lattice of this rank. Global variables cannot be read on the GPU, so the kernels
    capture a copy of this struct (see local_lattice() below) to find neighbors and boundary signs.
*/
struct LocalLattice {
    int width_x;  //number of lattice sites on x (no halo)
    int width_t;  //number of lattice sites on t (no halo)
    bool first_t; //the rank holds the first time slice of the global lattice (antiperiodic boundary in time)
    bool last_t;  //the rank holds the last time slice of the global lattice
};

inline LocalLattice local_lattice(){
    LocalLattice lat;
    lat.width_x = mpi::width_x;
    lat.width_t = mpi::width_t;
    lat.first_t = (mpi::rank2d % mpi::ranks_t == 0);
    lat.last_t  = ((mpi::rank2d+1) % mpi::ranks_t == 0);
    return lat;
}

/*
mu = 0 -> t
mu = 1 -> x
                    t
    0  +--------------------------+  Nt
       |                          |
       |            up    		    |
       |                          |
    x  |    left    n  right      |
       |                          |
       |           down           |
       |                          |
    Nx +--------------------------+ Nt

    n = x * (width_t+2) + t is the local index of the site (x,t), halo included
*/
HOST_DEVICE inline void get_neighbors(const LocalLattice& lat, const int n, int& right, int& down, int& left, int& up, double &rsign, double &lsign){
   const int t = n % (lat.width_t+2);

   //Periodic boundary is already considered in the halo exchange due to the rank topology
   //Neighbor coordinates
   right   = n+1;                //Right (x,t+1)
   down    = n+(lat.width_t+2);  //Down  (x+1,t)
   left    = n-1;                //Left  (x,t-1)
   up      = n-(lat.width_t+2);  //Up    (x-1,t)

   rsign = (lat.last_t  && t == lat.width_t) ? -1 : 1;  //sign for the right boundary in time
   lsign = (lat.first_t && t == 1) ? -1 : 1;            //sign for the left boundary in time
}

HOST_DEVICE inline void get_corners(const LocalLattice& lat, const int n, int& x1_t_1, int& x_1_t_1,int& x_1_t1, int& x1_t1){
   const int row = lat.width_t+2;
   x1_t_1  = n+row-1; //(x+1,t-1)
   x_1_t_1 = n-row-1; //(x-1,t-1)
   x_1_t1  = n-row+1; //(x-1,t+1)
   x1_t1   = n+row+1; //(x+1,t+1)
}

#endif
