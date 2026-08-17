#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "mpi_setup.h"


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
*/
inline void get_neighbors(const int x, const int t, int& right, int& down, int& left, int& up, double &rsign, double &lsign){
   int xp = x+1;
   int xm = x-1;
   int tp = t+1;
   int tm = t-1;
   int n = x * (mpi::width_t+2) + t;

   //Periodic boundary is already considered in the halo exchange due to the rank topology
   //Neighbor coordinates
   right   = x*(mpi::width_t+2)+tp;    //Right
   down    = xp*(mpi::width_t+2)+t;    //Down
   left    = x*(mpi::width_t+2)+tm;    //Left
   up      = xm*(mpi::width_t+2)+t;    //Up

   rsign = 1;
   lsign = 1;

	if ((mpi::rank2d+1) % mpi::ranks_t == 0){
		rsign = (t == mpi::width_t) ? -1 : 1;  //sign for the right boundary in time
	} 
	if (mpi::rank2d % mpi::ranks_t == 0){
		lsign = (t == 1) ? -1 : 1;             //sign for the left boundary in time	
	}
}

#endif