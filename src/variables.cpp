#include "variables.h"

namespace mpi{
    int rank = 0;
    int size = 1; 
    int maxSizeH = 2*(LV::Nx+2)*(LV::Nt+2); //maxSize with halos and spin included. For 1 rank this is the default value
    int sitesH = (mpi::width_x+2)*(mpi::width_t+2);
    int ranks_x = 1;
    int ranks_t = 1;
    int width_x = LV::Nx;
    int width_t = LV::Nt;
    int rank2d = 0; //linearize rank
    int coords[2] = {0,0}; //rank coords
    int top = 0; 
    int bot = 0; 
    int right = 0; 
    int left = 0;
    int bot_left = 0;
    int bot_right = 0;
    int top_left = 0;
    int top_right = 0;

    MPI_Comm cart_comm; //cartesian communicator
    //Datatypes for reading/writing gauge confs and rhs
    MPI_Datatype column_type;
    MPI_Datatype column_type_vec;
    MPI_Datatype global_conf_type;
    MPI_Datatype global_conf_resized;
    MPI_Datatype local_conf_type;
    MPI_Datatype local_conf_resized;

}

namespace CG{
	int max_iter = 10000; //Maximum number of iterations for the conjugate gradient method
	double tol = 1e-10;   //Tolerance for convergence
    bool print_convergence_message = false; //printing convergence message, useful for testing
}

namespace BiCG{
	int max_iter = 10000; //Maximum number of iterations for the bi-cgstab method
	double tol = 1e-10;   //Tolerance for convergence
    bool print_convergence_message = false; //printing convergence message, useful for testing
}

namespace sim_params {
    double beta=1;
    double m0=1;
    int MD_steps = 10;
    double trajectory_length = 1.0;
    int Ntherm = 10;
    int Nmeas = 10;
    int Nsteps = 1;
    std::string start_time_str = "";
    double tm = 0; //twisted mass term (default value)
    double csw = 0; //clover term (default value)
}
