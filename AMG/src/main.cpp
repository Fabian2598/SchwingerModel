#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "jackknife.h"
#include "tests.h"
#include "mpi.h"

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

    using namespace SAPV;
    MPI_Init(&argc, &argv);
    int rank, size; 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //srand(19);

    srand(time(0));
    
    Coordinates(); //Builds array with coordinates of the lattice points x * Nt + t
    MakeBlocks(); //Makes lattice blocks 
    periodic_boundary(); //Builds LeftPB and RightPB (periodic boundary for U_mu(n))
    Aggregates(); //build aggregates
    CheckAggregates();
    CheckBlocks(); //Check blocks dimensions
    
    m0 = -0.18840579710144945; //Globally declared
    AMGV::SAP_test_vectors_iterations = 2; //This parameter can change things quite a lot.
    FGMRESV::fgmres_restart_length = 25;
    AMGV::Nit = 3;
    //Parameters in variables.cpp
    if (rank == 0)
        print_parameters();
    
    GaugeConf GConf = GaugeConf(Nx, Nt);
    GConf.initialize(); //Initialize a random gauge configuration

    
    double beta = 2;
    int nconf = 20;
    if (LV::Nx == 64)
        nconf = 0;
    else if (LV::Nx == 128)
        nconf = 0;
    else if (LV::Nx == 256)
        nconf = 20;  
    
    //Reading Conf
    {
        std::ostringstream NameData;
        NameData << "../../confs/b" << beta << "_" << LV::Nx << "x" << LV::Nt << "/m-018/2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt << "_b" << 
        format(beta).c_str() << "_m" << format(m0).c_str() << "_" << nconf << ".ctxt";
        //GConf.read_conf(NameData.str());
    }

    sap.set_params(GConf.Conf, m0); //Setting gauge conf and m0 for SAP 

    spinor rhs(Ntot, c_vector(2, 0)); //right hand side
    spinor x0(Ntot, c_vector(2, 0)); //initial guess
    std::ostringstream FileName;
    FileName << "../../confs/rhs/rhs_conf" << nconf << "_" << Nx << "_Nt" << Nt << ".rhs";
    //read_rhs(rhs,FileName.str());
    //Random right hand side
    random_rhs(rhs,10);

    // Save rhs to a .txt file
    if (rank == 0){
        std::ostringstream FileName;
        FileName << "rhs_conf" << nconf << "_" << Nx << "_Nt" << Nt
                 << ".rhs";
        //save_rhs(rhs,FileName.str());
    }
    
    
    
    clock_t start, end;
    double elapsed_time;
    double startT, endT;

    spinor x_bi(Ntot, c_vector(2, 0));
    spinor x_cg(Ntot, c_vector(2, 0));
    spinor x_gmres(Ntot, c_vector(2, 0));
    spinor x_fsap(Ntot, c_vector(2, 0));
    spinor x_sap(Ntot, c_vector(2, 0));
    spinor x_famg(Ntot,c_vector(2,0));
    spinor x_amg(Ntot,c_vector(2,0));

    Tests tests(GConf,rhs,x0, m0);

    if (rank == 0){
        tests.BiCG(x_bi,10000,false,true); 
        //tests.GMRES(x_gmres,25, 100,false,true);
        //tests.CG(x_cg);
    }

    //tests.FGMRES_sap(x_fsap,false,true);
    //tests.SAP(x_sap,200,true);


    int Meas = 1;
    std::vector<double> iterations(Meas,0);
    if (rank == 0) std::cout << "--------------Flexible GMRES with AMG preconditioning--------------" << std::endl;

    for(int i = 0; i < Meas; i++){
        if (rank == 0) std::cout << "Meas " << i << std::endl;

        tests.FGMRES_2grid(x_famg,false, true);
    }
    if (rank == 0){
        std::cout << "Average iteration number over " << Meas << " runs: " << mean(iterations) << " +- " 
        << Jackknife_error(iterations, 5) << std::endl;
    
    }
            
    MPI_Barrier(MPI_COMM_WORLD);
    //tests.twoLevel(x_amg,false,true);
    //tests.check_solution(x_amg);
    
    MPI_Finalize();

    return 0;
}