#include <time.h> 
#include <ctime>
#include <chrono>
#include "mpi_setup.h"
#include "hmc.h"


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi::size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi::rank);
        
    srand((mpi::rank+1)*time(0));
    
    int Ntherm, Nmeas, Nsteps, Nm0; //Simulation parameters
    double beta; //Beta range
    double trajectory_length; //HMC parameters
    int MD_steps;
    double m0; //bare mass
	int saveconf = 0; //Save configurations

    //To call the sequential program one has to choose ranks_x = ranks_t = 1
    if (mpi::rank == 0){
         //---Input data---//
        std::cerr << "  -----------------------------" << std::endl;
        std::cerr << "|  Two-flavor Schwinger model   |" << std::endl;
        std::cerr << "| Hybrid Monte Carlo simulation |" << std::endl;
        std::cerr << "  -----------------------------" << std::endl;
        std::cerr << "Nx " << LV::Nx << " Nt " << LV::Nt << std::endl;
        std::cerr << "ranks_x: " << std::endl;
        std::cin >> mpi::ranks_x;
        std::cerr << "ranks_t: " << std::endl;
        std::cin >> mpi::ranks_t;
        std::cerr << "m0: " << std::endl;
        std::cin >> m0;
        std::cerr << "Molecular dynamics steps: " << std::endl;
        std::cin >> MD_steps;
        std::cerr << "Trajectory length: " << std::endl;
        std::cin >> trajectory_length; 
        std::cerr << "beta: " << std::endl;
        std::cin >> beta;
        std::cerr << "Thermalization: " << std::endl;
        std::cin >> Ntherm;
        std::cerr << "Measurements: " << std::endl;
        std::cin >> Nmeas;
        std::cerr << "Step (sweeps between measurements): " << std::endl;
        std::cin >> Nsteps;
        std::cerr << "Save configurations yes/no (1 or 0): " << std::endl;
        std::cin >> saveconf;
        std::cerr << std::endl;
    }
    
    MPI_Bcast(&mpi::ranks_x, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(&mpi::ranks_t, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(&m0, 1, MPI_DOUBLE,  0, MPI_COMM_WORLD);
    MPI_Bcast(&MD_steps, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(&trajectory_length, 1, MPI_DOUBLE,  0, MPI_COMM_WORLD);
    MPI_Bcast(&beta, 1, MPI_DOUBLE,  0, MPI_COMM_WORLD);
    MPI_Bcast(&Ntherm, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(&Nmeas, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(&Nsteps, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(&saveconf, 1, MPI_INT,  0, MPI_COMM_WORLD);

    sim_params::m0 = m0;
    sim_params::beta = beta;
    sim_params::MD_steps = MD_steps;
    sim_params::trajectory_length = trajectory_length;
    sim_params::Ntherm = Ntherm;
    sim_params::Nmeas = Nmeas;
    sim_params::Nsteps = Nsteps;
    
    initializeMPI(); //2D rank topology
        
    //Get time
    std::string start_time_str;
    if (mpi::rank == 0) {
        auto now = std::chrono::system_clock::now();
        std::time_t now_c = std::chrono::system_clock::to_time_t(now);
        std::ostringstream tss;
        // format: YYYY-MM-DD HH:MM:SS 
        tss << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S");
        start_time_str = tss.str();
    }
    // broadcast start_time_str length and content so other ranks could log if needed
    int tlen = static_cast<int>(start_time_str.size());
    MPI_Bcast(&tlen, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (mpi::rank != 0) start_time_str.resize(tlen);
    MPI_Bcast(start_time_str.data(), tlen, MPI_CHAR, 0, MPI_COMM_WORLD);
    sim_params::start_time_str = start_time_str;


    std::ostringstream NameData;
    NameData << "2D_U1_" << LV::Nx << "x" << LV::Nt << "_m0";
    {
        std::ostringstream m0_stream;
        m0_stream << std::setprecision(4) << sim_params::m0;
        NameData << m0_stream.str();
    }
    NameData << "_SimData.txt";
    //Metadata file with simulation parameters
    std::ofstream Datfile;
    if (mpi::rank == 0){
        Datfile.open(NameData.str());
        Datfile << "#Date and time\n";
        Datfile << sim_params::start_time_str << "\n";
        Datfile << "#Host\n";
        const char* hostname = std::getenv("HOSTNAME");
        Datfile << (hostname ? hostname : "unknown") << "\n";
        Datfile << "#Nx      #Nt\n";
        Datfile << LV::Nx << std::setw(10) << LV::Nt << "\n";
        Datfile << "#ranks_x     #ranks_t     #ranks\n";
        Datfile << mpi::ranks_x << std::setw(15) << mpi::ranks_t << std::setw(15) << mpi::size << "\n";
        Datfile << "#beta     #Ntherm     #Nmeas     #Nsteps\n";
        Datfile <<  std::setprecision(17) << sim_params::beta << std::setw(11) 
        << sim_params::Ntherm << std::setw(11) << sim_params::Nmeas << std::setw(11) << sim_params::Nsteps << "\n";
        Datfile << "#trajectory_length     #MD_steps\n";
        Datfile << std::setprecision(17) << sim_params::trajectory_length << std::setw(30)
         << sim_params::MD_steps << "\n";
        Datfile << "#CG max iterations     #CG relative tolerance\n";
        Datfile << CG::max_iter << std::setw(30) << std::setprecision(17) << CG::tol << "\n";
        Datfile << "#m0\n";
        Datfile << std::setprecision(17) << sim_params::m0 << "\n";
        Datfile.close();
    }
    print_parameters();

    GaugeConf GConf = GaugeConf();  //Initial gauge configuration         
    HMC hmc = HMC(GConf,MD_steps, trajectory_length, Ntherm, Nmeas, Nsteps, beta, LV::Nx, LV::Nt, m0,saveconf);   
    double begin = MPI_Wtime();
    hmc.HMC_algorithm();
    double end = MPI_Wtime();

    if (mpi::rank == 0){
        std::cout << "Average plaquette value / volume: Ep = " << hmc.getEp() << " dEp = " << hmc.getdEp() << std::endl;
        std::cout << "Average gauge action / volume: gS = " << hmc.getgS() << " dgS = " << hmc.getdgS() << std::endl;
        std::cout << "Acceptance rate: " << hmc.getacceptance_rate(Nmeas+Nsteps*Nmeas) << std::endl;
        double elapsed_secs = end - begin;
        std::cout << "Execution time = " << elapsed_secs << " s" << std::endl;
        std::cout << "-------------------------------" << std::endl;
        Datfile.open(NameData.str(),std::ios::app);
        Datfile << "#Ep                           #dEp\n";
        Datfile << std::setprecision(17) << hmc.getEp() << std::setw(30) << hmc.getdEp() << "\n";
        Datfile << "#gS                           #dgS\n";
        Datfile << std::setprecision(17) << hmc.getgS() << std::setw(30) << hmc.getdgS() << "\n";
        Datfile << "#Acceptance rate\n";
        Datfile << std::setprecision(17) << hmc.getacceptance_rate(Nmeas+Nsteps*(Nmeas-1)) << "\n";
        Datfile << "#Execution time\n";
        Datfile << std::setprecision(17) << elapsed_secs;
        Datfile.close();
    }
    
    MPI_Finalize();

	return 0;
}

