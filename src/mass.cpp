#include "gauge_conf.h"
#include "bi_cgstab.h"
#include <cstring>

//Read configurations from a list of files
//Confs is passed by reference so we can fill the caller's buffers.
void read_confs_from_list(const int nconf, std::vector<spinor*>& Confs, const std::vector<std::string>& filePaths){
    GaugeConf GConf;
    for(int conf=0; conf<nconf; conf++){
        GConf.ReadConf(filePaths[conf]);
        // Copy the read configuration into the pre-allocated spinor buffer
        *Confs[conf] = GConf.Conf;
    }
}

//Write pion correlator to disk
void write_correlator(const std::vector<double>& data,const std::vector<double>& error, std::string Name){
    std::ofstream Datfile;
    Datfile.open(Name);
    if (!Datfile.is_open()) {
        std::cerr << "Error opening file: " << Name << std::endl;
        return;
    }
    for (int i = 0; i < data.size(); i++) {
        Datfile << i << std::setw(30) << std::scientific << data[i] << std::setw(30) <<  std::scientific << error[i] <<"\n";
    }
    Datfile.close();
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi::size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi::rank);

    std::vector<spinor*> Confs; //Vector with the gauge configurations
    int nconf = 0; 
    double m0, beta; 
    std::string listFilePath;

    if (mpi::rank == 0){
         //---Input data---//
        std::cerr << "  -----------------------------" << std::endl;
        std::cerr << "|  Pion correlator computation  |" << std::endl;
        std::cerr << "  -----------------------------" << std::endl;
        std::cerr << "Nx " << LV::Nx << " Nt " << LV::Nt << std::endl;
        std::cerr << "ranks_x: " << std::endl;
        std::cin >> mpi::ranks_x;
        std::cerr << "ranks_t: " << std::endl;
        std::cin >> mpi::ranks_t;
        std::cerr << "m0: " << std::endl;
        std::cin >> m0;
        std::cerr << "beta: " << std::endl;
        std::cin >> beta;
        std::cerr << "File with list of confs (ls -1 *.ctxt > confFiles.txt): ";
        std::cin >> listFilePath;
        std::cerr << " " << std::endl;
    }
    MPI_Bcast(&mpi::ranks_x, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(&mpi::ranks_t, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(&m0, 1, MPI_DOUBLE,  0, MPI_COMM_WORLD);
    MPI_Bcast(&beta, 1, MPI_DOUBLE,  0, MPI_COMM_WORLD);
    sim_params::m0 = m0;
    sim_params::beta = beta;

    initializeMPI(); //initialize 2d communicator and MPI datatypes 
    
 
    std::vector<std::string> filePaths; //file path of configurations to read
    std::string filePath;   

    // Rank 0 opens the list file, reads paths and broadcasts them to all ranks.
    if (mpi::rank2d == 0) {
        std::ifstream listFile(listFilePath);
        if (!listFile.is_open()) {
            std::cerr << "Error: Could not open the file containing the list of paths: " << listFilePath << std::endl;
            return 1;
        }
        while (std::getline(listFile, filePath)) {
            // skip empty lines
            if (!filePath.empty()) filePaths.push_back(filePath);
        }
        listFile.close();
        nconf = static_cast<int>(filePaths.size());
    }

    // Broadcast number of configurations to all ranks
    MPI_Bcast(&nconf, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Broadcast each file path (length + chars)
    for (int i = 0; i < nconf; ++i) {
        int len = 0;
        if (mpi::rank2d == 0) len = static_cast<int>(filePaths[i].size()) + 1; // include null
        MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
        char* buf = new char[len];
        if (mpi::rank2d == 0) std::strcpy(buf, filePaths[i].c_str());
        MPI_Bcast(buf, len, MPI_CHAR, 0, MPI_COMM_WORLD);
        if (mpi::rank2d != 0) filePaths.push_back(std::string(buf));
        delete[] buf;
    }
    
    if (mpi::rank2d == 0)
        std::cout << "#" << nconf <<  " confs in " << listFilePath << std::endl;
    
    for(int confID = 0; confID<nconf; confID++){
        spinor* temp = new spinor(mpi::maxSizeH);
        Confs.push_back(temp);
    }


    
    std::vector<std::vector<double>> CorrMat; //Correlation function for each conf.
    if (mpi::rank2d == 0){
        CorrMat.resize(LV::Nt, std::vector<double>(nconf, 0)); // Resize CorrMat
        std::cout << "Reading configurations from list ...";
    }
    read_confs_from_list(nconf, Confs, filePaths); //Read configurations and store them in Confs
    if (mpi::rank2d == 0) std::cout << " Done!" << std::endl;
    
   

    //Buffers
    spinor source1(mpi::maxSizeH), source2(mpi::maxSizeH); //source vector
    spinor Dcol1(mpi::maxSizeH), Dcol2(mpi::maxSizeH); //D^-1 source 
    spinor x0(mpi::maxSizeH); //Initial solution

    for(int nx = 1; nx<=mpi::width_x; nx++)
        for(int nt = 1; nt<=mpi::width_t; nt++)
            for(int mu=0; mu<2; mu++)
                x0.val[idx(nx,nt,mu)]=1;

        
      
    if (mpi::rank2d == 0){
        int nx = 1, nt = 1;
        int n0 = nx*(mpi::width_t+2) + nt;
        source1.val[2*n0]   = 1; //First lattice site with mu=0
        source2.val[2*n0+1] = 1; //First lattice site with mu=1
    }

    
    //--------Compute c(nt) for the pion--------//
    
    for(int confID = 0; confID<nconf; confID++){
        if (confID % 100 == 0 && mpi::rank2d == 0)
            std::cout << "--------Computing c(nt) for conf " << confID << "--------" << std::endl; 
        //We only need two sources, equivalent to extracting the first two columns of D^-1
        exchange_halo(Confs[confID]->val);
        bi_cgstab(*Confs[confID], source1, x0, Dcol1); //D^-1 source = D^-1((nx,nt),0)
        bi_cgstab(*Confs[confID], source2, x0, Dcol2); //D^-1 source = D^-1((nx,nt),1)
        //Up to this part everything works fine in parallel

        //This particular part only works fine if I have one rank. For more than one rank 
        //I have to think on how to modify this to consider the contribution for the other ranks
        //Note: I will have to transfer from the coordinates in one rank to its global coordinates
        int n;
        for(int t=1; t<=mpi::width_t; t++){
            double local_contribution = 0;
            double global_contribution = 0;
            for(int x=1; x<=mpi::width_x; x++){
                n = x*(mpi::width_t+2)+t;
                local_contribution += std::real(Dcol1.val[2*n] * std::conj(Dcol1.val[2*n]))
                + std::real(Dcol1.val[2*n+1]    * std::conj(Dcol1.val[2*n+1]))  
                + std::real(Dcol2.val[2*n]      * std::conj(Dcol2.val[2*n]))
                + std::real(Dcol2.val[2*n+1]    * std::conj(Dcol2.val[2*n+1])); 
            }

            local_contribution *= 1.0/std::sqrt(LV::Nx); //Average over spatial coordinates
            MPI_Allreduce(&local_contribution, &global_contribution, 1, MPI_DOUBLE, MPI_SUM, mpi::cart_comm);
            if (mpi::rank2d == 0) CorrMat[t-1][confID] = global_contribution;
    
        } 
            
    }
  
    //Write c(nt) and its error into a file
    //Only the root rank performs this part
    if (mpi::rank2d == 0){
        std::vector<double> Corr(LV::Nt,0), dCorr(LV::Nt,0); //Correlation function averaged over configurations and its error
        for(int t=0; t<LV::Nt; t++){
            for(int confID=0; confID<nconf; confID++){
                Corr[t] += CorrMat[t][confID]; //Sum over configurations
            }
            Corr[t] /= nconf; //Average over configurations
            dCorr[t] = Jackknife_error(CorrMat[t], 20); //Jackknife error with 20 bins.
            std::cout << "c(" << t << ") = " << Corr[t] << " +/- " << dCorr[t] << std::endl;
        } 
        std::ostringstream Name;
        Name << "2D_U1_" << LV::Nx << "x" << LV::Nt << "_b" << beta << "_m" << format(m0) << "_" << "corr" << ".txt";
        write_correlator(Corr, dCorr, Name.str());
    }
    
    for (auto ptr : Confs) delete ptr;
    

    MPI_Finalize();
    return 0;
}