
#include <string>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "gauge_conf.h"
#include "conjugate_gradient.h"


int nconf; //This could be an input parameter from terminal
double m0, beta; //Mass and coupling constant

int max_iter = 100000;
int coord;

std::string listFilePath;

//Formats decimal numbers
static std::string format(const double& number) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << number;
    std::string str = oss.str();
    str.erase(str.find('.'), 1); //Removes decimal dot 
    return str;
}


void write(const std::vector<double>& data,const std::vector<double>& error, const std::string& Name){
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
   

int main(int argc, char **argv){
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi::size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi::rank);
    allocate_lattice_arrays();
    // Input parameters
    //beta = 1.0; m0=0.1; nconf=1000; listFilePath = "confs.txt"; 
    // dir/s/b *.txt > confs.txt
    beta = 2, m0 = 0.016129032258064502, nconf = 1000; listFilePath = "confFiles.txt";
    /*if (rank == 0){
        std::cout << "beta: ";
        std::cin >> beta;
        std::cout << "m0: ";
        std::cin >> m0;
        std::cout << "Number of configurations (nconf): ";
        std::cin >> nconf;
        std::cout << "List of confs: ";
        std::cin >> listFilePath;
        std::cout << " " << std::endl;
    }*/
    
    
    periodic_boundary(); //Compute right and left periodic boundary

    // Read the list of file paths from the specified file
    std::ifstream listFile(listFilePath);

    if (!listFile.is_open()) {
        std::cerr << "Error: Could not open the file containing the list of paths: " << listFilePath << std::endl;
        return 1;
    }

    std::vector<std::string> filePaths;
    std::string filePath;
    //Read each file path from the list
    while (std::getline(listFile, filePath)) {
        filePaths.push_back(filePath);
    }
    listFile.close();

    std::vector<spinor> Confs(nconf);
    GaugeConf GConf;
    std::cout << "Reading configurations...";
    for(int conf=0; conf<nconf; conf++){
        //GConf.readBinary(filePaths[conf])
        GConf.read_conf(filePaths[conf]);
        Confs[conf] = GConf.Conf;
    }
    std::cout << " Done!" << std::endl;


    double* CorrMat; double* CorrMat_loc;
    CorrMat = new double[LV::Nt*nconf]; //Resize CorrMat
    CorrMat_loc = new double[LV::Nt*nconf];

    
    
    //We only need two sources, equivalent to extracting the first two columns of D^-1
    spinor source1(mpi::maxSize);
    spinor source2(mpi::maxSize);
    spinor Dcol1(mpi::maxSize);
    spinor Dcol2(mpi::maxSize);
    spinor x0(mpi::maxSize);

    if (mpi::rank == 0){
        source1.mu0[0] = 1;
        source2.mu1[0] = 1;
    }
    int x_width = (mpi::rank != mpi::size-1) ? (LV::Nx/mpi::size) : (LV::Nx/mpi::size) + (LV::Nx%mpi::size);

    /*
    //--------Compute c(nt) for the pion--------//
    for(int conf = 0; conf<nconf; conf++){
        if (conf % 100 == 0) { std::cout << "--------Computing c(nt) for conf " << conf << "--------" << std::endl;} 
        Dcol1 = bi_cgstab(Confs[conf], source1, x0, m0, max_iter, 1e-10, false); //D^-1 source = D^-1((nx,nt),0)
        Dcol2 = bi_cgstab(Confs[conf], source2, x0, m0, max_iter, 1e-10, false); //D^-1 source = D^-1((nx,nt),1)   
        for(int t=0; t<LV::Nt; t++){
            CorrMat_loc[t*nconf+conf] = 0;
            for(int x=0; x<x_width; x++){
                coord = x*LV::Nt+t; //x*Nt + t
                CorrMat_loc[t*nconf+conf] += std::real(Dcol1.mu0[coord] * std::conj(Dcol1.mu0[coord]))
                + std::real(Dcol1.mu1[coord] * std::conj(Dcol1.mu1[coord]))  
                + std::real(Dcol2.mu0[coord] * std::conj(Dcol2.mu0[coord]))
                + std::real(Dcol2.mu1[coord] * std::conj(Dcol2.mu1[coord])); 
            }
            CorrMat_loc[t*nconf+conf] *= 1.0/std::sqrt(LV::Nx); //Average over spatial coordinates
        }            
    }

    MPI_Allreduce(&CorrMat_loc, &CorrMat, LV::Nt*nconf, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
        
    //Write c(nt) and its error into a file
    if (mpi::rank == 0) {
        std::vector<std::vector<double>> CorrelationMatrix(LV::Nt, std::vector<double>(nconf, 0)); //Correlation function for each conf.
        std::vector<double> Corr(LV::Nt,0), dCorr(LV::Nt,0); //Correlation function averaged over configurations and its error
        for(int t=0; t<LV::Nt; t++){
        for(int conf=0; conf<nconf; conf++){
            CorrelationMatrix[t][conf] = CorrMat[t*nconf+conf];
            Corr[t] += CorrMat[t*nconf+conf]; //Sum over configurations
        }
        Corr[t] /= nconf; //Average over configurations
        dCorr[t] = Jackknife_error(CorrelationMatrix[t], 20); //Jackknife error with 20 bins.
        std::cout << "c(" << t << ") = " << Corr[t] << " +/- " << dCorr[t] << std::endl;
        } 

        std::ostringstream Name;
        Name << "2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt << "_b" << beta << "_m" << format(m0) << "_" << "corr" << ".txt";
        write(Corr, dCorr, Name.str());
    }
    */
    
    delete[] CorrMat;
    delete[] CorrMat_loc;
    free_lattice_arrays();
    
    MPI_Finalize();

    return 0;
}

