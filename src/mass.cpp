
#include <string>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "gauge_conf.h"
#include "conjugate_gradient.h"

std::vector<c_matrix> Confs; //Confs[conf][coordinate][spin]

int nconf; //This could be an input parameter from terminal
double m0, beta; //Mass and coupling constant

int max_iter = 10000;
std::vector<std::vector<double>> CorrMat; //Correlation function for each conf.
std::vector<double> Corr(LV::Nt,0), dCorr(LV::Nt,0); //Correlation function averaged over configurations and its error

c_matrix source1, source2; //source vector
c_matrix Dcol1, Dcol2; //D^-1 source 
c_matrix x0(LV::Ntot, c_vector(2,1)); //Initial solution for inverting Dirac matrix
int coord;

std::string listFilePath;


c_matrix canonical_vector(const int& i, const int& N1, const int& N2) {
	c_matrix e_i(N1, c_vector (N2,0.0));
	int j = i / N2; //Lattice sitece
	int k = i % N2; //spin
	e_i[j][k] = 1.0;
	return e_i;
}

//Formats decimal numbers
static std::string format(const double& number) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << number;
    std::string str = oss.str();
    str.erase(str.find('.'), 1); //Removes decimal dot 
    return str;
}


void read_confs(const int& nconf,std::vector<std::string>& filePaths){
    int x, t, mu;
    double re, im;
    for(int conf=0; conf<nconf; conf++){
        std::ifstream file(filePaths[conf]);
        if (!file.is_open()) {
            std::cerr << "Error opening file " << filePaths[conf] << std::endl;
            return;
        }
        while (file >> x >> t >> mu >> re >> im) {
            Confs[conf][x * LV::Nt + t][mu] = c_double(re, im);
        }
        file.close();
    }

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
   

int main(){
    // Input parameters
    //beta = 1.0; m0=0.1; nconf=1000; listFilePath = "confs.txt"; 
    // dir/s/b *.txt > confs.txt
    
    std::cout << "beta: ";
    std::cin >> beta;
    std::cout << "m0: ";
    std::cin >> m0;
    std::cout << "Number of configurations (nconf): ";
    std::cin >> nconf;
    std::cout << "List of confs: ";
    std::cin >> listFilePath;
    std::cout << " " << std::endl;
    

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
        //std::cout << "File path: " << filePath << std::endl; // Print the file path
    }
    listFile.close();



    //------These calls are necessary to use D_phi------//
    //initialize_matrices(); //Intialize gamma matrices, identity and unit vectors
    Coordinates(); //Compute vectorized coordinates
    periodic_boundary(); //Compute right and left periodic boundary
    //--------------------------------------------------//

    Confs.resize(nconf, c_matrix(LV::Ntot, c_vector(2, 0))); //Resize Confs
    CorrMat.resize(LV::Nt, std::vector<double>(nconf, 0)); // Resize CorrMat

    std::cout << "Reading configurations...";
    read_confs(nconf,filePaths); //Read configurations and store them in Confs
    std::cout << " Done!" << std::endl;
    //--------Compute c(nt) for the pion--------//
    for(int conf = 0; conf<nconf; conf++){
        if (conf % 100 == 0) { std::cout << "--------Computing c(nt) for conf " << conf << "--------" << std::endl;} 
        //We only need two sources, equivalent to extracting the first two columns of D^-1
        source1 = canonical_vector(0, LV::Ntot, 2); 
        source2 = canonical_vector(1, LV::Ntot, 2); 
        Dcol1 = bi_cgstab(Confs[conf], source1, x0, m0, max_iter, 1e-10, false); //D^-1 source = D^-1((nx,nt),0)
        Dcol2 = bi_cgstab(Confs[conf], source2, x0, m0, max_iter, 1e-10, false); //D^-1 source = D^-1((nx,nt),1)
        
        for(int t=0; t<LV::Nt; t++){
            CorrMat[t][conf] = 0;
            for(int x=0; x<LV::Nx; x++){
                coord = Coords[x][t]; //x*Nt + t
                CorrMat[t][conf] += std::real(Dcol1[coord][0] * std::conj(Dcol1[coord][0]))
                + std::real(Dcol1[coord][1] * std::conj(Dcol1[coord][1]))  
                + std::real(Dcol2[coord][0] * std::conj(Dcol2[coord][0]))
                + std::real(Dcol2[coord][1] * std::conj(Dcol2[coord][1])); 
            }
            CorrMat[t][conf] *= 1.0/std::sqrt(LV::Nx); //Average over spatial coordinates
        }            
    }
    
        
    //Write c(nt) and its error into a file
    for(int t=0; t<LV::Nt; t++){
        for(int conf=0; conf<nconf; conf++){
            Corr[t] += CorrMat[t][conf]; //Sum over configurations
        }
        Corr[t] /= nconf; //Average over configurations
        dCorr[t] = Jackknife_error(CorrMat[t], 20); //Jackknife error with 20 bins.
        std::cout << "c(" << t << ") = " << Corr[t] << " +/- " << dCorr[t] << std::endl;
    } 
    
    std::ostringstream Name;
    Name << "2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt << "_b" << beta << "_m" << format(m0) << "_" << "corr" << ".txt";
    write(Corr, dCorr, Name.str());

    return 0;
}

