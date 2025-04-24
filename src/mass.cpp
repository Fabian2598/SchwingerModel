
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
std::vector<double> Corr(Nt,0), dCorr(Nt,0); //Correlation function averaged over configurations and its error

c_matrix source; //source vector
c_matrix Dcol; //D^-1 source 
c_matrix x0(Ntot, c_vector(2,1)); //Initial solution for inverting Dirac matrix
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
        //Check how I actually want to do this part in general ...
        /*
        std::ostringstream formattedName;
        formattedName << "../confs/b1_16x16/2D_U1_Ns" << Ns
                    << "_Nt" << Nt
                    << "_b" << format(beta).c_str()
                    << "_m" << format(m0).c_str()
                    << "_" << conf << ".txt";

        std::string fileName = formattedName.str();
        */

        std::ifstream file(filePaths[conf]);
        if (!file.is_open()) {
            std::cerr << "Error opening file " << filePaths[conf] << std::endl;
            return;
        }
        while (file >> x >> t >> mu >> re >> im) {
            Confs[conf][x * Nt + t][mu] = c_double(re, im);
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
    initialize_matrices(); //Intialize gamma matrices, identity and unit vectors
    Coordinates(); //Compute vectorized coordinates
    periodic_boundary(); //Compute right and left periodic boundary
    //--------------------------------------------------//

    Confs.resize(nconf, c_matrix(Ns * Nt, c_vector(2, 0))); //Resize Confs
    CorrMat.resize(Nt, std::vector<double>(nconf, 0)); // Resize CorrMat

    std::cout << "Reading configurations...";
    read_confs(nconf,filePaths); //Read configurations and store them in Confs
    std::cout << " Done!" << std::endl;
    //--------Compute c(nt) for the pion--------//
    for(int conf = 0; conf<nconf; conf++){
        if (conf % 100 == 0) { std::cout << "--------Computing c(nt) for conf " << conf << "--------" << std::endl;} 
        for (int bet=0; bet<2; bet++){
            source = canonical_vector(bet, Ntot, 2); //Source indexed by bet (spin component of (nx,nt=0,0) )   
            Dcol = bi_cgstab(Confs[conf], source, x0, m0, max_iter, 1e-10, false); //D^-1 source = Dcol((nx,nt),alpha)
            
            for(int t=0; t<Nt; t++){
                for(int alf=0; alf<2; alf++){  
                    for(int x=0; x<Ns; x++){
                        coord = Coords[x][t]; //x*Nt + t
                        CorrMat[t][conf] += std::real(Dcol[coord][alf] * std::conj(Dcol[coord][alf])); //D^-1 (nt,0)_{alpha,beta} Confs is a column of the dirac matrix
                    }
                
                }
                CorrMat[t][conf] *= 1.0/std::sqrt(Ns); //Average over spatial coordinates
            }
        }
    }    
    
    //Write c(nt) and its error into a file
    for(int t=0; t<Nt; t++){
        for(int conf=0; conf<nconf; conf++){
            Corr[t] += CorrMat[t][conf]; //Sum over configurations
        }
        Corr[t] /= nconf; //Average over configurations
        dCorr[t] = Jackknife_error(CorrMat[t], 20); //Jackknife error with 20 bins.
        std::cout << "c(" << t << ") = " << Corr[t] << " +/- " << dCorr[t] << std::endl;
    } 
    
    std::ostringstream Name;
    Name << "2D_U1_Ns" << Ns << "_Nt" << Nt << "_b" << beta << "_m" << format(m0) << "_" << "corr" << ".txt";
    write(Corr, dCorr, Name.str());


    return 0;
}

