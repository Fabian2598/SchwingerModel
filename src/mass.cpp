
#include <string>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "gauge_conf.h"
#include "conjugate_gradient.h"

std::vector<c_matrix> Confs; //Confs[conf][coordinate][spin]

int nconf; //This could be an input parameter from terminal
double m0;
int max_iter = 10000;
std::vector<std::vector<double>> CorrMat; //Correlation function for each conf.
std::vector<double> Corr(Nt,0), dCorr(Nt,0); //Correlation function for each conf.

c_matrix source;
c_matrix Dcol; //D^-1 source 
c_matrix x0(Ntot, c_vector(2,1)); //Initial solution for inverting Dirac matrix
int coord;
double beta;


c_matrix canonical_vector(const int& i, const int& N1, const int& N2) {
	c_matrix e_i(N1, c_vector (N2,0.0));
	int j = i / N2;
	int k = i % N2;
	e_i[j][k] = 1.0;
	return e_i;
}

//Formats decimal numbers
static std::string format(const double& number) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << number;
    std::string str = oss.str();
    str.erase(str.find('.'), 1); //Removes decimal dot 
    return str;
}


void read_confs(const int& nconf){
    int x, t, mu;
    double re, im;
    for(int conf=0; conf<nconf; conf++){
        //Check how I actually want to do this part in general ...
        std::ostringstream formattedName;
        formattedName << "../confs/b1_16x16/2D_U1_Ns" << Ns
                    << "_Nt" << Nt
                    << "_b" << format(beta).c_str()
                    << "_m" << format(m0).c_str()
                    << "_" << conf << ".txt";

        std::string fileName = formattedName.str();
        std::ifstream file(fileName);
        if (!file.is_open()) {
            std::cerr << "Error opening file " << fileName << std::endl;
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
    /*
    std::cout << "beta: ";
    std::cin >> beta;
    std::cout << "m0: ";
    std::cin >> m0;
    std::cout << "Enter number of configurations (nconf): ";
    std::cin >> nconf;
    std::cout << " " << std::endl;
    */

    beta = 1.0; m0=0.1; nconf=1000;

    //These calls are necessary to use D_phi
    initialize_matrices(); //Intialize gamma matrices, identity and unit vectors
    Coordinates(); //Compute vectorized coordinates
    periodic_boundary(); //Compute right and left periodic boundary
    //---------------------------------------------//

    Confs.resize(nconf, c_matrix(Ns * Nt, c_vector(2, 0))); //Resize Confs
    CorrMat.resize(Nt, std::vector<double>(nconf, 0)); // Resize CorrMat

    std::cout << "Reading configurations...";
    read_confs(nconf); //Read configurations and store them in Confs
    std::cout << " Done!" << std::endl;
    //--------Compute c(nt) for the pion--------//
    for(int conf = 0; conf<nconf; conf++){
        for (int bet=0; bet<2; bet++){
            source = canonical_vector(bet, Ntot, 2); 
            Dcol = bi_cgstab(Confs[conf], source, x0, m0, max_iter, 1e-10, false); //D^-1 source
            for(int t=0; t<Nt; t++){
                for(int alf=0; alf<2; alf++){  
                    CorrMat[t][conf] -= std::real(Dcol[t][alf] * std::conj(Dcol[t][alf])); //|D^-1 (nt,0)_{alpha,beta}|^2  Dcol is a column of the dirac matrix
                }
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
    Name << "2D_U1_Ns" << Ns << "_Nt" << Nt << "_b" << beta << "_m" << m0 << "_" << "corr" << ".txt";
    write(Corr, dCorr, Name.str());
    
    

    return 0;
}

