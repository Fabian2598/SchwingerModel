#include <iostream>
#include <vector>
#include <complex>
#include <iomanip>
#include <fstream>
#include <string>

/*
Program for converting gauge configurations in a binary format to human readable text
The lattice dimensions are fixed on lines 23 and 24
Inputs: 
-file path to the binary configuration
-name for the human readable text

The output file has the following layout
x, t, mu, Real part, Im part
*/

typedef std::complex<double> c_double;
//------------Lattice parameters--------------//
namespace LV {
    //Lattice dimensions//
    constexpr int Nx= 64; 
    constexpr int Nt= 64; 
    constexpr int Ntot = Nx*Nt; //Total number of lattice points
}

struct spinor {
    c_double* mu0;
    c_double* mu1;
    int size;
    //Constructor
    spinor(int N = LV::Ntot) : size(N) {
        mu0 = new c_double[N]();
        mu1 = new c_double[N]();
    }

    // Copy constructor (deep copy)
    spinor(const spinor& other) : size(other.size) {
        mu0 = new c_double[size];
        mu1 = new c_double[size];
        std::copy(other.mu0, other.mu0 + size, mu0);
        std::copy(other.mu1, other.mu1 + size, mu1);
    }

    // Assignment operator (deep copy)
    spinor& operator=(const spinor& other) {
        if (this != &other) {
            if (size != other.size) {
                delete[] mu0;
                delete[] mu1;
                size = other.size;
                mu0 = new c_double[size];
                mu1 = new c_double[size];
            }
            std::copy(other.mu0, other.mu0 + size, mu0);
            std::copy(other.mu1, other.mu1 + size, mu1);
        }
        return *this;
    }

    // Destructor
    ~spinor() {
        delete[] mu0;
        delete[] mu1;
    }

    inline void clearBuffer(){
        for(int n = 0; n<size; n++){
            mu0[n] = 0;
            mu1[n] = 0;
        }
    }
};

void read_binary(spinor& GlobalConf, const std::string& name){
    using namespace LV;
    std::ifstream infile(name, std::ios::binary);
    if (!infile) {
       std::cerr << "File " << name << " not found " << std::endl;
        exit(1);
    } 
    for (int x = 0; x < Nx; x++) {
    for (int t = 0; t < Nt; t++) {
        int n = x * Nt + t;
        for (int mu = 0; mu < 2; mu++) {
            int x_read, t_read, mu_read;
            double re, im;
            infile.read(reinterpret_cast<char*>(&x_read), sizeof(int));
            infile.read(reinterpret_cast<char*>(&t_read), sizeof(int));
            infile.read(reinterpret_cast<char*>(&mu_read), sizeof(int));
            infile.read(reinterpret_cast<char*>(&re), sizeof(double));
            infile.read(reinterpret_cast<char*>(&im), sizeof(double));
            if (mu == 0)
                GlobalConf.mu0[n] = c_double(re, im); //n --> for testing
            else
                GlobalConf.mu1[n] = c_double(re, im);
    }
    }
    }
    infile.close(); 
}

/*
Save Gauge configuration
x, t, mu, Real part, Im part
*/ 
void SaveConf(spinor& Conf, const std::string& Name) {
    std::ofstream Datfile(Name);
    if (!Datfile.is_open()) {
        std::cerr << "Error opening file: " << Name << std::endl;
        return;
    }
    using namespace LV;
    for (int x = 0; x < Nx; x++) {
        for (int t = 0; t < Nt; t++) {
            int n = x * Nt + t;
            for (int mu = 0; mu < 2; mu++) {
                const double& re = std::real(mu == 0 ? Conf.mu0[n] : Conf.mu1[n]);
                const double& im = std::imag(mu == 0 ? Conf.mu0[n] : Conf.mu1[n]); 
                Datfile << x
                        << std::setw(10) << t
                        << std::setw(10) << mu
                        << std::setw(30) << std::setprecision(17) << std::scientific << re
                        << std::setw(30) << std::setprecision(17) << std::scientific << im
                        << "\n";
            }
        }
    }
    Datfile.close();
}

int main(){
    using namespace LV;
    std::cout << "Nx " << LV::Nx << "  Nt" << LV::Nt << std::endl;
    std::string FilePath, Name;
    int nconf; //This could be an input parameter from terminal
    std::cout << "Conf path ";
    std::cin >> FilePath;
    std::cout << "Name for text file ";
    std::cin >> Name;
    std::cout << " " << std::endl;
    spinor Conf;
    read_binary(Conf, FilePath);
    SaveConf(Conf, Name);

    //For printing the confs
    /*
    for (int x = 0; x < Nx; x++) {
        for (int t = 0; t < Nt; t++) {
            int n = x * Nt + t;
            for (int mu = 0; mu < 2; mu++) {
                const double& re = std::real(mu == 0 ? Conf.mu0[n] : Conf.mu1[n]);
                const double& im = std::imag(mu == 0 ? Conf.mu0[n] : Conf.mu1[n]); 
                std::cout << x
                        << std::setw(5) << t
                        << std::setw(5) << mu
                        << std::setw(15) << std::setprecision(4) << std::scientific << re
                        << std::setw(15) << std::setprecision(4) << std::scientific << im
                        << "\n";
            }
        }
    }
    */

}