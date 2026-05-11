#include <iostream>
#include <vector>
#include <complex>
#include <iomanip>
#include <fstream>
#include <string>
typedef std::complex<double> c_double;


//------------Lattice parameters--------------//
namespace LV {
    //Lattice dimensions//
    constexpr int Nx= 256; 
    constexpr int Nt= 256; 
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

void read_conf(spinor& GlobalConf ,const std::string& name){
    std::ifstream infile(name);
    if (!infile) {
        std::cerr << "File " << name << " not found " << std::endl;
        exit(1);
    }
    int x, t, mu;
    double re, im; 
    while (infile >> x >> t >> mu >> re >> im) {
        if (mu == 0)
            GlobalConf.mu0[x*LV::Nt+t] = c_double(re, im); //x*LV::Nt+t;
        else
            GlobalConf.mu1[x*LV::Nt+t] = c_double(re, im); //x*LV::Nt+t;
    }
    infile.close();
    //std::cout << "Conf read from " << name << std::endl;                
}

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
    //std::cout << "Binary conf read from " << name << std::endl;     
}

/*
Save Gauge configuration
*/   
void SaveConf(spinor& GlobalConf, const std::string& Name) {
    using namespace LV;
    std::ofstream Datfile(Name,std::ios::binary);
    if (!Datfile.is_open()) {
        std::cerr << "Error opening file: " << Name << std::endl;
        return;
    }
    for (int x = 0; x < Nx; x++) {
    for (int t = 0; t < Nt; t++) {
    int n = x * Nt + t;
    for (int mu = 0; mu < 2; mu++) {
        const double& re = std::real(mu == 0 ? GlobalConf.mu0[n] : GlobalConf.mu1[n]);
        const double& im = std::imag(mu == 0 ? GlobalConf.mu0[n] : GlobalConf.mu1[n]);                          
        Datfile.write(reinterpret_cast<const char*>(&x), sizeof(int));
        Datfile.write(reinterpret_cast<const char*>(&t), sizeof(int));
        Datfile.write(reinterpret_cast<const char*>(&mu), sizeof(int));
        Datfile.write(reinterpret_cast<const char*>(&re), sizeof(double));
        Datfile.write(reinterpret_cast<const char*>(&im), sizeof(double)); 
    }
    }
    }
    Datfile.close();
} 

c_double TrivialAverage(std::vector<spinor> Confs,int nconf) {
    c_double result=0;
    for (int conf_id = 0; conf_id <nconf; conf_id++){
        result += Confs[conf_id].mu0[15];
    }
    return result/(nconf*1.0);

}


void read_rhs(std::vector<std::vector<c_double>>& vec,const std::string& name){
    std::ifstream infile(name);
    if (!infile) {
        std::cerr << "File " << name << " not found" << std::endl;
    }
    int x, t, mu;
    double re, im;
    //x, t, mu, real part, imaginary part
    while (infile >> x >> t >> mu >> re >> im) {
        vec[x*LV::Nt + t][mu] = c_double(re, im); 
    }
    infile.close();
  
}

void save_rhs(std::vector<std::vector<c_double>>& rhs,const std::string& name){
    std::ofstream rhsfile(name,std::ios::binary);
    if (!rhsfile.is_open()) {
        std::cerr << "Error opening rhs.txt for writing." << std::endl;
    } 
    else {
        //x, t, mu, real part, imaginary part
    for (int x = 0; x < LV::Nx; x++) {
    for (int t = 0; t < LV::Nt; t++) {
    int n = x * LV::Nt + t;
    for (int mu = 0; mu < 2; mu++) {
        const double& re = std::real(mu == 0 ? rhs[n][0] :  rhs[n][1]);
        const double& im = std::imag(mu == 0 ?  rhs[n][0] :  rhs[n][1]);          
        rhsfile.write(reinterpret_cast<const char*>(&x), sizeof(int));
        rhsfile.write(reinterpret_cast<const char*>(&t), sizeof(int));
        rhsfile.write(reinterpret_cast<const char*>(&mu), sizeof(int));
        rhsfile.write(reinterpret_cast<const char*>(&re), sizeof(double));
        rhsfile.write(reinterpret_cast<const char*>(&im), sizeof(double)); 
    }
    }
    }
    }

}

void read_rhs_binary(std::vector<std::vector<c_double>>& rhs, const std::string& name){
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
            rhs[n][mu] = c_double(re, im); 
           
    }
    }
    }
    infile.close(); 
}



    
//Formats decimal numbers
static std::string format(const double& number) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << number;
    std::string str = oss.str();
    str.erase(str.find('.'), 1); //Removes decimal dot 
    return str;
}

//For converting right-hand-sides
int main(){
    std::cout << "Nx " << LV::Nx << "  Nt" << LV::Nt << std::endl;
    std::string file;
    std::cout << "file: ";
    std::cin >> file;
    std::cout << " " << std::endl;
    std::vector<std::vector<c_double>> rhs(LV::Ntot, std::vector<c_double>(2, 0)); //right hand side
    read_rhs(rhs,file);

    c_double trivial1 = 0;
    c_double trivial2 = 0;
    for(int n = 0; n < LV::Ntot; n++){
        trivial1 += rhs[n][0]/(1.0*LV::Ntot);
    }
    std::cout << "trivial computation 1 " << trivial1 << std::endl;
    save_rhs(rhs,file);
    read_rhs_binary(rhs,file);
    for(int n = 0; n < LV::Ntot; n++){
        trivial2 += rhs[n][0]/(1.0*LV::Ntot);
    }
    std::cout << "trivial computation 2 " << trivial2 << std::endl;
    if (std::abs(trivial1-trivial2) > 1e-8)
        std::cout << "something wrong" << std::endl;
    else
        std::cout << "all good" << std::endl;



    return 0;
}
//For converting gauge confs to binary format
/*
int main(){
    std::cout << "Nx " << LV::Nx << "  Nt" << LV::Nt << std::endl;
    std::string listFilePath;
    int nconf; //This could be an input parameter from terminal
    std::cout << "Number of configurations (nconf): ";
    std::cin >> nconf;
    std::cout << "List of confs: ";
    std::cin >> listFilePath;
    std::cout << " " << std::endl;
    std::vector<spinor> Confs(nconf); //Confs[conf][coordinate][spin]

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

    std::cout << "Reading configurations...";
    for(int conf_id = 0; conf_id<nconf; conf_id++){
        read_conf(Confs[conf_id],filePaths[conf_id]);
    }
    std::cout << " Confs. read" << std::endl;

    c_double trivial1, trivial2;
    trivial1 = TrivialAverage(Confs,nconf);
    std::cout << "Trivial computation " << trivial1 << std::endl;

    std::cout << "Writing binary confs..." << std::endl;
    for(int conf_id = 0; conf_id<nconf; conf_id++){
        SaveConf(Confs[conf_id],filePaths[conf_id]);
    }
    std::cout << " done!" << std::endl;

    std::cout << "Checking that Trivial computation is the same with the new file format" << std::endl;
    std::cout << "Reading binary confs...";
    for(int conf_id = 0; conf_id<nconf; conf_id++){
        read_binary(Confs[conf_id],filePaths[conf_id]);
    }
    std::cout << " Confs. read" << std::endl;

    trivial2 = TrivialAverage(Confs,nconf);
    std::cout << "Trivial computation " << trivial2 << std::endl;
    if (std::abs(trivial1 - trivial2) > 1e-8)
        std::cout << "Something is wrong " << std::endl;
    else
        std::cout << "All good" << std::endl;

    
    return 0;
}
*/



