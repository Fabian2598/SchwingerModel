#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <complex>  
#include <vector>

//-------Lattice dimensions-------//
constexpr int Ns = 16, Nt = 16;
constexpr int Ntot = Ns * Nt;
//--------------------------------//

double pi = 3.14159265359;
typedef std::complex<double> c_double;
typedef std::vector<c_double> c_vector;
typedef std::vector<c_vector> c_matrix;


c_matrix Conf(Ntot, c_vector(2, 0)); //Gauge configuration
std::vector<std::vector<int>>Coords(Ns, std::vector<int>(Nt, 0)); //Vectorized coordinates
//Periodic boundary
std::vector<std::vector<std::vector<int>>>LeftPB(Ns, std::vector<std::vector<int>>(Nt, std::vector<int>(2, 0))); 
std::vector<std::vector<std::vector<int>>>RightPB(Ns, std::vector<std::vector<int>>(Nt, std::vector<int>(2, 0))); 

std::vector<c_matrix> gamma_mat(2,c_matrix(2, c_vector(2, 0)));  //Pauli matrices
c_double I_number(0, 1); //imaginary number
c_matrix Identity(Ntot, c_vector(2, 0));
std::vector<std::vector<int>>hat_mu(Ntot, std::vector<int>(2, 0)); //hat_mu[mu][2] = {hat_mu_t, hat_mu_x}

//Vectorized coordinates
void Coordinates() {
    for (int x = 0; x < Ns; x++) {
        for (int t = 0; t < Nt; t++) {
            Coords[x][t] = x * Nt + t;
        }
    }
}

inline int mod(int a, int b) {
    int r = a % b;
    return r < 0 ? r + b : r;
}

//Intialize gamma matrices, identity and unit vectors
void initialize_matrices() {
    //sigma_0 --> time
    gamma_mat[0][0][0] = 0; gamma_mat[0][0][1] = 1;
    gamma_mat[0][1][0] = 1; gamma_mat[0][1][1] = 0;
    //simgma_1 --> space
    gamma_mat[1][0][0] = 0; gamma_mat[1][0][1] = -I_number;
    gamma_mat[1][1][0] = I_number; gamma_mat[1][1][1] = 0;
    //2d identity
    Identity[0][0] = 1; Identity[0][1] = 0;
    Identity[1][0] = 0; Identity[1][1] = 1;
    hat_mu[0] = { 1, 0 }; //hat_t
    hat_mu[1] = { 0, 1 }; //hat_x
}

//right periodic boundary x+hat{mu}
//left periodic boundary x-hat{mu}
//hat_mu[0] = { 1, 0 }; //hat_t
//hat_mu[1] = { 0, 1 }; //hat_x
void periodic_boundary() {
    for (int x = 0; x < Ns; x++) {
        for (int t = 0; t < Nt; t++) {
            for (int mu = 0; mu < 2; mu++) {
                RightPB[x][t][mu] = Coords[mod(x + hat_mu[mu][1], Ns)][mod(t + hat_mu[mu][0], Nt)];
                LeftPB[x][t][mu] = Coords[mod(x - hat_mu[mu][1], Ns)][mod(t - hat_mu[mu][0], Nt)];

            }
        }
    }
}

//right fermionic boundary (antiperiodic in time) x+hat{mu}
inline c_double rfb(const c_matrix& phi, const int& x, const int& t, const int& mu, const int& bet) {
    //time
    if (mu == 0) {
        if (t == Nt - 1) {
            return -phi[Coords[x][0]][bet];
            //return phi[Coords[x][0]][bet];
        }
        else {
            return phi[Coords[x][t + 1]][bet];
        }
    }
    else {
        //periodic
        return phi[Coords[mod(x + 1, Ns)][t]][bet];
    }
}

//left fermionic boundary (antiperiodic in time) x-hat{mu}
inline c_double lfb(const c_matrix& phi, const int& x, const int& t, const int& mu, const int& bet) {
    //time
    if (mu == 0) {
        if (t == 0) {
            return -phi[Coords[x][Nt - 1]][bet];
            //return phi[Coords[x][Nt - 1]][bet];
        }
        else {
            return phi[Coords[x][t - 1]][bet];
        }
    }
    else {
        //periodic
        return phi[Coords[mod(x - 1, Ns)][t]][bet];
    }
}


//Random U1 variable
c_double RandomU1() {
    double cociente = ((double)rand() / (RAND_MAX));
    double theta = 2.0 * pi * cociente;
    c_double z(cos(theta), sin(theta));
    return z;
}


//D phi
c_matrix D_phi(const c_matrix& U, const c_matrix& phi, const double& m0) {
    int Ntot = Ns * Nt;
    c_matrix Dphi(Ntot,c_vector(2, 0)); //Dphi[Ntot][2]
    for (int x = 0; x < Ns; x++) {
        for (int t = 0; t < Nt; t++) {
            int n = Coords[x][t];
            for (int alf = 0; alf < 2; alf++) {
                Dphi[n][alf] = (m0 + 2) * phi[n][alf];
                for (int bet = 0; bet < 2; bet++) {
                    for (int mu = 0; mu < 2; mu++) {
                        Dphi[n][alf] += -0.5 * (
                            (Identity[alf][bet] - gamma_mat[mu][alf][bet]) * U[n][mu] * rfb(phi, x, t, mu, bet)
                            + (Identity[alf][bet] + gamma_mat[mu][alf][bet]) * std::conj(U[LeftPB[x][t][mu]][mu]) * lfb(phi, x, t, mu, bet)
                            );
                    }
                }
            }
        }
    }
    return Dphi;
}

c_matrix canonical_vector(const int& i, const int& N1, const int& N2) {
    c_matrix e_i(N1, c_vector(N2, 0.0));
    int j = i / N2;
    int k = i % N2;
    e_i[j][k] = 1.0;
    return e_i;
}

//For writing the Dirac matrix
void save_matrix(c_matrix& Matrix, char* Name) {
    char NameData[500], Data_str[500];
    sprintf(NameData, Name);
    std::ofstream Datfile;
    Datfile.open(NameData);
    for (int i = 0; i < Matrix.size(); i++) {
        for (int j = 0; j < Matrix[i].size(); j++) {
            sprintf(Data_str, "%-30d%-30d%-30.17g%-30.17g\n", i, j, std::real(Matrix[i][j]), std::imag(Matrix[i][j]));
            Datfile << Data_str;
        }
    }
    Datfile.close();
}

//For writing the right hand side
void save_vector(c_matrix& Matrix, char* Name) {
    char NameData[500], Data_str[500];
    sprintf(NameData, Name);
    std::ofstream Datfile;
    Datfile.open(NameData);
    for (int i = 0; i < Matrix.size(); i++) {
        sprintf(Data_str, "%-30d%-30.17g%-30.17g\n", 2*i, std::real(Matrix[i][0]), std::imag(Matrix[i][0]));
        Datfile << Data_str;
        sprintf(Data_str, "%-30d%-30.17g%-30.17g\n", 2*i+1, std::real(Matrix[i][1]), std::imag(Matrix[i][1]));
        Datfile << Data_str;
    }
    
    Datfile.close();
}


//Formats decimal numbers
std::string format(const double& number) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << number;
    std::string str = oss.str();
    str.erase(str.find('.'), 1); //Removes decimal dot
    return str;
}


int main() {
    srand(time(0));//srand(0); //srand(time(0));
    initialize_matrices(); //Initialize gamma matrices, identity and unit vectors
    Coordinates(); //Vectorized coordinates
    periodic_boundary(); //Builds LeftPB and RightPB (periodic boundary for U_mu(n))
    std::cout << "Ns = " << Ns << " Nt = " << Nt << std::endl;
    std::cout << "Dirac matrix dimension = " << (2 * Ntot) << " X " << (2 * Ntot) << std::endl;
    //These parameter have to be same as in the simulation
    double m0 = -0.81; 
    double beta = 2;
    std::cout << "m0 " << m0 << " beta " << beta << std::endl;
    int n;// = 1; //Configuration number (this has to match the data file name)
    
    char PATH[500];
    //Path to the configuration

    for(n=0; n<1; n++) {
    std::cout << "Configuration number  " << n << std::endl;
    sprintf(PATH, "b2_16x16/m-016/2D_U1_Ns%d_Nt%d_b%s_m%s_%d.txt", Ns, Nt, format(beta).c_str(), format(-0.16).c_str(), n);
    std::cout << "Reading configuration from " << PATH << std::endl;
    std::ifstream infile(PATH);
    if (!infile) {
        std::cerr << "File not found" << std::endl;
        return 1;
    }
    int x, t, mu;
    double re, im;
    while (infile >> x >> t >> mu >> re >> im) {
        Conf[Coords[x][t]][mu] = c_double(re, im);
    }
    infile.close();

    //Random right hand side//
    c_matrix Phi(Ntot, c_vector(2, 0));
    for (int i = 0; i < Ntot; i++) {
        for (int j = 0; j < 2; j++) {
            Phi[i][j] = 1.0 * RandomU1();
        }
    }

    //Store the matrix//
    c_matrix D(2 * Ntot, c_vector(2 * Ntot, 0));
    for (int col = 0; col < 2 * Ntot; col++) {
        c_matrix v = canonical_vector(col, Ntot, 2);  
        c_matrix Dv = D_phi(Conf, v, m0); //column
        int count = 0;
        for (int i = 0; i < Dv.size(); i++) {
            for (int j = 0; j < Dv[i].size(); j++) {
                D[count][col] = Dv[i][j];
                count += 1;
            }
        }
    }
    char NameD[500], NamePhi[500];
    //Save right hand side and dirac matrix
    sprintf(NameD, "D%d.dat", n);  //sprintf(NamePhi, "Phi%d.dat", n);
    save_matrix(D, NameD);
    //save_vector(Phi, NamePhi);
    }
    return 0;
}