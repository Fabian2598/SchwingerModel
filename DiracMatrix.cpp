#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <complex>  
#include <vector>

/*
    Given a gauge configuration in binary format, this program assembles the Dirac matrix for the Schwinger model.
    Only the non-zero elements of the Dirac matrix are stored in a file.
    The file format is:
    row_index col_index real_part imag_part

    The path to the configuration is specified in the code.
    The values of m0 and beta are read from the user. They have to coincide with the values used in the simulation.
*/

//-------Lattice dimensions-------//
constexpr int Nx = 256, Nt = 256;
constexpr int Ntot = Nx * Nt;
//--------------------------------//

double pi = 3.14159265359;
double m0, beta;
typedef std::complex<double> c_double;
typedef std::vector<c_double> c_vector;
typedef std::vector<c_vector> c_matrix;
typedef std::vector<c_vector> spinor;


c_matrix Conf(Ntot, c_vector(2, 0)); //Gauge configuration
std::vector<std::vector<int>>Coords(Nx, std::vector<int>(Nt, 0)); //Vectorized coordinates

//---Neighbor coordinates---//
std::vector<std::vector<int>>LeftPB = std::vector<std::vector<int>>(Ntot, std::vector<int>(2,0)); 
std::vector<std::vector<int>>RightPB = std::vector<std::vector<int>>(Ntot, std::vector<int>(2,0)); 
std::vector<std::vector<c_double>>SignL =std::vector<std::vector<c_double>>(Ntot,std::vector<c_double>(2,0)); 
std::vector<std::vector<c_double>>SignR = std::vector<std::vector<c_double>>(Ntot,std::vector<c_double>(2,0)); 

c_double I_number(0, 1); //imaginary number


/*
    Random U(1) number --> For storing a random right hand side vector
*/
c_double RandomU1() {
    double cociente = ((double)rand() / (RAND_MAX));
    double theta = 2.0 * pi * cociente;
    c_double z(cos(theta), sin(theta));
    return z;
}

/*
    Vectorized coordinate
*/
void Coordinates() {
    for (int x = 0; x < Nx; x++) {
        for (int t = 0; t < Nt; t++) {
            Coords[x][t] = x * Nt + t;
        }
    }
}

/*
	Modulo operation
*/
inline int mod(int a, int b) {
	int r = a % b;
	return r < 0 ? r + b : r;
}

/*
	Periodic boundary conditions used for the link variables U_mu(n).
	This function builds the RightPB and LeftPB, which
	store the neighbor coordinates.
	This prevents recalculation every time we call the operator D.
	The function is only called once at the beginning of the program.

	right periodic boundary x+hat{mu}
	left periodic boundary x-hat{mu}
	hat_mu[0] = { 1, 0 } --> hat_t
	hat_mu[1] = { 0, 1 } --> hat_x
*/
inline void periodic_boundary() {
	std::vector<std::vector<int>>hat_mu(2, std::vector<int>(2, 0));
	hat_mu[0] = { 1, 0 }; //hat_t
	hat_mu[1] = { 0, 1 }; //hat_x
	for (int x = 0; x < Nx; x++) {
		for (int t = 0; t < Nt; t++) {
			for (int mu = 0; mu < 2; mu++) {
				RightPB[Coords[x][t]][mu] = Coords[mod(x + hat_mu[mu][1], Nx)][mod(t + hat_mu[mu][0], Nt)]; 
				LeftPB[Coords[x][t]][mu] = Coords[mod(x - hat_mu[mu][1], Nx)][mod(t - hat_mu[mu][0], Nt)];
				SignR[Coords[x][t]][mu] = (mu == 0 && t == Nt - 1) ? -1 : 1; //sign for the right boundary in time
				SignL[Coords[x][t]][mu] = (mu == 0 && t == 0) ? -1 : 1; //sign for the left boundary in time
			}
		}
	}
}

/*
    Dirac operator
    U: gauge configuration
    phi: spinor field
    Dphi: result of the Dirac operator applied to phi
    m0: mass parameter
*/
void D_phi(const c_matrix& U, const spinor& phi, spinor &Dphi,const double& m0) {
	for (int n = 0; n < Ntot; n++) {
		//n = x * Nt + t
		Dphi[n][0] = (m0 + 2) * phi[n][0] -0.5 * ( 
			U[n][0] * SignR[n][0] * (phi[RightPB[n][0]][0] - phi[RightPB[n][0]][1]) 
		+   U[n][1] * SignR[n][1] * (phi[RightPB[n][1]][0] + I_number * phi[RightPB[n][1]][1])
		+   std::conj(U[LeftPB[n][0]][0]) * SignL[n][0] * (phi[LeftPB[n][0]][0] + phi[LeftPB[n][0]][1])
		+   std::conj(U[LeftPB[n][1]][1]) * SignL[n][1] * (phi[LeftPB[n][1]][0] - I_number * phi[LeftPB[n][1]][1])
		);

		Dphi[n][1] = (m0 + 2) * phi[n][1] -0.5 * ( 
			U[n][0] * SignR[n][0] * (-phi[RightPB[n][0]][0] + phi[RightPB[n][0]][1]) 
		+   U[n][1] * SignR[n][1] * (-I_number*phi[RightPB[n][1]][0] + phi[RightPB[n][1]][1])
		+   std::conj(U[LeftPB[n][0]][0]) * SignL[n][0] * (phi[LeftPB[n][0]][0] + phi[LeftPB[n][0]][1])
		+   std::conj(U[LeftPB[n][1]][1]) * SignL[n][1] * (I_number*phi[LeftPB[n][1]][0] + phi[LeftPB[n][1]][1])
		);
			
	}
	
}


void read_nonbinary(c_matrix& GlobalConf, const std::string& name){
    std::ifstream infile(name);
    if (!infile) {
        std::cerr << "File not found" << std::endl;
        exit(1);
    }
    int x, t, mu;
    double re, im;

    while (infile >> x >> t >> mu >> re >> im) {
        GlobalConf[Coords[x][t]][mu] = c_double(re, im);
    }
    infile.close();
}

void read_binary(c_matrix& GlobalConf, const std::string& name){
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
            GlobalConf[n][mu] = c_double(re, im); //n --> for testing      
    }
    }
    }
    infile.close();
    //std::cout << "Binary conf read from " << name << std::endl;     
}


/*
    Canonical vector e_i used for extracting the columns of the Dirac matrix.
*/
spinor canonical_vector(const int& i, const int& N1, const int& N2) {
    spinor e_i(N1, c_vector(N2, 0.0));
    int j = i / N2;
    int k = i % N2;
    e_i[j][k] = 1.0;
    return e_i;
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
    srand(time(0));
    Coordinates(); //Vectorized coordinates
    periodic_boundary(); //Builds LeftPB and RightPB (periodic boundary for U_mu(n))
    std::string PATH;
    int n = 0; //Configuration number

    std::cout << "Assembling Dirac matrix for the Schwinger model" << std::endl;
    std::cout << "Nx = " << Nx << " Nt = " << Nt << std::endl;
    std::cout << "Dirac matrix dimension = " << (2 * Ntot) << " X " << (2 * Ntot) << std::endl;

    std::cout << "m0 (the same as in the simulation) = ";
    std::cin >> m0; //Mass parameter
    std::cout << "beta (the same as in the simulation) = ";
    std::cin >> beta; //Beta parameter
    std::cout << "Conf. path: ";
    std::cin >> PATH; //Beta parameter
    std::cout << "Conf ID: ";
    std::cin >> n;
    std::cout << " " << std::endl;

    std::cout << "m0 = " << m0 << " beta = " << beta << std::endl;
    std::cout << "Conf. path " << PATH << std::endl;


    //------Path to the configuration--------//
    /*n= 0;
    PATH = "confs/b2_" + std::to_string(Nx) +
                   "x" + std::to_string(Nt) + "/m-018/2D_U1_Ns" + std::to_string(Nx) +
                   "_Nt" + std::to_string(Nt) +
                   "_b" + format(beta) +
                   "_m" + format(m0) +
                   "_" + std::to_string(n) + ".ctxt";
    */
    //--------------------------------------//
   
    
    std::cout << "Reading configuration from " << PATH << std::endl;
    read_nonbinary(Conf,PATH);
    //read_binary(Conf,PATH);

    //------------Store the matrix---------------//
    std::string Name = "DiracMatrix_" + std::to_string(Nx) +
                   "x" + std::to_string(Nt) +
                   "_b" + format(beta) +
                   "_m" + format(m0) +
                   "_" + std::to_string(n) + ".bin";
    std::ofstream Datfile(Name,std::ios::binary);
    if (!Datfile.is_open()) {
        std::cerr << "Error opening file: " << Name << std::endl;
        return 0;
    }
    
    spinor Dv(Ntot, c_vector(2, 0)); //Result of the Dirac operator applied to the canonical vector, i.e. a column of D
    spinor v(Ntot, c_vector(2, 0)); //Canonical vector
    int row, col;
    std::cout << "Writing Dirac matrix " << Name << std::endl;
    int nonzero = 0;
    //----Loop over the columns of the Dirac matrix----//
    //The Dirac matrix is 2*Ntot x 2*Ntot, where
    for (col = 0; col < 2 * Ntot; col++) {
        v = canonical_vector(col, Ntot, 2);  
        D_phi(Conf, v, Dv, m0); //Extracting the column
        for (int n = 0; n < Ntot; n++) {
            for (int mu = 0; mu < 2; mu++) {
                row = 2 * n + mu; //Row index in the Dirac matrix
                //We only store the non-zero elements of the Dirac matrix
                if ( std::abs(std::real(Dv[n][mu])) > 1e-10 || std::abs(std::imag(Dv[n][mu])) > 1e-10 ){  
                    nonzero += 1;
                    const double& re = std::real(Dv[n][mu]);
                    const double& im = std::imag(Dv[n][mu]);          
                    Datfile.write(reinterpret_cast<const char*>(&row), sizeof(int));
                    Datfile.write(reinterpret_cast<const char*>(&col), sizeof(int));
                    Datfile.write(reinterpret_cast<const char*>(&re), sizeof(double));
                    Datfile.write(reinterpret_cast<const char*>(&im), sizeof(double));
                }       
            }
        }
    }
    std::cout << "Number of non-zero entries " << nonzero << std::endl;

    return 0;
}