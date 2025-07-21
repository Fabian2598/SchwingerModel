#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <complex>  
#include <vector>

/*
    Given a gauge configuration, this program assembles the Dirac matrix for the Schwinger model.
    Only the non-zero elements of the Dirac matrix are stored in a file.
    The file format is:
    row_index col_index real_part imag_part

    The path to the configuration is specified in the code.
    Lattice dimensions are set to 16x16.
    The values of m0 and beta are read from the user. They have to coincide with the values used in the simulation.
*/

//-------Lattice dimensions-------//
constexpr int Nx = 16, Nt = 16;
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

/*
    Saves a spinor vector to a file.
    The file format is:
    row_index real_part imag_part
    where row_index = 2*i for the first component and 2*i+1 for the second component.
*/
void save_vector(const spinor& psi, const std::string& Name) {
    std::ofstream Datfile(Name);
    if (!Datfile.is_open()) {
        std::cerr << "Error opening file: " << Name << std::endl;
        return;
    }
    for (int i = 0; i < Ntot; i++) {
        Datfile << 2*i
                << std::setw(30) << std::setprecision(17) << std::real(psi[i][0])
                << std::setw(30) << std::setprecision(17) << std::imag(psi[i][0])
                << "\n";

        Datfile << 2*i+1
                << std::setw(30) << std::setprecision(17) << std::real(psi[i][1])
                << std::setw(30) << std::setprecision(17) << std::imag(psi[i][1])
                << "\n";
    }

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

    std::cout << "Assembling Dirac matrix for the Schwinger model" << std::endl;
    std::cout << "Nx = " << Nx << " Nt = " << Nt << std::endl;
    std::cout << "Dirac matrix dimension = " << (2 * Ntot) << " X " << (2 * Ntot) << std::endl;

    std::cout << "m0 (the same as in the simulation) = ";
    std::cin >> m0; //Mass parameter
    std::cout << "beta (the same as in the simulation) = ";
    std::cin >> beta; //Beta parameter
    std::cout << " " << std::endl;
    std::cout << "m0 = " << m0 << " beta = " << beta << std::endl;


    //------Path to the configuration--------//
    int n = 0; //Configuration number
    std::string PATH = "2D_U1_Ns" + std::to_string(Nx) +
                   "_Nt" + std::to_string(Nt) +
                   "_b" + format(beta) +
                   "_m" + format(m0) +
                   "_" + std::to_string(n) + ".ctxt";
    //--------------------------------------//
   
    
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

    //--------Random right hand side--------//
    spinor Phi(Ntot, c_vector(2, 0));
    for (int i = 0; i < Ntot; i++) {
        for (int j = 0; j < 2; j++) {
            Phi[i][j] = 1.0 * RandomU1();
        }
    }
    std::string rhs_path = "RHS.txt";
    save_vector(Phi, rhs_path); //Save the random right hand side vector
    std::cout << "Random right hand side vector saved in " << rhs_path << std::endl;
    //--------------------------------------//

    //------------Store the matrix---------------//
    std::string Name = "DiracMatrix.txt";
    std::ofstream Datfile(Name);
    if (!Datfile.is_open()) {
        std::cerr << "Error opening file: " << Name << std::endl;
        return 0;
    }
    
    spinor Dv(Ntot, c_vector(2, 0)); //Result of the Dirac operator applied to the canonical vector, i.e. a column of D
    spinor v(Ntot, c_vector(2, 0)); //Canonical vector
    int row, col;
    std::cout << "Writing Dirac matrix " << Name << std::endl;

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
                    Datfile << row
                    << std::setw(30) << col
                    << std::setw(30) << std::setprecision(17) << std::real(Dv[n][mu])
                    << std::setw(30) << std::setprecision(17) << std::imag(Dv[n][mu])
                    << "\n";
                }       
            }
        }
    }

    return 0;
}