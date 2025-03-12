#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "amg.h"
#include "statistics.h"

//Formats decimal numbers
static std::string format(const double& number) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << number;
    std::string str = oss.str();
    str.erase(str.find('.'), 1); //Removes decimal dot
    return str;
}


int main() {
    srand(0);//srand(time(0));
    initialize_matrices(); //Initialize gamma matrices, identity and unit vectors
    Coordinates(); //Vectorized coordinates
    periodic_boundary(); //Builds LeftPB and RightPB (periodic boundary for U_mu(n))
    Aggregates(); //Aggregates
    std::cout << "Ns = " << Ns << " Nt = " << Nt << std::endl;
    std::cout << "block_x = " << block_x << " block_t = " << block_t << std::endl;
    std::cout << "x_elements = " << x_elements << " t_elements = " << t_elements << std::endl;
    std::cout << "Number of test vectors = " << Ntest << std::endl;
	std::cout << "Dirac matrix dimension = " << (2*Ntot) << std::endl;
    std::cout << "Dc dimension = " << Ntest*Nagg << std::endl;

    GaugeConf GConf = GaugeConf(Ns, Nt);
    GConf.initialization(); //Initialize the gauge configuration
    PrintAggregates();
    double m0 = -1;

    AMG amg = AMG(GConf, Ns, Nt, Ntest, m0);
    amg.tv_init(1, 3); //test vectors intialization
    for (int i = 0; i < Ntest; i++) {
		std::cout << "------tv" << i << "------" << std::endl;
        PrintVector(amg.test_vectors[i]);
    }
	
	c_matrix P(2*Ntot,c_vector(Ntest*Nagg,0));
    for (int col = 0; col < Ntest*Nagg; col++) {
        c_matrix v = amg.P_v(canonical_vector(col, Ntest, Nagg));  //column
        int count = 0;
        for (int i = 0; i < v.size(); i++) {
            for (int j = 0; j < v[i].size(); j++) {
                P[count][col] = v[i][j];
                count += 1;
            }
        }
    }
	std::cout << "------P------" << std::endl;
	//PrintVector(P);
	save_matrix(P, "P.dat");
    c_matrix Pt(Ntest*Nagg, c_vector(2*Ntot,0));
    for (int col = 0; col < 2 * Ntot; col++) {
        c_matrix v = amg.Pt_v(canonical_vector(col, Ntot, 2));  //column
        int count = 0;
        for (int i = 0; i < v.size(); i++) {
            for (int j = 0; j < v[i].size(); j++) {
                Pt[count][col] = v[i][j];
                count += 1;
            }
        }
    }
    std::cout << "------P^T------" << std::endl;
    //PrintVector(Pt);
	save_matrix(Pt, "Pt.dat");

	c_matrix D(2*Ntot, c_vector(2 * Ntot, 0));
    for (int col = 0; col < 2 * Ntot; col++) {
		c_matrix v = canonical_vector(col, Ntot, 2);  //column
		c_matrix Dv = D_phi(GConf.Conf, v, m0);
		int count = 0;
		for (int i = 0; i < Dv.size(); i++) {
			for (int j = 0; j < Dv[i].size(); j++) {
				D[count][col] = Dv[i][j];
				count += 1;
			}
		}
    }
    save_matrix(D, "D.dat");

    c_matrix Dc(Ntest * Nagg, c_vector(Ntest * Nagg, 0));
    for (int col = 0; col < Ntest * Nagg; col++) {
        c_matrix v = amg.Pt_D_P(canonical_vector(col, Ntest, Nagg)); //column
        int count = 0;
        for (int i = 0; i < v.size(); i++) {
            for (int j = 0; j < v[i].size(); j++) {
                Dc[count][col] = v[i][j];
                count += 1;
            }
        }
    }
    save_matrix(Dc, "Dc.dat");
	PrintVector(Dc);
   
    return 0;
}
