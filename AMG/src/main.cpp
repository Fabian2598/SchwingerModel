#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include "amg.h"


int main() {
    srand(time(0));
    initialize_matrices(); //Initialize gamma matrices, identity and unit vectors
    Coordinates(); //Vectorized coordinates
    periodic_boundary(); //Builds LeftPB and RightPB (periodic boundary for U_mu(n))
    Aggregates(); //Aggregates
    
    std::cout << "Ns = " << Ns <<  " Nt = " << Nt << std::endl;
    std::cout << "block_x = " << block_x << " block_t = " << block_t << std::endl;
    std::cout << "x_elements = " << x_elements << " t_elements = " << t_elements << std::endl;
    std::cout << "Number of test vectors = " << Ntest << std::endl;
    GaugeConf GConf = GaugeConf(Ns,Nt);
    GConf.initialization(); //Initialize the gauge configuration
    PrintAggregates();
    double m0 = 1;
    AMG amg = AMG(GConf,Ns,Nt,Ntest,m0);
    amg.tv_init(0.5); //test vectors intialization
	//Print P e_i where e_i are the basis vectors. The result should return the test vectors chopped over the aggregates.
    std::vector<std::vector<std::vector<std::complex<double>>>> tv = amg.test_vectors;
    for(int i =0; i<Ntest; i++){
        std::cout << "------tv" << i << "------" << std::endl;
        PrintVector(tv[i]);
    }
    int count = 0;
    for(int i = 0; i < Ntest; i++){
        for(int j =0; j<block_x*block_t; j++) {
            
            std::cout << "------e_" << count << "------" << std::endl;
            std::vector<std::vector<std::complex<double>>> e_i(Ntest, std::vector<std::complex<double> >(block_x*block_t,0.0));
            e_i[i][j] = 1.0;
            std::vector<std::vector<std::complex<double>>> Pe = amg.P_v(e_i);
            
            PrintVector(Pe);
            count += 1;
        }
    }


	return 0;
}
