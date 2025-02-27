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
    amg.tv_init(1,0); //test vectors intialization

    c_matrix phi(Ntot, c_vector(2, 0)); //For multiplying Pt phi
	for (int i = 0; i < Ntot; i++) {
		for (int j = 0; j < 2; j++) {
			phi[i][j] = 1.0 * RandomU1();
		}
	}
    std::cout << "------phi------" << std::endl;
    PrintVector(phi);
	c_matrix Pt(Nagg*Ntest, c_vector(2*Ntot, 0));
    c_matrix P(2*Ntot, c_vector(Nagg*Ntest, 0));

	for(int col = 0; col < 2*Ntot; col++){
		c_matrix v = amg.Pt_v(canonical_vector(col,Ntot,2)); //column
		int count = 0;
		for(int i = 0; i < v.size(); i++){
			for(int j = 0; j < v[i].size(); j++){
				Pt[count][col] = v[i][j];
				count += 1;
			}
		}	
	} 

    for(int col = 0; col < Ntest*Nagg; col++){
		c_matrix v = amg.P_v(canonical_vector(col,Ntest,Nagg)); //column
		int count = 0;
		for(int i = 0; i < v.size(); i++){
			for(int j = 0; j < v[i].size(); j++){
				P[count][col] = v[i][j];
				count += 1;
			}
		}	
	} 
	std::cout << "------P------" << std::endl;
	for(int i = 0; i < P.size(); i++){
		for(int j = 0; j < P[i].size(); j++){
			std::cout << P[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "------P^T------" << std::endl;
	for(int i = 0; i < Pt.size(); i++){
		for(int j = 0; j < Pt[i].size(); j++){
			std::cout << Pt[i][j] << " ";
		}
		std::cout << std::endl;
	}
	//Flatten phi
	c_vector phi_flat(2*Ntot, 0);
	int count = 0;
	for(int i = 0; i < phi.size(); i++){
		for(int j = 0; j < phi[i].size(); j++){
			phi_flat[count] = phi[i][j];
			count += 1;
		}
	}
	std::cout << "------phi_flat------" << std::endl;
	for(int i = 0; i < phi_flat.size(); i++){
		std::cout << phi_flat[i] << " ";
	}
	std::cout << std::endl;
    std::cout << "------P^T phi_flat------" << std::endl;
    PPrintVector(Pt*phi_flat);

    c_matrix Pt_phi = amg.Pt_v(phi);
    std::cout << "------P^T phi------" << std::endl;
    PrintVector(Pt_phi);
    
	return 0;
}

/*
    AMG amg = AMG(GConf,Ns,Nt,Ntest,m0);
    amg.tv_init(0.5); //test vectors intialization
	//Print P e_i where e_i are the basis vectors. The result should return the test vectors chopped over the aggregates.
    std::vector<std::vector<std::vector<std::complex<double>>>> tv = amg.test_vectors;
    for(int i =0; i<Ntest; i++){
        std::cout << "------tv" << i << "------" << std::endl;
        PrintVector(tv[i]);
    }
    for(int i = 0; i < Ntest*block_x*block_t; i++){
        std::cout << "------ P e_" << i << "------" << std::endl;
        std::vector<std::vector<std::complex<double>>> e_i = canonical_vector(i, Ntest,block_x*block_t);
        PrintVector(e_i);
        std::vector<std::vector<std::complex<double>>> Pe = amg.P_v(e_i); 
        PrintVector(Pe);    
    }
    for(int i = 0; i < 2*Ntot;i++){
        std::cout << "------P^T e_" << i << "------" << std::endl;
        std::vector<std::vector<std::complex<double>>> e_i = canonical_vector(i, Ntot, 2);
        //PrintVector(e_i);
        std::vector<std::vector<std::complex<double>>> Pe = amg.Pt_v(e_i); 
        PrintVector(Pe);
    }
    
*/