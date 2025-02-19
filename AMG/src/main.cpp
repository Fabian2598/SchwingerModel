#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include "gauge_conf.h"


int main() {
    std::cout << "Ns = " << Ns <<  " Nt = " << Nt << std::endl;
    std::cout << "block_x = " << block_x << " block_t = " << block_t << std::endl;
    std::cout << "x_elements = " << x_elements << " t_elements = " << t_elements << std::endl;
    GaugeConf GConf = GaugeConf(Ns,Nt);
    Coordinates(); //Vectorized coordinates
    Aggregates(); //Aggregates
    PrintAggregates();


	return 0;
}
