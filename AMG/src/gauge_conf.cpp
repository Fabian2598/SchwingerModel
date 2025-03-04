#include "gauge_conf.h"

//Vectorized coordinates
void Coordinates() {
	for (int x = 0; x < Ns; x++) {
		for (int t = 0; t < Nt; t++) {
			Coords[x][t] = x * Ns + t;
			//Coords[t][x] = t * Nt + x;
		}
	}
}

//Aggregates A_j = L_j x {0,1}
void Aggregates(){
	for (int x = 0; x < block_x; x++) {
		for (int t = 0; t < block_t; t++) {
			int x0 = x * x_elements, t0 = t * t_elements;
			int x1 = (x + 1) * x_elements, t1 = (t + 1) * t_elements;
			int aggregate = x * block_x + t;
			int count = 0;
			//x and t are redefined in the following loop
			for (int x = x0; x < x1; x++) {
				for (int t = t0; t < t1; t++) {
					int i = x * Ns + t;
					Agg[aggregate][count] = i;
					count++;
				}
			}
			if(count != x_elements*t_elements) {
				std::cout << "Aggregate " << aggregate << " has " << count << " elements" << std::endl;
			}
			//Once the loops are finished count should be x_elements*t_elements
		}
	}
}

//Print aggregates. Useful for debugging
void PrintAggregates() {
	for (int i = 0; i < block_x * block_t; i++) {
		std::cout << "-------Aggregate-----" << i << std::endl;
		for (int j = 0; j < x_elements * t_elements; j++) {
			std::cout << Agg[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

//Aggregates A_j_0 = L_j x {0}, A_j_1 = L_j x {1}
void AggregatesV2() {
	for (int x = 0; x < block_x; x++) {
		for (int t = 0; t < block_t; t++) {
			int x0 = x * x_elements, t0 = t * t_elements;
			int x1 = (x + 1) * x_elements, t1 = (t + 1) * t_elements;
			int aggregate = x * block_x + t;
			int count = 0;
			//x and t are redefined in the following loop
			for (int x = x0; x < x1; x++) {
				for (int t = t0; t < t1; t++) {
					int i = x * Ns + t;
					Agg[aggregate][count] = i;
					count++;
				}
			}
			if (count != x_elements * t_elements) {
				std::cout << "Aggregate " << aggregate << " has " << count << " elements" << std::endl;
			}
			//Once the loops are finished count should be x_elements*t_elements
		}
	}
}

//Print aggregates. Useful for debugging
void PrintAggregatesV2() {
	for (int i = 0; i < block_x * block_t; i++) {
		std::cout << "-------Aggregate-----" << i << std::endl;
		for (int j = 0; j < x_elements * t_elements; j++) {
			std::cout << Agg[i][j] << " ";
		}
		std::cout << std::endl;
	}
}



//Random U1 variable
std::complex<double> RandomU1() {
	double cociente = ((double) rand() / (RAND_MAX));
    double theta = 2.0*pi * cociente;
	std::complex<double> z(cos(theta), sin(theta));
	return z;
}

//Initialize a random conf
void GaugeConf::initialization() {
	for (int i = 0; i < Ntot; i++) {
		for (int mu = 0; mu < 2; mu++) {
			Conf[i][mu] = RandomU1(); //Conf[Ns x Nt][mu in {0,1}]
		}
	}
}

void GaugeConf::Compute_Plaquette01() {
	//U_mv(x) = U_m(x) U_v(x+m) U*_m(x+v) U*_v(x)
	//mu = 0 time direction, mu = 1 space direction
	int Ntot = Ns * Nt;
	for (int x = 0; x < Ns; x++) {
		for (int t = 0; t < Nt; t++) {
			//int Coord0 = Coords[x][t], Coord1 = Coords[x][modulo(t - 1, Nt)], Coord2 = Coords[modulo(x - 1, Ns)][t];
			int Coord0 = Coords[x][t], Coord1 = LeftPB[x][t][0], Coord2 = LeftPB[x][t][1];
			Plaquette01[Coord0] = Conf[Coord0][0] * Conf[Coord1][1] * std::conj(Conf[Coord2][0]) * std::conj(Conf[Coord0][1]);
		}
	}
}

//Compute staple at coordinate (x,t) in the mu-direction
void GaugeConf::Compute_Staple() {
    // WARNING: Some references define the staple as the conjugate of this:
    //U_v(x) U_m(x+v) U*_v(x+m) + U*_v(x-v) U_m(x-v) U_v(x+m-v)
    //mu = 0 time direction, mu = 1 space direction
    for (int x = 0; x < Ns; x++) {
        for (int t = 0; t < Nt; t++) {
            //These coordinates could change depending on the conventions 
			int x1 = LeftPB[x][t][1]; //Coords[modulo(x - 1, Ns) ,t]
			int x_1 = RightPB[x][t][1]; //Coords[modulo(x + 1, Ns) ,t]
			int t1 = LeftPB[x][t][0]; //Coords[x, modulo(t - 1, Nt)]
			int t_1 = RightPB[x][t][0]; //Coords[x, modulo(t + 1, Nt)]
            int i = Coords[x][t];
            for (int mu = 0; mu < 2; mu++) {
                if (mu == 0) {
                    const std::complex<double>& conf1 = Conf[i][1];
                    const std::complex<double>& conf2 = Conf[x1][0];
                    const std::complex<double>& conf3 = Conf[t1][1];
                    const std::complex<double>& conf4 = Conf[x_1][1];
                    const std::complex<double>& conf5 = Conf[x_1][0];
                    const std::complex<double>& conf6 = Conf[ x_1_t1[x][t] ][1]; //Coords[mod(x + 1, Ns)][mod(t - 1, Nt)];
                    Staples[i][mu] = conf1 * conf2 * std::conj(conf3) +
                        std::conj(conf4) * conf5 * conf6;
                }
                else {
                    const std::complex<double>& conf1 = Conf[i][0];
                    const std::complex<double>& conf2 = Conf[t1][1];
                    const std::complex<double>& conf3 = Conf[x1][0];
                    const std::complex<double>& conf4 = Conf[t_1][0];
                    const std::complex<double>& conf5 = Conf[t_1][1];
                    const std::complex<double>& conf6 = Conf[ x1_t_1[x][t] ][0]; //Coords[mod(x - 1, Ns)][mod(t + 1, Nt)];
                    Staples[i][mu] = conf1 * conf2 * std::conj(conf3) +
                        std::conj(conf4) * conf5 * conf6;
                }
            }
        }
    }
}

//Save Gauge configuration
void SaveConf(std::vector<std::vector<std::complex<double>>>& Conf,char* Name) {
	char NameData[500], Data_str[500];
	sprintf(NameData, Name);
	std::ofstream Datfile;
	Datfile.open(NameData);
	for (int x = 0; x < Ns; x++) {
		for (int t = 0; t < Nt; t++) {
			int i = x * Ns + t;
			for (int mu = 0; mu < 2; mu++) {
				sprintf(Data_str, "%-30d%-30d%-30d%-30.17g%-30.17g\n", x, t, mu, std::real(Conf[i][mu]), std::imag(Conf[i][mu]));
				Datfile << Data_str;
			}
		}
	}
	Datfile.close();
}



double GaugeConf::MeasureSp_HMC() {
	//Plaquettes have to be computed at the HMC update
	double Sp = 0.0;
	for (int i = 0; i < Ntot; i++) {
		Sp += std::real(Plaquette01[i]);	
	}
	return Sp;
}


