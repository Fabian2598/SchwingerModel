#include "gauge_conf.h"



//Random U1 variable
c_double RandomU1() {
	double cociente = ((double) rand() / (RAND_MAX));
    double theta = 2.0*pi * cociente;
	c_double z(cos(theta), sin(theta));
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
	for (int x = 0; x < Nx; x++) {
		for (int t = 0; t < Nt; t++) {
			//int Coord0 = Coords[x][t], Coord1 = Coords[x][modulo(t + 1, Nt)], Coord2 = Coords[modulo(x + 1, Ns)][t];
			int Coord0 = Coords[x][t], Coord1 = RightPB[x][t][0], Coord2 = RightPB[x][t][1];
			Plaquette01[Coord0] = Conf[Coord0][0] * Conf[Coord1][1] * std::conj(Conf[Coord2][0]) * std::conj(Conf[Coord0][1]);
		}
	}
}

//Compute staple at coordinate (x,t) in the mu-direction
void GaugeConf::Compute_Staple() {
    // WARNING: Some references define the staple as the conjugate of this:
    //U_v(x) U_m(x+v) U*_v(x+m) + U*_v(x-v) U_m(x-v) U_v(x+m-v)
    //mu = 0 time direction, mu = 1 space direction
    for (int x = 0; x < Nx; x++) {
        for (int t = 0; t < Nt; t++) {
            //These coordinates could change depending on the conventions 
			int x1 = RightPB[x][t][1]; //Coords[modulo(x + 1, Ns) ,t]
			int x_1 = LeftPB[x][t][1]; //Coords[modulo(x - 1, Ns) ,t]
			int t1 = RightPB[x][t][0]; //Coords[x, modulo(t + 1, Nt)]
			int t_1 = LeftPB[x][t][0]; //Coords[x, modulo(t - 1, Nt)]
            int i = Coords[x][t];
            for (int mu = 0; mu < 2; mu++) {
                if (mu == 0) {
                    const c_double& conf1 = Conf[i][1];
                    const c_double& conf2 = Conf[x1][0];
                    const c_double& conf3 = Conf[t1][1];
                    const c_double& conf4 = Conf[x_1][1];
                    const c_double& conf5 = Conf[x_1][0];
                    const c_double& conf6 = Conf[ x_1_t1[x][t] ][1]; //Coords[mod(x - 1, Ns)][mod(t + 1, Nt)];
                    Staples[i][mu] = conf1 * conf2 * std::conj(conf3) +
                        std::conj(conf4) * conf5 * conf6;
                }
                else {
                    const c_double& conf1 = Conf[i][0];
                    const c_double& conf2 = Conf[t1][1];
                    const c_double& conf3 = Conf[x1][0];
                    const c_double& conf4 = Conf[t_1][0];
                    const c_double& conf5 = Conf[t_1][1];
                    const c_double& conf6 = Conf[ x1_t_1[x][t] ][0]; //Coords[mod(x + 1, Ns)][mod(t - 1, Nt)];
                    Staples[i][mu] = conf1 * conf2 * std::conj(conf3) +
                        std::conj(conf4) * conf5 * conf6;
                }
            }
        }
    }
}

//Save Gauge configuration //I could make this a member function ...
void SaveConf(c_matrix& Conf,char* Name) {
	char NameData[500], Data_str[500];
	sprintf(NameData, Name);
	std::ofstream Datfile;
	Datfile.open(NameData);
	using namespace LV;
	for (int x = 0; x < Nx; x++) {
		for (int t = 0; t < Nt; t++) {
			int i = x * Nx + t;
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


double GaugeConf::Compute_gaugeAction(const double& beta) {
	double action = 0;
	for (int i = 0; i < Ntot; i++) {
        action += beta * std::real(1.0-Plaquette01[i]);
	}
	return action;
}
