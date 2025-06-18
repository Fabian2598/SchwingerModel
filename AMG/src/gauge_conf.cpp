#include "gauge_conf.h"

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


void GaugeConf::PrintConf() {
	for (int x = 0; x < Nx; x++) {
		for (int t = 0; t < Nt; t++) {
			for (int mu = 0; mu < 2; mu++) {
				std::cout << "x " << x << " t " << t << " mu " << mu << "   " << Conf[Coords[x][t]][mu] << std::endl;
			}
		}
	}
}

//Save a two-dimensional complex vector to a file
void GaugeConf::SaveConf(char* Name) {
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

