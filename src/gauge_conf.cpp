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


void GaugeConf::printPlaquette() {
	for (int i = 0; i < Ntot; i++) {
		std::cout << Plaquette01[i] << " ";
	}
}

/*
void GaugeConf::Compute_Plaquette01() {
	for (int x = 0; x < Ns; x++) {
		for (int t = 0; t < Nt; t++) {
			int Coord0 = Coords[x][t], Coord1 = Coords[modulo(x + 1, Ns)][t], Coord2 = Coords[x][modulo(t + 1, Nt)];
			Plaquette01[Coords[x][t]] = Conf[Coord0][0] * Conf[Coord1][1] * std::conj(Conf[Coord2][0]) * std::conj(Conf[Coord0][1]);
			Plaquette01_prime[Coords[x][t]] = Conf_copy[Coord0][0] * Conf_copy[Coord1][1] * std::conj(Conf_copy[Coord2][0]) * std::conj(Conf_copy[Coord0][1]);
		}
	}
}
*/

/*
double GaugeConf::MeasureSp_HMC() {
	//Plaquettes have to be computed at the HMC update
	double Sp = 0.0;
	for (int x = 0; x < Ns; x++) {
		for (int t = 0; t < Nt; t++) {
			Sp += std::real(Plaquette01[Coords[x][t]]);
		}
	}
	return Sp;

}
*/


//Print conf
/*
std::ostream& operator<<(std::ostream& out, const GaugeConf& GConf) {
	std::vector<std::vector<std::complex<double>>> Conf = GConf.getConf();
	for (int i = 0; i < GConf.getNs(); i++) {
		for (int j = 0; j < GConf.getNt(); j++) {
			int N0 = GConf.Coordinate(i, j);
			out << "[" << Conf[N0][0] << ", " << Conf[N0][1] << "] "; //Conf[Ns x Nt][mu in {0,1}]
		}
		out << "\n";
	}
	return out;
}
*/
