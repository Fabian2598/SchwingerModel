#include "gauge_conf.h"


c_double RandomU1() {
	//Random angle in (0,2*pi) with uniform distribution 
	double cociente = ((double) rand() / (RAND_MAX));
    double theta = 2.0*pi * cociente;
	c_double z(cos(theta), sin(theta));
	return z;
}

void GaugeConf::initialization() {
	for (int n = 0; n < Ntot; n++) {
		Conf.mu0[n] = RandomU1(); 
		Conf.mu1[n] = RandomU1(); 
	}
}

void GaugeConf::Compute_Plaquette01() {
	//U_mv(x) = U_m(x) U_v(x+m) U*_m(x+v) U*_v(x)
	//mu = 0 time direction, mu = 1 space direction
    #pragma omp parallel for
    for (int n = 0; n<Ntot; n++){
        //int Coord0 = Coords[x][t], Coord1 = Coords[x][modulo(t + 1, Nt)], Coord2 = Coords[modulo(x + 1, Ns)][t];
		Plaquette01[n] = Conf.mu0[n] * Conf.mu1[RightPB[2*n]] * std::conj(Conf.mu0[RightPB[2*n+1]]) * std::conj(Conf.mu1[n]);
    }		
}

//Compute staple at coordinate (x,t) in the mu-direction
void GaugeConf::Compute_Staple() {
    // WARNING: Some references define the staple as the conjugate of this:
    //U_v(x) U_m(x+v) U*_v(x+m) + U*_v(x-v) U_m(x-v) U_v(x+m-v)
    //mu = 0 time direction, mu = 1 space direction
    #pragma omp parallel for
    for (int n = 0; n < Ntot; n++) {
        //These coordinates could change depending on the conventions 
		int x1 = RightPB[2*n+1];  //Coords[modulo(x + 1, Ns) ,t]
		int x_1 = LeftPB[2*n+1];  //Coords[modulo(x - 1, Ns) ,t]
		int t1 = RightPB[2*n];  //Coords[x, modulo(t + 1, Nt)]
		int t_1 = LeftPB[2*n];  //Coords[x, modulo(t - 1, Nt)]
        for (int mu = 0; mu < 2; mu++) {
            if (mu == 0) {
                const c_double& conf1 = Conf.mu1[n];
                const c_double& conf2 = Conf.mu0[x1]; 
                const c_double& conf3 = Conf.mu1[t1];
                const c_double& conf4 = Conf.mu1[x_1];
                const c_double& conf5 = Conf.mu0[x_1];
                const c_double& conf6 = Conf.mu1[x_1_t1[n]]; //Coords[mod(x - 1, Ns)][mod(t + 1, Nt)];
                Staples.mu0[n] = conf1 * conf2 * std::conj(conf3) +
                    std::conj(conf4) * conf5 * conf6;
            }
            else {
                const c_double& conf1 = Conf.mu0[n];
                const c_double& conf2 = Conf.mu1[t1];
                const c_double& conf3 = Conf.mu0[x1];
                const c_double& conf4 = Conf.mu0[t_1];
                const c_double& conf5 = Conf.mu1[t_1];
                const c_double& conf6 = Conf.mu0[x1_t_1[n]]; //Coords[mod(x + 1, Ns)][mod(t - 1, Nt)];
                Staples.mu1[n] = conf1 * conf2 * std::conj(conf3) +
                    std::conj(conf4) * conf5 * conf6;
            }
        }
        
    }
}

/*
Save Gauge configuration
*/ 
void SaveConf(const GaugeConf& GConf, const std::string& Name) {
    std::ofstream Datfile(Name);
    if (!Datfile.is_open()) {
        std::cerr << "Error opening file: " << Name << std::endl;
        return;
    }
    using namespace LV;
    
    for (int x = 0; x < Nx; x++) {
        for (int t = 0; t < Nt; t++) {
            int n = x * Nx + t;
            for (int mu = 0; mu < 2; mu++) {
                const double& re = std::real(mu == 0 ? GConf.Conf.mu0[n] : GConf.Conf.mu1[n]);
                const double& im = std::imag(mu == 0 ? GConf.Conf.mu0[n] : GConf.Conf.mu1[n]);
                Datfile << x
                        << std::setw(30) << t
                        << std::setw(30) << mu
                        << std::setw(30) << std::setprecision(17) << std::scientific << re
                        << std::setw(30) << std::setprecision(17) << std::scientific << im
                        << "\n";
            }
        }
    }
}


double GaugeConf::MeasureSp_HMC() {
	//Plaquettes have to be computed during the HMC update
	double Sp = 0.0;
	for (int n = 0; n < Ntot; n++) {
		Sp += std::real(Plaquette01[n]);	
	}
	return Sp;
}


double GaugeConf::Compute_gaugeAction(const double& beta) {
	double action = 0;
    #pragma omp parallel for reduction(+:action)
	for (int n = 0; n < Ntot; n++) {
        action += beta * std::real(1.0-Plaquette01[n]);
	}
	return action;
}

void GaugeConf::read_conf(const std::string& name){
    std::ifstream infile(name);
    if (!infile) {
        std::cerr << "File " << name << " not found " << std::endl;
        exit(1);
    }
    int x, t, mu;
    double re, im; 
    while (infile >> x >> t >> mu >> re >> im) {
        if (mu == 0)
            Conf.mu0[x*Nt+t] = c_double(re, im); 
        else
            Conf.mu1[x*Nt+t] = c_double(re, im); 
    }
    infile.close();
    std::cout << "Conf read from " << name << std::endl;
    
}