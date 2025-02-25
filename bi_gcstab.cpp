#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <ctime>

//This program solves a linear problem using the conjugate gradient method.
//I just have to modify it to solve the Dirac equation and connect it with the rest of the Schwinger model.

typedef std::complex<double> c_double;
typedef std::vector<c_double> c_vector;
typedef std::vector<c_vector> c_matrix;

c_vector A_d(const c_matrix& A,const c_vector& x){
    //Matrix-vector multiplication
    c_vector y(A.size(),0);
    for(int i=0;i<A.size();i++){
        for(int j=0;j<A[i].size();j++){
            y[i] += A[i][j]*x[j];
        }
    }
    return y;
}

c_double dot(const c_vector& x,const c_vector& y){
    //Dot product of two vectors
    c_double z = 0;
    for(int i=0;i<x.size();i++){
        z += x[i]*std::conj(y[i]); //Complex dot product
    }
    return z;
}

//Overload + - and * operators
template <typename T>
c_vector operator*(const T& v1, const c_vector& v2){
    c_vector y(v2.size(),0);
    for(int i=0;i<v2.size();i++){
        y[i] = v1*v2[i];
    }
    return y;
}

c_vector operator+(const c_vector& v1, const c_vector& v2){
    c_vector y(v1.size(),0);
    for(int i=0;i<v1.size();i++){
        y[i] = v1[i] + v2[i];
    }
    return y;
}

c_vector operator-(const c_vector& v1, const c_vector& v2){
    c_vector y(v1.size(),0);
    for(int i=0;i<v1.size();i++){
        y[i] = v1[i] - v2[i];
    }
    return y;
}

int main (){
    //srand(time(0));
    srand(0);
    int N = 8;
    c_matrix A(N,c_vector(N,0));
    c_vector b(N,0);
    for (int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            A[i][j] = c_double(1.0*(std::rand() % 10), 1.0 * (std::rand() % 10));
        }
    }

    for (int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){        
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "--------------------" << std::endl;  
    std::cout << "b = " << std::endl;  
    for (int i = 0; i<N; i++){
        b[i] = c_double(1.0*(std::rand() % 10), 1.0 * (std::rand() % 10));
        std::cout << b[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "--------------------" << std::endl;  
    c_vector x(N,0);
    for (int i = 0; i<N; i++){x[i] = 0.0;}
    c_vector r(N,0);
    c_vector r_tilde(N, 0);
    c_vector d(N,0);
    c_vector Ad(N,0);
	c_vector s(N, 0);
	c_vector t(N, 0);
    c_double alpha,beta, rho_i, omega, rho_i_2;
    int k = 0;
    int max_iter = 100;
    double tol = 1e-12;
    double err = 1;
    r = b - A_d(A,x);
    r_tilde = r;
    while(k<max_iter && err>tol){
        rho_i = dot(r, r_tilde); //r . r_dagger
		if (k == 0) {
			d = r;
		}
        else {
            beta = alpha * rho_i / (omega * rho_i_2); //c_double
			d = r + beta * (d - omega * Ad); //c_vector
        }
        Ad = A_d(A,d); //c_vector
        alpha = rho_i/dot(Ad,r_tilde); //c_double
		s = r - alpha * Ad; //c_vector
        err = std::real(dot(s, s));
		if (err < tol){
			x = x + alpha * d;
			break;
		}
		t = A_d(A, s); //c_vector
		omega = dot(s,t) / dot(t, t); //c_double
		r = s - omega * t; //c_vector
        x = x + alpha*d + omega*s; //c_vector
		rho_i_2 = rho_i;
        std::cout << "Error: " << err << std::endl;
        k++;
    }


    std::cout << "Number of iterations: " << k << std::endl;
    std::cout << "Final Error: " << err << std::endl;
    std::cout << "Solution: ";
    for(int i = 0; i<N; i++){
        std::cout << x[i] << " ";
    }
    std::cout << std::endl;

    x = A_d(A, x);
        //this should print b
    std::cout << "A A^-1 b = b: ";
    for (int i = 0; i < N; i++) {
        std::cout << x[i] << " ";
    }
    std::cout << std::endl;
    

return 0;
}

