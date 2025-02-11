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
    int N = 4;
    c_matrix A(N,c_vector(N,0));
    c_vector b(N,0);
    for (int i = 0; i<N; i++){
        for(int j = i+1; j<N; j++){
            A[i][j] = c_double(1.0*(std::rand() % 10), 1.0 * (std::rand() % 10));
            A[j][i] = std::conj(A[i][j]);
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
    c_vector d(N,0);
    c_vector Ad(N,0);
    c_double alpha,beta;
    int k = 0;
    int max_iter = 100;
    double tol = 1e-10;
    double err = 1;
    r = b - A_d(A,x);
    d = r;
    c_double r_norm2 = dot(r,r);
    while(k<max_iter && err>tol){
        Ad = A_d(A,d); //c_vector
        alpha = r_norm2/dot(d,Ad); //c_double
        x = x + alpha*d; //c_vector
        r = r - alpha*Ad; //c_vector
		err = std::real(dot(r, r)); //c_double
        beta = err/r_norm2; //c_double
        d = r + beta*d; //c_vector
        std::cout << "Error: " << err << std::endl;
		r_norm2 = err;
        k++;
    }
    std::cout << "Number of iterations: " << k << std::endl;
    std::cout << "Solution: ";
    for(int i = 0; i<N; i++){
        std::cout << x[i] << " ";
    }
    std::cout << std::endl;

    x = A_d(A, x);
        //this should print b
    for (int i = 0; i < N; i++) {
        std::cout << x[i] << " ";
    }
    std::cout << std::endl;
    

return 0;
}

