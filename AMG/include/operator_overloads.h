#ifndef OPERATOR_OVERLOADS_H
#define OPERATOR_OVERLOADS_H
#include <vector>
#include <complex>
//#include <cblas.h>
//#include <omp.h>

typedef std::complex<double> c_double;
typedef std::vector<c_double> c_vector;
typedef std::vector<c_vector> c_matrix; 
typedef std::vector<c_vector> spinor;


//Complex dot product between two vectors 
inline c_double dot(const c_vector& x, const c_vector& y) {
    c_double z = 0;
    for (int i = 0; i < x.size(); i++) {
        z += x[i] * std::conj(y[i]);
    }
    return z;
}


//Complex dot product between spinors psi[ntot][2]
// A.B = sum_i A_i conj(B_i) 
inline c_double dot(const spinor& x, const spinor& y) {
    c_double z = 0;
    for (int i = 0; i < x.size(); i++) {
        for (int j = 0; j < x[i].size(); j++) {
            z += x[i][j] * std::conj(y[i][j]);
        }
    }
    return z;
}


//I have to check if this is actually faster than the above version.
/*
inline c_double dot(const spinor& x, const spinor& y) {
    // Flatten both spinors to 1D arrays (column-major order)
    int rows = x.size();
    int cols = x[0].size();
    std::vector<c_double> flat_x(rows * cols), flat_y(rows * cols);
    for (int j = 0; j < cols; ++j)
        for (int i = 0; i < rows; ++i) {
            flat_x[j * rows + i] = x[i][j];
            flat_y[j * rows + i] = y[i][j];
        }

    c_double result;
    cblas_zdotc_sub(rows * cols,
        reinterpret_cast<const void*>(flat_y.data()), 1,
        reinterpret_cast<const void*>(flat_x.data()), 1,
        reinterpret_cast<void*>(&result)
    );
    return result;
}
*/

//Scalar times complex vector
template <typename T>
inline c_vector operator*(const T& lambda, const c_vector& A) {
    c_vector B(A.size(), 0);
    for (int i = 0; i < A.size(); i++) {
        B[i] = lambda * A[i];
    }
    return B;
}
//Vector addition
inline c_vector operator+(const c_vector& A, const c_vector& B) {
    c_vector C(A.size(), 0);
    for (int i = 0; i < A.size(); i++) {
        C[i] = A[i] + B[i];
    }
    return C;
}

//Vector subtraction
inline c_vector operator-(const c_vector& A, const c_vector& B) {
    c_vector C(A.size(), 0);
    for (int i = 0; i < A.size(); i++) {
        C[i] = A[i] - B[i];
    }
    return C;
}

//Matrix vector multiplication
inline c_vector operator*(const c_matrix& A, const c_vector& v) {
    c_vector w(A.size(), 0);
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            w[i] += A[i][j] * v[j];
        }
    }
    return w;
}


//scalar multiplication of a matrix
template <typename T>
inline c_matrix operator*(const T& lambda, const c_matrix& A) {
    c_matrix B(A.size(), c_vector(A[0].size(), 0));
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            B[i][j] = lambda * A[i][j];
        }
    }
    return B;
}

//Matrix addition
//Also works for spinors
inline c_matrix operator+(const c_matrix& A, const c_matrix& B) {
    c_matrix C(A.size(), c_vector(A[0].size(), 0));
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    return C;
}

//Matrix subtraction
//Also works for spinors
inline c_matrix operator-(const c_matrix& A, const c_matrix& B) {
    c_matrix C(A.size(), c_vector(A[0].size(), 0));
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

//This also works for spinors.
inline void PrintComplexMatrix(const c_matrix& v ){
    for(int i = 0; i < v.size(); i++){
        for(int j = 0; j < v[i].size(); j++){
            std::cout << v[i][j] << " ";
        }
		std::cout << std::endl;
    }
    std::cout << std::endl;
}

inline void PrintComplexVector(const c_vector& v ){
	//for c_vector
    for(int i = 0; i < v.size(); i++){
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;
}

#endif 
