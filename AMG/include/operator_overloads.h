#ifndef OPERATOR_OVERLOADS_H
#define OPERATOR_OVERLOADS_H
#include <vector>
#include <complex>

typedef std::complex<double> c_double;
typedef std::vector<c_double> c_vector;
typedef std::vector<c_vector> c_matrix; 

//Complex dot product 
inline c_double dot(const c_vector& x, const c_vector& y) {
    c_double z = 0;
    for (int i = 0; i < x.size(); i++) {
        z += x[i] * std::conj(y[i]);
    }
    return z;
}

//Complex dot product between vectors arranged like matrices (not matrix multiplication)
// A.B = sum_i A_i conj(B_i) 
inline c_double dot(const c_matrix& x, const c_matrix& y) {
    c_double z = 0;
    for (int i = 0; i < x.size(); i++) {
        for (int j = 0; j < x[i].size(); j++) {
            z += x[i][j] * std::conj(y[i][j]);
        }
    }
    return z;
}

//Dot product between two vectors arranged like matrices (not matrix multiplication)
// A.B = sum_i sum_j A_ij B_ij (not conjugate)
inline c_double dot_v2(const c_matrix& x, const c_matrix& y) {
    c_double z = 0;
    for (int i = 0; i < x.size(); i++) {
        for (int j = 0; j < x[i].size(); j++) {
            z += x[i][j] * y[i][j];
        }
    }
    return z;
}

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
inline c_matrix operator-(const c_matrix& A, const c_matrix& B) {
    c_matrix C(A.size(), c_vector(A[0].size(), 0));
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

#endif 
