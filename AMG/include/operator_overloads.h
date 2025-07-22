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


/*
    dot product between two complex vectors 
*/
inline c_double dot(const c_vector& x, const c_vector& y) {
    c_double z = 0;
    for (int i = 0; i < x.size(); i++) {
        z += x[i] * std::conj(y[i]);
    }
    return z;
}


/*
    dot product between two spinors of the form psi[ntot][2]
    A.B = sum_i A_i conj(B_i) 
*/
inline c_double dot(const spinor& x, const spinor& y) {
    c_double z = 0;
    for (int i = 0; i < x.size(); i++) {
        for (int j = 0; j < x[i].size(); j++) {
            z += x[i][j] * std::conj(y[i][j]);
        }
    }
    return z;
}

/*
    Scalar times a complex vector
*/
template <typename T>
inline c_vector operator*(const T& lambda, const c_vector& A) {
    c_vector B(A.size(), 0);
    for (int i = 0; i < A.size(); i++) {
        B[i] = lambda * A[i];
    }
    return B;
}
/*
    Complex vector addition
*/
inline c_vector operator+(const c_vector& A, const c_vector& B) {
    c_vector C(A.size(), 0);
    for (int i = 0; i < A.size(); i++) {
        C[i] = A[i] + B[i];
    }
    return C;
}

/*
    Complex vector subtraction
*/
inline c_vector operator-(const c_vector& A, const c_vector& B) {
    c_vector C(A.size(), 0);
    for (int i = 0; i < A.size(); i++) {
        C[i] = A[i] - B[i];
    }
    return C;
}

/*
    Matrix-vector multiplication 
*/
inline c_vector operator*(const c_matrix& A, const c_vector& v) {
    c_vector w(A.size(), 0);
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            w[i] += A[i][j] * v[j];
        }
    }
    return w;
}


/*
    Scalar times a complex matrix.
    Also works for spinors, since they are just matrices with 2 columns.
*/
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

/*
    Complex matrix addition
    Also works for spinors, since they are just matrices with 2 columns.
*/
inline c_matrix operator+(const c_matrix& A, const c_matrix& B) {
    c_matrix C(A.size(), c_vector(A[0].size(), 0));
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    return C;
}

/*
    Complex matrix subtraction
    Also works for spinors, since they are just matrices with 2 columns.
*/
inline c_matrix operator-(const c_matrix& A, const c_matrix& B) {
    c_matrix C(A.size(), c_vector(A[0].size(), 0));
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

/*
    Print a complex matrix
    Also works for spinors
*/
inline void PrintComplexMatrix(const c_matrix& v ){
    for(int i = 0; i < v.size(); i++){
        for(int j = 0; j < v[i].size(); j++){
            std::cout << v[i][j] << " ";
        }
		std::cout << std::endl;
    }
    std::cout << std::endl;
}

/*
    Print a complex vector
*/
inline void PrintComplexVector(const c_vector& v ){
    for(int i = 0; i < v.size(); i++){
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;
}

#endif 


/*
 //Operator overloads --> They create too many copies and are not efficient
    //Scalar times a complex vector

template <typename T>
inline c_vector operator*(const T& lambda, const c_vector& A) {
    c_vector B(A.size(), 0);
    for (int i = 0; i < A.size(); i++) {
        B[i] = lambda * A[i];
    }
    return B;
}

    //Complex vector addition

inline c_vector operator+(const c_vector& A, const c_vector& B) {
    c_vector C(A.size(), 0);
    for (int i = 0; i < A.size(); i++) {
        C[i] = A[i] + B[i];
    }
    return C;
}


    //Complex vector subtraction

inline c_vector operator-(const c_vector& A, const c_vector& B) {
    c_vector C(A.size(), 0);
    for (int i = 0; i < A.size(); i++) {
        C[i] = A[i] - B[i];
    }
    return C;
}


    //Matrix-vector multiplication 

inline c_vector operator*(const c_matrix& A, const c_vector& v) {
    c_vector w(A.size(), 0);
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            w[i] += A[i][j] * v[j];
        }
    }
    return w;
}



    //Scalar times a complex matrix.
    //Also works for spinors, since they are just matrices with 2 columns.

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


    //Complex matrix addition
    //Also works for spinors, since they are just matrices with 2 columns.

inline c_matrix operator+(const c_matrix& A, const c_matrix& B) {
    c_matrix C(A.size(), c_vector(A[0].size(), 0));
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    return C;
}


    //Complex matrix subtraction
    //Also works for spinors, since they are just matrices with 2 columns.

inline c_matrix operator-(const c_matrix& A, const c_matrix& B) {
    c_matrix C(A.size(), c_vector(A[0].size(), 0));
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}


*/