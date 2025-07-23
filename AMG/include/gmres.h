#ifndef GMRES_H_INCLUDED
#define GMRES_H_INCLUDED

#include "dirac_operator.h"
#include "operator_overloads.h"
#include "variables.h"
#include <cmath>
#include <iostream>


/*
Implementation of GMRES for inverting D
    dim1: dimension of the first index of the spinor 
  	dim2: dimension of the second index of the spinor 	
        ->The solution is of the form x[dim1][dim2]
 	m0: mass parameter
	m: restart length
  	restarts: number of restarts of length m
 	tol: tolerance for the solver 
   
    The convergence criterion is ||r|| < ||phi|| * tol
*/
class GMRES{
    public:
    GMRES(const int& dim1, const int& dim2, const int& m, const int& restarts, const double& tol) : 
    dim1(dim1), dim2(dim2), m(m), restarts(restarts), tol(tol) {
        VmT = std::vector<c_matrix>(m+1, c_matrix(dim1, c_vector(dim2, 0))); 
        Hm = c_matrix(m+1 , c_vector(m, 0)); 
        gm = c_vector (m + 1, 0); 
        sn = c_vector (m, 0);
        cn = c_vector (m, 0);
        eta = c_vector (m, 0);
        r = spinor(dim1, c_vector(dim2, 0));
        w = spinor (dim1, c_vector(dim2, 0));
        Dx = spinor (dim1, c_vector(dim2, 0));  
    }	

	~GMRES() { };

    int gmres(const spinor&phi, const spinor& x0, spinor& x, const bool& print_message = false); //1 --> converged, 0 --> not converged

    private:

    const int dim1; //dimension of the first index of the spinor
    const int dim2; //dimension of the second index of the spinor
    const int m; //Restart length
    const int restarts; //Number of restarts
    const double tol; //Tolerance for the solver

    spinor r;  //r[coordinate][spin] residual
    //VmT[column vector index][vector arrange in matrix form]
    std::vector<c_matrix> VmT; //V matrix transpose-->dimensions exchanged
    c_matrix Hm; //H matrix (Hessenberg matrix)
    c_vector gm; 

    //Elements of rotation matrix |sn[i]|^2 + |cn[i]|^2 = 1
    c_vector sn;
    c_vector cn;
    //Solution to the triangular system
    c_vector eta;
    spinor w;
    spinor Dx; //auxiliary spinor
    c_double beta; //not 1/g^2 from simulations


    /*
    Matrix-vector operation
    This is defined in the derived classes
    */
    virtual void func(const spinor& in, spinor& out) = 0; 

    /*
    Rotations to transform Hessenberg matrix to upper triangular form
    cn: cosine components of the rotation
    sn: sine components of the rotation
    H: Hessenberg matrix
    j: index of the column being processed
    */
    void rotation(const int& j); 

    /*
    Solves an upper triangular system Ax = b, where A is an upper triangular matrix of dimension n
    A: upper triangular matrix
    b: right-hand side vector
    n: dimension of the matrix
    out: output vector where the solution will be stored
    */
    void solve_upper_triangular(const c_matrix& A, const c_vector& b, const int& n, c_vector& out);

    void setZeros(){
        //We set all these to zero to avoid memory issues
        for(int i = 0; i < m + 1; i++) {
            gm[i] = 0.0; //gm vector
            for(int j = 0; j < dim1; j++) {
                for(int k = 0; k < dim2; k++) {
                    VmT[i][j][k] = 0.0;
                }
            }
            for(int j = 0; j < m; j++) {
                Hm[i][j] = 0.0; //Hm matrix
            }
        }
        for(int i = 0; i < m; i++) {
            sn[i] = 0.0; //sn vector
            cn[i] = 0.0; //cn vector
            eta[i] = 0.0; //eta vector
        }
        //r[0][0] = 0.0;
        for(int i = 0; i < dim1; i++) {
            for(int j = 0; j < dim2; j++) {
                r[i][j] = 0.0; //r vector
                w[i][j] = 0.0; //w vector
                Dx[i][j] = 0.0; //Dx vector
            }
        }
    }
};


/*
    GMRES for the fine level
    This class inherits from GMRES and implements the matrix-vector product for the fine level
*/
class GMRES_fine_level : public GMRES {

public:
    GMRES_fine_level(const int& dim1, const int& dim2, const int& m, const int& restarts, const double& tol,
    const c_matrix& U, const double& m0) : GMRES(dim1, dim2, m, restarts, tol), U(U), m0(m0){
    };
    ~GMRES_fine_level() { };
    
private:
    const c_matrix& U; //reference to Gauge configuration. This is to avoid copying the matrix
    const double& m0; //reference to mass parameter
    /*
    Implementation of the function that computes the matrix-vector product for the fine level
    */
    void func(const spinor& in, spinor& out) override {
        D_phi(U, in, out, m0);
    }
};

//PREVIOUS IMPLEMENTATION//

/*
    Implementation of GMRES for inverting D
    func: matrix-vector operation 
    dim1: dimension of the first index of the spinor 
  	dim2: dimension of the second index of the spinor 	
        ->The solution is of the form x[dim1][dim2]
    U: gauge configuration
 	phi: right hand side
 	x0: initial guess
 	m0: mass parameter
	m: restart length
  	restarts: number of restarts of length m
 	tol: tolerance for the solver 
 	print_message: if true, print the convergence message
        
    The convergence criterion is ||r|| < ||phi|| * tol
*/
spinor gmres(spinor (*func)(const c_matrix&, const c_matrix&, const double&),const int& dim1, const int& dim2,
const c_matrix& U,const spinor& phi, const spinor& x0, const double& m0,const int& m, const int& restarts, const double& tol, const bool& print_message);


/*
    Rotations to transform Hessenberg matrix to upper triangular form
    cn: cosine components of the rotation
    sn: sine components of the rotation
    H: Hessenberg matrix
    j: index of the column being processed
*/
void rotation(c_vector& cn, c_vector& sn, c_matrix& H, const int& j);

/*
    Solves an upper triangular system Ax = b, where A is an upper triangular matrix of dimension n
    A: upper triangular matrix
    b: right-hand side vector
    n: dimension of the matrix
    out: output vector where the solution will be stored
*/
void solve_upper_triangular(const c_matrix& A, const c_vector& b, const int& n, c_vector& out);


#endif 