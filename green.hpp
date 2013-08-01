#ifndef GREEN_HPP
#define GREEN_HPP

#include <iostream>
#include <cmath>
#include <armadillo>

/* !!-- General heat kernel and Green's functions  --!! */

/*
Computes the heat kernel H, for a matrix A at time t.
*/
void heat_kernel(arma::mat& H, double t, arma::mat A);

/*
Computes the discrete Green's function G for A.  Note that
this code assumes that A is non-singular.
*/
void green_function(arma::mat& G, arma::mat A);

/*
Computes the variant of the Green's function, G_alpha, for A,
as defined in "Discrete Green's Functions", Chung&Yau.
*/
void green_function_alpha(arma::mat& G, arma::mat A, double alpha);

/* !!-- Chebyshev polynomial functions --!! */
void init_cheby();
double chebyshev_2(int,double);

/* !!-- Path & lattice functions --!! */

/*
Computes the discrete green's function of the lattice P_m x P_n
for entry (x1,x2),(y1,y2).
*/
double lattice_green_val(int m, int n, int x1, int x2, int y1, int y2);

/*
Computes the discrete green's function of the lattice P_m x P_n
as a matrix.
*/
void lattice_green(arma::mat& G,int m,int n);

/*
Computes the normalized Laplacian for the path P_n.
*/
void laplace_path(arma::mat& M,int n);

/*
Computes the normalized Laplacian for the lattice P_m x P_n.
*/
void laplace_path_2d(arma::mat& M,int m,int n);


#endif
