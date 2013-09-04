#ifndef SANDPILE_HPP
#define SANDPILE_HPP

#include <armadillo>

/*
Computes the Z^2 sandpile with initially num_chips at the origin and 
zeros everywhere else.
*/

void sandpile_lattice(arma::mat& C, arma::mat& F, int num_chips, int dim_m, int dim_n);

/*
Computes the Z^2 sandpile with intial condition given by Ci.
*/

void sandpile_lattice(arma::mat& Ci, arma::mat& Ce, arma::mat& F, int dim_m, int dim_n);

/*
Computes the k-regular tree sandpile, for a given initial condition.
depth represents the radius of the regular tree.
*/

void regular_tree_lattice(arma::mat& Ci, arma::mat& Ce, arma::mat& F, int num_chips, int k, int depth);

#endif
