/*
Sanity check for the lattice_green code.  Computes the discrete Green's function
by lattice_green and two other methods, and checks that they agree.

Requires the armadillo library.

To compile:
g++ test.cpp green.cpp -larmadillo
*/

#include <iostream>

#include "green.hpp"

using namespace std;
using namespace arma;

double max_val(double a, double b, double c)
{
  a = fabs(a);
  b = fabs(b);
  c = fabs(c);
  
  if(a > b && a > c)
    return a;
  if(b > c)
    return b;
  return c;
}

int main()
{
  mat G1,G2,G3,P,X;
  
  init_cheby();
  lattice_green(G1,2,2);
  laplace_path_2d(P,2,2);
  green_function(G2,P);
  G3 = inv(P);
  
  mat D1 = G1 - G2;
  mat D2 = G2 - G3;
  mat D3 = G1 - G3;
  
  double max_pos = max_val(D1.max(), D2.max(), D3.max()); 
  double max_neg = max_val(D1.min(), D2.min(), D3.min());

  double max_diff = ((max_pos > max_neg) ? max_pos : max_neg);
  
  if(max_diff > 0.000001)
  {
    cout << "The green's functions did not agree." << endl;
  }
  else
  {
    cout << "Green's functions agree." << endl;
  }  
  return 0;
}
