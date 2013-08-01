#include "green.hpp"

using namespace arma;
using namespace std;

void heat_kernel(mat& H, double t, mat A )
{
  if(A.n_rows != A.n_cols)
    return;
  int N = A.n_rows;

  vec eigval;
  mat eigvec;
  
  eig_sym(eigval, eigvec, A);

  H.set_size(N,N);
  H.fill(0);
  
  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
      for(int k=0; k<N; k++)
        H(i,j) += exp(-t * eigval[k]) * eigvec(i,k) * eigvec(j,k);
}

void green_function(mat& G, mat A)
{
  if(A.n_rows != A.n_cols)
    return;
  int N = A.n_rows;

  vec eigval;
  mat eigvec;
  
  eig_sym(eigval, eigvec, A);
    
  G.set_size(N,N);
  G.fill(0);
  
  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
      for(int k=0; k<N; k++)
        G(i,j) += (1. /  eigval[k]) * eigvec(i,k) * eigvec(j,k);
}

void green_function_alpha(mat& G, mat A, double alpha)
{
  if(A.n_rows != A.n_cols)
    return;
  int N = A.n_rows;

  vec eigval;
  mat eigvec;
  
  eig_sym(eigval, eigvec, A);
    
  G.set_size(N,N);
  G.fill(0);
  
  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
      for(int k=0; k<N; k++)
        G(i,j) += (1. /  (eigval[k] + alpha)) * eigvec(i,k) * eigvec(j,k);
}

int cheby_coeff[500][500];

void init_cheby()
{
  static bool init = false;

  if(init)
    return;

  for(int i=0; i<500; i++)
    for(int j=0; j<500; j++)
      cheby_coeff[i][j] = 0;
  cheby_coeff[0][0] = 1;
  cheby_coeff[1][1] = 2;
  for(int i=2; i<500; i++)
  {
    cheby_coeff[i][0] = - cheby_coeff[i-2][0];
    for(int j=1; j<500; j++)
      cheby_coeff[i][j] = 2*cheby_coeff[i-1][j-1] - cheby_coeff[i-2][j];
  } 
  
  init = true;
}

double chebyshev_2(int n, double x)
{
  double ret=0.;
  for(int i=0; i<=n; i++)
    ret += cheby_coeff[n][i] * pow(x,(double) i);
  return ret;
}

double lattice_green_val(int m, int n, int x1, int x2, int y1, int y2)
{
  if(x1 > y1)
  {
    int z1 = x1, z2 = x2;
    x1 = y1;
    x2 = y2;
    y1 = z1;
    y2 = z2;
  }
  
  double result = 0.;
  for(int i=1; i<=n; i++)
  {
    double c = M_PI * (double) i / (double) (n+1);
    double tmp = 8. * sin(c * x2) * sin(c*y2);
    tmp *= chebyshev_2(x1-1, 2.-cos(c));
    tmp *= chebyshev_2(m-y1, 2.-cos(c));
    tmp *= (i%2) ? 1 : -1;
    
    tmp /= (double) (n+1);
    tmp /= chebyshev_2(m, 2.-cos(c));
    
    result += tmp;
  }
  
  return result;
}

void lattice_green(mat& G, int m, int n)
{
  G.set_size(m*n,m*n);
 
  for(int x1=0; x1<m; x1++)
    for(int x2=0; x2<n; x2++)
      for(int y1=0; y1<m; y1++)
        for(int y2=0; y2<n; y2++)
        {
          //this ordering is strange, but seems to be what is required to get
          //the matrix to look right... TODO justify this
          int i = m-1-x2 + x1*n;
          int j = y2 + y1*n;
          if(i > j)
            continue;
          G(i,j) = lattice_green_val(m,n,x1+1,x2+1,y1+1,y2+1);
          G(j,i) = G(i,j); 
        }
}

void laplace_path(mat& M, int n)
{
  M.set_size(n,n);
  for(int x=0; x<n; x++)
    for(int y=0; y<n; y++)
    {
      if(x==y)
        M(x,y) = 1.;
      else if(x == y-1  ||  x == y+1)
        M(x,y) = -.5;
      M(x,y) = 0.;
    }
}

void laplace_path_2d(mat& M, int m, int n)
{
  M.set_size(m*n,m*n);
 
  for(int x1=1; x1<=m; x1++)
    for(int x2=1; x2<=n; x2++)
      for(int y1=1; y1<=m; y1++)
        for(int y2=1; y2<=n; y2++)
        {
          int i = (x2-1) + (x1-1)*n;
          int j = (y2-1) + (y1-1)*n;
          if(i > j)
            continue;
           
          int d1 = (x1-y1+2*(m+1))%(m+1);
          int d2 = (x2-y2+2*(n+1))%(n+1);
            
          if(x1==y1 && x2==y2)
            M(i,j) = 1;
          else if((d1==1 || d1==m) && x2==y2)
            M(i,j) = -1. / 4.;
          else if((d2==1 || d2==n) && x1==y1)
            M(i,j) = -1. / 4.;
          else
            M(i,j) = 0.;
          M(j,i) = M(i,j); 
        }
}


/*int main()
{
  init_cheby();
 
  int N=20;
  int M=20;
  //pngd_image* img = pngd_make_image("green.png", N*M, N*M);

  double** g = new double*[N*M];
  for(int i=0; i<N*M; i++)
    g[i] = new double[N*M];

  double max = 0.;
  for(int i=0; i<N*M; i++)
  {
    for(int j=0; j<N*M; j++)
    {
      g[i][j] = green(N,N, (i/N)+1,(i%N)+1, (j/M)+1,(j%M)+1 );
      if(g[i][j] > max)
        max = g[i][j];
    }
    //cout << endl;
  }
  
  char name[100];
  
  for(int i=0; i<N*M; i++)
  {
    sprintf(name, "green/mat_%d.png", i);
    pngd_image* img = pngd_make_image(name,N,M);
    for(int x=0; x<N; x++)
      for(int y=0; y<M; y++)
      {
        int shade = (255. * (g[i][x+y*N] / max));
        pngd_draw_pixel(img,x,y,shade,shade,shade); 
      }
    pngd_finalise(img);
  }
  
  //for(int i=0; i<N*M; i++)
  //{
  //  for(int j=0; j<N*M; j++)
  //  {
  //    int shade = (255. * (g[i][j] / max));
  //    pngd_draw_pixel(img,i,j,shade,shade,shade);
  //  }
  //}
  //pngd_finalise(img);
  return 0;
}*/
