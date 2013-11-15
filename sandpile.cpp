#include <deque>

#include "sandpile.hpp"

using namespace std;
using namespace arma;

void sandpile_lattice(mat& Ci, mat& Ce, mat& F, int dim_m, int dim_n)
{
  unsigned int** grid   = new unsigned int*[2000];
  unsigned int** firing = new unsigned int*[2000]; 

  int or_x,or_y;
  or_x = or_y = 1000;
  int top_x = or_x - dim_m/2;
  int top_y = or_x - dim_n/2;

  for(int i=0; i<2000; i++)
  {
    grid[i]   = new unsigned int[2000];
    firing[i] = new unsigned int[2000];
  }

  for(int i=0; i<2000; i++)
    for(int j=0; j<2000; j++)
      grid[i][j] = firing[i][j] = 0;

  deque<int> topple;

  for(int x=0; x<dim_m; x++)
    for(int y=0; y<dim_n; y++)
    {
      grid[top_x+x][top_y+y] = Ci(x+y*dim_m);
      if(grid[top_x+x][top_y+y] >= 4)
      {
	topple.push_back(top_x+x);
	topple.push_back(top_y+y);
      }
    }
  
  while(topple.size() > 0)
  {
    int t_x = topple.front();
    topple.pop_front();
    int t_y = topple.front();
    topple.pop_front();
    
    int t_amount = grid[t_x][t_y] / 4;
    
    if(t_amount != 0)
    {
      for(int i=-1; i<=1; i++)
        for(int j=-1; j<=1; j++)
        {
          // do not topple diagonally
          if((i!=0 && j !=0) ||  (i==0 && j==0))
            continue;
          
          if(t_x+i < top_x || t_y+j < top_y || t_x + i >= dim_m+top_x || t_y + j >= dim_n+top_y)
          {
            //cout << "skipping: " << t_x + i << ' ' << t_y + j << endl;
            continue;
          }
          grid[t_x + i][t_y + j]+=t_amount;
          if(grid[t_x + i][t_y + j] >= 4)
          {
            // was already going to topple, don't add it to the queue again
            if(grid[t_x+i][t_y+j]-t_amount >= 4)
              continue;
            topple.push_back(t_x+i);
            topple.push_back(t_y+j); 
          }
        }
    } 
    grid[t_x][t_y] -= t_amount*4;
    firing[t_x][t_y] += t_amount;
  }
    
  Ce.set_size(dim_m*dim_n,1);
  F.set_size(dim_m*dim_n,1);

  for(int x=0; x<dim_m; x++)
    for(int y=0; y<dim_n; y++)
    {
      Ce(x+y*dim_m) = grid[top_x+x][top_y+y];
      F(x+y*dim_m) = firing[top_x+x][top_y+y];
    }

  for(int i=0; i<2000; i++)
  {
    delete[] grid[i];
    delete[] firing[i];
  }  
  delete[] grid;
  delete[] firing;
}

void sandpile_lattice(mat& C, mat& F, int N, int dim_m, int dim_n)
{
  unsigned int** grid   = new unsigned int*[2000];
  unsigned int** firing = new unsigned int*[2000]; 

  for(int i=0; i<2000; i++)
  {
    grid[i]   = new unsigned int[2000];
    firing[i] = new unsigned int[2000];
  }

  for(int i=0; i<2000; i++)
    for(int j=0; j<2000; j++)
      grid[i][j] = firing[i][j] = 0;
      
  deque<int> topple;
  
  int or_x = 1000, or_y = 1000;
  grid[or_x][or_y] = N;
  topple.push_back(or_x);
  topple.push_back(or_y);
  
  while(topple.size() > 0)
  {
    int t_x = topple.front();
    topple.pop_front();
    int t_y = topple.front();
    topple.pop_front();
    
    int t_amount = grid[t_x][t_y] / 4;
    
    if(t_amount != 0)
    {
      for(int i=-1; i<=1; i++)
        for(int j=-1; j<=1; j++)
        {
          // do not topple diagonally
          if((i!=0 && j !=0) ||  (i==0 && j==0))
            continue;
          
          grid[t_x + i][t_y + j]+=t_amount;
          if(grid[t_x + i][t_y + j] >= 4)
          {
            // was already going to topple, don't add it to the queue again
            if(grid[t_x+i][t_y+j]-t_amount >= 4)
              continue;
            topple.push_back(t_x+i);
            topple.push_back(t_y+j); 
          }
        }
    } 
    grid[t_x][t_y] -= t_amount*4;
    firing[t_x][t_y] += t_amount;
  }
    
  int top_x = or_x - dim_m/2;
  int top_y = or_x - dim_n/2;

  C.set_size(dim_m*dim_n,1);
  F.set_size(dim_m*dim_n,1);

  for(int x=0; x<dim_m; x++)
    for(int y=0; y<dim_n; y++)
    {
      C(x+y*dim_m) = grid[top_x+x][top_y+y];
      F(x+y*dim_m) = firing[top_x+x][top_y+y];
    }

  for(int i=0; i<2000; i++)
  {
    delete[] grid[i];
    delete[] firing[i];
  }  
  delete[] grid;
  delete[] firing;
}

void regular_tree_lattice(mat& Ci, mat& Ce, mat& F, int num_chips, int k, int depth)
{
  int size = k * (pow(k-1,depth) - 1.) / (k - 2.) + 1. + 0.0000001;
  Ce.set_size(size,1);
  F.set_size(size,1);

  Ce.fill(0);
  F.fill(0);

  deque<int> topple;
  Ce(0,0) = num_chips;
  topple.push_back(0);
  topple.push_back(0);

  while(topple.size())
  {
    int r,a,c;
    r = topple.front();
    topple.pop_front();
    a = topple.front();
    topple.pop_front();
    
    if(r==0)
      c = 0;
    else
      c = k * (pow(k-1,r-1) - 1.) / (k - 2.) + 1. + 0.0000001 + a;

    int n_amount = Ce(c,0) / k;

    cerr << "toppling: " << r << ' ' << a << ' ' << c << ' ' << n_amount << endl;

    for(int i=0; i<k; i++)
    {
      int r2,a2;

      if(c==0)
      {
        r2 = 1;
        a2 = i;
      }
      else
      {
        if(i==k-1)
        {
          r2 = r-1;
          a2 = a / (k-1) + 0.000001;
        }
        else
        {
          r2 = r+1;
          a2 = a * (k-1) + i;
        }
      }

      int c2 = k * (pow(k-1,r2-1) - 1.) / (k - 2.) + 1. + 0.0000001 + a2;

      Ce(c2,0) += n_amount;
      if(Ce(c2,0) > k - 0.0000001 && Ce(c2,0) - n_amount < k + 0.000001)
      {
        topple.push_back(r2);
        topple.push_back(a2);
      }
    }
    Ce(c,0) -= k * n_amount;
    F(c,0) += n_amount;
  }  
}
