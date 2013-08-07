/*
Draws final configurations for the Z^2 sandpile, over a range of initial
confiigurations.  In each case, the initial configuration consists of
some number of chips placed at the origin.

The output is a directory sandpile/ which is assumed to exist prior to
calling the program.

Requires libpng.

To compile:
g++ sandpile-draw.cpp png_draw.cpp -lpng
 */

#include <iostream>
#include <deque>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "png_draw.h"

using namespace std;

unsigned int** grid;

int distance(int x1, int y1, int x2, int y2)
{
  return abs(x1-x2) + abs(y1-y2);
}

void graph(int N)
{
  grid = new unsigned int*[2000];
  
  for(int i=0; i<2000; i++)
    grid[i] = new unsigned int[2000];
  
  for(int i=0; i<2000; i++)
    for(int j=0; j<2000; j++)
      grid[i][j] = 0;
      
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
  }
  
  char fpath[100];
  sprintf(fpath, "sandpile_images/spile%d.png", N);
  
  int width = 800, height = 800;
  
  pngd_image* img = pngd_make_image(fpath, width+1, height+1);
  
  int color[5][3];
  color[0][0] =   0; color[0][1] = 159; color[0][2] = 107;
  color[1][0] =   0; color[1][1] = 135; color[1][2] = 189;
  color[2][0] = 255; color[2][1] = 211; color[2][2] =   0;
  color[3][0] = 196; color[3][1] =   2; color[3][2] =  51;

  
  int cid;
  for(int x = or_x-width/2; x<=or_x+width/2; x++)
  {
    for(int y = or_y-height/2; y<=or_y+height/2; y++)
    {
      cid = grid[x][y];
      pngd_draw_pixel(img, x-(or_x - width/2), y-(or_y - height/2), color[cid][0], color[cid][1], color[cid][2]);
    }
  }
  pngd_finalise(img);
  
  for(int i=0; i<2000; i++)
    delete[] grid[i];
  delete[] grid;
}

int main()
{
  for(int N=100; N<100000; N+=1000)
    graph(N);

  return 0;
}
