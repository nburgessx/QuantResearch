/*  A program to calculate the Lennard-Jones Forces between a set of water molecules */

#include <stdio.h>
#include <math.h>
#include <iostream>

#include <vector>
#include <chrono>
#include <random>

/* constants for water molecule */
#define sigma 0.3165555f /* nanometers */
#define epsilon 0.6501696f /* kJ/mol */

/* numbers of molecules */
#define n 10000000 


using namespace std;
int main() {

  float dx, dy, dz, r, mag, t1, t2;
  
  vector<float>  x (n,0.0);
  vector<float>  y (n,0.0);
  vector<float>  z (n,0.0);
  vector<float>  fx(n,0.0);
  vector<float>  fy(n,0.0);
  vector<float>  fz(n,0.0);


  std::default_random_engine generator;
  std::uniform_real_distribution<float> distribution(0.0,1.0);
#pragma omp parallel for  
  for (int i = 0; i < n; i++) {
    /* x, y, z are positions of water molecules in space (nanometers) */ 
    x[i] = 10.0*distribution(generator);
    y[i] = 10.0*distribution(generator);
    z[i] = 10.0*distribution(generator);
    fx[i] = 0.0;
    fy[i] = 0.0;
    fz[i] = 0.0;
  }

  auto start = std::chrono::steady_clock::now();
  
#define BLOCK 16
#define MIN(x,y) ((x)<(y)?(x):(y))
  
#pragma omp parallel for schedule(dynamic) proc_bind(spread) private(dx,dy,dz,r,mag)
  for (int ii = 0; ii < n; ii+=BLOCK) {
    int iend = MIN(ii+BLOCK,n);
    for (int j = 0; j < n; j++) {
#pragma omp simd
      for (int i=ii;i<iend;i++) {
	
	if (i != j) {
	  /* r_j - r_i  vector components */ 
	  dx = x[j] - x[i];
	  dy = y[j] - y[i];
	  dz = z[j] - z[i];
	  
	  /* calculate distrance retween particles i and j */  
	  r = (dx*dx + dy*dy + dz*dz);
	  r = 1.f/sqrtf(r);
	  
	  /* force magnitude */
	  mag = -24.f*epsilon*((2.f*powf(sigma,12.f)*powf(r,13.f)) - (powf(sigma,6.f)*powf(r,7.f)));
	  /* force components */
	  fx[i] = fx[i] + mag*dx;
	  fy[i] = fy[i] + mag*dy;
	  fz[i] = fz[i] + mag*dy;
	}
	
      }
    }
  }
    
  auto finish = std::chrono::steady_clock::now();
  cout << "TIME TO COMPUTE LJ FORCES " << std::chrono::duration_cast<std::chrono::duration<double> >(finish-start).count() << "s\n";
  cout << "Gflops = " << 60.*1.e-9*n*n/std::chrono::duration_cast<std::chrono::duration<double> >(finish-start).count() << "\n";

  
}
