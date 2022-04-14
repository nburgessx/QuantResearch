/*  A program to calculate the Lennard-Jones Forces between a set of water molecules */

#include <stdio.h>
#include <math.h>
#include <iostream>

#include <vector>
#include <chrono>
#include <random>

/* constants for water molecule */
#define sigma 0.3165555 /* nanometers */
#define epsilon 0.6501696 /* kJ/mol */

/* numbers of molecules */
#define n 100000 


using namespace std;
int main() {

  double dx, dy, dz, r, mag, t1, t2;
  
  vector<double>  x (n,0.0);
  vector<double>  y (n,0.0);
  vector<double>  z (n,0.0);
  vector<double>  fx(n,0.0);
  vector<double>  fy(n,0.0);
  vector<double>  fz(n,0.0);


  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);
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

#pragma omp parallel for schedule(dynamic) proc_bind(spread) 
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {

      if (i != j) {
	/* r_j - r_i  vector components */ 
	dx = x[j] - x[i];
	dy = y[j] - y[i];
	dz = z[j] - z[i];

	/* calculate distrance retween particles i and j */  
	r = (dx*dx + dy*dy + dz*dz);
	r = sqrt(r);

	/* force magnitude */
	mag = -24.*epsilon*((2.*pow(sigma,12)/pow(r,13)) - (pow(sigma,6)/pow(r,7)));
	/* force components */
	fx[i] = fx[i] + mag*dx;
	fy[i] = fy[i] + mag*dy;
	fz[i] = fz[i] + mag*dz;
      }
      
    }
  }

  auto finish = std::chrono::steady_clock::now();
  cout << "TIME TO COMPUTE LJ FORCES " << std::chrono::duration_cast<std::chrono::duration<double> >(finish-start).count() << "s\n";
  cout << "Gflops = " << 60.*1.e-9*n*n/std::chrono::duration_cast<std::chrono::duration<double> >(finish-start).count() << "\n";

  
}
