#include <iostream>
#include <vector>
#include <chrono>


int main(int argc,char**argv){
  int n=500000000;
  std::vector<float> x(n,1.0);
  std::vector<float> y(n,1.0);
  std::vector<float> z(n,0.0);


  auto start = std::chrono::high_resolution_clock::now();
#pragma omp parallel for num_threads(4)
  for(int i=0;i<n;i++){
    z[i]=x[i]+y[i];
  }
  auto finish = std::chrono::high_resolution_clock::now();

  std::cout << std::chrono::duration_cast<std::chrono::duration<double> >(finish-start).count() << "s\n";

  return 0;
}
