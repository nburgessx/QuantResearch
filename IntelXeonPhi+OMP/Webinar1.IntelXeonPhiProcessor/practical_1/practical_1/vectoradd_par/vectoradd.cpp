#include <iostream>
#include <vector>
#include <chrono>
int main(int argc,char**argv){
  int n=500000000;
  std::vector<float> x(n,1.0),y(n,1.0),z(n,0.0);
  auto start = std::chrono::steady_clock::now();
  for(int i=0;i<n;i++){
    z[i]=x[i]+y[i];
  }
  auto finish = std::chrono::steady_clock::now();
  std::cout << std::chrono::duration_cast<std::chrono::duration<double> >(finish-start).count() << "s\n";
  return 0;
}
