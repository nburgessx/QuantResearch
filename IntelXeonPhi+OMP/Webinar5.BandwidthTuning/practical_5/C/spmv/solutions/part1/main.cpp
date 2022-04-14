#include <iostream>
#include <vector>
#include <chrono>

#include "sparse.h"



int main(int argc,char** argv){





  int64_t m=100000;
  int64_t n=100000;
  sparse_t<double> A(m,n);
  std::vector<double> b(m,1.0);
  auto start = std::chrono::steady_clock::now();
  auto x = cg(A,b);
  auto finish = std::chrono::steady_clock::now();
  std::cout << "Time to solve:" << std::chrono::duration_cast<std::chrono::duration<double> >(finish-start).count() << "s\n";
  

  return 0;
}
