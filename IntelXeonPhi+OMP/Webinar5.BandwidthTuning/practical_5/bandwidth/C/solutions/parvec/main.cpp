#include <iostream>
#include <chrono>
#include <vector>
#include <stdint.h>


#define NTIMES 20
#define NTHREADS 68

typedef double real_t;
void add (std::vector<real_t>& out,const std::vector<real_t>& in1,const std::vector<real_t>& in2);
void fill(std::vector<real_t>& out,      std::vector<real_t>& in1,      std::vector<real_t>& in2);


int main(int argc,char** argv){
  int64_t n=500000000;
  std::vector<real_t> z(n),x(n),y(n);
  fill(z,x,y);
  /*Exclude first run from timing.*/
  add(z,x,y);

  /*Begin timing now.*/
  auto start=std::chrono::steady_clock::now();
  for(int i=0;i<NTIMES;i++){ 
#pragma noinline
    add(z,x,y);
  }
  auto end =std::chrono::steady_clock::now();


  /*Compute effective bandwidth.*/
  int64_t nloads=2;
  int64_t nwrites=1;
  auto nbytes =(nloads+nwrites)*n*sizeof(real_t);
  double nGbytes=nbytes*1e-09*NTIMES;
  std::cout << "Effective GB: " << nGbytes << "\n";
  std::cout << "Bandwidth   : " << nGbytes / std::chrono::duration_cast<std::chrono::duration<double> >(end-start).count() << "GB/s\n";
  return 0;
}


void add(std::vector<real_t>& out,const std::vector<real_t>& in1,const std::vector<real_t>& in2){
  auto n=out.size();
#pragma omp parallel for simd num_threads(NTHREADS)
  for(auto i=0;i<n;i++) {
    out[i]=in1[i]+in2[i];
  }
}

void fill(std::vector<real_t>& out,      std::vector<real_t>& in1,      std::vector<real_t>& in2){
  auto n=out.size();
#pragma omp parallel for simd num_threads(NTHREADS)
  for(auto i=0;i<n;i++) {
    out[i]=0.0;
    in1[i]=1.0;
    in2[i]=1.0;
  }
}





