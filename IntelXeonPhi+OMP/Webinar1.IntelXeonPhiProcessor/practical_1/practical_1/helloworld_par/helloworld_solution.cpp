#include <iostream>
#include <omp.h>
int main(int argc,char**argv){

#pragma omp parallel num_threads(4)
  {
    int id=omp_get_thread_num();
    std::cout<< "Hello from thread "<<id<<"\n";
  }

  return 0;
}
