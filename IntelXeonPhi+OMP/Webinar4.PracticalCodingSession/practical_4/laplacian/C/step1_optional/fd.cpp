#include <stdint.h>

#include <cmath>
#include <iostream>
#include <vector>
#include <chrono>

#include "coeffs12.h"


#define HALF 6
#define ORDER 12


typedef float float32_t;
typedef double float64_t;


int64_t get_stencil(int64_t n,int64_t i);
template <typename Real>

void fd12(std::vector<Real>& in,std::vector<Real>& out,
    Real    dx,Real    dy,Real    dz,
    int64_t nx,int64_t ny,int64_t nz);


std::vector<float32_t> input(
    float32_t x0,float32_t y0,float32_t z0,
    float32_t dx,float32_t dy,float32_t dz,
    int64_t   nx,int64_t   ny,int64_t   nz);
std::vector<float32_t> exact_laplacian(
    float32_t x0,float32_t y0,float32_t z0,
    float32_t dx,float32_t dy,float32_t dz,
    int64_t   nx,int64_t   ny,int64_t   nz);
template <typename Real>
void fd12(std::vector<Real>& in,std::vector<Real>& out,int64_t nx,int64_t ny,int64_t nz);
float64_t compute_err(std::vector<float> exact,std::vector<float> approx);



int main(int argc,char** argv){




  int64_t nx=512;
  int64_t ny=512;
  int64_t nz=512;
  float32_t dx=0.5;
  float32_t dy=0.5;
  float32_t dz=0.5;
  float32_t x0=0.0;
  float32_t y0=0.0;
  float32_t z0=0.0;
  


  auto in   = input     (x0,y0,z0,dx,dy,dz,nx,ny,nz);
  auto exact=exact_laplacian(x0,y0,z0,dx,dy,dz,nx,ny,nz);
  std::vector<float32_t> out(nx*ny*nz,0.0);

  auto start = std::chrono::steady_clock::now();
  fd12(in,out,dx,dy,dz,nx,ny,nz);
  auto finish = std::chrono::steady_clock::now();
  std::cout << (nx*ny*nz) / std::chrono::duration_cast<std::chrono::duration<double> >(finish-start).count() << " nodes/s\n";

  std::cout << "error = " << compute_err(exact,out) << std::endl;



  return 0;
}


#define id(ix,iy,iz) ((ix) + nx*((iy)+(ny*(iz))))

template <typename Real>
void fd12(std::vector<Real>& in,std::vector<Real>& out,
    Real    dx,Real    dy,Real    dz,
    int64_t nx,int64_t ny,int64_t nz){


#pragma omp parallel for \
  num_threads(68) \
  proc_bind(spread) \
  schedule(static,1)
  for(int64_t ix=0;ix<nx;ix++){
#pragma omp parallel for \
    num_threads(2) \
    proc_bind(close) \
    schedule(static,1)
    for(int64_t iy=0;iy<ny;iy++){
      for(int64_t iz=0;iz<nz;iz++){
        int64_t jx=get_stencil(nx,ix);
        int64_t jy=get_stencil(ny,iy);
        int64_t jz=get_stencil(nz,iz);
        Real d2fdx2=0.0;
        Real d2fdy2=0.0;
        Real d2fdz2=0.0;


        /*Second order x-derivative.*/
        for(int64_t i=-jx;i<=ORDER-jx ;i++){
          d2fdx2 += coeffs[jx][i+jx]*in[id(ix+i,iy,iz)]/(dx*dx);
        }
        /*Second order y-derivative.*/
        for(int64_t i=-jy;i<=ORDER-jy ;i++){
          d2fdy2 += coeffs[jy][i+jy]*in[id(ix,iy+i,iz)]/(dy*dy);
        }
        /*Second order z-derivative.*/
        for(int64_t i=-jz;i<=ORDER-jz ;i++){
          d2fdz2 += coeffs[jz][i+jz]*in[id(ix,iy,iz+i)]/(dz*dz);
        }

        /*Approximate Laplacian.*/
        out[id(ix,iy,iz)]=d2fdx2+d2fdy2+d2fdz2;
      }
    }
  }

}

int64_t get_stencil(int64_t n,int64_t i){
  if(i<HALF){
    return i;
  }
  if(i>=n-HALF){
    auto m=n-HALF;
    return HALF + (i-m) + 1;
  }
  return HALF;
}


/*Compute error between approximate and exact.*/
float64_t compute_err(std::vector<float> exact,std::vector<float> approx){
  float64_t err=100.0;
  float64_t nrm=1e-6;


  err=0.0;
  nrm=0.0;
  for(int64_t i=0;i<exact.size();i++){
    float64_t e=exact[i];
    float64_t a=approx[i];
    err += (e-a)*(e-a);
    nrm += e*e;
  }
  err=sqrt(err);
  nrm=sqrt(nrm);





  return err/nrm;
}



/*Exact laplacian of test input.*/
std::vector<float32_t> exact_laplacian(
    float32_t x0,float32_t y0,float32_t z0,
    float32_t dx,float32_t dy,float32_t dz,
    int64_t   nx,int64_t   ny,int64_t   nz){
  std::vector<float32_t> out(nx*ny*nz,0.0);
  for(int64_t ix=0;ix<nx;ix++)
    for(int64_t iy=0;iy<ny;iy++)
      for(int64_t iz=0;iz<nz;iz++){
        float32_t x=x0+ix*dx;
        float32_t y=y0+iy*dy;
        float32_t z=z0+iz*dz;
        float32_t d2fdx2=-sin(x+y+z);
        float32_t d2fdy2=-sin(x+y+z);
        float32_t d2fdz2=-sin(x+y+z);
        out[id(ix,iy,iz)]=d2fdx2+d2fdy2+d2fdz2;
      }
  return out;
}

/*To be used as testing input for 
 * finite difference stencil.*/
std::vector<float32_t> input(
    float32_t x0,float32_t y0,float32_t z0,
    float32_t dx,float32_t dy,float32_t dz,
    int64_t   nx,int64_t   ny,int64_t   nz){

  std::vector<float32_t> out(nx*ny*nz,0.0);
  for(int64_t ix=0;ix<nx;ix++)
    for(int64_t iy=0;iy<ny;iy++)
      for(int64_t iz=0;iz<nz;iz++){
        float32_t x=x0+ix*dx;
        float32_t y=y0+iy*dy;
        float32_t z=z0+iz*dz;
        out[id(ix,iy,iz)]=sin(x+y+z);
      }
  return out;
}








