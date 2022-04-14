# Finite Difference Practical Exercise


This exercise extracts a highly common computational kernel seen
across a wide variety of industries, a Laplacian calculation of a
function sampled at regularly spaced intervals in 3D. This
kernel can be the foundation for both a wave propagation solver
(such as acoustic wave equation) or a diffusion equation solver
(such as the Laplace equation, or heat equation). As such
it is a commonly optimized routine.

This exercise will lead you on the first steps to
optimizing this routine. So far you have learned
the basics of parallelizing and vectorizing code
using OpenMP. You have also learned some ways
these tools can be misued to produce an 
incorrect answer. Here you will put these
tools to use as well as learn some new
features.


## Suggestions before proceeding

In a production environment this work should be done under
careful version control, but here to keep things simple
it is recommended you keep a directory for experimenting
and new directories for every new major change. This is
so that you can keep track of performance gains or losses 
as you make changes to the code.

Another suggestion is to make sure anything that is parameterized
be done so in a way that it can be easily changed. For example
all parameters could be defined in a header file which is 
included in "fd.cpp." This is useful especially if you ever
do work which must perform well across architectures.

Finally: Make sure the test passes after each change to
the code. The computed relative error should not exceed 1e-04
(roughly corresponding to four digits of accuracy).


## What the starting code comes with
The starting code comes with the straightforward serial implementation
of the stencil for you to optimize with parallelism and vectorization.
It also contains a timer for this code, and finally it has a test
routine to make sure the error is what it should. Note that error
will be a little different with each change to code because
floating point arithmetic is not associative. The important
thing to guarantee here is that the error is sufficiently small -
which for this application is around 1e-04 or 1e-05 (which
due to the way I normalize the error, corresponds to about
4 or 5 digits of accuracy).



## Step 1: Parallelizing

(solution found in "step1" directory)

### First steps
The code can be parallelized normally (by adding the appropriate pragma),
but it might be worth exploring different options such as thread affinity
and hyperthreading. To utilize hyperthreading, just set the
num_threads() parameter to a larger number than the number of
available cores. Generally speaking it is not usually worthwhile
to try anything other than 1X threads per core, 2X threads per core, or 4X threads per core.
Thread affinity can interact in interesting ways with hyperthreading as well, so I recommend
trying many combinations.


### Nested Parallelism With Hot Teams(Optional)
(solution found in "step1_optional directory")

Instead of simply overprescribing the number of threads to achieve hyperthreading,
you can also explicitly specify how the hyperthreads parallelize the code. For example
you could do something like:


```c
#pragma omp parallel for num_threads(64) 
for(int i=0;i<big;i++){
  #pragma omp parallel for num_threads(4)
  for(int j=0;j<small;j++){
  // code here
  }
}
```

This can potentially give big gains. If you try this
I recommend enabling Intel's "hot teams" mode which significantly lowers the overhead of the inner
parallel loop. This can be be done in a shell script which executes the code
```bash
 #Enables Intel hot teams mode for efficient parallel inner loops.
 export KMP_HOT_TEAMS_MODE=1
 #Sets the maximum nested-ness of parallelism to 2.
 export KMP_HOT_TEAMS_MAX_LEVEL=2
 ./fd
```


Another advantage of this approach is it gives much simpler control over affinity of the threads, since the
outer threads can be forced to occupy whole cores while the inner threads forced to
occupy only the associated hyperthreads.

## Vectorization



### Step 2: Vectorizing inner loop as-is

(solution found in "step2" directory)

Vectorization of this code is where the interesting work happens. As a first approach
to vectorizing it never hurts to just try the inner loop, which in this case is the
loop over the stencil derivative calculation. This can be achieved by
adding the "omp simd" pragma on the three inner loops.




Unfortunately this does not provide any performance gain. One problem is that
the inner loop is just too complicated to easily vectorize, you can confirm this
by outputting the vectorization report and seeing what sort of memory access patterns
were assumed - in this case it should be gather/scatter because of the complex indexing,


Finally on the KNL we already mentioned it has an effective 32-float wide vector lane, but the
stencil length (what we just vectorized over) is only 13 wide. We could try to combine
all three derivative calculations into one, but that won't work either because their respective indexing
are all different (unit stride in x direction, gather/scatter in y and z). This leads into 
a very common optimization tactic



### Loop reordering

(solution found in "step3" directory)

You may have noticed from the optimization report that the x-direction has unit stride access.
It is very common in stencil calculations as well as n-dimensional array calculations
to re-order for-loops so that the inner-most loop is looping over the fastest dimension, and then vectorize
over this. This will solve the problem of complex indexing and should also fully occupy the vector lanes.
Here is a simple example in C of this strategy:


```c


int m=32;
float A[m][m][m];

for(int k=0;k<m;k++)
  for(int i=0;i<m;i++)
    for(int j=0;j<m;j++)
      A[i][j][k]=1.0;

```



which can be translated to the equivalent


```c


int m=32;
float A[m][m][m];

for(int i=0;i<m;i++)
  for(int j=0;j<m;j++)
    for(int k=0;k<m;k++)
      A[i][j][k]=1.0;

```

and observe that because of C's storage order for arrays, the inner-most loop is always
accessing contiguous data (which for vectorization translates to unit-stride access).

The code as provided for finite difference however is a little more
"realistic" because it is also incorporating boundary conditions (or
one possible way of doing so). When boundary conditions like this are
included in the stencil computations, you usually have to separate out
the conditions before reordering loops (or else possibly apply the
wrong boundary conditions, which in a physics context usually results
in wildly incorrect answers).

One way to achieve this would be to do what is called "loop fission," this is where
the loop is split into three identical loops each containing part of the original loop
body. In this case we could split out the x,y, and z derivatives and this will allow
loop reordering in the y and z dimensions. The x dimension is a little more complicated
but it is already the fastest dimension so let's leave that for now and aim at the bigger
cause of trouble.


Applying these optimizations the code should end up looking like the following

```c

    for(int64_t iy=0;iy<ny;iy++){
      int64_t jy=get_stencil(ny,iy);
      for(int64_t iz=0;iz<nz;iz++){
        int64_t jz=get_stencil(nz,iz);
        /*Second order x-derivative.*/
        for(int64_t ix=0;ix<nx;ix++){ 
        int64_t jx=get_stencil(nx,ix);
          for(int64_t i=-jx;i<=ORDER-jx ;i++){
              out[id(ix,iy,iz)] += coeffs[jx][i+jx]*in[id(ix+i,iy,iz)]/(dx*dx);
            }
        }
        /*Second order y-derivative.*/
        for(int64_t i=-jy;i<=ORDER-jy ;i++){
          for(int64_t ix=0;ix<nx;ix++){ 
              out[id(ix,iy,iz)] += coeffs[jy][i+jy]*in[id(ix,iy+i,iz)]/(dy*dy);
            }
        }
        /*Second order z-derivative.*/
        for(int64_t i=-jz;i<=ORDER-jz ;i++){
          for(int64_t ix=0;ix<nx;ix++){ 
            out[id(ix,iy,iz)] += coeffs[jz][i+jz]*in[id(ix,iy,iz+i)]/(dz*dz);
          }
        }
    }
  }

```

This code combined several steps into one:

* Reordered loops to make the x direction the inner-most loop
  * Guaranteeing unit-stride indexing
* Applied loop fission to separate the boundary logic for different dimensions

This optimization is very portable, therefore it should yield the best performance on both
KNL and non-KNL architectures, and we will see later that this holds true even without vectorization.





## Summary and Look Ahead


In this practical exercise we took a common computational kernel,
a laplacian calculation for regularly spaced data, and parallelized
and vectorized it. The parallelism was more or less straightforward
except that some tunable parameters can be changed to yield better results.
The vectorization required us to take a closer look at the code and modify it
to be more amenable to using faster vector code.

These techniques are very important across a wide variety of situations,
but significantly more performance is possible by further optimizing to
make better use of memory - in a later practical exercise we will do this.


