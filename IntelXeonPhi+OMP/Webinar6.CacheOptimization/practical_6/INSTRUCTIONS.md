# Optimizing a Ray Tracing code


Ray tracing in inhomogeneous media is a high-frequency
approximation to the acoustic wave equation. Moving
to high frequencies effectively eliminates complex
wave features such as diffraction, which
from the point of view of programming means
that each of the rays may be treated independently
and therefore are easy to parallelize. For this reason
it is a popular approximation in areas requiring
the simulation of waves, but where full PDE simulation
can be prohibitively expensive.



Ray tracing exercises many of the optimizations seen from
the accompanying lecture. A ray for example can be described
using its three position coordinates (i.e. three floats) and its momentum
vector (which is itself three floats). The simplest way of implementing
ray tracing is to define a struct representing a single ray and
then to make an array of these, or the "array of structs"
approach.


The array-of-structs approach may be found in the "original" directory
which contains an unoptimized ray-tracer. The first step will be to transform
this to struct-of-arrays for better data access patterns.



## Step 1: Intel SDLT

(solution found in "step1" directory)

For step 1, take the provided code and use Intel SDLT
to abstract away the data layout. Try both struct-of-array
as well as array-of-structs and determine their relative performance.


An example usage of intel SDLT can be as follows

```cpp
#include <sdlt/sdlt.h>
class point_t{
  public:
    float x,y;
};
SDLT_PRIMITIVE(point_t,x,y);
//typedef sdlt::aos1d_container<point_t> Container;
typedef sdlt::soa1d_container<point_t> Container;
int n=50000;
Container points(n);
#pragma omp simd
for(int i=0;i<n;i++){
  #pragma forceinline recursive
  do_something(points[i]);
}
```



note the "pragma forceinline recursive." This was not described in
the lecture, but it is an intel-specific pragma that inlines function calls
and any function calls inside them. This is not a strictly necessary
step, but because of the way the sdlt library works (using templates)
the maximum benefit comes from inlining functions.

Test the performance of this code for varying numbers of rays
compared to the unoptimized version.


## Step 2: Blocking for Cache
(solution found in "step2" directory)


Another useful optimization for this code is to 
limit the number of rays that are stepped
forward in time. That means find an optimal
block size of rays that will fit in cache,
and loop through the rays on increments of
that block size. 


For example

```cpp
for(int i=0;i<n;i++){
  for(int j=0;j<n;j++){
    do_something(i,j);
  }
}
```

could become something like

```cpp
for(int i=0;i<n;i+=BLOCKI){
  int iend=MIN(n,i+BLOCKI);
  for(int j=0;j<n;j++){
    for(int ii=i;ii<iend;ii++)
      do_something(ii,j);
  }
}
```


The trick with blocking is understanding
which loops need to be blocked and how
to permute the resulting nested loops.
You have to think about how the data
is being accessed to get this right
and it may take some trial and error.
