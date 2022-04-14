# Vector Add

This example will show how to vectorize a simple vector add
using OpenMP.


## Instructions
This example shows how to instruct the compiler
to generate vector instructions. A relatively
recent addition to the OpenMP standard
was a set of directives designed to target
this functionality of modern CPUs.

While the functionality of this new feature 
is extensive, this example will use it in the
simplest possible way just for illustration. 
A later lecture will discuss the details and subtleties of
vectorization.

To enable vectorization and parallelism for the
vector add example, simply use the following
OpenMP directive.


```cpp
#pragma omp simd
{
//For-loop to paralellize here
}
```

and in Fortran

```fortran
!$OMP SIMD
//Do-loop to paralellize here
!$OMP END SIMD
```



# (Optional)
Compare the performance of the unoptimized vector add
with the multicore vector add and the single-core
vectorized vector add.
