# Vector Add

This example will show how to vectorize a simple vector add
using OpenMP.


## Instructions
This example shows how to instruct the compiler to
generate parallel and vector instructions at the same time.  A
relatively recent addition to the OpenMP standard was a set of
directives designed to target this functionality of modern CPUs.

While the functionality of this new feature 
is extensive, this example will use it in the
simplest possible way just for illustration. 
A later lecture will discuss the details and subtleties of
vectorization combined with parallelization.

To enable vectorization and parallelism for the
vector add example, simply use the following
OpenMP directive.


```cpp
#pragma omp parallel for simd num_threads(4)
{
//For-loop to paralellize here
}
```

```fortran
!$OMP PARALLEL DO SIMD num_threads(4)
//Do-loop to paralellize here
!$OMP PARALLEL DO SIMD 
```



# (Optional)
Compare the performance of the unoptimized vector add
with the multicore vector add and the single-core
vectorized vector add.

Note that the performance improvement for this example is not too
significant.  We will show examples in later examples that get
significant performance improvements by adding both vecotirization and
parallelization.
