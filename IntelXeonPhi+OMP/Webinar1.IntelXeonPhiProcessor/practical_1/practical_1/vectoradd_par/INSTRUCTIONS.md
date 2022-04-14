# Vector Add

This example will show how to parallelize a simple vector add
using OpenMP.


## Instructions

In this example, unlike the hello world example, we will make
use of some extra features of OpenMP. In this case
we will make use of the fact that OpenMP can automatically
parallelize (certain) for-loops, and also make smart
choices about how which iterates the different threads execute,
this is known as "scheduling," and this course
will show later how to exercise more control over scheduling
in order to improve performance. For now we let OpenMP
make the scheduling choices


To achieve this, take the existing vector add code in this
directory and add the following OpenMP directive

```cpp
#pragma omp parallel for num_threads(4)
{
//For-loop to paralellize here
}
```

and in Fortran

```fortran
!$OMP PARALLEL DO num_hreads(4)
//Do-loop to paralellize here
!$OMP END PARALLEL DO
```

Note that this code uses some C++11 features, namely type deduction
and time measurement using the "chrono" standard library. Most C++ 
compilers will not assume C++11 or beyond, but you can usually
specify it with the command line flag "--std=c++11."

A makefile is supplied but you should familiarize yourself
with the correct compiler flags before using it.


## Optional

Change the value of num_threads to investigate the impact of different values.
Often the ideal num_threads would be equal to the number of cores on your machine,
but if there is hyperthreading that may not be the ideal parameter. It is
always best to parameterize options like this and investigate objectively
which produces the best result.
