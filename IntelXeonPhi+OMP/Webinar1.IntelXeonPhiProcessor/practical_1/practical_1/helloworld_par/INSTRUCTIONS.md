# Hello World

This simple example will establish that your environment is set up to build
a minimal running example using the Intel C++ or Fortran compiler and OpenMP
for parallelism.

## Instructions


### Step 1
Modify the hello world example to run
in parallel with four threads using OpenMP and print out each thread's
number.

The OpenMP directive for this simple kind of parallelism
is


```cpp
#pragma omp parallel num_threads(4)
{
  //Code goes here
}
```

```fortran
integer :: id
!$OMP PARALLEL num_threads(4), private(id)
  //Code goes here
!$OMP END PARALLEL
```

NOTE: The private clause ensures that each thread has is own copy of the variable, 'id'.
This is necessary in fortran since it cannot be declared within the parallel region.

Inside an OpenMP block the current thread id can be discovered with

``` cpp
 int id=omp_get_thread_num();
```

``` fortran
 id=omp_get_thread_num()
```


Don't forget to include the openmp header in C

```cpp
#include <omp.h>
```

or the openmp module in fortran

```fortran
use omp_lib
```



and you also need to link the openmp runtime, which can be done with the command line flag

```shell
icpc -qopenmp helloworld.cpp -o helloworld
```
or

```shell
ifort -qopenmp helloworld.cpp -o helloworld
```


### Step 2 (Optional)

On the KNL the output from should be nondterministic. The reason
for this is that all threads are writing to the same output.
This can be resolved by using a critical region which ensures
that only one thread can execute the specified region at a time.

This can be achieved with

```cpp
#pragma omp critical
{
//Code to declare as "critical" goes here
}
```

or

```fortran
!$OMP CRITICAL
//Code to declare as "critical" goes here
!$OMP END CRITICAL
```


rebuilding and running this code should give the expected result.
