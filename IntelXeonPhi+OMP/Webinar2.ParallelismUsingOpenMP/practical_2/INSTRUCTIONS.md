# Using OpenMP for parallelism


# Instructions

In this directory you will find
a function which computes a Chebyshev
expansion for an input function f.

The expansion is computed such that
the resulting error is within
numerical roundoff, which is
confirmed by computing the norm of
the error.


There is a C and fortran in the
C and fortran directories
respectively. You can build the unoptimized code
by running 

```
>>> make
>>> ./cheby.exe 5000
```


The runtime of this example is quadratic, so 
increasing the input argument can significantly
increase runtime.

In order to make the most meaningful comparisons
and see if your code changes indeed improved performance,
it is highly advised you make any changes in a 
fresh directory with all the original files in it.
The original files are contained in the
"original." There are accompanying
solutions as well


The overview of steps are as follows:

  1. Compile and run the code on the host and record the times reported for find_Tcoeffs and evaluate_Tpoly.
  2. Parallelise the loop in gen_data using OpenMP. Use an omp for pragma with a static schedule.
  3. Parallelise the loop in the main function which is used to evaluate the norm of the errors. Use a reduction to get the final value of var.
  4. Parallelise the loops in find_Tcoeffs and evaluate_Tpoly. How does each function scale, going from 1 to 288 threads in suitable steps?
    1. You can change number of threads with environment variable OMP_NUM_THREADS, which you can set in a shell script that executes the code
    2. (Optional) Try changing the schedules between dynamic, static, auto, to see if it makes a difference
 
# Step 5: More Practice 


In order to practise more of the OpenMP introduced in the talk.) So far we have created a new parallel region for each parallel loop. Instead, starting from the original serial code, create a single parallel region at the start of the program as follows:
```cpp
#pragma omp parallel default(none) \
private(i,j,dif,t1,t2) \
shared(argc,argv,n1,xbar) \
shared(fr,coeffs,fe,var) \
shared(macheps,nrm2,tol,stderr)
for (i=1; i<argc; i++)
//…
```
  The existing omp for directives should become orphaned directives. 
Ensure that the master thread is the only thread capable of printing.
 1.   Ensure that only a single thread allocates and frees the pointers.
 2.  Note that t1 and t2 are private. This time we want to print the cumulative time taken by all threads. Sum the time taken by each thread into a new shared variable tall. Make sure you protect access to tall. Reset tall to zero after it has been used to print the cumulative time taken by find_Tcoeffs and reuse it to time evaluate_Tpoly in the same way.
 3.  You could also declare t1 and t2 as shared pointers and allocate them at the start of the program, outside of the parallel region (i.e. use t1[i] to record the time for thread i). Which OpenMP runtime function would you use to ensure that t1 and t2 each point to an area of memory large enough to record the time for all threads in the team?



