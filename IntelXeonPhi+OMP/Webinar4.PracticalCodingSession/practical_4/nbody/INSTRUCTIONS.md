# Nbody Practical Exercise

This exercise extracts a highly another common computational kernel
seen across a wide variety applications, an "nbody" force computation
where a number of particles interact with each and every other
particles through a long-range force.  Typically this is shown using
the gravitational force, but here we use compute an appoximate force
between molecules using what is known as the Lennard-Jones force.

Our purpose is to show several optimizations that other practicals
will not cover, but also to use this as a test case to show what raw
performace we can expext to get out the Xeon Phi processor (Knight's
Landing).  The peak performance we should expect is 3 Tflops.  Let's
see what we can get on a "real" code.


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


## What the starting code comes with 

The starting code comes with the straightforward serial implementation
of Lennard-Jones force calculation.  It is characterized by two loops
over the number of molecules in our model system, an inner and outer
loop.


## Step 1: Parallelizing

The code can be parallelized normally (by adding the appropriate
pragma), but it might be worth exploring different options such as
thread affinity and hyperthreading. To utilize hyperthreading, just
set the num_threads() parameter to a larger number than the number of
available cores. Generally speaking it is not usually worthwhile to
try anything other than 1X threads per core, 2X threads per core, or
4X threads per core.  Thread affinity can interact in interesting ways
with hyperthreading as well, so I recommend trying many combinations.


### Step 2: Vectorizing inner loop as-is

Vectorization of this code is where the interesting work happens. As a
first approach to vectorizing it never hurts to just try the inner
loop, which in this case is the loop over the stencil derivative
calculation. This can be achieved by adding the "omp simd" pragma on
the inner loop.


### Step 3: Strength Reduction and Precision Control.

This code compute the sqrt and the pow functions on the inner loop.
The can be computed in single precision using sqrtf and powf.
Secondly, we can add the compiler option -fp-model fast=2 will force
the code to use less accurate implementations of transcendental
functions and make more aggresive compiler options.


### Step 4: Remove the masked vectorization

Removing the if test in the inner loop will increase the performance
of this vectorized loop.  But, we need to change the way the distance
between molecules is calculated by added a (very) small value to the
variable "r" in the code (1.e-20) so that the force between the molecules
is zero when the loop indices i and j are equal.



