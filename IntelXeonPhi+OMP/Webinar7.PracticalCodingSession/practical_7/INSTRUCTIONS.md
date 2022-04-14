# Further Optimizing Finite Differences

In practical exercise 4 we looked at a finite difference code
and how to parallelize and vectorize it.

Here we investigate advanced optimizations designed to not only
utilize all vector lanes and available cores, but also
make efficient use of other important resources which 
impact performance. These include the various memory subsystems of
the target architecture such as high bandwidth memory (MCDRAM),
cache. Also investigated are hardware threads and thread affinity

Step 1:

First run the code as is to get a basline.

Step 2:

Increase the problem size to be 1024 x 1024 x 1024 and run it again.
Increase the problem again to 2048 x 2048 x 2048.

Step 3: 

To explore why the performance changes lets use Intel vtune to
examine what happened.

First, to get access to vtune on TACC use the command:
module load vtune

Now to perform an analysis:
amplxe-cl -collect general-exploration ./fd
and to view the results
amplxe-cl -report summary -r r000ge

CPI is the cycles per instruction.  Ideally this should be 0.5, but a
value around 1 is usually pretty good.  Look at L1 hit rate and L2 hit
rate (measured in percentages).  VPU untilization gives the efficiency
with with the VPU's are untilized.

Step 4: 

Block the three loops in the laplacian computations as described in
yesterday's presentation to improve cashe usage for the larger
problems. Experiment with different block sizes in each dimension.

Step 5 (optional):

Repeat the above steps with the nbody code