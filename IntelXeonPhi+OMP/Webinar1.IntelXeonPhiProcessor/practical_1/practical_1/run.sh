#!/bin/bash
#----------------------------------------------------
# Example SLURM job script to run OpenMP applications
# on TACC's Stampede system.
#----------------------------------------------------
#SBATCH -J webinar_exercise           # Job name
#SBATCH -o webinar_exercise.o%j   # Name of stdout output file(%j expands to jobId)
#SBATCH -e webinar_exercise.o%j   # Name of stderr output file(%j expands to jobId)
#SBATCH -p normal             # Submit to the 'normal' or 'development' queue
#SBATCH -N 1                  # Total number of nodes requested (16 cores/node)
#SBATCH -n 1                  # Total number of mpi tasks requested
#SBATCH -t 01:30:00           # Run time (hh:mm:ss) - 1.5 hours

# SBATCH -A A-NAG-and-Intel-Webina #Use webinar allocation


# Sets number of threads to number of cores in this case, which is 68
export OMP_NUM_THREADS=68
# Run the OpenMP application
./main


#An alternative partition which has KNLs under different settings.
#Safe to ignore this for now, but it will come up in 
#later webinars.
#SBATCH -p Flat-Quadrant    # Submit to the 'normal' or 'development' queue
