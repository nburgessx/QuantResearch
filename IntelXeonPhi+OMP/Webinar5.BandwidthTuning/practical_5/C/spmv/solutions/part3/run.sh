export KMP_AFFINITY=granularity=fine,balanced
export OMP_NUM_THREADS=68
numactl -m1 ./main
