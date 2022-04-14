export OMP_PLACES=threads
numactl -m1 ./ray 5000000 1000
