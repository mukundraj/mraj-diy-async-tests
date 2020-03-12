#!/bin/bash
#COBALT -t 15 
#COBALT -n 2
#COBALT -q debug-cache-quad 
#COBALT --attrs mcdram=cache:numa=quad 
#COBALT -A CSC249ADCD01 
echo "Starting Cobalt job script on 2 nodes with 64 ranks on each node" 
# aprun -n 128 -N 64 --env OMP_NUM_THREADS=4 -cc depth -d 4 -j 4 myprogram.exe

# executable

exe=./ptrace-iexchange

# inout file
infile="/home/mraj/data/nek5000.nc"

# max number of advection steps
max_steps=1024
# max_steps=128
# max_steps=8

# seed rate (seed particles every this many grid points in each dimension)
sr=4

# domain extents
mins="0 0 0"
maxs="511 511 511"

# without prediction
prediction=0

num_procs=128
opts="--blocks 128 --max-rounds 9999 --synthetic 0 --check 0"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
aprun -n $num_procs -N 64 --env -d=2 -j 2 $exe $args
# mpiexec -n $num_procs $exe $args