#!/bin/bash

#SBATCH --job-name=p1
#SBATCH --account=pedal
#SBATCH --partition=knlall
#SBATCH --constraint knl,quad,cache
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=32
#SBATCH --output=p0.%j.%N.out
#SBATCH --error=p0.%j.%N.error
#SBATCH --mail-user=mraj@lcrc.anl.gov
#SBATCH --time=00:30:00

# number of procs

# executable
exe=./ptrace-iexchange

# inout file
infile="/home/mraj/datasets/nek/nek5000.nc"

# max number of advection steps
max_steps=1024

# seed rate (seed particles every this many grid points in each dimension)
sr=64

# domain extents
mins="0 0 0"
maxs="511 511 511"

# without prediction
prediction=0

num_procs=8
opts="--blocks 8 --max-rounds 9999 --synthetic 0 --check 0"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
mpiexec -n $num_procs $exe $args