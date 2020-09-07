#!/bin/bash

#SBATCH --job-name=p1
#SBATCH --account=pedal
#SBATCH --partition=knlall
#SBATCH --constraint knl,quad,cache
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --output=p0.%j.%N.out
#SBATCH --error=p0.%j.%N.error
#SBATCH --mail-user=mraj@anl.gov
#SBATCH --time=00:30:00



# executable
exe=./btrace

# inout file
# infile="/nfs/proj-tpeterka/hguo/data/nek5000.nc"
infile="/home/mraj/datasets/nek/nek5000.nc"
# infile="/nfs/proj-tpeterka/jiang/data/isabel0.nc"
# infile="/home/mraj/datasets/isabel/isabel0.nc"

# max number of advection steps
max_steps=512

# seed rate (seed particles every this many grid points in each dimension)
sr=64
#sr=3
# sr=128

# domain extents
mins="0 0 0"
maxs="511 511 511"
# maxs="499 499 99"

# without prediction
prediction=0

num_procs=8
opts="--blocks 8 --max-rounds 9999 --synthetic 0 --check 0"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
mpiexec -n $num_procs $exe $args

num_procs=16
opts="--blocks 16 --max-rounds 9999 --synthetic 0 --check 0"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
mpiexec -n $num_procs $exe $args

num_procs=32
opts="--blocks 32 --max-rounds 9999 --synthetic 0 --check 0"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
mpiexec -n $num_procs $exe $args

num_procs=64
opts="--blocks 64 --max-rounds 9999 --synthetic 0 --check 0"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
mpiexec -n $num_procs $exe $args
