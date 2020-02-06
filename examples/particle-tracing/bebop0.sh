#!/bin/bash

#SBATCH --job-name=p1
#SBATCH --account=pedal
#SBATCH --partition=knlall
#SBATCH --constraint knl,quad,cache
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --output=p0.%j.%N.out
#SBATCH --error=p0.%j.%N.error
#SBATCH --mail-user=mraj@anl.gov
#SBATCH --time=02:00:00

# number of procs

# executable
exe=./ptrace-iexchange

# inout file
infile="/home/mraj/datasets/nek/nek5000.nc"

# max number of advection steps
max_steps=1024
# max_steps=128
# max_steps=8

# seed rate (seed particles every this many grid points in each dimension)
sr=8

# domain extents
mins="0 0 0"
maxs="511 511 511"

# without prediction
prediction=0

num_procs=2
opts="--blocks 2 --max-rounds 9999 --synthetic 0 --check 0"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
mpiexec -n $num_procs $exe $args

num_procs=4
opts="--blocks 4 --max-rounds 9999 --synthetic 0 --check 0"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
mpiexec -n $num_procs $exe $args

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

# num_procs=64
# opts="--blocks 64 --max-rounds 9999 --synthetic 0 --check 0"
# args="$opts $infile $max_steps $sr $mins $maxs $prediction"
# mpiexec -n $num_procs $exe $args
# num_procs=128
# opts="--blocks 128 --max-rounds 9999 --synthetic 0 --check 0"
# args="$opts $infile $max_steps $sr $mins $maxs $prediction"
# mpiexec -n $num_procs $exe $args
# num_procs=256
# opts="--blocks 256 --max-rounds 9999 --synthetic 0 --check 0"
# args="$opts $infile $max_steps $sr $mins $maxs $prediction"
# mpiexec -n $num_procs $exe $args

# with prediction
prediction=1

num_procs=2
opts="--blocks 2 --max-rounds 9999 --synthetic 0 --check 0"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
mpiexec -n $num_procs $exe $args

num_procs=4
opts="--blocks 4 --max-rounds 9999 --synthetic 0 --check 0"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
mpiexec -n $num_procs $exe $args

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

# num_procs=64
# opts="--blocks 64 --max-rounds 9999 --synthetic 0 --check 0"
# args="$opts $infile $max_steps $sr $mins $maxs $prediction"
# mpiexec -n $num_procs $exe $args

# num_procs=128
# opts="--blocks 128 --max-rounds 9999 --synthetic 0 --check 0"
# args="$opts $infile $max_steps $sr $mins $maxs $prediction"
# mpiexec -n $num_procs $exe $args

# num_procs=256
# opts="--blocks 256 --max-rounds 9999 --synthetic 0 --check 0"
# args="$opts $infile $max_steps $sr $mins $maxs $prediction"
# mpiexec -n $num_procs $exe $args
