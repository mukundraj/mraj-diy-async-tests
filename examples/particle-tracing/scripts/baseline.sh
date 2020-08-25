#!/bin/bash

# number of procs
num_procs=4

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
opts="--blocks 8 --max-rounds 9999 --synthetic 0 --check 1"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
mpiexec -n $num_procs $exe $args


# mpiexec -n $num_procs $exe