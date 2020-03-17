#!/bin/bash
#COBALT -t 30
#COBALT -n 128
#COBALT -q default
#COBALT --attrs mcdram=cache:numa=quad 
#COBALT -A CSC249ADCD01 

# executable
exe=./ptrace-iexchange

# inout file
infile="/home/mraj/data/nek5000.nc"

# max number of advection steps
max_steps=1024

# seed rate (seed particles every this many grid points in each dimension)
sr=4

# domain extents
mins="0 0 0"
maxs="511 511 511"

# without prediction
prediction=0

num_procs=1024
opts="--blocks 1024 --max-rounds 9999 --synthetic 0 --check 0"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
aprun -n $num_procs -N 8 $exe $args

num_procs=2048
opts="--blocks 2048 --max-rounds 9999 --synthetic 0 --check 0"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
aprun -n $num_procs -N 16 $exe $args

num_procs=4096
opts="--blocks 4096 --max-rounds 9999 --synthetic 0 --check 0"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
aprun -n $num_procs -N 32 $exe $args

num_procs=8192
opts="--blocks 8192 --max-rounds 9999 --synthetic 0 --check 0"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
aprun -n $num_procs -N 64 $exe $args


# without prediction
prediction=1

num_procs=1024
opts="--blocks 1024 --max-rounds 9999 --synthetic 0 --check 0"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
aprun -n $num_procs -N 8 $exe $args

num_procs=2048
opts="--blocks 2048 --max-rounds 9999 --synthetic 0 --check 0"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
aprun -n $num_procs -N 16 $exe $args

num_procs=4096
opts="--blocks 4096 --max-rounds 9999 --synthetic 0 --check 0"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
aprun -n $num_procs -N 32 $exe $args

num_procs=8192
opts="--blocks 8192 --max-rounds 9999 --synthetic 0 --check 0"
args="$opts $infile $max_steps $sr $mins $maxs $prediction"
aprun -n $num_procs -N 64 $exe $args