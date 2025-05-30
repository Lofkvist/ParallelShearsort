#!/bin/bash -l

#SBATCH -A uppmax2025-2-247

#SBATCH -p node

#SBATCH -N 2

# Request 30 minutes of runtime for the job
#SBATCH -t 00:30:00

# Set the names for the error and output files
#SBATCH -o job.%J.out


# Add modules
ml gcc openmpi

#Run makefile and add acess
make
chmod +x quicksort


# Set the number of runs
N_RUNS=5
PIVOTING_STRAT=1
INPUT=/proj/uppmax2025-2-247/nobackup/A3/inputs/input125000000.txt

# Run with 1 proc
echo "Pivoting Strategy: "$PIVOTING_STRAT

echo "1 procs:"
for i in $(seq 1 $N_RUNS); do
    mpirun -np 1 ./quicksort $INPUT output.txt $PIVOTING_STRAT
done


echo "2 procs:"
for i in $(seq 1 $N_RUNS); do
    mpirun -np 2 ./quicksort $INPUT output.txt $PIVOTING_STRAT
done

# Run with 4 procs
echo "4 procs:"
for i in $(seq 1 $N_RUNS); do
    mpirun -np 4 ./quicksort $INPUT output.txt $PIVOTING_STRAT
done


# Run with 8 procs
echo "8 procs:"
for i in $(seq 1 $N_RUNS); do
    mpirun -np 8 ./quicksort $INPUT output.txt $PIVOTING_STRAT
done


# Run with 16 procs
echo "16 procs:"
for i in $(seq 1 $N_RUNS); do
    mpirun -np 16 ./quicksort $INPUT output.txt $PIVOTING_STRAT
done


# Run with 32 procs
echo "32 procs:"
for i in $(seq 1 $N_RUNS); do
    mpirun -np 32 ./quicksort $INPUT output.txt $PIVOTING_STRAT
done
