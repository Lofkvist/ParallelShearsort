#!/bin/bash

#SBATCH -A uppmax2025-2-247

#SBATCH -p node

#SBATCH -N 2

# Request 30 minutes of runtime for the job
#SBATCH -t 00:30:00

ml gcc openmpi
echo "=== STRONG SCALABILITY TEST ==="

STRONG_NS=(1000 2000 3000 4000)
CORES=(1 2 4 8 16)
OUTPUT_FILE="strong_new_results.txt"

> "$OUTPUT_FILE"
echo "Strong Scalability Results" >> "$OUTPUT_FILE"

for N in "${STRONG_NS[@]}"; do
  elements=$((N * N))
  input_file="../inputs/input${elements}.txt"
  if [ ! -f "$input_file" ]; then
    echo "Warning: Input file $input_file not found, skipping N=$N"
    continue
  fi

  echo "=== N=$N (Elements=$elements) ===" >> "$OUTPUT_FILE"
  for cores in "${CORES[@]}"; do
    echo "Running with $cores cores for N=$N"
    echo "== $cores cores ==" >> "$OUTPUT_FILE"
    mpirun -n "$cores" ../shearsort "$input_file" >> "$OUTPUT_FILE"
  done
done

echo ""
echo "=== WEAK SCALABILITY TEST ==="

declare -A WEAK_CONFIG=(
  [1]=1000
  [4]=2000
  [9]=3000
  [16]=4000
)

OUTPUT_FILE_WEAK="weak_new_results.txt"
> "$OUTPUT_FILE_WEAK"
echo "Weak Scalability Results - Constant work per core" >> "$OUTPUT_FILE_WEAK"

for cores in "${!WEAK_CONFIG[@]}"; do
  N=${WEAK_CONFIG[$cores]}
  elements=$((N * N))
  input_file="../inputs/input${elements}.txt"

  if [ ! -f "$input_file" ]; then
    echo "Warning: Input file $input_file not found, skipping cores=$cores"
    continue
  fi

  echo "Running weak scaling with $cores cores, N=$N"
  echo "== $cores cores, N=$N ==" >> "$OUTPUT_FILE_WEAK"
  mpirun -n "$cores" ../shearsort "$input_file" >> "$OUTPUT_FILE_WEAK"
done

echo "=== SCALABILITY TESTING COMPLETE ==="

