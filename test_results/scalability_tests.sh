#!/bin/bash

# Scalability Test Script (Strong and Weak) up to 32 cores

echo "=== STRONG SCALABILITY TEST ==="

# Configuration for Strong Scaling
FIXED_N=1000 # Matrix side length (N*N elements)
CORES=(1 2 4 8 16)

INPUT_FILE="../inputs/input$((FIXED_N * FIXED_N)).txt"
OUTPUT_FILE="strong_scalability_results.txt"

# Clear previous results
> "$OUTPUT_FILE"

echo "Strong Scalability Results - Fixed N=$FIXED_N (${FIXED_N}*${FIXED_N} elements)" >> "$OUTPUT_FILE"
echo "Cores,ExecutionTime(seconds)" >> "$OUTPUT_FILE"

if [ ! -f "$INPUT_FILE" ]; then
  echo "Error: Input file $INPUT_FILE not found. Please generate it before running."
  exit 1
fi

for cores in "${CORES[@]}"; do
  echo "Testing strong scaling with $cores cores..."
  execution_time=$(mpirun -n "$cores" ../shearsort "$INPUT_FILE")
  if [ $? -ne 0 ]; then
    echo "Run failed for $cores cores during strong scaling. Skipping..."
    continue
  fi
  echo "$cores,$execution_time" >> "$OUTPUT_FILE"
  echo " $cores cores: $execution_time seconds"
done

echo "Strong scalability results saved to $OUTPUT_FILE"
echo ""

#############################################

echo "=== WEAK SCALABILITY TEST ==="

# Configuration for Weak Scaling (constant work per core)
declare -A WEAK_CONFIG=(
  [1]=1000
  [4]=2000
  [9]=3000
  [16]=4000
)

OUTPUT_FILE_WEAK="weak_scalability_results.txt"
> "$OUTPUT_FILE_WEAK"

echo "Weak Scalability Results - Constant work per core" >> "$OUTPUT_FILE_WEAK"
echo "Cores,N,Elements,ExecutionTime(seconds)" >> "$OUTPUT_FILE_WEAK"

for cores in "${!WEAK_CONFIG[@]}"; do
  N=${WEAK_CONFIG[$cores]}
  elements=$((N * N))
  input_file="../inputs/input$elements.txt"

  echo "Testing weak scaling with $cores cores, N=$N (${elements} elements)..."

  if [ ! -f "$input_file" ]; then
    echo "Warning: Input file $input_file not found. Skipping..."
    continue
  fi

  execution_time=$(mpirun -n "$cores" ../shearsort "$input_file")
  if [ $? -ne 0 ]; then
    echo "Run failed for $cores cores during weak scaling. Skipping..."
    continue
  fi

  echo "$cores,$N,$elements,$execution_time" >> "$OUTPUT_FILE_WEAK"
  echo " $cores cores, N=$N: $execution_time seconds"
done

echo "Weak scalability results saved to $OUTPUT_FILE_WEAK"
echo ""
echo "=== SCALABILITY TESTING COMPLETE ==="
