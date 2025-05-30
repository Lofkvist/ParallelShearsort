# ParallelShearsort

A parallel implementation of the **Shear Sort** algorithm using MPI for distributed memory systems.

## Overview

Shear Sort is a parallel sorting algorithm designed for sorting elements arranged in a 2D square matrix. The matrix undergoes alternating phases of row-wise and column-wise sorting, with row sorting alternating between ascending and descending order on even and odd rows respectively, and column sorting always ascending. The process repeats for $\lceil \log_2(N) \rceil + 1$ iterations, where $N$ is the total number of elements.

Example matrix input:

```
16 2  3  13
5  11 10 8
9  7  6  12
4  14 15 1
```

and resulting output:

```
1  2  3  4
8  7  6  5
9  10 11 12
16 15 14 13
```

This project distributes the matrix rows among MPI processes, allowing the sorting phases and matrix transpositions to be performed in parallel with MPI communication.

## Features

- Efficient row distribution among MPI processes with load balancing.
- Parallel row-wise sorting: even rows ascending, odd rows descending.
- Parallel column-wise sorting via distributed matrix transpose.
- Verification of sorted matrix correctness.
- Timing of the parallel execution.

## Usage

Build the project using `make` and run the executable:

```bash
mpirun -n <num_procs> ./shearsort <input_file>
```

- `<num_procs>` is the number of processes used
- `<input_file>` is a text file with one line, where the first number `n` denotes the side length of the matrix, followed by `n^2` integers to be sorted (inputs are provided in directory `inputs` with varying `n`)

## Example input file

```
4 16 2 3 13 5 11 10 8 9 7 6 12 4 14 15 1
```
The program writes the execution time of the slowest process to standard output as a `float`
## Requirements

- MPI (e.g., OpenMPI or MPICH)
- C compiler (e.g., `gcc`)
- `make`

Author: Carl Löfkvist