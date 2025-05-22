#ifndef SHEADSORT_H
#define SHEADSORT_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

typedef struct {
    int* rows;
    int row_count;
} Rows_t;



// Read matrix from file
void read_file(char *file_name, int **matrix, int *side_len);

// Shearsort algorithm
void shearsort(int *matrix, int side_len);

// Shearsort algorithm
void shearsort_par(Rows_t* rows, int side_len);

int distribute_from_root(int* global_array, int side_len, Rows_t** my_rows);

// Transpose the distributed matrix
void transpose_distributed_matrix(Rows_t *local_rows, int side_len);

// Gather all data to root
void gather_to_root(int *global_array, int side_len, Rows_t *my_rows);

// Sort a single column in ascending order
void sort_column(int *matrix, int side_len, int col_idx);

// Print the matrix to stdout
void print_matrix(int *matrix, int side_len);

// Comparison functions for qsort
int compare_asce(const void *a, const void *b);
int compare_desc(const void *a, const void *b);

#endif // SHEADSORT_H
