#ifndef SHEADSORT_H
#define SHEADSORT_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

typedef struct {
   int* rows;
} Rows_t;

typedef struct {
   int rank;
   int size;
} mpi_info_t;

typedef struct {
   int start_row;
   int row_count;
   int base;
   int remainder;
} row_dist_t;

/**
* Read matrix from file
* @param file_name Name of the input file containing matrix data
* @param matrix Pointer to matrix array (allocated by function)
* @param side_len Pointer to store the side length of the square matrix
*/
void read_file(char *file_name, int **matrix, int *side_len);

/**
* Get MPI rank and communicator size information
* @return mpi_info_t structure containing rank and size
*/
mpi_info_t get_mpi_info();

/**
* Distribute rows from global matrix to local processes
* @param global_array The complete matrix array (only valid on root process)
* @param side_len Side length of the square matrix
* @param my_rows Pointer to local rows structure (allocated by function)
* @return Total number of elements in local rows
*/
int distribute_rows(int *global_array, int side_len, Rows_t **my_rows);

/**
* Sequential shearsort algorithm
* @param local_rows Local rows structure containing process's portion of matrix
* @param side_len Side length of the square matrix
*/
void shearsort(Rows_t *local_rows, int side_len);

/**
* Parallel shearsort algorithm
* @param rows Local rows structure containing process's portion of matrix
* @param side_len Side length of the square matrix
*/
void shearsort_par(Rows_t* rows, int side_len);

/**
* Transpose the distributed matrix across all processes
* @param local_rows Local rows structure containing process's portion of matrix
* @param side_len Side length of the square matrix
*/
void transpose_square_matrix(Rows_t *local_rows, int side_len);

/**
* Gather all distributed rows to root process
* @param global_matrix Buffer to store complete matrix (only used on root)
* @param side_len Side length of the square matrix
* @param my_rows Local rows structure containing process's portion of matrix
*/
void gather_rows(int *global_matrix, int side_len, Rows_t *my_rows);

/**
* Sort a single column in ascending order
* @param matrix The matrix array
* @param side_len Side length of the square matrix
* @param col_idx Index of the column to sort
*/
void sort_column(int *matrix, int side_len, int col_idx);

/**
* Print the matrix to stdout
* @param matrix The matrix array to print
* @param side_len Side length of the square matrix
*/
void print_matrix(int *matrix, int side_len);

/**
* Comparison function for ascending order sorting
* @param a Pointer to first element
* @param b Pointer to second element
* @return Negative if a < b, positive if a > b, zero if equal
*/
int compare_asce(const void *a, const void *b);

/**
* Comparison function for descending order sorting
* @param a Pointer to first element
* @param b Pointer to second element
* @return Positive if a < b, negative if a > b, zero if equal
*/
int compare_desc(const void *a, const void *b);

/**
* Check if matrix is sorted according to shearsort pattern
* @param matrix The matrix array to check
* @param side_len Side length of the square matrix
* @return 1 if matrix is shearsorted, 0 otherwise
*/
int is_shearsorted(int *matrix, int side_len);

#endif // SHEADSORT_H