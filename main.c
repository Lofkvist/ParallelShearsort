#include "header_files/shearsort.h"
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int compare_matrices(const int *a, const int *b, int side_len);
int *deep_copy_matrix(const int *src, int side_len);
void serial_transpose(int *out, const int *in, int side_len);

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv); // Initialize MPI early
	int comm_size, my_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if (argc != 2) {
		if (my_rank == 0)
			printf("USAGE: ./shearsort <file>\n");
		MPI_Finalize();
		return 1;
	}

	char *file = argv[1];
	int *global_matrix = NULL;
	int *orig_copy = NULL;
	int *expected_transpose = NULL;
	int side_len = 0;

	// Root process reads the file
	if (my_rank == 0) {
		read_file(file, &global_matrix, &side_len);
		if (comm_size > side_len) {
			printf("Number of processes must not exceed side length of matrix: ");
			printf("%d > %d\n", comm_size, side_len);
			free(global_matrix);
			MPI_Abort(MPI_COMM_WORLD, 1);
		} else if (side_len == 0) {
			printf("Error reading file or file is empty\n");
			free(global_matrix);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		// Create copies for verification
		orig_copy = deep_copy_matrix(global_matrix, side_len);
		expected_transpose = malloc(side_len * side_len * sizeof(int));
		serial_transpose(expected_transpose, global_matrix, side_len);
	}

	// Broadcast side_len to all processes
	MPI_Bcast(&side_len, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Allocate global_matrix on non-root processes
	if (my_rank != 0) {
		global_matrix = malloc(side_len * side_len * sizeof(int));
	}

	// Each proc is assigned an array of Row_t
	Rows_t *my_rows;
	distribute_from_root(global_matrix, side_len, &my_rows);

	// Perform distributed transpose
	transpose_distributed_matrix(my_rows, side_len);

	// Gather results back to root
	gather_to_root(global_matrix, side_len, my_rows);
    MPI_Barrier(MPI_COMM_WORLD);

	if (my_rank == 0) {
		printf("Original matrix:\n");
		print_matrix(orig_copy, side_len);
		printf("\n");
	}
    MPI_Barrier(MPI_COMM_WORLD);

	// Verify results on root process
	if (my_rank == 0) {
		int same = compare_matrices(global_matrix, expected_transpose, side_len);
		printf("Transposed matrix:\n");
		print_matrix(global_matrix, side_len);
		printf("\nTranspose correct? %s\n", same ? "YES" : "NO");

		if (!same) {
			printf("\nExpected transpose:\n");
			print_matrix(expected_transpose, side_len);
		}

		// Clean up
		free(orig_copy);
		free(expected_transpose);
	}

	// Clean up
	if (my_rows) {
		if (my_rows->rows) {
			free(my_rows->rows);
		}
		free(my_rows);
	}
	free(global_matrix);

	MPI_Finalize();
	return 0;
}

// Deep copy a 2D matrix (flattened as 1D array)
int *deep_copy_matrix(const int *src, int side_len) {
	if (!src)
		return NULL;

	int *copy = malloc(side_len * side_len * sizeof(int));
	if (!copy)
		return NULL;

	for (int i = 0; i < side_len * side_len; i++) {
		copy[i] = src[i];
	}
	return copy;
}

// Compare two matrices, return 1 if equal, 0 if not
int compare_matrices(const int *a, const int *b, int side_len) {
	if (!a || !b)
		return 0;

	for (int i = 0; i < side_len * side_len; i++) {
		if (a[i] != b[i])
			return 0;
	}
	return 1;
}

// Serial transpose: out = transpose of in
void serial_transpose(int *out, const int *in, int side_len) {
	for (int i = 0; i < side_len; i++) {
		for (int j = 0; j < side_len; j++) {
			out[j * side_len + i] = in[i * side_len + j];
		}
	}
}