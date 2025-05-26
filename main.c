// main.c
// Main function for parallel shear sort using MPI
// Author: Carl Löfkvist
// Date: 2024-05-26

#include "header_files/shearsort.h"
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int compare_matrices(const int *a, const int *b, int side_len);
int *deep_copy_matrix(const int *src, int side_len);
void serial_transpose(int *out, const int *in, int side_len);
int is_shearsorted(int *matrix, int side_len);

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    mpi_info_t mpi = get_mpi_info();

	if (argc != 2) {
		if (mpi.rank == 0)
			printf("USAGE: ./shearsort <file>\n");
		MPI_Finalize();
		return 1;
	}

	char *file = argv[1];
	int *global_matrix = NULL;
	int side_len;

	// Root process reads the file
	if (mpi.rank == 0) {
		read_file(file, &global_matrix, &side_len);
		if (mpi.size > side_len) {
            printf("============================================================\n");
			printf("Number of processes exceeded side length of matrix: %d > %d\n", mpi.size, side_len);
            printf("============================================================\n");
			free(global_matrix);
			MPI_Abort(MPI_COMM_WORLD, 1);
		} else if (side_len == 0) {
            printf("====================================\n");
			printf("Error reading file or file is empty\n");
            printf("====================================\n");
			free(global_matrix);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}

	// Broadcast side_len to all processes
	MPI_Bcast(&side_len, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Allocate global_matrix on non-root processes
	if (mpi.rank != 0) {
		global_matrix = malloc(side_len * side_len * sizeof(int));
	}

	double start, end, elapsed, max_elapsed;

	// Local 
	Rows_t *local_rows;

	distribute_rows(global_matrix, side_len, &local_rows);

	// Perform distributed transpose
	start = MPI_Wtime();
	shearsort(local_rows, side_len);
	end = MPI_Wtime();
	elapsed = end - start;
	// Gather results back to root
	gather_rows(global_matrix, side_len, local_rows);

	MPI_Reduce(&elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (mpi.rank == 0) {
		int sheadsorted = is_shearsorted(global_matrix, side_len);
		if (sheadsorted) {
            printf("%f\n", max_elapsed);
		} else {
			printf("Matrix not sorted!\n");
            return 1;
        }
	}


	// Clean up
	if (local_rows) {
		if (local_rows->rows) {
			free(local_rows->rows);
		}
		free(local_rows);
	}
	free(global_matrix);

	MPI_Finalize();
	return 0;
}
