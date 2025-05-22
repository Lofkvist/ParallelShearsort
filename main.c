#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "header_files/shearsort.h"

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
    int side_len = 0;

	if (my_rank == 0) {
	    read_file(file, &global_matrix, &side_len);
		if (comm_size > side_len) {
		    printf("Number of processes must not exceed side length "
		           "of matrix: ");
		    printf("%d > %d\n", comm_size, side_len);
		    free(global_matrix);
		    MPI_Abort(MPI_COMM_WORLD, 1);
		} else if (side_len == 0) {
		    printf("Error reading file or file is empty\n");
		    free(global_matrix);
		    MPI_Abort(MPI_COMM_WORLD, 1);
	    }
    }
    MPI_Bcast(&side_len, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Each proc is assigned an array of Row_t, sometimes of len 0, sometimes
    // more
    Rows_t *my_rows;
    distribute_from_root(global_matrix, side_len, &my_rows);
    if (my_rank == 0) {
        printf("Original matrix:\n");
        print_matrix(global_matrix, side_len);
    }

    printf("Proc %d:\n", my_rank);
    for (int i = 0; i < my_rows->row_count * side_len; i++) {
        printf("%d ", my_rows->rows[i]);
        if ((i+1) % side_len == 0) {
            printf("\n");
        }
    }
    printf("\n");
    transpose_distributed_matrix(my_rows, side_len);

    gather_to_root(global_matrix, side_len, my_rows);
    if (my_rank == 0) {
        printf("Transposed and gatheredmatrix:\n");
        print_matrix(global_matrix, side_len);
    }

    free(global_matrix);
    MPI_Finalize();
    return 0;
}