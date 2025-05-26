// shearsort.c
// Parallel shear sort and helper functions using MPI
// Author: Carl Löfkvist
// Date: 2024-05-26

#include "../header_files/shearsort.h"

static row_dist_t get_row_distribution(int side_len, int rank, int comm_size) {
	row_dist_t dist;
	dist.base = side_len / comm_size;
	dist.remainder = side_len % comm_size;
	dist.start_row = rank * dist.base + (rank < dist.remainder ? rank : dist.remainder);
	dist.row_count = dist.base + (rank < dist.remainder ? 1 : 0);
	return dist;
}

mpi_info_t get_mpi_info() {
	mpi_info_t info;
	MPI_Comm_rank(MPI_COMM_WORLD, &info.rank);
	MPI_Comm_size(MPI_COMM_WORLD, &info.size);
	return info;
}

void shearsort(Rows_t *local_rows, int side_len) {
	mpi_info_t mpi = get_mpi_info();

	// Instead of 4 calculations, just:
	row_dist_t dist = get_row_distribution(side_len, mpi.rank, mpi.size);

	int k, i, global_row;
	int *row_start;
	const int d = (int)ceil(log2(side_len * side_len));
	for (k = 0; k < d + 1; k++) {

		// ====== Sort Rows ======
		for (i = 0; i < dist.row_count; i++) {
			global_row = dist.start_row + i;
			row_start = &local_rows->rows[i * side_len];
			if (global_row % 2 == 0) {
				qsort(row_start, side_len, sizeof(int), compare_asce);
			} else {
				qsort(row_start, side_len, sizeof(int), compare_desc);
			}
		}

		// ====== Sort Columns ======
		// Transpose (rows => cols)
		transpose_square_matrix(local_rows, side_len);

		// Sort Columns
		for (i = 0; i < dist.row_count; i++) {
			row_start = &local_rows->rows[i * side_len];
			qsort(row_start, side_len, sizeof(int), compare_asce);
		}

		// Transpose (cols => rows)
		transpose_square_matrix(local_rows, side_len);
	}
}

// Populate matrix with contents from file
void read_file(char *file_name, int **matrix, int *side_len) {
	FILE *file = fopen(file_name, "r");
	if (file == NULL) {
		perror("Error opening file");
		fclose(file);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (fscanf(file, "%d", side_len) != 1) {
		fprintf(stderr, "Error reading element count\n");
		fclose(file);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	*matrix = malloc((*side_len) * (*side_len) * sizeof(int));
	if (*matrix == NULL) {
		fprintf(stderr, "Memory allocation failed\n");
		fclose(file);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	for (int i = 0; i < (*side_len) * (*side_len); i++) {
		if (fscanf(file, "%d", &(*matrix)[i]) != 1) {
			fprintf(stderr, "Error reading matrix element %d\n", i);
			fclose(file);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}

	fclose(file);
}

// Gather rows of all processes to root process
void gather_rows(int *global_array, int side_len, Rows_t *my_rows) {
    mpi_info_t mpi = get_mpi_info();
    row_dist_t dist = get_row_distribution(side_len, mpi.rank, mpi.size);

    if (mpi.rank == 0) {
        // Copy own rows into global array
        memcpy(&global_array[dist.start_row * side_len], my_rows->rows,
               dist.row_count * side_len * sizeof(int));

        // Receive rows from other processes
        for (int rank = 1; rank < mpi.size; rank++) {
            int r_start = rank * dist.base + (rank < dist.remainder ? rank : dist.remainder);
            int recv_count = ((rank < dist.remainder) ? dist.base + 1 : dist.base) * side_len;

            MPI_Status status;
            MPI_Recv(&global_array[r_start * side_len], recv_count, MPI_INT, rank, 0,
                     MPI_COMM_WORLD, &status);

            int actual_count;
            MPI_Get_count(&status, MPI_INT, &actual_count);
            if (actual_count != recv_count) {
                fprintf(stderr, "Rank 0: MPI_Recv mismatch from rank %d! Expected %d, got %d\n",
                        rank, recv_count, actual_count);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    } else {
        // Non-root processes send their block
        int send_count = dist.row_count * side_len;
        MPI_Send(my_rows->rows, send_count, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}


// Distribute rows from root in blocks
int distribute_rows(int *global_array, int side_len, Rows_t **my_rows) {
    mpi_info_t mpi = get_mpi_info();

    // Determine rows per process (block distribution)
    row_dist_t dist = get_row_distribution(side_len, mpi.rank, mpi.size);

    int count = dist.row_count * side_len;

    // Allocate Rows_t
    *my_rows = malloc(sizeof(Rows_t));
    if (*my_rows == NULL) {
        fprintf(stderr, "Rank %d: Failed to allocate Rows_t\n", mpi.rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int *recv_buffer = malloc(count * sizeof(int));
    if (recv_buffer == NULL) {
        fprintf(stderr, "Rank %d: Failed to allocate recv_buffer\n", mpi.rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    (*my_rows)->rows = recv_buffer;

    if (mpi.rank == 0) {
        for (int rank = 0; rank < mpi.size; rank++) {
            int r_start = rank * dist.base + (rank < dist.remainder ? rank : dist.remainder);
            int send_rows = (rank < dist.remainder) ? dist.base + 1 : dist.base;
            int send_count = send_rows * side_len;

            if (rank == 0) {
                memcpy(recv_buffer, &global_array[r_start * side_len],
                       send_count * sizeof(int));
            } else {
                MPI_Send(&global_array[r_start * side_len], send_count, MPI_INT,
                         rank, 0, MPI_COMM_WORLD);
            }
        }
    } else {
        MPI_Status status;

        MPI_Recv(recv_buffer, count, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

        int actual_count;
        MPI_Get_count(&status, MPI_INT, &actual_count);
        if (actual_count != count) {
            fprintf(stderr, "Rank %d: MPI_Recv mismatch! Expected %d, got %d\n",
                    mpi.rank, count, actual_count);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    return count;
}



void transpose_square_matrix(Rows_t *local_rows, int side_len) {
	mpi_info_t mpi = get_mpi_info();
	row_dist_t dist = get_row_distribution(side_len, mpi.rank, mpi.size);

	// Calculate send/recv counts
	int send_counts[mpi.size], send_displs[mpi.size];
	send_displs[0] = 0;
	int total_send = 0;
	int i, rank, col, row;

    for (rank = 0; rank < mpi.size; rank++) {
        row_dist_t rank_dist = get_row_distribution(side_len, rank, mpi.size);
        send_counts[rank] = dist.row_count * rank_dist.row_count;

        if (rank > 0) {  // Add this check
            send_displs[rank] = send_displs[rank - 1] + send_counts[rank - 1];
        }
        total_send += send_counts[rank];
    }

	// Send/Recieve buffers
	int *send_buffer = malloc(total_send * sizeof(int));
	int *recv_buffer = malloc(total_send * sizeof(int));

    assert(send_buffer != NULL && recv_buffer != NULL);


	// send_idx[] keeps track of the current index to be populated in
	// each
	int send_idx[mpi.size];
	memcpy(send_idx, send_displs, mpi.size * sizeof(int));

	// Fill send buffer
	for (col = 0; col < side_len; col++) {
		int dest = (col < (dist.base + 1) * dist.remainder) ? col / (dist.base + 1)
		                                                  : (col - dist.remainder) / dist.base;
		if (dest != mpi.rank) {
			for (row = 0; row < dist.row_count; row++) {
				send_buffer[send_idx[dest]++] = local_rows->rows[row * side_len + col];
			}
		}
	}

	// Exchange values
	MPI_Alltoallv(send_buffer, send_counts, send_displs, MPI_INT, recv_buffer, send_counts,
	              send_displs, MPI_INT, MPI_COMM_WORLD);

	// Interleave received data by row
	int *temp_buffer = malloc(total_send * sizeof(int));
	int temp_idx = 0;
	int elems_per_row, offset;
	for (row = 0; row < dist.row_count; row++) {
		for (rank = 0; rank < mpi.size; rank++) {
			if (rank != mpi.rank) {
				elems_per_row = send_counts[rank] / dist.row_count;
				offset = send_displs[rank] + row * elems_per_row;
				memcpy(&temp_buffer[temp_idx], &recv_buffer[offset], elems_per_row * sizeof(int));
				temp_idx += elems_per_row;
			}
		}
	}

	// Perform transposition
	int counter = 0;
	int global_i, global_j, idx1, idx2, temp;
	for (row = 0; row < dist.row_count; row++) {
		// === Exchange elements locally ===
        global_i = dist.start_row + row;
		for (i = row + 1; i < dist.row_count; i++) {
			global_j = dist.start_row + i;
			idx1 = row * side_len + global_j;
			idx2 = i * side_len + global_i;

			temp = local_rows->rows[idx1];
			local_rows->rows[idx1] = local_rows->rows[idx2];
			local_rows->rows[idx2] = temp;
		}

        // === Place recieved elements ===
		for (i = 0; i < side_len; i++) {
			if (i >= dist.start_row && i < dist.start_row + dist.row_count) {
				// Skip these, already transposed within proc
				continue;
			}
			local_rows->rows[row * side_len + i] = temp_buffer[counter++];
		}
	}

	free(send_buffer);
	free(temp_buffer);
	free(recv_buffer);
}

void print_matrix(int *matrix, int side_len) {
	int i;
	for (i = 0; i < side_len * side_len; i++) {
		printf("%3d ", matrix[i]);
		if ((i + 1) % side_len == 0) {
			printf("\n");
		}
	}
	printf("\n");
}

int compare_asce(const void *a, const void *b) { return (*(int *)a - *(int *)b); }
int compare_desc(const void *a, const void *b) { return (*(int *)b - *(int *)a); }

int is_shearsorted(int *matrix, int side_len) {
	int i, j;

	// Check each row
	for (i = 0; i < side_len; i++) {
		for (j = 0; j < side_len - 1; j++) {
			int a = matrix[i * side_len + j];
			int b = matrix[i * side_len + j + 1];
			if (i % 2 == 0) { // even-indexed row: left to right
				if (a > b)
					return 0;
			} else { // odd-indexed row: right to left
				if (a < b)
					return 0;
			}
		}
	}

	// Check each column (top to bottom)
	for (j = 0; j < side_len; j++) {
		for (i = 0; i < side_len - 1; i++) {
			int a = matrix[i * side_len + j];
			int b = matrix[(i + 1) * side_len + j];
			if (a > b)
				return 0;
		}
	}

	return 1; // Passed all checks
}