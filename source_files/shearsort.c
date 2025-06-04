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

void shearsort(int *local_rows, int side_len) {
	mpi_info_t mpi = get_mpi_info();
	row_dist_t dist = get_row_distribution(side_len, mpi.rank, mpi.size);

	int k, i, global_row;
	int *row_start;
	const int d = (int)ceil(log2(side_len));
	for (k = 0; k < d + 1; k++) {
		// ====== Row Sort Phase ======
		for (i = 0; i < dist.row_count; i++) {
			global_row = dist.start_row + i;
			row_start = &local_rows[i * side_len];
			if (global_row % 2 == 0) {
				qsort(row_start, side_len, sizeof(int), compare_asce);
			} else {
				qsort(row_start, side_len, sizeof(int), compare_desc);
			}
		}
		// ====== Column Sort Phase ======
		// Transpose (rows => cols)
		transpose_square_matrix(local_rows, side_len);

		// Sort Columns
		for (i = 0; i < dist.row_count; i++) {
			row_start = &local_rows[i * side_len];
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

void gather_rows(int *global_array, int side_len, int *my_rows) {
    mpi_info_t mpi = get_mpi_info();
    row_dist_t dist = get_row_distribution(side_len, mpi.rank, mpi.size);

    if (mpi.rank == 0) {
        // Copy own rows into global array
        memcpy(&global_array[dist.start_row * side_len], my_rows,
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
        MPI_Send(my_rows, send_count, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}


int distribute_rows(int *global_array, int side_len, int **my_rows) {
    mpi_info_t mpi = get_mpi_info();
    row_dist_t dist = get_row_distribution(side_len, mpi.rank, mpi.size);

    int count = dist.row_count * side_len;

    int *recv_buffer = malloc(count * sizeof(int));
    if (recv_buffer == NULL) {
        fprintf(stderr, "Rank %d: Failed to allocate recv_buffer\n", mpi.rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    *my_rows = recv_buffer;

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



void transpose_square_matrix(int *local_rows, int side_len) {
	mpi_info_t mpi = get_mpi_info();
	row_dist_t dist = get_row_distribution(side_len, mpi.rank, mpi.size);

	int send_counts[mpi.size], send_displs[mpi.size];
	send_displs[0] = 0;
	int total_send = 0;

	for (int rank = 0; rank < mpi.size; rank++) {
		row_dist_t rank_dist = get_row_distribution(side_len, rank, mpi.size);
		send_counts[rank] = dist.row_count * rank_dist.row_count;
		if (rank > 0) {
			send_displs[rank] = send_displs[rank - 1] + send_counts[rank - 1];
		}
		total_send += send_counts[rank];
	}

	int *send_buffer = malloc(total_send * sizeof(int));
	int *recv_buffer = malloc(total_send * sizeof(int));
	assert(send_buffer != NULL && recv_buffer != NULL);

	int send_idx[mpi.size];
	memcpy(send_idx, send_displs, mpi.size * sizeof(int));

	for (int col = 0; col < side_len; col++) {
		int dest = (col < (dist.base + 1) * dist.remainder) ? col / (dist.base + 1)
		                                                    : (col - dist.remainder) / dist.base;
		if (dest != mpi.rank) {
			for (int row = 0; row < dist.row_count; row++) {
				send_buffer[send_idx[dest]++] = local_rows[row * side_len + col];
			}
		}
	}

	MPI_Alltoallv(send_buffer, send_counts, send_displs, MPI_INT,
	              recv_buffer, send_counts, send_displs, MPI_INT,
	              MPI_COMM_WORLD);

	int *temp_buffer = malloc(total_send * sizeof(int));
	int temp_idx = 0;
	for (int row = 0; row < dist.row_count; row++) {
		for (int rank = 0; rank < mpi.size; rank++) {
			if (rank != mpi.rank) {
				int elems_per_row = send_counts[rank] / dist.row_count;
				int offset = send_displs[rank] + row * elems_per_row;
				memcpy(&temp_buffer[temp_idx], &recv_buffer[offset], elems_per_row * sizeof(int));
				temp_idx += elems_per_row;
			}
		}
	}

	int counter = 0;
	for (int row = 0; row < dist.row_count; row++) {
		int global_i = dist.start_row + row;
		for (int i = row + 1; i < dist.row_count; i++) {
			int global_j = dist.start_row + i;
			int idx1 = row * side_len + global_j;
			int idx2 = i * side_len + global_i;
			int temp = local_rows[idx1];
			local_rows[idx1] = local_rows[idx2];
			local_rows[idx2] = temp;
		}

		for (int i = 0; i < side_len; i++) {
			if (i >= dist.start_row && i < dist.start_row + dist.row_count) continue;
			local_rows[row * side_len + i] = temp_buffer[counter++];
		}
	}

	free(send_buffer);
	free(recv_buffer);
	free(temp_buffer);
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