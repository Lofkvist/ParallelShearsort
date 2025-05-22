#include "../header_files/shearsort.h"

void read_file(char *file_name, int **matrix, int *side_len) {
	FILE *file = fopen(file_name, "r");
	if (file == NULL) {
		perror("Error opening file");
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

void gather_to_root(int *global_array, int side_len, Rows_t *my_rows) {
	int my_rank, comm_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	int *send_buffer = my_rows->rows; // Use rows as the buffer

	if (my_rank == 0) {
		// Root process receives rows from all processes (including
		// itself)
		for (int i = 0, count = 0; i < side_len; i++) {
			int source = i % comm_size;
			if (source == 0) {
				// Copy directly from own buffer
				memcpy(&global_array[i * side_len], &send_buffer[(count++) * side_len],
				       side_len * sizeof(int));
			} else {
				// Receive row from other process
				MPI_Recv(&global_array[i * side_len], side_len, MPI_INT, source, i, MPI_COMM_WORLD,
				         MPI_STATUS_IGNORE);
			}
		}
	} else {
		// Non-root processes send their rows
		for (int i = 0, count = 0; i < side_len; i++) {
			if (i % comm_size == my_rank) {
				MPI_Send(&send_buffer[count * side_len], side_len, MPI_INT, 0, i, MPI_COMM_WORLD);
				count++;
			}
		}
	}
}

int distribute_from_root(int *global_array, int side_len, Rows_t **my_rows) {
	int my_rank, comm_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	// Count how many rows each process will get (cyclic distribution)
	int local_row_count = 0;
	for (int i = 0; i < side_len; i++) {
		if (i % comm_size == my_rank)
			local_row_count++;
	}

	// Allocate Rows_t struct
	*my_rows = malloc(sizeof(Rows_t));
	(*my_rows)->row_count = local_row_count;

	// Allocate buffer to hold all local rows in a flat 1D array
	int *recv_buffer = malloc(local_row_count * side_len * sizeof(int));
	(*my_rows)->rows = recv_buffer; // point rows to buffer

	if (my_rank == 0) {
		// Root process sends rows to others, keeps own rows
		for (int i = 0; i < side_len; i++) {
			int target = i % comm_size;
			if (target == 0) {
				// Copy directly into own buffer
				memcpy(&recv_buffer[(i / comm_size) * side_len], &global_array[i * side_len],
				       side_len * sizeof(int));
			} else {
				// Send row i to process target
				MPI_Send(&global_array[i * side_len], side_len, MPI_INT, target, i, MPI_COMM_WORLD);
			}
		}
	} else {
		// Other processes receive their rows
		for (int i = 0, count = 0; i < side_len; i++) {
			if (i % comm_size == my_rank) {
				MPI_Recv(&recv_buffer[count * side_len], side_len, MPI_INT, 0, i, MPI_COMM_WORLD,
				         MPI_STATUS_IGNORE);
				count++;
			}
		}
	}

	return local_row_count * side_len;
}

void shearsort(int *matrix, int side_len) {
	int d = (int)ceil(log2(side_len * side_len));
	int k, i;

	for (k = 0; k < d + 1; k++) {
		// Sort rows
		for (i = 0; i < side_len; i++) {
			int *row_start = &matrix[i * side_len];
			if (i % 2 == 0) {
				qsort(row_start, side_len, sizeof(int), compare_asce);
			} else {
				qsort(row_start, side_len, sizeof(int), compare_desc);
			}
		}

		for (i = 0; i < side_len; i++) {
			sort_column(matrix, side_len, i);
		}
	}
}

void transpose_distributed_matrix(Rows_t *local_rows, int side_len) {
	int my_rank, comm_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	int local_row_count = local_rows->row_count;
	int send_count = side_len * local_row_count - local_row_count * local_row_count;
	int recv_count = send_count;
	int *send_buffer = malloc(send_count * sizeof(int));
	int *recv_buffer = malloc(recv_count * sizeof(int));

	int send_counts[comm_size], recv_counts[comm_size];
	int send_displs[comm_size], recv_displs[comm_size];
	int base = side_len / comm_size;
	int procs_with_extra = side_len % comm_size;
	int i, j;

	for (i = 0; i < comm_size; i++) {
		if (i == my_rank) {
			send_counts[i] = recv_counts[i] = 0;
		} else {
			int col_count = base + (i < procs_with_extra ? 1 : 0);
			send_counts[i] = recv_counts[i] = local_row_count * col_count;
		}
	}

	send_displs[0] = 0;
	for (i = 1; i < comm_size; i++) {
		send_displs[i] = send_displs[i - 1] + send_counts[i - 1];
	}
	memcpy(recv_displs, send_displs, sizeof(send_displs)); // Alltoallv symmetry

	int current_idx[comm_size];
	memcpy(current_idx, send_displs, sizeof(send_displs));

	// Fill send buffer
	for (i = 0; i < local_row_count; i++) {
		for (j = 0; j < side_len; j++) {
			int col_owner = j % comm_size;
			if (col_owner == my_rank)
				continue;

			int val = local_rows->rows[i * side_len + j];
			int target_index = current_idx[col_owner]++;
			send_buffer[target_index] = val;

			printf("[Rank %d] Will send value %d at row %d, col %d (global) to rank %d\n", my_rank,
			       val, my_rank + i * comm_size, j, col_owner);
		}
	}

	// Exchange values
	MPI_Alltoallv(send_buffer, send_counts, send_displs, MPI_INT, recv_buffer, recv_counts,
	              recv_displs, MPI_INT, MPI_COMM_WORLD);

	// Print what was received
	printf("[Rank %d] Received values:", my_rank);
	for (i = 0; i < recv_count; i++) {
		printf(" %d", recv_buffer[i]);
	}
	printf("\n");

	for (int i = 0; i < local_row_count; i++) {
		for (int j = 0; j < side_len; j++) {
			int global_row = my_rank + i * comm_size;
			int col_owner = j % comm_size;
			int j_local = (j - col_owner) / comm_size;

			if (col_owner != my_rank)
				continue; // Only handle internal elements

			if ((j - my_rank) % comm_size == 0 && i == j_local) {
				// Diagonal element (i == j_local): no action needed or could print/verify
				continue;
			}

			if ((j - my_rank) % comm_size != 0 || i >= j_local)
				continue; // Avoid duplicate lower-triangle swaps

			int idx1 = i * side_len + j;
			int idx2 = j_local * side_len + global_row;

			printf("[Rank %d] Swapping (%d,%d)<->(%d,%d): %d <-> %d\n", my_rank, global_row, j,
			       my_rank + j_local * comm_size, global_row, local_rows->rows[idx1],
			       local_rows->rows[idx2]);

			int tmp = local_rows->rows[idx1];
			local_rows->rows[idx1] = local_rows->rows[idx2];
			local_rows->rows[idx2] = tmp;
		}
	}

	// Handle the first my_rank items
    int rev = 0;
    int col = 1 + my_rank;
    int rec_idx = 0;

    // Pupulate first my_rank items in rows
    for (i = 0; i < my_rank; i++) {
        for (j = 0; j < local_row_count; j++) {
            int idx = j * side_len + col;
            local_rows->rows[idx] = recv_buffer[rec_idx++];
        }

        col += comm_size;

        if (col >= side_len) {
            rev++;
            col = 1 + my_rank + rev;
        }
    }

    // Rest of row updates
    col = my_rank+1;
    for (i = 0; i < recv_count/local_row_count && rec_idx < recv_count; i++) {
        for (j = 0; j < local_row_count && rec_idx < recv_count; j++) {
            int idx = j * side_len + col;
            local_rows->rows[idx] = recv_buffer[rec_idx++];
        }

        col += comm_size;

        if (col >= side_len) {
            rev++;
            col = 1 + my_rank + rev;
        }
    }

	free(send_buffer);
	free(recv_buffer);
}

void sort_column(int *matrix, int side_len, int col_idx) {
	int *col = malloc(side_len * sizeof(int));
	int i;
	for (i = 0; i < side_len; i++) {
		col[i] = matrix[col_idx + side_len * i];
	}

	qsort(col, side_len, sizeof(int), compare_asce);

	for (i = 0; i < side_len; i++) {
		matrix[col_idx + side_len * i] = col[i];
	}
	free(col);
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
