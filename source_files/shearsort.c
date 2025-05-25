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

	// Determine rows per process (block distribution)
	int base_rows = side_len / comm_size;
	int remainder = side_len % comm_size;
	int start_row = my_rank * base_rows + (my_rank < remainder ? my_rank : remainder);
	int local_row_count = base_rows + (my_rank < remainder ? 1 : 0);

	if (my_rank == 0) {
		// Root process receives blocks from all ranks
		memcpy(&global_array[start_row * side_len], my_rows->rows,
		       local_row_count * side_len * sizeof(int));

		for (int rank = 1; rank < comm_size; rank++) {
			int r_base = side_len / comm_size;
			int r_extra = rank < remainder ? 1 : 0;
			int r_count = r_base + r_extra;
			int r_start = rank * r_base + (rank < remainder ? rank : remainder);

			MPI_Recv(&global_array[r_start * side_len], r_count * side_len, MPI_INT, rank, 0,
			         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	} else {
		// Non-root sends its entire block
		MPI_Send(my_rows->rows, local_row_count * side_len, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
}

int distribute_from_root_block(int *global_array, int side_len, Rows_t **my_rows) {
	int my_rank, comm_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	// Determine rows per process (block distribution)
	int base_rows = side_len / comm_size;
	int remainder = side_len % comm_size;
	int start_row = my_rank * base_rows + (my_rank < remainder ? my_rank : remainder);
	int local_row_count = base_rows + (my_rank < remainder ? 1 : 0);

	// Allocate Rows_t
	*my_rows = malloc(sizeof(Rows_t));
	(*my_rows)->row_count = local_row_count;
	int *recv_buffer = malloc(local_row_count * side_len * sizeof(int));
	(*my_rows)->rows = recv_buffer;

	if (my_rank == 0) {
		for (int rank = 0; rank < comm_size; rank++) {
			int r_base = side_len / comm_size;
			int r_extra = rank < remainder ? 1 : 0;
			int r_count = r_base + r_extra;
			int r_start = rank * r_base + (rank < remainder ? rank : remainder);

			if (rank == 0) {
				memcpy(recv_buffer, &global_array[r_start * side_len],
				       r_count * side_len * sizeof(int));
			} else {
				MPI_Send(&global_array[r_start * side_len], r_count * side_len, MPI_INT, rank, 0,
				         MPI_COMM_WORLD);
			}
		}
	} else {
		MPI_Recv(recv_buffer, local_row_count * side_len, MPI_INT, 0, 0, MPI_COMM_WORLD,
		         MPI_STATUS_IGNORE);
	}

	return local_row_count * side_len;
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

void shearsort(Rows_t *local_rows, int side_len) {
    int d = (int)ceil(log2(side_len * side_len));
    int my_rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    
    // Calculate global starting row index for this process
    int base = side_len / comm_size;
    int remainder = side_len % comm_size;
    int start_row = my_rank * base + (my_rank < remainder ? my_rank : remainder);
    
    for (int k = 0; k < d + 1; k++) {
        // Sort rows with snake pattern
        for (int i = 0; i < local_rows->row_count; i++) {
            int global_row = start_row + i;
            int *row_start = &local_rows->rows[i * side_len];
            if (global_row % 2 == 0) {
                qsort(row_start, side_len, sizeof(int), compare_asce);
            } else {
                qsort(row_start, side_len, sizeof(int), compare_desc);
            }
        }
        
        // Transpose distributed matrix
        transpose_distributed_matrix(local_rows, side_len);
        
        // After transpose, check what the actual row width should be
        // This might not be side_len anymore!
        int current_row_width = side_len; // You may need to calculate this differently
        
        for (int i = 0; i < local_rows->row_count; i++) {
            int *row_start = &local_rows->rows[i * current_row_width];
            qsort(row_start, current_row_width, sizeof(int), compare_asce);
        }
        
        // Transpose back
        transpose_distributed_matrix(local_rows, side_len);
        
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

	// Fill send buffer - pack data COLUMN-WISE by destination process
	int send_idx[comm_size];
	memcpy(send_idx, send_displs, sizeof(send_displs));

	int remainder = side_len % comm_size;

	for (int col = 0; col < side_len; col++) {
		int dest;
		if (col < (base + 1) * remainder)
			dest = col / (base + 1);
		else
			dest = (col - remainder) / base;

		if (dest == my_rank)
			continue;

		for (int row = 0; row < local_row_count; row++) {
			int value = local_rows->rows[row * side_len + col];
			send_buffer[send_idx[dest]++] = value;
		}
	}

	// Exchange values
	MPI_Alltoallv(send_buffer, send_counts, send_displs, MPI_INT, recv_buffer, recv_counts,
	              recv_displs, MPI_INT, MPI_COMM_WORLD);

	// Interleave received data by row
	int *temp_buffer = malloc(recv_count * sizeof(int));
	int temp_idx = 0;

	for (int row = 0; row < local_row_count; row++) {
		for (int sender = 0; sender < comm_size; sender++) {
			if (sender == my_rank)
				continue;

			int elems_per_row = recv_counts[sender] / local_row_count;
			int offset = recv_displs[sender] + row * elems_per_row;

			memcpy(&temp_buffer[temp_idx], &recv_buffer[offset], elems_per_row * sizeof(int));
			temp_idx += elems_per_row;
		}
	}

	memcpy(recv_buffer, temp_buffer, recv_count * sizeof(int));
	free(temp_buffer);

	// Transpose within the local block (swap elements across the diagonal)
	// Only swap within the square block owned by this process
	int base_rows = side_len / comm_size;
	int start_row = my_rank * base_rows + (my_rank < remainder ? my_rank : remainder);

	// Transpose within the local block, considering the global positions
	for (int i = 0; i < local_row_count; i++) {
		for (int j = i + 1; j < local_row_count; j++) {
			int global_i = start_row + i;
			int global_j = start_row + j;

			// Transpose only if the full (i,j) and (j,i) are in the local block
			int idx1 = i * side_len + global_j;
			int idx2 = j * side_len + global_i;

			int temp = local_rows->rows[idx1];
			local_rows->rows[idx1] = local_rows->rows[idx2];
			local_rows->rows[idx2] = temp;
		}
	}

	// Add the recieved elements in the row
	int counter = 0;
	for (int i = 0; i < local_row_count; i++) {
		int diag_col = start_row + i;
		int local_diag_idx = i * side_len + diag_col;

		for (int j = 0; j < side_len; j++) {
			// Skip elements before and after the diagonal in the local square block
			if (j >= start_row && j < start_row + local_row_count) {
				if (j >= diag_col - i && j <= diag_col + (local_row_count - 1 - i)) {
					continue; // skip elements surrounding the diagonal
				}
			}

			int local_idx = i * side_len + j;
			local_rows->rows[local_idx] = recv_buffer[counter++];
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
