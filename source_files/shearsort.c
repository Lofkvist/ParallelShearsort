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
    int* row_start;
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

	int rank, r_start;
	if (mpi.rank == 0) {
		memcpy(&global_array[dist.start_row * side_len], my_rows->rows,
		       dist.row_count * side_len * sizeof(int));

		for (rank = 1; rank < mpi.size; rank++) {
			r_start = rank * dist.base + (rank < dist.remainder ? rank : dist.remainder);

			MPI_Recv(&global_array[r_start * side_len], dist.row_count * side_len, MPI_INT, rank, 0,
			         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	} else {
		// Non-root sends its entire block
		MPI_Send(my_rows->rows, dist.row_count * side_len, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
}

// Distrubute rows from root in blocks
int distribute_rows(int *global_array, int side_len, Rows_t **my_rows) {
	mpi_info_t mpi = get_mpi_info();

	// Determine rows per process (block distribution)
	row_dist_t dist = get_row_distribution(side_len, mpi.rank, mpi.size);

	// Allocate Rows_t
	*my_rows = malloc(sizeof(Rows_t));
	int *recv_buffer = malloc(dist.row_count * side_len * sizeof(int));
	(*my_rows)->rows = recv_buffer;

	int rank, r_start;
	if (mpi.rank == 0) {
		for (rank = 0; rank < mpi.size; rank++) {
			r_start = rank * dist.base + (rank < dist.remainder ? rank : dist.remainder);

			if (rank == 0) {
				memcpy(recv_buffer, &global_array[r_start * side_len],
				       dist.row_count * side_len * sizeof(int));
			} else {
				MPI_Send(&global_array[r_start * side_len], dist.row_count * side_len, MPI_INT,
				         rank, 0, MPI_COMM_WORLD);
			}
		}
	} else {
		MPI_Recv(recv_buffer, dist.row_count * side_len, MPI_INT, 0, 0, MPI_COMM_WORLD,
		         MPI_STATUS_IGNORE);
	}

	return dist.row_count * side_len;
}

void transpose_square_matrix(Rows_t *local_rows, int side_len) {
	mpi_info_t mpi = get_mpi_info();
	row_dist_t dist = get_row_distribution(side_len, mpi.rank, mpi.size);

	// Calculate send/recv counts
	int send_counts[mpi.size], send_displs[mpi.size];
	send_displs[0] = 0;
	int total_send = 0;
	int i, j;

	for (i = 0; i < mpi.size; i++) {
		send_counts[i] = dist.row_count * (dist.base + dist.remainder);
		send_displs[i] = send_displs[i - 1] + send_counts[i - 1];
		total_send += send_counts[i];
	}
	send_counts[i] = 0;

    // Send/Recieve buffers
	int *send_buffer = malloc(total_send * sizeof(int));
	int *recv_buffer = malloc(total_send * sizeof(int));

	// send_idx[] keeps track of the current index to be populated in
	// each
	int send_idx[mpi.size];
	memcpy(send_idx, send_displs, mpi.size * sizeof(int));

	// Fill send buffer
	for (i = 0; i < side_len; i++) {
		int dest = (i < (dist.base + 1) * dist.remainder) ? i / (dist.base + 1)
		                                                  : (i - dist.remainder) / dist.base;
		if (dest != mpi.rank) {
			for (j = 0; j < dist.row_count; j++) {
				send_buffer[send_idx[dest]++] = local_rows->rows[j * side_len + i];
			}
		}
	}

	// Exchange values
	MPI_Alltoallv(send_buffer, send_counts, send_displs, MPI_INT, recv_buffer, send_counts,
	              send_displs, MPI_INT, MPI_COMM_WORLD);

	// Interleave received data by row
	int *temp_buffer = malloc(total_send * sizeof(int));
	int temp_idx = 0;
	int row, sender, elems_per_row, offset;
	for (row = 0; row < dist.row_count; row++) {
		for (sender = 0; sender < mpi.size; sender++) {
			if (sender != mpi.rank) {
				elems_per_row = send_counts[sender] / dist.row_count;
				offset = send_displs[sender] + row * elems_per_row;
				memcpy(&temp_buffer[temp_idx], &recv_buffer[offset], elems_per_row * sizeof(int));
				temp_idx += elems_per_row;
			}
		}
	}

	// Exchange local element within proc
	int global_i, global_j, idx1, idx2, temp;
	for (i = 0; i < dist.row_count; i++) {
		global_i = dist.start_row + i;
		for (j = i + 1; j < dist.row_count; j++) {
			global_j = dist.start_row + j;
			idx1 = i * side_len + global_j;
			idx2 = j * side_len + global_i;

			temp = local_rows->rows[idx1];
			local_rows->rows[idx1] = local_rows->rows[idx2];
			local_rows->rows[idx2] = temp;
		}
	}

	// Place received elements
	int counter = 0;
	for (i = 0; i < dist.row_count; i++) {
		for (j = 0; j < side_len; j++) {
			if (j >= dist.start_row && j < dist.start_row + dist.row_count) {
                // Skip these, already transposed within proc
				continue;
			}
			local_rows->rows[i * side_len + j] = temp_buffer[counter++];
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