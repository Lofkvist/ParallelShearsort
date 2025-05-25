#include <stdio.h>

int main() {
    int my_rank = 0;
    int comm_size = 3;
    int side_len = 7;
    int local_row_count = 3;
    int recv_count = 12;

    int i, j, k;
    // Turn this into
    int vals[] = {10, 38, 13, 14, 11, 45, 18, 35, 48, 22, 42, 41};
    // This
    //int vals[] = {10, 13, 38, 14, 11, 18, 45, 35, 48, 42, 22, 41};
    int new_vals[recv_count];
    int rev = 0;
    int col = 1 + my_rank;
    
    int elems_per_row = recv_count / local_row_count;
    int idx = 0;
    for (i = 0; i < local_row_count; i++) {
        
        int offset = i * elems_per_row;
        new_vals[offset] = vals[idx++];
    }

    return 0;
}
