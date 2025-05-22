#include <stdio.h>

int main() {
    int my_rank = ;
    int comm_size = 3;
    int side_len = 7;
    int local_row_count = 3;
    int recv_count = 12;

    int i, j;
    
    int rev = 0;
    int col = 1 + my_rank;
    for (i = 0; i < recv_count/local_row_count; i++) {
        for (j = 0; j < local_row_count; j++) {
            int idx = (j) * side_len + col;
            printf("%d\n", idx);
        }

        col += comm_size;

        if (col >= side_len) {
            rev++;
            col = 1 + my_rank + rev;
        }
    }

    return 0;
}
