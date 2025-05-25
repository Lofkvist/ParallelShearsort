#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void shuffle(int *array, int size) {
    for (int i = size - 1; i > 0; --i) {
        int j = rand() % (i + 1);
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s N output_file\n", argv[0]);
        return 1;
    }

    int N = atoi(argv[1]);
    const char *filename = argv[2];
    int total = N * N;

    int *data = malloc(total * sizeof(int));
    if (!data) {
        perror("malloc");
        return 1;
    }

    // Fill the array with numbers from 1 to N^2
    for (int i = 0; i < total; ++i) {
        data[i] = i + 1;
    }

    srand(time(NULL));
    shuffle(data, total);

    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("fopen");
        free(data);
        return 1;
    }

    fprintf(f, "%d\n", N);
    for (int i = 0; i < total; ++i) {
        fprintf(f, "%d ", data[i]);
    }
    fprintf(f, "\n");

    fclose(f);
    free(data);
    return 0;
}
