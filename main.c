#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>

const int lower_bound = 0;
const int upper_bound = 2;

extern inline double get_random_double(double a, double b) {
    return a + (rand()/((double)RAND_MAX / (b-a)));
}

double** allocate_memory_for_matrix(double** matrix, int n, int m) {
    matrix = calloc(n, sizeof(double*));
    for (int i = 0; i < n; ++i) {
        matrix[i] = calloc(m, sizeof(double));
    }
    return matrix;
}

void deallocate_memory_for_matrix(double** matrix, int rows) {
    for (int i = 0; i < rows; ++i) {
        free(matrix[i]);
    }
    free(matrix);
}

void init_random_matrix(double** matrix, int n, int m) {
    if (!matrix) {
        printf("Received nullptr\nPointer to matrix is awaited\n");
        exit(1);
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            matrix[i][j] = get_random_double(lower_bound, upper_bound);
        }
    }
}

void print_matrix(double** matrix, int n, int m) {
    // print first matrix
    printf("%d x %d matrix\n", n, m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            printf("%lf\t", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void multiply_matrix_by_row(double** first_matrix, 
                                double** second_matrix,
                                double** result_matrix,
                                int n, int m, int k) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            for (int l = 0; l < k; ++l) {
                result_matrix[i][l] += first_matrix[i][j] * second_matrix[j][l];
            }
        }
    }
}

void* multiply_matrix_by_row_parallel(void* thread) {
    
}

int main() {
    srand(time(NULL));
    double** first_matrix;
    double** second_matrix;
    double** result_matrix;
    int n = 1, m = 2, k = 1;

    first_matrix = allocate_memory_for_matrix(first_matrix, n, m);
    second_matrix = allocate_memory_for_matrix(second_matrix, m, k);
    result_matrix = allocate_memory_for_matrix(result_matrix, n, k);
    
    init_random_matrix(first_matrix, n, m);
    init_random_matrix(second_matrix, m, k);
    print_matrix(first_matrix, n, m);
    print_matrix(second_matrix, m, k);

    multiply_matrix_by_row(first_matrix, second_matrix,
                            result_matrix, n, m, k);
    print_matrix(result_matrix, n, k);

    deallocate_memory_for_matrix(first_matrix, n);
    deallocate_memory_for_matrix(second_matrix, m);
    deallocate_memory_for_matrix(result_matrix, n);

    return 0;
}