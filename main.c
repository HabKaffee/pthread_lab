#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

const int lower_bound = 0;
const int upper_bound = 2;
const double eps = 1e-6;
long thread_count = 0;

pthread_mutex_t mutex;

struct args_parallel {
    long thread_num;
    double** first;
    double** second;
    double** result;
    int n;
    int m;
    int k;
} typedef ARGS;

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
        for (int l = 0; l < k; ++l) {
            result_matrix[i][l] = 0.0;
            for (int j = 0; j < m; ++j) {
                result_matrix[i][l] += first_matrix[i][j] * second_matrix[j][l];
            }
        }
    }
}

void* multiply_matrix_by_column_parallel(void* arguments) {
    ARGS* args = (ARGS*) arguments;
    // printf("%x\n", args);
    long thread_num = args->thread_num;
    // printf("Thread num %ld\n", thread_num);
    long local_m = (args->m) / thread_count; // column division
    // printf("local_m %ld\n", local_m);
    long start_ind = thread_num * local_m;
    // printf("start_ind %ld\n", start_ind);
    long end_ind = (thread_num + 1) * local_m - 1;
    // printf("end_ind %ld\n", end_ind);
    if (thread_num == (thread_count - 1)) {
        end_ind = (end_ind < (args->m - 1)) ? (args->m - 1) : end_ind;
    }
    for (int column = start_ind; column <= end_ind; ++column) {
        for (int k = 0; k < args->k; ++k) {
            for (int i = 0; i < args->n; ++i) {
                pthread_mutex_lock(&mutex);
                args->result[i][k] += args->first[i][column] * args->second[column][k];
                pthread_mutex_unlock(&mutex);
            }
        }
    }
}

void* multiply_matrix_by_row_parallel(void* arguments) {
    ARGS* args = (ARGS*) arguments;
    long thread_num = args->thread_num;
    long local_n = (args->n) / thread_count;
    long start_ind = thread_num * local_n;
    long end_ind = (thread_num + 1) * local_n - 1;
    if (thread_num == (thread_count - 1)) {
        end_ind = (end_ind < (args->n - 1)) ? (args->n - 1) : end_ind;
    }
    for (int i = start_ind; i <= end_ind; ++i) {
        for (int k = 0; k < args->k; ++k) {
            for (int j = 0; j < args->m; ++j) {
                args->result[i][k] += args->first[i][j] * args->second[j][k];
            }
        }
    }
}

bool validate_result(double** first, double** second, int n, int m, double eps) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (fabs(first[i][j] - second[i][j]) > eps) {
                return false;
            }
        }
    }
    return true;
}

ARGS* initialise_args_for_pthread(ARGS* args_array, double** first_mat, 
                                  double** second_mat, double** result_mat, 
                                  int n, int m, int k) {
    args_array = malloc(sizeof(ARGS) * thread_count);

    for (int i = 0; i < thread_count; ++i) {
        args_array[i].thread_num = i;
        args_array[i].first = first_mat;
        args_array[i].second = second_mat;
        args_array[i].result = &(&result_mat)[0][0];
        args_array[i].n = n;
        args_array[i].m = m;
        args_array[i].k = k;
    }
    return args_array;
}

void uninitialise_args_for_pthread(ARGS* args_array) {
    for (int i = 0; i < thread_count; ++i) {
        args_array[i].first  = NULL;
        args_array[i].second = NULL;
        args_array[i].result = NULL;
    }
    free(args_array);
}

int main(int argc, char** argv) {
    // Valid imput guards
    if (!argv[1]) {
        printf("Number of thread is needed\nNull is provided\n");
        exit(2);
    } else if (strtol(argv[1], NULL, 10) < 1) {
        printf("Needed more than 0 threads, but %ld is provided\n", strtol(argv[1], NULL, 10));
        exit(3);
    }
    srand(time(NULL));

    thread_count = strtol(argv[1], NULL, 10);
    pthread_t* thread_handles = malloc(thread_count * sizeof(pthread_t));

    double** first_matrix;
    double** second_matrix;
    double** result_row_col_matrix;
    double** result_col_row_matrix;
    double** sequential_matrix;

    int n = 10, m = 5, k = 15;

    first_matrix = allocate_memory_for_matrix(first_matrix, n, m);
    second_matrix = allocate_memory_for_matrix(second_matrix, m, k);
    result_row_col_matrix = allocate_memory_for_matrix(result_row_col_matrix, n, k);
    result_col_row_matrix = allocate_memory_for_matrix(result_col_row_matrix, n, k);
    sequential_matrix = allocate_memory_for_matrix(sequential_matrix, n, k);

    init_random_matrix(first_matrix, n, m);
    init_random_matrix(second_matrix, m, k);
    print_matrix(first_matrix, n, m);
    print_matrix(second_matrix, m, k);
    
    // Row by column multiplication in parallel
    ARGS* args_array;
    
    args_array = initialise_args_for_pthread(args_array, first_matrix, second_matrix, result_row_col_matrix, n, m, k);

    pthread_mutex_init(&mutex, NULL);

    for (int i = 0; i < thread_count; ++i) {
        pthread_create(&thread_handles[i], NULL, multiply_matrix_by_row_parallel, (void*) &args_array[i]);
    }
    for (int i = 0; i < thread_count; ++i) {
        pthread_join(thread_handles[i], NULL);
    }
    
    print_matrix(result_row_col_matrix, n, k);

    // Column by row multiplication in parallel
    args_array = initialise_args_for_pthread(args_array, first_matrix, second_matrix, result_col_row_matrix, n, m, k);

    for (int i = 0; i < thread_count; ++i) {
        pthread_create(&thread_handles[i], NULL, multiply_matrix_by_column_parallel, (void*) &args_array[i]);
    }
    for (int i = 0; i < thread_count; ++i) {
        pthread_join(thread_handles[i], NULL);
    }
    
    pthread_mutex_destroy(&mutex);
    print_matrix(result_col_row_matrix, n, k);

    // Sequential row by column multiplication
    multiply_matrix_by_row(first_matrix, second_matrix, sequential_matrix, n, m, k);
    print_matrix(sequential_matrix, n, k);
    
    printf("Row * col parallel = Sequential %s\n", validate_result(result_row_col_matrix, sequential_matrix, n, k, eps) ? "true" : "false");
    printf("Col * row parallel = Sequential %s\n", validate_result(result_col_row_matrix, sequential_matrix, n, k, eps) ? "true" : "false");
    uninitialise_args_for_pthread(args_array);
    
    deallocate_memory_for_matrix(first_matrix, n);
    deallocate_memory_for_matrix(second_matrix, m);
    deallocate_memory_for_matrix(result_col_row_matrix, n);
    deallocate_memory_for_matrix(result_row_col_matrix, n);
    deallocate_memory_for_matrix(sequential_matrix, n);

    free(thread_handles);
    return 0;
}