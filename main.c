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
time_t time_start_seq;
time_t time_end_seq;
time_t time_start_row;
time_t time_end_row;
time_t time_start_col;
time_t time_end_col;


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
            for (int j = 0; j < m; ++j) {
                result_matrix[i][l] += first_matrix[i][j] * second_matrix[j][l];
            }
        }
    }
}

void* multiply_matrix_by_column_parallel(void* arguments) {
    ARGS* args = (ARGS*) arguments;
    long thread_num = args->thread_num;
    long local_m = (args->m) / thread_count; // column division
    long start_ind = thread_num * local_m;
    long end_ind = (thread_num + 1) * local_m - 1;
    if (thread_num == (thread_count - 1)) {
        end_ind = (end_ind < (args->m - 1)) ? (args->m - 1) : end_ind;
    }
    double** temp_matrix;
    temp_matrix = allocate_memory_for_matrix(temp_matrix, args->n, args->k);
    for (int column = start_ind; column <= end_ind; ++column) {
        for (int k = 0; k < args->k; ++k) {
            for (int i = 0; i < args->n; ++i) {
                temp_matrix[i][k] += args->first[i][column] * args->second[column][k];
            }
        }
    }
    for (int i = 0; i < args->n; ++i) {
        for (int j = 0; j < args->k; ++j) {
            pthread_mutex_lock(&mutex);
            args->result[i][j] += temp_matrix[i][j];
            pthread_mutex_unlock(&mutex);
        }
    }
    deallocate_memory_for_matrix(temp_matrix, args->n);
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

void* multiply_matrix_by_block_parallel(void* arguments) {
    ARGS* args = (ARGS*) arguments;
    long thread_num = args->thread_num;
    long local_n = (args->n)/floor(sqrt(thread_count)+1);
    long local_m = (args->m)/floor(sqrt(thread_count)+1);
    long local_k = (args->k)/floor(sqrt(thread_count)+1);
    
    long start_i = (thread_num * args->n) / args->n;//thread_num*local_n;
    long end_i = start_i + local_n; // (thread_num + 1) * local_n - 1;
    long start_j = (thread_num * args->m) % args->m; //thread_num*local_n;
    long end_j = start_j + local_m;//(thread_num + 1) * local_m - 1;
    long start_k = (thread_num * args->k) / args->k;//thread_num*local_n;
    long end_k = start_k + local_k;//(thread_num + 1) * local_k - 1;
    // if (thread_num == (thread_count - 1)) {
    //     end_i = (end_i < (args->n - 1)) ? (args->n - 1) : end_i;
    //     end_j = (end_j < (args->m - 1)) ? (args->m - 1) : end_j;
    //     end_k = (end_k < (args->k - 1)) ? (args->k - 1) : end_k;
    // }
    printf("%ld %ld %ld i (%ld %ld) j (%ld %ld) k (%ld %ld) \n", local_n, local_m, local_k, start_i, end_i, start_j, end_j, start_k, end_k);
    double** temp_matrix;
    temp_matrix = allocate_memory_for_matrix(temp_matrix, args->n, args->k);
    for (int i = start_i; i <= end_i; ++i) {
        for (int k = start_k; k <= end_k; ++k) {
            for (int j = start_j; j <= end_j; ++j) {
                temp_matrix[i][k] += args->first[i][j] * args->second[j][k];
            }
        }
    }
    for (int i = start_i; i <= end_i; ++i) {
        for (int k = start_k; k <= end_k; ++k) {
            pthread_mutex_lock(&mutex);
            args->result[i][k] += temp_matrix[i][k];
            pthread_mutex_unlock(&mutex);
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
        args_array[i].first  = &(&first_mat)[0][0];
        args_array[i].second = &(&second_mat)[0][0];
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
    double** result_block_matrix;
    double** sequential_matrix;

    int n = 2048, m = 2048, k = 2048;
    // int n = 4096, m = 4096, k = 4096;
    first_matrix = allocate_memory_for_matrix(first_matrix, n, m);
    second_matrix = allocate_memory_for_matrix(second_matrix, m, k);
    result_row_col_matrix = allocate_memory_for_matrix(result_row_col_matrix, n, k);
    result_col_row_matrix = allocate_memory_for_matrix(result_col_row_matrix, n, k);
    result_block_matrix = allocate_memory_for_matrix(result_block_matrix, n, k);
    sequential_matrix = allocate_memory_for_matrix(sequential_matrix, n, k);

    init_random_matrix(first_matrix, n, m);
    init_random_matrix(second_matrix, m, k);
    // print_matrix(first_matrix, n, m);
    // print_matrix(second_matrix, m, k);

    // first_matrix[0][0] = 1;
    // first_matrix[0][1] = 2;
    // first_matrix[0][2] = 3;
    // first_matrix[0][3] = 4;
    // first_matrix[1][0] = 4;
    // first_matrix[1][1] = 3;
    // first_matrix[1][2] = 2;
    // first_matrix[1][3] = 1;
    // first_matrix[2][0] = 0;
    // first_matrix[2][1] = 1;
    // first_matrix[2][2] = 0;
    // first_matrix[2][3] = 1;

    // second_matrix[0][0] = 1;
    // second_matrix[0][1] = 5;
    // second_matrix[1][0] = 2;
    // second_matrix[1][1] = 6;
    // second_matrix[2][0] = 3;
    // second_matrix[2][1] = 7;
    // second_matrix[3][0] = 4;
    // second_matrix[3][1] = 8;
    // Row by column multiplication in parallel
    ARGS* args_array;
    
    args_array = initialise_args_for_pthread(args_array, first_matrix, second_matrix, result_row_col_matrix, n, m, k);

    pthread_mutex_init(&mutex, NULL);
    time(&time_start_row);
    for (int i = 0; i < thread_count; ++i) {
        pthread_create(&thread_handles[i], NULL, multiply_matrix_by_row_parallel, (void*) &args_array[i]);
    }
    for (int i = 0; i < thread_count; ++i) {
        pthread_join(thread_handles[i], NULL);
    }
    time(&time_end_row);
    printf("Parallel row by column done\n");
    
    // print_matrix(result_row_col_matrix, n, k);

    // Column by row multiplication in parallel
    args_array = initialise_args_for_pthread(args_array, first_matrix, second_matrix, result_col_row_matrix, n, m, k);
    time(&time_start_col);
    for (int i = 0; i < thread_count; ++i) {
        pthread_create(&thread_handles[i], NULL, multiply_matrix_by_column_parallel, (void*) &args_array[i]);
    }
    for (int i = 0; i < thread_count; ++i) {
        pthread_join(thread_handles[i], NULL);
    }
    time(&time_end_col);
    printf("Parallel column by row done\n");

    // args_array = initialise_args_for_pthread(args_array, first_matrix, second_matrix, result_block_matrix, n, m, k);
    // for (int i = 0; i < thread_count; ++i) {
    //     pthread_create(&thread_handles[i], NULL, multiply_matrix_by_block_parallel, (void*) &args_array[i]);
    // }
    // for (int i = 0; i < thread_count; ++i) {
    //     pthread_join(thread_handles[i], NULL);
    // }
    // printf("Parallel block multiplication done\n");

    // print_matrix(result_block_matrix, n, k);

    pthread_mutex_destroy(&mutex);
    // print_matrix(result_col_row_matrix, n, k);

    // Sequential row by column multiplication
    time(&time_start_seq);
    multiply_matrix_by_row(first_matrix, second_matrix, sequential_matrix, n, m, k);
    time(&time_end_seq);
    // print_matrix(sequential_matrix, n, k);
    
    printf("Row * col parallel = Sequential | %s\n", validate_result(result_row_col_matrix, sequential_matrix, n, k, eps) ? "true" : "false");
    printf("Col * row parallel = Sequential | %s\n", validate_result(result_col_row_matrix, sequential_matrix, n, k, eps) ? "true" : "false");

    printf("Seq = %ld\nRow parallel = %ld\nCol parallel = %ld\n", time_end_seq - time_start_seq, time_end_row - time_start_row, time_end_col - time_start_col);
    // printf("  Block parallel   = Sequential | %s\n", validate_result(result_block_matrix, sequential_matrix, n, k, eps)   ? "true" : "false");
    uninitialise_args_for_pthread(args_array);
    
    deallocate_memory_for_matrix(first_matrix, n);
    deallocate_memory_for_matrix(second_matrix, m);
    deallocate_memory_for_matrix(result_col_row_matrix, n);
    deallocate_memory_for_matrix(result_row_col_matrix, n);
    deallocate_memory_for_matrix(result_block_matrix, n);
    deallocate_memory_for_matrix(sequential_matrix, n);

    free(thread_handles);
    return 0;
}