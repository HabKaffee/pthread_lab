#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

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
    long thread_rank = args->thread_num;
    long block_size = sqrt((double)thread_count);
    long local_n = fmax(args->n / block_size, 1);
    long local_m = fmax(args->m / block_size, 1);

    long block_row_st = 0;
    long block_col_st = 0;

    long inner_i = 0;
    while (thread_rank > 0) {
        inner_i += local_m;
        if ((inner_i < args->m) && (inner_i + local_m <= args->m)) {
            block_col_st += local_m;
            --thread_rank;
        } else {
            block_col_st = 0;
            inner_i = 0;
            block_row_st += local_n;
            --thread_rank;
        }
    }
    printf("block_row_st = %ld, block_cos_st = %ld\nlocal_n = %ld, local_m = %ld\n", block_row_st, block_col_st, local_n, local_m);
    int to_add = fmin(local_n, local_m);
    pthread_mutex_lock(&mutex);
    for (int block_i = 1; (block_i < args->n) && (block_i < args->m); block_i += to_add) {
        printf("block_i = %d\n", block_i);
        for (int i = block_row_st; i < args->n; i += block_i) {
            // printf("i = %d\n", i);
            for (int k = 0; k < args->k; k += block_i) {
                for (int j = block_col_st; j < args->m; j += block_i) {
                    args->result[i][k] += args->first[i][j] * args->second[j][k];
                    // printf("res[%d][%d] = %lf\n", i, k, args->result[i][k]);
                }
            }
        }
    }
    pthread_mutex_unlock(&mutex);
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

double get_wall_time() {
    struct timeval time;
    gettimeofday(&time, NULL);
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
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
    long n = 100, m = 100, k = 100;
    if (argv[2] && argv[3] && argv[4]) {
        n = strtol(argv[2], NULL, 10);
        m = strtol(argv[3], NULL, 10);
        k = strtol(argv[4], NULL, 10);
    }
    pthread_t* thread_handles = malloc(thread_count * sizeof(pthread_t));

    double** first_matrix;
    double** second_matrix;
    double** result_row_col_matrix;
    double** result_col_row_matrix;
    double** result_block_matrix;
    double** sequential_matrix;

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

    // Row by column multiplication in parallel
    ARGS* args_array;
    
    args_array = initialise_args_for_pthread(args_array, first_matrix, second_matrix, result_row_col_matrix, n, m, k);

    pthread_mutex_init(&mutex, NULL);
    double time_st_row = get_wall_time();
    for (int i = 0; i < thread_count; ++i) {
        pthread_create(&thread_handles[i], NULL, multiply_matrix_by_row_parallel, (void*) &args_array[i]);
    }
    for (int i = 0; i < thread_count; ++i) {
        pthread_join(thread_handles[i], NULL);
    }
    double time_end_row = get_wall_time();
    // printf("Parallel row by column done\n");
    
    // print_matrix(result_row_col_matrix, n, k);

    // Column by row multiplication in parallel
    args_array = initialise_args_for_pthread(args_array, first_matrix, second_matrix, result_col_row_matrix, n, m, k);
    double time_st_col = get_wall_time();
    for (int i = 0; i < thread_count; ++i) {
        pthread_create(&thread_handles[i], NULL, multiply_matrix_by_column_parallel, (void*) &args_array[i]);
    }
    for (int i = 0; i < thread_count; ++i) {
        pthread_join(thread_handles[i], NULL);
    }
    double time_end_col = get_wall_time();
    printf("Parallel column by row done\n");

    args_array = initialise_args_for_pthread(args_array, first_matrix, second_matrix, result_block_matrix, n, m, k);
    for (int i = 0; i < thread_count; ++i) {
        pthread_create(&thread_handles[i], NULL, multiply_matrix_by_block_parallel, (void*) &args_array[i]);
    }
    for (int i = 0; i < thread_count; ++i) {
        pthread_join(thread_handles[i], NULL);
    }
    printf("Parallel block multiplication done\n");

    print_matrix(result_block_matrix, n, k);

    pthread_mutex_destroy(&mutex);
    // print_matrix(result_col_row_matrix, n, k);

    // Sequential row by column multiplication
    double time_st_seq = get_wall_time();
    multiply_matrix_by_row(first_matrix, second_matrix, sequential_matrix, n, m, k);
    double time_end_seq = get_wall_time();
    print_matrix(sequential_matrix, n, k);
    
    printf("Row * col parallel = Sequential | %s\n", validate_result(result_row_col_matrix, sequential_matrix, n, k, eps) ? "true" : "false");
    printf("Col * row parallel = Sequential | %s\n", validate_result(result_col_row_matrix, sequential_matrix, n, k, eps) ? "true" : "false");
    printf("  Block parallel   = Sequential | %s\n", validate_result(result_block_matrix, sequential_matrix, n, k, eps)   ? "true" : "false");
    // printf("%d x %d and %d x %d matrix multiplication with %ld threads\n", n, m, m, k, thread_count);
    // printf("Seq = %lf\nRow parallel = %lf\nCol parallel = %lf\n", difftime(time_end_seq, time_start_seq), 
    //                                                               difftime(time_end_row, time_start_row), 
    //                                                               difftime(time_end_col, time_start_col));
    // printf("Seq = %lf\nRow parallel = %lf\nCol parallel = %lf\n", time_end_seq-time_st_seq, time_end_row - time_st_row, time_end_col - time_st_col);
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