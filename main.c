#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>

const int lower_bound = 0;
const int upper_bound = 2;

pthread_mutex_t mutex;

long thread_count = 0;

struct args_parallel {
    long thread_num;
    double** result;
    double** first;
    double** second;
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
        for (int l = 0; l < k; ++l) {
            result_matrix[i][l] = 0.0;
            for (int j = 0; j < m; ++j) {
                result_matrix[i][l] += first_matrix[i][j] * second_matrix[j][l];
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
            args->result[i][k] = 0.0;
            for (int j = 0; j < args->m; ++j) {
                pthread_mutex_lock(&mutex);
                args->result[i][k] += args->first[i][j] * args->second[j][k];
                pthread_mutex_unlock(&mutex);
            }
        }
    }
}

int main(int argc, char** argv) {
    if (!argv[1]) {
        printf("Number of thread is needed\nNull is provided\n");
        exit(2);
    } else if (strtol(argv[1], NULL, 10) < 1) {
        printf("Needed more than 0 threads, but %ld is provided\n", strtol(argv[1], NULL, 10));
        exit(3);
    }
    thread_count = strtol(argv[1], NULL, 10);
    pthread_t* thread_handles = malloc(thread_count * sizeof(pthread_t));

    srand(time(NULL));
    double** first_matrix;
    double** second_matrix;
    double** result_matrix;
    
    int n = 10, m = 5, k = 15;

    first_matrix = allocate_memory_for_matrix(first_matrix, n, m);
    second_matrix = allocate_memory_for_matrix(second_matrix, m, k);
    result_matrix = allocate_memory_for_matrix(result_matrix, n, k);

    init_random_matrix(first_matrix, n, m);
    init_random_matrix(second_matrix, m, k);
    print_matrix(first_matrix, n, m);
    print_matrix(second_matrix, m, k);
    
    ARGS* args_array = malloc(sizeof(ARGS) * thread_count);

    for (int i = 0; i < thread_count; ++i) {
        args_array[i].thread_num = i;
        args_array[i].first = first_matrix;
        args_array[i].second = second_matrix;
        args_array[i].result = &*result_matrix;
        args_array[i].n = n;
        args_array[i].m = m;
        args_array[i].k = k;
    }

    pthread_mutex_init(&mutex, NULL);

    for (int i = 0; i < thread_count; ++i) {
        pthread_create(&thread_handles[i], NULL, multiply_matrix_by_row_parallel, (void*) &args_array[i]);
    }
    for (int i = 0; i < thread_count; ++i) {
        pthread_join(thread_handles[i], NULL);
    }
    pthread_mutex_destroy(&mutex);
    print_matrix(result_matrix, n, k);

    multiply_matrix_by_row(first_matrix, second_matrix, result_matrix, n, m, k);
    print_matrix(result_matrix, n, k);
    
    deallocate_memory_for_matrix(first_matrix, n);
    deallocate_memory_for_matrix(second_matrix, m);
    deallocate_memory_for_matrix(result_matrix, n);

    free(args_array);
    free(thread_handles);
    return 0;
}