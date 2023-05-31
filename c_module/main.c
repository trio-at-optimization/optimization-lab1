#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

typedef double F_TYPE;
const F_TYPE min_delta = 1e-8f;

//typedef F_TYPE (* FP)(const F_TYPE*, const int, const int, const int);

F_TYPE f(const F_TYPE* x, const int n, F_TYPE k, const int pos) {
    F_TYPE res = 0.0f;

//    #pragma omp parallel for shared(x) reduction(+:res) num_threads(omp_get_max_threads())
    for (int i = 0; i < n; ++i) {
        if(i == pos)
            res += k * x[i] * x[i];
        else
            res += x[i] * x[i];
    }

    return res;
}

F_TYPE* grad(F_TYPE* ans, F_TYPE* x, int n, F_TYPE k, int pos) {
    F_TYPE xd[1000];

//    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        xd[i] = x[i];
    }

//    #pragma omp for simd
    for (int i = 0; i < n; ++i) {
        F_TYPE temp;
        temp = xd[i];
        xd[i] += min_delta;
        ans[i] = (f(xd, n, k, pos) - f(x, n, k, pos)) / min_delta;
        xd[i] = temp;
    }
    return ans;
}

F_TYPE fd(F_TYPE alpha, const F_TYPE* x, const F_TYPE* d, int n, F_TYPE k, int pos) {
    F_TYPE x_alpha[1000];

//    #pragma omp for simd
    for (int i = 0; i < n; ++i) {
        x_alpha[i] = x[i] - alpha * d[i];
    }
    return f(x_alpha, n, k, pos);
}

F_TYPE dichotomy_search(F_TYPE a, F_TYPE b, F_TYPE eps, F_TYPE* x, F_TYPE* d, int n, F_TYPE k, int pos) {
    while (b - a > eps) {
        F_TYPE c = (a + b) / 2;
        if (fd(c - eps, x, d, n, k, pos) < fd(c + eps, x, d, n, k, pos)) {
            b = c;
        } else {
            a = c;
        }
    }
    return (a + b) / 2;
}

int gradient_descent_dichotomy_with_end_condition(const F_TYPE* x0, F_TYPE step_size, F_TYPE eps, int max_iter, F_TYPE minimum, int n, F_TYPE k, int pos) {
    F_TYPE x[n];

    for (int i = 0; i < n; ++i) {
        x[i] = x0[i];
    }

    int steps = 0;
    for (int i = 0; i < max_iter; ++i) {
        if (fabs(f(x, n, k, pos) - minimum) < eps) {
            break;
        }

        F_TYPE grad_x[1000];
        grad(grad_x, x, n, k, pos);
        F_TYPE alpha = dichotomy_search(0, 100, min_delta, x, grad_x, n, k, pos);

//        #pragma omp for simd
        for (int j = 0; j < n; ++j) {
            x[j] -= alpha * step_size * grad_x[j];
        }

        ++steps;
    }

    return steps;
}

void do_research(const int n_st, const int n_end, const int n_step, const F_TYPE k_st, const F_TYPE k_end, const F_TYPE k_step, const int number_of_repeats, const F_TYPE eps, const int max_iter, const F_TYPE minimum, const char* filename, double *table)
{
    int total_iterations = (int)(((double)((double)(n_end - n_st) / (double)(n_step + 1))) * ((k_end - k_st) / k_step));
    int completed_iterations = 0;
    int last_completed_iterations = 0;
    int j_count = (int)((k_end-k_st)/k_step);

    double elapsed, eta;
    double start_time = omp_get_wtime();
    double last_update_time = start_time;

#pragma omp parallel for num_threads((int)((double)omp_get_max_threads()*0.9))
    for (int i = n_st; i < n_end; i += n_step) {
//        #pragma omp parallel for num_threads((int)((double)omp_get_max_threads()*0.9))
        for (int j = 0; j < j_count; ++j) {
            double sum = 0;

            for (int d = 0; d < number_of_repeats; ++d)
            {
                srand(time(NULL));
                int pos = rand() % i;
                for (int c = 0; c < number_of_repeats; ++c) {
                    F_TYPE x0[1000];

                    for (int p = 0; p < i; ++p)
                        x0[p] = rand() % 60001 - 30000;

                    sum += gradient_descent_dichotomy_with_end_condition(x0, 1.0f, eps, max_iter, minimum, i, j*k_step + k_st, pos);
                }
            }

            sum /= (number_of_repeats*number_of_repeats);

#pragma omp critical
            {
                //printf("%5d %5d\n", i, j);
                table[i*1001 + j] = sum;

#pragma omp atomic
                ++completed_iterations;

               printf("n=%-5d, k=%-5.0f r=%7.3f | %d %d\n", i, j*k_step + k_st, table[i*1001 + j], total_iterations, completed_iterations);
            }



            // Update progress bar



//            double current_time = omp_get_wtime();
//            if (current_time - last_update_time >= 1.0) {
//                int progress = (int) (((double) completed_iterations / (double) total_iterations) * 100);
//                elapsed = current_time - start_time;
//                eta = ((double)(total_iterations - completed_iterations)) / ((double)(completed_iterations - last_completed_iterations))/(current_time - last_update_time);
//                last_completed_iterations = completed_iterations;
//                last_update_time = current_time;
//
//                printf("%5d %5.0f ", i, j*k_step + k_st);
//                printf("[");
//                for (int k = 0; k < progress / 2; ++k) {
//                    printf("#");
//                }
//                for (int k = progress / 2; k < 50; ++k) {
//                    printf("-");
//                }
//                printf("] %d%%", progress);
//                printf(" Ela: %.1f s  ETA: %.1f s\r\r", elapsed, eta);
//                fflush(stdout);
//            }
        }

#pragma omp critical
        {
            FILE* out = fopen(filename, "a");
            for(int j = 0; j < j_count; ++j)
            {
                fprintf(out, "%5d %7.2f - %g\n", i, j*k_step + k_st, table[i*1001 + j]);
            }
            fclose(out);
        }
    }
}


void research(int n_st, int n_end, int n_step, F_TYPE k_st, F_TYPE k_end, F_TYPE k_step, int number_of_repeats, char* filename, double *table)
{
    F_TYPE eps = 1e-3f;
    int max_iter = 10000;
    F_TYPE minimum = 0.0f;

//    FILE* out = fopen("result.txt", "r");
//    if(out != NULL)
//    {
//        fscanf(out, "%d%d%d%d%d%d", &n_st, &n_end, &n_step, &k_st, &k_end, &k_step);
//        fscanf(out, "%d%lf%d%lf", &number_of_repeats, &eps, &max_iter, &minimum);
//
//        int k = 0, i, j;
//        F_TYPE value;
//        while((k = fscanf(out, "%5d %5d = %g", &i, &j, &value)))
//        {
//            table[i*
//        }
//
//
//        fclose(out)
//    }


    do_research(n_st, n_end, n_step, k_st, k_end, k_step, number_of_repeats, eps, max_iter, minimum, filename, table);
}


int main() {
    printf("Hello\n");
    clock_t start_time, end_time;
    double total_time;
    start_time = clock();


    int n_st = 1, n_end = 1001, n_step = 80;
    F_TYPE k_st = 1, k_end = 1001, k_step = 80;


    double* table = malloc(1001 * 10001 * sizeof(double));
    research(n_st, n_end, n_step, k_st, k_end, k_step, 3, "result_1_1000_n2_step_0_01_num_research_999.txt", table);
    end_time = clock();
    total_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    printf("Program execution time: %f seconds\n", total_time);

//    start_time = clock();
//    printf("SAVE...");
//    FILE* out = fopen("result.txt", "a");
//    for(int i = n_st; i < n_end; i += n_step)
//    {
//        for(int j = k_st; j < k_end; j += k_step)
//        {
//            fprintf(out, "%5d %5d - %g\n", i, j, table[i*1001 + j]);
//        }
//    }
//    fclose(out);
    free(table);
//    printf(" OK\n");

//    end_time = clock();
//    total_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
//    printf("Program save time: %f seconds\n", total_time);
    return 0;
}
