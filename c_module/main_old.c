#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

typedef double F_TYPE;
const F_TYPE min_delta = 1e-8f;

typedef F_TYPE (* FP)(const F_TYPE*, const int, const int);

F_TYPE* grad(FP f, F_TYPE* ans, F_TYPE* x, int n, int k) {
    F_TYPE xd[1000];

    for (int i = 0; i < n; ++i) {
        xd[i] = x[i];
    }

    F_TYPE temp;
    for (int i = 0; i < n; ++i) {
        temp = xd[i];
        xd[i] += min_delta;
//        printf("%f %f\n", f(xd, n, k), f(x, n, k));
        ans[i] = (f(xd, n, k) - f(x, n, k)) / min_delta;
        xd[i] = temp;
    }

//    printf("grd x:[");
//    for(int j = 0; j < n - 1; ++j)
//        printf("%10.3f, ", x[j]);
//    printf("%10.3f] - [", x[n-1]);
//    for(int j = 0; j < n - 1; ++j)
//        printf("%10.3f, ", ans[j]);
//    printf("%10.3f]\n", ans[n-1]);

    return ans;
}

//int gradient_descent_constant_with_end_condition(FP f, const F_TYPE* x0, F_TYPE lr, F_TYPE eps, int max_iter, F_TYPE minimum, int n, int k) {
//    F_TYPE x[1000];
//    int steps = 0;
//
//    for (int i = 0; i < n; ++i) {
//        x[i] = x0[i];
//    }
//
//    for (int i = 0; i < max_iter; ++i) {
//        if (fabsf(f(x, n, k) - minimum) < eps) {
//            break;
//        }
//
//        ++steps;
//        F_TYPE grad_res[1000];
//        grad(f, grad_res, x, n, k);
//
//        for (int j = 0; j < n; ++j) {
//            x[j] -= lr * grad_res[j];
//        }
//    }
//
//    return steps;
//}


F_TYPE f(const F_TYPE* x, const int n, const int k) {
    F_TYPE res = 0.0f;

    res += (F_TYPE)k * x[0] * x[0];

    #pragma omp parallel for num_threads(omp_get_max_threads())
    for (int i = 1; i < n; ++i) {
        res += x[i] * x[i];
    }

    return res;
}

F_TYPE fd(FP f, F_TYPE alpha, const F_TYPE* x, const F_TYPE* d, int n, int k) {
    F_TYPE x_alpha[1000];

    for (int i = 0; i < n; ++i) {
        x_alpha[i] = x[i] - alpha * d[i];
    }
    return f(x_alpha, n, k);
}

F_TYPE dichotomy_search(FP f, F_TYPE a, F_TYPE b, F_TYPE eps, F_TYPE* x, F_TYPE* d, int n, int k) {
    while (b - a > eps) {
        F_TYPE c = (a + b) / 2;
        if (fd(f, c - eps, x, d, n, k) < fd(f, c + eps, x, d, n, k)) {
            b = c;
        } else {
            a = c;
        }
    }
    return (a + b) / 2;
}

int gradient_descent_dichotomy_with_end_condition(FP f, const F_TYPE* x0, F_TYPE step_size, F_TYPE eps, int max_iter, F_TYPE minimum, int n, int k) {
    F_TYPE x[n];

    for (int i = 0; i < n; ++i) {
        x[i] = x0[i];
    }

    int steps = 0;
    for (int i = 0; i < max_iter; ++i) {
//        printf("%5d [", i);
//        for(int j = 0; j < n - 1; ++j)
//            printf("%10.3f, ", x[j]);
//        printf("%10.3f] - %10.3f\n", x[n-1], f(x, n, k));


        if (fabsf(f(x, n, k) - minimum) < eps) {
            break;
        }

        F_TYPE grad_x[1000];
        grad(f, grad_x, x, n, k);
        F_TYPE alpha = dichotomy_search(f, 0, 1, min_delta, x, grad_x, n, k);

        for (int j = 0; j < n; ++j) {
            x[j] -= alpha * step_size * grad_x[j];
        }

        steps += 1;
    }

    return steps;
}

void research(int n_st, int n_end, int n_step, int k_st, int k_end, int k_step, int number_of_repeats, double *table)
{
    F_TYPE eps = 1e-3f;
    int max_iter = 10000;
    F_TYPE minimum = 0.0f;

    #pragma omp parallel for num_threads(omp_get_max_threads())
    for(int i = n_st; i <= n_end; i += n_step)
    {
        #pragma omp parallel for num_threads(omp_get_max_threads())
        for(int j = k_st; j <= k_end; j += k_step)
        {
            double sum = 0;
            for(int c = 0; c < number_of_repeats; ++c)
            {
                srand(time(NULL));
                F_TYPE x0[1000];

                for(int p = 0; p < i; ++p)
                    x0[p] = rand()%60001 - 30000;

                sum += gradient_descent_dichotomy_with_end_condition(f, x0, 1.0f, eps, max_iter, minimum, i, j);
            }

            sum /= number_of_repeats;
            printf("%5d %5d\n", i, j);
            table[i*1001 + j] = sum;
        }
    }


    for(int i = n_st; i < n_end; i += n_step)
    {
        for(int j = k_st; j < k_end; j += k_step)
        {
            printf("%5d %5d - %g\n", i, j, table[i*1001 + j]);
        }
    }
}

int main() {
    clock_t start_time, end_time;
    double total_time;
    start_time = clock(); // Запоминаем время начала выполнения программы

    double* table = malloc(1001 * 1001 * sizeof(double));
    research(1, 32, 1, 1, 32, 1, 10, table);
    free(table);
    end_time = clock(); // Запоминаем время окончания выполнения программы
    total_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    printf("Program execution time: %f seconds\n", total_time);
    return 0;
}
