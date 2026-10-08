/**
 * Solves ∇u = f on [a, b] x [a, b] for a function u that satisfies
 * u(a, y) = u(b, y)
 * u(x, a) = u(x, b)
 *
 * for all x, y.
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../util.h"
#include <time.h>
#include <string.h>

/* Number of timed repetitions of the solve */
#define RUNS 10

double u(double x, double y)
{
    return sin(x + y);
}

double f(double x, double y)
{
    return -2.0 * sin(x + y);
}

void set_zero(size_t m, size_t n, double *x)
{
    memset(x, 0, m * n * sizeof(double));
}

/* Turns 2D index (i, j) of an ... x n array with ghost-cells into
 * the offset from the start of the array. */
#define GHOST(i, j, n) ((i) + 1) * ((n) + 2) + ((j) + 1)

void update_ghost_cells(size_t m, size_t n, double *x)
{
    /* Left and right ghost columns */
    for (size_t i = 0; i < m; i++) {
        x[GHOST(i, -1, n)] = x[GHOST(i, n - 1, n)];
        x[GHOST(i, n, n)]  = x[GHOST(i, 0, n)];
    }

    /* Top and bottom ghost rows, including the corners. Must run after the
     * column update, because it copies the ghost columns of rows 0 and m-1. */
    for (ssize_t k = -1; k <= (ssize_t)n; k++) {
        x[GHOST(-1, k, n)] = x[GHOST(m - 1, k, n)];
        x[GHOST(m, k, n)]  = x[GHOST(0, k, n)];
    }
}

double *discretize_spooky(size_t m, size_t n, double h, double a,
                          double (*fun)(double, double))
{
    double *dis = malloc((m + 2) * (n + 2) * sizeof(double));
    const size_t tile = 64;

    for (size_t bi = 0; bi < m; bi += tile) {
        for (size_t bj = 0; bj < n; bj += tile) {
            size_t bi_end = min(bi + tile, m);
            size_t bj_end = min(bj + tile, n);
            for (size_t i = bi; i < bi_end; i++) {
                for (size_t j = bj; j < bj_end; j++) {
                    double x = a + i * h;
                    double y = a + j * h;
                    dis[GHOST(i, j, n)] = fun(x, y);
                }
            }
        }
    }

    update_ghost_cells(m, n, dis);
    return dis;
}

/**
 * Loop body of sor function
 **/
#define SOR_BODY(n, u_dis, f_dis, i, j, omega, h)                              \
{                                                                              \
    u_dis[GHOST(i, j, n)] =                                                    \
       (1 - omega) *  u_dis[GHOST(i    , j    , n)] +                          \
       omega / 4   * (u_dis[GHOST(i - 1, j    , n)] +                          \
                      u_dis[GHOST(i + 1, j    , n)] +                          \
                      u_dis[GHOST(i    , j - 1, n)] +                          \
                      u_dis[GHOST(i    , j + 1, n)] -                          \
                      h * h * f_dis[GHOST(i, j, n)]);                          \
}

void sor(double *u_dis, double *f_dis,
         size_t m, size_t n, double omega, double h)
{
    const size_t block_size = 64;

    // Black update with blocking
    for (size_t bi = 0; bi < m; bi += block_size) {
        size_t bi_end = min(bi + block_size, m);
        for (size_t bj = 0; bj < n; bj += block_size) {
            size_t bj_end = min(bj + block_size, n);

            // Process black cells in current tile. Rows are handled in pairs
            // (i, i + 1); a trailing unpaired row is covered by the first
            // inner loop on its own.
            for (size_t i = bi; i < bi_end; i += 2) {
                for (size_t j = bj; j < bj_end; j += 2) {
                    SOR_BODY(n, u_dis, f_dis, i, j, omega, h);
                }
                for (size_t j = bj + 1; j < bj_end && (i + 1) < m; j += 2) {
                    SOR_BODY(n, u_dis, f_dis, i + 1, j, omega, h);
                }
            }
        }
    }
    update_ghost_cells(m, n, u_dis);

    // Red update with blocking
    for (size_t bi = 0; bi < m; bi += block_size) {
        size_t bi_end = min(bi + block_size, m);
        for (size_t bj = 0; bj < n; bj += block_size) {
            size_t bj_end = min(bj + block_size, n);

            for (size_t i = bi; i < bi_end; i += 2) {
                for (size_t j = bj + 1; j < bj_end; j += 2) {
                    SOR_BODY(n, u_dis, f_dis, i, j, omega, h);
                }
                for (size_t j = bj; j < bj_end && (i + 1) < m; j += 2) {
                    SOR_BODY(n, u_dis, f_dis, i + 1, j, omega, h);
                }
            }
        }
    }
    update_ghost_cells(m, n, u_dis);
}

int sor_solve(double *u_dis, double *f_dis,
              size_t m, size_t n, double h, int max_iter)
{
    /**
     * Math black magic computing the number of iterations necessary for
     * convergence.
     **/
    double log_spectral_radius = log1p(-2.0 * sin(M_PI * h) /
                                    (1.0 + sin(M_PI * h)));
    int iter = 200.0 * log(h) / log_spectral_radius;
    iter = min(iter, max_iter);

    double omega = 2.0 / (1.0 + sin(M_PI * h));
    for (int t = 0; t < iter; t++) {
        sor(u_dis, f_dis, m, n, omega, h);
    }

    return iter;
}

double Linf(double *u_dis, size_t m, size_t n)
{
    double maximum = -1;
    const size_t tile = 64;

    for (size_t bi = 0; bi < m; bi += tile) {
        for (size_t bj = 0; bj < n; bj += tile) {
            double local_max = -1;
            size_t bi_end = min(bi + tile, m);
            size_t bj_end = min(bj + tile, n);
            for (size_t i = bi; i < bi_end; i++) {
                for (size_t j = bj; j < bj_end; j++) {
                    local_max = fmax(fabs(u_dis[GHOST(i, j, n)]), local_max);
                }
            }
            maximum = fmax(local_max, maximum);
        }
    }
    return maximum;
}

int main(int argc, char **argv)
{
    if (argc != 3) {
        printf("Usage: %s <N> <MAX_ITER>\n"
               "\tN       : Discretize f on an N x N grid.\n"
               "\tMAX_ITER: Run for at most this many iterations\n",
               argv[0]);
        return EXIT_FAILURE;
    }

    size_t n     = atol(argv[1]);
    int max_iter = atoi(argv[2]);

    double a = 0;
    double b = 2 * M_PI;
    double h = (b - a) / n;

    double *f_dis = discretize_spooky(n, n, h, a, f);
    double *u_dis = malloc((n + 2) * (n + 2) * sizeof(double));

    double times[RUNS];
    int iter = 0;

    for (int run = 0; run < RUNS; run++) {
        set_zero(n + 2, n + 2, u_dis);
        double start = wtime();
        iter = sor_solve(u_dis, f_dis, n, n, h, max_iter);
        double stop  = wtime();
        times[run] = stop - start;
    }

    /* Mean and sample standard deviation of the time per run */
    double avg_time = 0.0;
    for (int run = 0; run < RUNS; run++) {
        avg_time += times[run];
    }
    avg_time /= RUNS;

    double variance = 0.0;
    for (int run = 0; run < RUNS; run++) {
        variance += (times[run] - avg_time) * (times[run] - avg_time);
    }
    double stddev_time = sqrt(variance / (RUNS - 1));

    size_t size_of_array = n * n * sizeof(double);
    double bandwidth = 3.0 * size_of_array * iter / 1e9 / avg_time;
    double gflops = 11.0 * n * n * iter / 1e9 / avg_time;


    double *real_answer = discretize_spooky(n, n, h, a, u);
    double norm_real = Linf(real_answer, n, n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            u_dis[GHOST(i, j, n)] -= real_answer[GHOST(i, j, n)];
        }
    }
    double norm_diff = Linf(u_dis, n, n);
    double relative_error = norm_diff / norm_real;

    fprintf(stderr, "Bandwidth: %.3lf GB/s\n", bandwidth);
    fprintf(stderr, "Compute rate: %.3lf GFLOPS\n", gflops);
    fprintf(stderr, "Relative error: %.3e\n", relative_error);
    fprintf(stderr, "Iterations per run: %d\n", iter);
    printf("%.3lf\n", bandwidth);
    printf("Time per run (%d runs): %.6lf s avg, %.6lf s stddev\n",
           RUNS, avg_time, stddev_time);

    free(u_dis);
    free(f_dis);
    free(real_answer);

    return EXIT_SUCCESS;
}
