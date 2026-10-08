# Red-Black SOR Poisson Solver with OpenMP

A red-black successive over-relaxation (SOR) solver for the Poisson problem with periodic boundary conditions, written in C in a sequential and an OpenMP version, together with a performance study of both.

This was a project for the Parallel Computing course at Radboud University. The solver builds on the starter code provided by the course (`util.h` and the original sequential solver); my contributions are the tiling of the loops, the OpenMP parallelization and the measurements.

## Repository Contents

- `src/sor_seq.c`: sequential solver, with the grid traversed in 64 x 64 tiles
- `src/sor_omp.c`: OpenMP solver
- `util.h`: timing and helper functions from the course's starter code
- `results/`: raw output of all measurements
- `docs/`: full project report

A brief overview is provided below. The full report at `docs/report.pdf` contains more detail on the methodology, the convergence behavior and the analysis.

## What I Wanted to Find Out

- How much faster does the solver get with OpenMP, and for which grid sizes?
- What limits the speedup: thread overhead, the number of cores, or memory?
- How does the cost per grid cell change when the grid no longer fits in cache?
- Does tiling the loops for cache use pay off?

## Approach

The grid is colored like a chessboard. A black cell only depends on its four red neighbors and the other way around, so all cells of one color can be updated at the same time without data races.

One iteration is:

```
update all black cells   (parallel over rows of 64 x 64 tiles)
    ↓
copy the periodic boundary into the ghost cells
    ↓
update all red cells     (parallel over rows of 64 x 64 tiles)
    ↓
copy the periodic boundary into the ghost cells
```

The iteration is repeated until the number of iterations given by the convergence formula is reached, or the `MAX_ITER` argument, whichever is smaller. At the end the result is compared with the exact solution `sin(x + y)` and the relative error is printed.

Both programs solve the problem 10 times and report the mean and standard deviation of the time per solve.

## Build and Run

```bash
cd src
gcc -O3 -march=native          -o sor_seq sor_seq.c -lm
gcc -O3 -march=native -fopenmp -o sor_omp sor_omp.c -lm

./sor_seq 512 1000                      # grid size, maximum iterations
OMP_NUM_THREADS=4 ./sor_omp 512 1000
```

## Setup

- Intel Core i7-1165G7 laptop: 4 physical cores, 8 hardware threads, 16 GB DDR4-3200
- GCC 13.3.0 on Ubuntu, flags `-O3 -march=native`
- Grid sizes 64 to 1024, 1 to 8 threads, at most 1000 iterations

## Results

Speedup of the OpenMP version over the sequential version:

| Grid | 1 thread | 2 threads | 4 threads | 8 threads |
|---|---|---|---|---|
| 64 x 64 | 1.19 | 0.70 | 0.69 | 0.16 |
| 128 x 128 | 1.07 | 1.02 | 0.83 | 0.41 |
| 256 x 256 | 0.83 | 1.08 | 1.19 | 0.54 |
| 512 x 512 | 0.98 | 1.49 | 1.78 | 1.38 |
| 1024 x 1024 | 0.99 | 1.61 | 2.85 | 2.33 |

- **Large grids benefit.** The best result is a speedup of 2.85 on 4 threads for the 1024 x 1024 grid, an efficiency of 71%.
- **Small grids get slower.** For 64 x 64 and 128 x 128, starting and coordinating the threads costs more than the work itself. The parallel loop also runs over rows of tiles, so a 64 x 64 grid has only one unit of work per sweep.
- **Eight threads are slower than four.** The machine has 4 physical cores, so the extra threads share cores.
- **Cache matters.** One cell update costs about 1 ns while both arrays fit in cache and 2.4 ns for the 1024 x 1024 grid, which does not fit in the 12 MB L3 cache.

### Tiling

Tiling did not pay off with these compiler flags. On the 512 x 512 grid:

| Version | Time per solve | Bandwidth |
|---|---|---|
| Starter code, untiled | 237 ± 27 ms | 26.5 GB/s |
| `sor_seq.c`, tiled | 328 ± 3 ms | 19.2 GB/s |

The speedups above are relative to the tiled sequential version. Relative to the untiled starter code, the 4-thread speedup on the 512 x 512 grid is about 1.5.

## Limitations

- All measurements were taken on a laptop with a desktop session running, not on a dedicated compute node.
- The sequential timing for the 64 x 64 grid was noisy, so the 64 x 64 row is only reliable qualitatively.
- The 8-thread runs varied a lot between repetitions.
- The comparison with the untiled starter code was only made for the 512 x 512 grid.
