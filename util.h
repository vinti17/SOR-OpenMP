#include <sys/time.h>

double wtime(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return 1e-6 * (double)tv.tv_usec  + (double)tv.tv_sec;
}

size_t CeilDiv(size_t a, size_t b)
{
    return (a + b - 1) / b;
}

int min(int a, int b)
{
    return (a < b) ? a : b;
}

#define MPI_SAFE(fncall) do {                                                 \
   int retcode = (fncall);                                                    \
   if (retcode != MPI_SUCCESS) {                                              \
       char msg[1024];                                                        \
       int resultlen;                                                         \
       MPI_Error_string(retcode, msg, &resultlen);                            \
       fprintf(stderr, "Call %s failed with %s(%d)\n",                        \
                    #fncall, msg, retcode);                                   \
       MPI_Abort(MPI_COMM_WORLD, retcode);                                    \
   }                                                                          \
} while (0)

#define MPI_CHECK_STAT(stat) do {                                             \
   if (stat != MPI_SUCCESS) {                                                 \
       char msg[1024];                                                        \
       int resultlen;                                                         \
       MPI_Error_string(stat, msg, &resultlen);                               \
       fprintf(stderr, "Failure %s(%d)\n",                                    \
                      msg, retcode);                                          \
       MPI_Abort(MPI_COMM_WORLD, stat);                                       \
   }                                                                          \
} while (0)
