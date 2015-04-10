#include <R.h>

double *real_vector(int n)
{
    return (double*) R_alloc(n, sizeof(double));
}

int *int_vector(int n)
{
    return (int*) R_alloc(n, sizeof(int));
}
