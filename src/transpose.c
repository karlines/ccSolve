void transpose(double *x, int nrx, int ncx, double *y)
{
    int i, j;

    for (i = 0; i < ncx; i++)
        for (j = 0; j < nrx; j++)
            y[j*ncx + i] = x[i*nrx + j];
}
