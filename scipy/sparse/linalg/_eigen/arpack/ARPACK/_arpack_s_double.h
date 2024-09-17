#include "_arpack.h"

// BLAS Routines used
void dswap_(int* n, double* x, int* incx, double* y, int* incy);
void dcopy_(int* n, double* x, int* incx, double* y, int* incy);

// LAPACK Routines used