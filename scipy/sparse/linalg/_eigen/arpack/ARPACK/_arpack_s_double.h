#include "_arpack.h"

// BLAS Routines used
void daxpy_(int* n, double* alpha, double* x, int* incx, double* y, int* incy);
void dswap_(int* n, double* x, int* incx, double* y, int* incy);
void dcopy_(int* n, double* x, int* incx, double* y, int* incy);
void dscal_(int* n, double* alpha, double* x, int* incx);
void dgemv_(char* trans, int* m, int* n, double* alpha, double* a, int* lda, double* x, int* incx, double* beta, double* y, int* incy);

// LAPACK Routines used
void dlartg_(double* f, double* g, double* c, double* s, double* r);