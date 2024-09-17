#include "_arpack.h"

// BLAS Routines used
void saxpy_(int* n, double* alpha, double* x, int* incx, double* y, int* incy);
void sswap_(int* n, double* x, int* incx, double* y, int* incy);
void scopy_(int* n, double* x, int* incx, double* y, int* incy);
void sscal_(int* n, double* alpha, double* x, int* incx);
void sgemv_(char* trans, int* m, int* n, double* alpha, double* a, int* lda, double* x, int* incx, double* beta, double* y, int* incy);

// LAPACK Routines used
void slartg_(double* f, double* g, double* c, double* s, double* r);