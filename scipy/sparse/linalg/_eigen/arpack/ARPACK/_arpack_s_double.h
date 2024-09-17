#include "_arpack.h"

// BLAS Routines used
void daxpy_(int* n, double* alpha, double* x, int* incx, double* y, int* incy);
void dswap_(int* n, double* x, int* incx, double* y, int* incy);
void dger_(int* m, int* n, double* alpha, double* x, int* incx, double* y, int* incy, double* a, int* lda);
void dcopy_(int* n, double* x, int* incx, double* y, int* incy);
void dscal_(int* n, double* alpha, double* x, int* incx);
double dnrm2_(int* n, double* x, int* incx);
void dgemv_(char* trans, int* m, int* n, double* alpha, double* a, int* lda, double* x, int* incx, double* beta, double* y, int* incy);
void dtrmm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n, double* alpha, double* a, int* lda, double* b, int* ldb);



void dgemv_(char* trans, int* m, int* n, double* alpha, double* a, int* lda, double* x, int* incx, double* beta, double* y, int* incy);

// LAPACK Routines used
void dlartg_(double* f, double* g, double* c, double* s, double* r);