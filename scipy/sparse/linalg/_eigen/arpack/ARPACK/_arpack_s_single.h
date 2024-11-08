#include "_arpack.h"

// BLAS Routines used
void saxpy_(int* n, float* alpha, float* x, int* incx, float* y, int* incy);
void scopy_(int* n, float* x, int* incx, float* y, int* incy);
void sgemv_(char* trans, int* m, int* n, float* alpha, float* a, int* lda, float* x, int* incx, float* beta, float* y, int* incy);
void sger_(int* m, int* n, float* alpha, float* x, int* incx, float* y, int* incy, float* a, int* lda);
float snrm2_(int* n, float* x, int* incx);
void sscal_(int* n, float* alpha, float* x, int* incx);
void sswap_(int* n, float* x, int* incx, float* y, int* incy);
void strmm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n, float* alpha, float* a, int* lda, float* b, int* ldb);


// LAPACK Routines used
void sgeqr2_(int* m, int* n, float* a, int* lda, float* tau, float* work, int* info);
void slacpy_(char* uplo, int* m, int* n, float* a, int* lda, float* b, int* ldb);
void slartg_(float* f, float* g, float* c, float* s, float* r);
void sorm2r_(char* side, char* trans, int* m, int* n, int* k, float* a, int* lda, float* tau, float* c, int* ldc, float* work, int* info);
void ssteqr_(char* compz, int* n, float* d, float* e, float* z, int* ldz, float* work, int* info);