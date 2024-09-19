#ifndef _MATFUNCS_SQRTM_H
#define _MATFUNCS_SQRTM_H

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <math.h>
#include "numpy/arrayobject.h"

#if defined(_MSC_VER)
    // MSVC nonsense
    #include <complex.h>
    #define EXPM_Z _Dcomplex
    #define EXPM_C _Fcomplex
#else
    // C99-compliant compilers
    #include <complex.h>
    #define EXPM_Z double complex
    #define EXPM_C float complex
#endif

typedef int user_select(const double, const double);

// LAPACK functions
void sgees_(char* jobvs, char* sort, user_select* select, int* n, float* a,
            int* lda, int* sdim, float* wr, float* wi, float* vs, int* ldvs,
            float* work, int* lwork, int* bwork, int* info);
void dgees_(char* jobvs, char* sort, user_select* select, int* n, double* a,
            int* lda, int* sdim, double* wr, double* wi, double* vs, int* ldvs,
            double* work, int* lwork, int* bwork, int* info);
void cgees_(char* jobvs, char* sort, user_select* select, int* n, EXPM_C* a,
            int* lda, int* sdim, EXPM_C* w, EXPM_C vs, int* ldvs, EXPM_C* work,
            int* lwork, float* rwork, int* bwork, int* info);
void zgees_(char* jobvs, char* sort, user_select* select, int* n, EXPM_Z* a,
            int* lda, int* sdim, EXPM_Z* w, EXPM_Z vs, int* ldvs, EXPM_Z* work,
            int* lwork, double* rwork, int* bwork, int* info);



#undef EXPM_C
#undef EXPM_Z

#endif // _MATFUNCS_EXPM_H
