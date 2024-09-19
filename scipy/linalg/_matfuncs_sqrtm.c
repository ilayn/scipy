#include "_matfuncs_sqrtm.h"
#include <math.h>

// For some reason these defs are not picked up from the header.
#if defined(_MSC_VER)
    // MSVC nonsense
    #include <complex.h>
    #define EXPM_Z _Dcomplex
    #define EXPM_C _Fcomplex
    #define CPLX_Z(real, imag) (_Cbuild(real, imag))
    #define CPLX_C(real, imag) (_FCbuild(real, imag))
#else
    // C99-compliant compilers
    #include <complex.h>
    #define EXPM_Z double complex
    #define EXPM_C float complex
    #define CPLX_Z(real, imag) (real + imag*I)
    #define CPLX_C(real, imag) (real + imag*I)
#endif

void
sqrtm_2x2_s(float* t1, float* t2, float* t3, float* t4)
{
    // Solves 2x2 matrix square root by the formulas given in
    // N. Higham, "Computing Real Square Roots of a Real Matrix", 1987,
    // https://doi.org/10.1016/0024-3795(87)90118-2

    float alpha, mu;
    if (*t3 != 0.0)
    {
        // Quasi-triangular block
        mu = sqrtf(-(*t2)*(*t3)); // Guaranteed to be positive number due to sgees
        if (*t1 > 0.0)
        {
            alpha = sqrtf((*t1 + hypotf(*t1, mu)) / 2.0);
        } else {
            alpha = mu / sqrtf(2.0*(-(*t1) + hypotf(*t1, mu)));
        }
        *t1 = alpha;
        *t4 = alpha;
        *t3 /= 2.0*alpha;
        *t2 /= 2.0*alpha;
    } else {
        *t1 = sqrtf(*t1);
        *t4 = sqrtf(*t4);
        *t2 /= *t1 + *t4;
    }
}


void
sqrtm_2x2_d(double* t1, double* t2, double* t3, double* t4)
{
    // Solves 2x2 matrix square root by the formulas given in
    // N. Higham, "Computing Real Square Roots of a Real Matrix", 1987,
    // https://doi.org/10.1016/0024-3795(87)90118-2

    double alpha, mu;
    if (*t3 != 0.0)
    {
        // Quasi-triangular block
        mu = sqrt(-(*t2)*(*t3)); // Guaranteed to be positive number due to dgees
        if (*t1 > 0.0)
        {
            alpha = sqrt((*t1 + hypot(*t1, mu)) / 2.0);
        } else {
            alpha = mu / sqrt(2.0*(-(*t1) + hypot(*t1, mu)));
        }
        *t1 = alpha;
        *t4 = alpha;
        *t3 /= 2.0*alpha;
        *t2 /= 2.0*alpha;
    } else {
        *t1 = sqrt(*t1);
        *t4 = sqrt(*t4);
        *t2 /= *t1 + *t4;
    }
}


check_spectrum_s(double* a, int n, int* singular)
{
    int tmp_int;
    // Input "a" is a 2D C-contiguous array.
    // First a Schur decomposition is performed and if there is no eigenvalue
    // on the negative real axis, the real sqrt is computed. If "a" is singular
    // a warning is generated as the results can be wrong.


    sgees_("V", "N", NULL, &n, a, &n, &tmp_int, wr, wi, vs, &n, work, &lwork, bwork, &info);


}










#undef EXPM_Z
#undef EXPM_C
#undef CPLX_Z
#undef CPLX_C
