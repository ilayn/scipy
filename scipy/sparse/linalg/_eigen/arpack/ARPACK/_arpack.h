#ifndef _ARPACK_H
#define _ARPACK_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "numpy/arrayobject.h"
#include <math.h>

#if defined(_MSC_VER)
    // MSVC definitions
    #include <complex.h>  // MSVC C++ header
    typedef _Dcomplex ARPACK_CPLX_TYPE;
    typedef _Fcomplex ARPACK_CPLXF_TYPE;
    #define ARPACK_cplx(real, imag) ((_Dcomplex){real, imag})
    #define ARPACK_cplxf(real, imag) ((_Fcomplex){real, imag})

#else
    // C99 compliant compilers
    #include <complex.h>
    typedef double complex ARPACK_CPLX_TYPE;
    typedef float complex ARPACK_CPLXF_TYPE;
    #define ARPACK_cplx(real, imag) ((real) + (imag)*I)
    #define ARPACK_cplxf(real, imag) ((real) + (imag)*I)

#endif



/*
 * ARPACK uses the so-called reverse-communication style that typically exits
 * the program with its in/out arguments to signal, in what stage the algorithm
 * is and what it needs. Then user modifies the arguments and calls again with
 * the necessary information. Thus the state of the whole program is sent back
 * and forth through in-place modified arguments. On top of this, ARPACK also
 * uses lots of variables through the Fortran's dreadful SAVE attribute that
 * persists the variable values across runs. Instead we move all those variables
 * into the reverse communication layer.
 *
 * For scalar arguments we use a C struct <-> Python dict bridge and for array
 * arguments we use NumPy array pointers that are originated from the user side
 * to modify things in-place and not to deal with ref counting or alloc/free.
 *
 * To generate random vectors that are used in ARPACK algorithm, we also expect
 * a NumPy generator object from the user for reproducible runs. However this
 * can be replaced with a different number generation routine.
*/

enum ARPACK_which {
    which_LM = 0,    // want the NEV eigenvalues of largest magnitude.
    which_SM = 1,    // want the NEV eigenvalues of smallest magnitude.
    which_LR = 2,    // want the NEV eigenvalues of largest real part.
    which_SR = 3,    // want the NEV eigenvalues of smallest real part.
    which_LI = 4,    // want the NEV eigenvalues of largest imaginary part.
    which_SI = 5,    // want the NEV eigenvalues of smallest imaginary part.
    which_LA = 6,    // compute the NEV largest (algebraic) eigenvalues. (sym)
    which_SA = 7,    // compute the NEV smallest (algebraic) eigenvalues. (sym)
    which_BE = 8     // compute NEV eigenvalues, half from each end of the spectrum. (sym)
};

enum ARPACK_ido {
    ido_OPX_V0 = -1,
    ido_FIRST = 0,
    ido_OPX = 1,
    ido_BX = 2,
    ido_USER_SHIFT = 3,
    ido_DONE = 99
};

/**
 * With the following structs, we collect all "SAVE"d Fortran variables to track
 * the problem and avoid reentry issues. It is not clean and is laborious but
 * otherwise a full rewrite of ARPACK is needed. There are additional variables
 * in the original Fortran code which are also "SAVE"d however upon inspection,
 * they are assigned and then used in the same call and thus used without saving.
**/

struct ARPACK_arnoldi_update_vars_s {
    int bmat;                // bmat = 0, standard, bmat = 1 generalized problem
    int info;                // problem outcome
    int ishift;              // problem parameter
    int maxiter;             // problem parameter
    int mode;                // problem parameter
    int n;                   // problem parameter
    int nconv;               // problem outcome
    int ncv;                 // problem parameter
    int nev;                 // problem parameter
    int np;                  // problem intermediate
    int numop;               // problem intermediate
    int numpb;               // problem intermediate
    int numreo;              // problem intermediate
    enum ARPACK_ido ido;     // naupd flow control
    enum ARPACK_which which; // naupd flow control
    int getv0_first;         // getv0 flow control
    int getv0_iter;          // getv0 flow control
    int getv0_itry;          // getv0 flow control
    int getv0_orth;          // getv0 flow control
    int aitr_iter;          // naitr flow control
    int aitr_j;             // naitr flow control
    int aitr_orth1;         // naitr flow control
    int aitr_orth2;         // naitr flow control
    int aitr_restart;       // naitr flow control
    int aitr_step3;         // naitr flow control
    int aitr_step4;         // naitr flow control
    int aitr_ierr;          // naitr flow control
    int aup2_initv;         // naupd2 flow control
    int aup2_iter;          // naupd2 flow control
    int aup2_getv0;         // naupd2 flow control
    int aup2_cnorm;         // naupd2 flow control
    int aup2_kplusp;        // naupd2 flow control
    int aup2_nev0;          // naupd2 internal compute
    int aup2_np0;           // naupd2 internal compute
    int aup2_numcnv;        // naupd2 internal compute
    int aup2_update;        // naupd2 flow control
    int aup2_ushift;        // naupd2 flow control
    float tol;               // problem parameter
    float getv0_rnorm0;      // getv0 internal compute
    float aitr_betaj;       // naitr internal compute
    float aitr_rnorm1;      // naitr internal compute
    float aitr_wnorm;       // naitr internal compute
    float aup2_rnorm;       // naup2 internal compute
};


struct ARPACK_arnoldi_update_vars_d {
    int bmat;                // bmat = 0, standard, bmat = 1 generalized problem
    int info;                // problem outcome
    int ishift;              // problem parameter
    int maxiter;             // problem parameter
    int mode;                // problem parameter
    int n;                   // problem parameter
    int nconv;               // problem outcome
    int ncv;                 // problem parameter
    int nev;                 // problem parameter
    int np;                  // problem intermediate
    int numop;               // problem intermediate
    int numpb;               // problem intermediate
    int numreo;              // problem intermediate
    enum ARPACK_ido ido;     // naupd flow control
    enum ARPACK_which which; // naupd flow control
    int getv0_first;         // getv0 flow control
    int getv0_iter;          // getv0 flow control
    int getv0_itry;          // getv0 flow control
    int getv0_orth;          // getv0 flow control
    int aitr_ierr;          // naitr flow control
    int aitr_iter;          // naitr flow control
    int aitr_j;             // naitr flow control
    int aitr_orth1;         // naitr flow control
    int aitr_orth2;         // naitr flow control
    int aitr_restart;       // naitr flow control
    int aitr_step3;         // naitr flow control
    int aitr_step4;         // naitr flow control
    int aup2_initv;         // naupd2 flow control
    int aup2_iter;          // naupd2 flow control
    int aup2_getv0;         // naupd2 flow control
    int aup2_cnorm;         // naupd2 flow control
    int aup2_kplusp;        // naupd2 flow control
    int aup2_nev0;          // naupd2 internal compute
    int aup2_np0;           // naupd2 internal compute
    int aup2_numcnv;        // naupd2 internal compute
    int aup2_update;        // naupd2 flow control
    int aup2_ushift;        // naupd2 flow control
    double tol;              // problem parameter
    double getv0_rnorm0;     // getv0 internal compute
    double aitr_betaj;      // naitr internal compute
    double aitr_rnorm1;     // naitr internal compute
    double aitr_wnorm;      // naitr internal compute
    double aup2_rnorm;      // naup2 internal compute
};


void snaupd(struct ARPACK_arnoldi_update_vars_s *V, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);
void dnaupd(struct ARPACK_arnoldi_update_vars_d *V, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);

/*
 * The following function is extracted out from classical ARPACK code in order
 * to replace the random number generation with NumPy Random C-API and more
 * importantly to be able to control the seed which is not possible with the
 * classical Fortran code. It can be replaced with some other random generator
 * implementation following the prototype below.
*/
void generate_random_vector_s(const int n, double* vec);
void generate_random_vector_d(const int n, double* vec);
void generate_random_vector_c(const int n, double* vec);
void generate_random_vector_z(const int n, double* vec);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* ifndef */