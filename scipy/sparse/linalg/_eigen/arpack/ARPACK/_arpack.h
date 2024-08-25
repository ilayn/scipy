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
 *
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
    ido_OPX_V0 = -1,  // compute  Y = OP * X to force the starting vector into the range of OP
    ido_FIRST = 0,    // compute  Y = OP * X
    ido_OPX = 1,      // compute  Y = B * X
    ido_BX = 2,       // compute the real and imaginary parts of the shifts
    ido_DONE = 99     // done
};

/*
 * With the following structs, we collect all "SAVE"d Fortran variables to track
 * the problem and avoid reentry issues. It is not clean and is laborious but
 * otherwise a full rewrite of ARPACK is needed. There are additional variables
 * in the original Fortran code which are also "SAVE"d however upon inspection,
 * they are assigned and then used in the same call and thus used without saving.
*/

struct ARPACK_snaupd_variables {
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
    int naitr_iter;          // naitr flow control
    int naitr_j;             // naitr flow control
    int naitr_orth1;         // naitr flow control
    int naitr_orth2;         // naitr flow control
    int naitr_restart;       // naitr flow control
    int naitr_step3;         // naitr flow control
    int naitr_step4;         // naitr flow control
    int naitr_ierr;          // naitr flow control
    int naup2_initv;         // naupd2 flow control
    int naup2_iter;          // naupd2 flow control
    int naup2_getv0;         // naupd2 flow control
    int naup2_cnorm;         // naupd2 flow control
    int naup2_kplusp;        // naupd2 flow control
    int naup2_nev0;          // naupd2 internal compute
    int naup2_np0;           // naupd2 internal compute
    int naup2_numcnv;        // naupd2 internal compute
    int naup2_update;        // naupd2 flow control
    int naup2_ushift;        // naupd2 flow control
    float tol;               // problem parameter
    float getv0_rnorm0;      // getv0 internal compute
    float naitr_betaj;       // naitr internal compute
    float naitr_rnorm1;      // naitr internal compute
    float naitr_wnorm;       // naitr internal compute
    float naup2_rnorm;       // naup2 internal compute
};


struct ARPACK_dnaupd_variables {
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
    int naitr_ierr;          // naitr flow control
    int naitr_iter;          // naitr flow control
    int naitr_j;             // naitr flow control
    int naitr_orth1;         // naitr flow control
    int naitr_orth2;         // naitr flow control
    int naitr_restart;       // naitr flow control
    int naitr_step3;         // naitr flow control
    int naitr_step4;         // naitr flow control
    int naup2_initv;         // naupd2 flow control
    int naup2_iter;          // naupd2 flow control
    int naup2_getv0;         // naupd2 flow control
    int naup2_cnorm;         // naupd2 flow control
    int naup2_kplusp;        // naupd2 flow control
    int naup2_nev0;          // naupd2 internal compute
    int naup2_np0;           // naupd2 internal compute
    int naup2_numcnv;        // naupd2 internal compute
    int naup2_update;        // naupd2 flow control
    int naup2_ushift;        // naupd2 flow control
    double tol;              // problem parameter
    double getv0_rnorm0;     // getv0 internal compute
    double naitr_betaj;      // naitr internal compute
    double naitr_rnorm1;     // naitr internal compute
    double naitr_wnorm;      // naitr internal compute
    double naup2_rnorm;      // naup2 internal compute
};


struct ARPACK_cnaupd_variables {
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
    int naitr_ierr;          // naitr flow control
    int naitr_iter;          // naitr flow control
    int naitr_j;             // naitr flow control
    int naitr_orth1;         // naitr flow control
    int naitr_orth2;         // naitr flow control
    int naitr_restart;       // naitr flow control
    int naitr_step3;         // naitr flow control
    int naitr_step4;         // naitr flow control
    int naup2_initv;         // naupd2 flow control
    int naup2_iter;          // naupd2 flow control
    int naup2_getv0;         // naupd2 flow control
    int naup2_cnorm;         // naupd2 flow control
    int naup2_kplusp;        // naupd2 flow control
    int naup2_nev0;          // naupd2 internal compute
    int naup2_np0;           // naupd2 internal compute
    int naup2_numcnv;        // naupd2 internal compute
    int naup2_update;        // naupd2 flow control
    int naup2_ushift;        // naupd2 flow control
    float tol;              // problem parameter
    float getv0_rnorm0;     // getv0 internal compute
    float naitr_betaj;      // naitr internal compute
    float naitr_rnorm1;     // naitr internal compute
    float naitr_wnorm;      // naitr internal compute
    float naup2_rnorm;      // naup2 internal compute
};


struct ARPACK_znaupd_variables {
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
    int naitr_ierr;          // naitr flow control
    int naitr_iter;          // naitr flow control
    int naitr_j;             // naitr flow control
    int naitr_orth1;         // naitr flow control
    int naitr_orth2;         // naitr flow control
    int naitr_restart;       // naitr flow control
    int naitr_step3;         // naitr flow control
    int naitr_step4;         // naitr flow control
    int naup2_initv;         // naupd2 flow control
    int naup2_iter;          // naupd2 flow control
    int naup2_getv0;         // naupd2 flow control
    int naup2_cnorm;         // naupd2 flow control
    int naup2_kplusp;        // naupd2 flow control
    int naup2_nev0;          // naupd2 internal compute
    int naup2_np0;           // naupd2 internal compute
    int naup2_numcnv;        // naupd2 internal compute
    int naup2_update;        // naupd2 flow control
    int naup2_ushift;        // naupd2 flow control
    double tol;              // problem parameter
    double getv0_rnorm0;     // getv0 internal compute
    double naitr_betaj;      // naitr internal compute
    double naitr_rnorm1;     // naitr internal compute
    double naitr_wnorm;      // naitr internal compute
    double naup2_rnorm;      // naup2 internal compute
};


struct ARPACK_dsaupd_variables {
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
    int naitr_ierr;          // naitr flow control
    int naitr_iter;          // naitr flow control
    int naitr_j;             // naitr flow control
    int naitr_orth1;         // naitr flow control
    int naitr_orth2;         // naitr flow control
    int naitr_restart;       // naitr flow control
    int naitr_step3;         // naitr flow control
    int naitr_step4;         // naitr flow control
    int naup2_initv;         // naupd2 flow control
    int naup2_iter;          // naupd2 flow control
    int naup2_getv0;         // naupd2 flow control
    int naup2_cnorm;         // naupd2 flow control
    int naup2_kplusp;        // naupd2 flow control
    int naup2_nev0;          // naupd2 internal compute
    int naup2_np0;           // naupd2 internal compute
    int naup2_numcnv;        // naupd2 internal compute
    int naup2_update;        // naupd2 flow control
    int naup2_ushift;        // naupd2 flow control
    double tol;              // problem parameter
    double getv0_rnorm0;     // getv0 internal compute
    double naitr_betaj;      // naitr internal compute
    double naitr_rnorm1;     // naitr internal compute
    double naitr_wnorm;      // naitr internal compute
    double naup2_rnorm;      // naup2 internal compute
};

/*
 * The following function is extracted out from classical ARPACK code in order
 * to replace the random number generation with NumPy Random C-API and more
 * importantly to be able to control the seed which is not possible with the
 * classical Fortran code. It can replaced with some other random generator
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