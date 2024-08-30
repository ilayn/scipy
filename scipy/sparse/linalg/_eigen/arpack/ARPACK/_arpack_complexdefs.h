#ifndef _ARPACK_COMPLEXDEFS_H
#define _ARPACK_COMPLEXDEFS_H

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

#endif // _ARPACK_COMPLEXDEFS_H