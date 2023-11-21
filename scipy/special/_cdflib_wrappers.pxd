from . cimport sf_error

from libc.math cimport NAN

from ._cdflib cimport (
    cdfbet_which3,
    cdfbet_which4,
)


cdef inline double get_result(
        char *name,
        char **argnames,
        double result,
        int status,
        double bound,
        int return_bound
) noexcept nogil:
    cdef char *arg
    """Get result and perform error handling from cdflib output."""
    if status < 0:
        arg = argnames[-(status + 1)]
        sf_error.error(name, sf_error.ARG,
                       "Input parameter %s is out of range", arg)
        return NAN
    if status == 0:
        return result
    if status == 1:
        sf_error.error(name, sf_error.OTHER,
                       "Answer appears to be lower than lowest search bound (%g)", bound)
        return bound if return_bound else NAN
    if status == 2:
        sf_error.error(name, sf_error.OTHER,
                       "Answer appears to be higher than highest search bound (%g)", bound)
        return bound if return_bound else NAN
    if status == 3 or status == 4:
        sf_error.error(name, sf_error.OTHER,
                       "Two internal parameters that should sum to 1.0 do not.")
        return NAN
    if status == 10:
        sf_error.error(name, sf_error.OTHER, "Computational error")
        return NAN
    sf_error.error(name, sf_error.OTHER, "Unknown error.")
    return NAN
    


cdef inline double btdtria(double p, double b, double x) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double y = 1.0 - x
        double result, bound
        int status
        char *argnames[5]

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "x"
    argnames[3] = "y"
    argnames[4] = "b"

    result, status, bound = cdfbet_which3(p, q, x, y, b)
    return get_result("btdtria", argnames, result, status, bound, 1)



cdef inline double btdtrib(double p, double a, double x) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double y = 1.0 - x
        double result, bound
        int status
        char *argnames[5]

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "x"
    argnames[3] = "y"
    argnames[4] = "a"

    result, status, bound = cdfbet_which4(p, q, x, y, a)
    return get_result("btdtrib", argnames, result, status, bound, 1)
