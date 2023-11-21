from . cimport sf_error

from libc.math cimport NAN

from ._cdflib cimport cdfbet_which3


cdef inline double btdtria(double p, double b, double x) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double y = 1.0 - x
        double result, bound
        int status
        char arg;

    result, status, bound = cdfbet_which3(p, q, x, y, b)
    if status < 0:
        arg = b'p' if status == -1 else b'x' if status == -2 else b'b'
        sf_error.error("btdtria", sf_error.ARG,
                       "Input parameter %c is out of range", arg)
        return NAN
    if status == 0:
        return result
    if status == 1:
        sf_error.error("btdtria", sf_error.OTHER,
                       "Answer appears to be lower than lowest search bound (%g)", bound)
        return bound
    if status == 2:
        sf_error.error("btdtria", sf_error.OTHER,
                       "Answer appears to be higher than highest search bound (%g)", bound)
        return bound
    if status == 3 or status == 4:
        sf_error.error("btdtria", sf_error.OTHER,
                       "Two parameters that should sum to 1.0 do not.")
        return NAN
    if status == 10:
        sf_error.error("btdtria", sf_error.OTHER, "Computational error")
        return NAN
    sf_error.error("btdtria", sf_error.OTHER, "Unknown error.")
    return NAN
