from . cimport sf_error

from libc.math cimport NAN

from ._cdflib cimport (
    cdfbet_which3,
    cdfbet_which4,
    cdfbin_which2,
    cdfbin_which3,
    cdfchi_which3,
    cdfchn_which1,
    cdfchn_which2,
    cdfchn_which3,
    cdfchn_which4,
    cdff_which4,
    cdfgam_which2,
    cdfgam_which3,
    cdfgam_which4,
    cdfnbn_which2,
    cdfnbn_which3,
    cdffnc_which1,
    cdffnc_which2,
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
        int status = 10
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
        int status = 10
        char *argnames[5]

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "x"
    argnames[3] = "y"
    argnames[4] = "a"

    result, status, bound = cdfbet_which4(p, q, x, y, a)
    return get_result("btdtrib", argnames, result, status, bound, 1)


cdef inline double bdtrik(double p, double xn, double pr) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double ompr = 1.0 - pr
        double result, bound
        int status
        char *argnames[5]

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "xn"
    argnames[3] = "pr"
    argnames[4] = "ompr"

    result, status, bound = cdfbin_which2(p, q, xn, pr, ompr)
    return get_result("btdtrik", argnames, result, status, bound, 1)


cdef inline double bdtrin(double s, double p, double pr) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double ompr = 1.0 - pr
        double result, bound
        int status = 10
        char *argnames[5]

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "s"
    argnames[3] = "pr"
    argnames[4] = "ompr"

    result, status, bound = cdfbin_which3(p, q, s, pr, ompr)
    return get_result("btdtrin", argnames, result, status, bound, 1)


cdef inline double chdtriv(double p, double x) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[3]

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "x"

    result, status, bound = cdfchi_which3(p, q, x)
    return get_result("chdtriv", argnames, result, status, bound, 1)


cdef inline double chndtr(double x, double df, double nc) noexcept nogil:
    cdef:
        double result, _, bound
        int status = 10
        char *argnames[3]

    argnames[0] = "x"
    argnames[1] = "df"
    argnames[2] = "nc"

    result, _, status, bound = cdfchn_which1(x, df, nc)
    return get_result("chndtr", argnames, result, status, bound, 1)


cdef inline double chndtridf(double x, double p, double nc) noexcept nogil:
    cdef:
        double result, bound
        int status = 10
        char *argnames[3]

    argnames[0] = "p"
    argnames[1] = "x"
    argnames[2] = "nc"

    result, status, bound = cdfchn_which3(p, x, nc)
    return get_result("chndtridf", argnames, result, status, bound, 1)


cdef inline double chndtrinc(double x, double df, double p) noexcept nogil:
    cdef:
        double result, bound
        int status = 10
        char *argnames[3]

    argnames[0] = "p"
    argnames[1] = "x"
    argnames[2] = "df"

    result, status, bound = cdfchn_which4(p, x, df)
    return get_result("chndtrinc", argnames, result, status, bound, 1)


cdef inline double chndtrix(double p, double df, double nc) noexcept nogil:
    cdef:
        double result, bound
        int status = 10
        char *argnames[3]

    argnames[0] = "p"
    argnames[1] = "df"
    argnames[2] = "nc"

    result, status, bound = cdfchn_which2(p, df, nc)
    return get_result("chndtrix", argnames, result, status, bound, 0)


cdef inline double fdtridfd(double dfn, double p, double f) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "f"
    argnames[3] = "dfn"

    result, status, bound = cdff_which4(p, q, f, dfn)
    return get_result("fdtridfd", argnames, result, status, bound, 1)


cdef inline double gdtria(double p, double shp, double x) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "x"
    argnames[3] = "shp"

    result, status, bound = cdfgam_which4(p, q, x, shp)
    return get_result("gdtria", argnames, result, status, bound, 1)


cdef inline double gdtrix(double scl, double shp, double p) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "shp"
    argnames[3] = "scl"

    result, status, bound = cdfgam_which2(p, q, shp, scl)
    return get_result("gdtrix", argnames, result, status, bound, 1)


cdef inline double nbdtrik(double p, double xn, double pr) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double ompr = 1.0 - pr
        double result, bound
        int status = 10
        char *argnames[5]

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "xn"
    argnames[3] = "pr"
    argnames[4] = "ompr"

    result, status, bound = cdfnbn_which2(p, q, xn, pr, ompr)
    return get_result("nbdtrik", argnames, result, status, bound, 1)


cdef inline double nbdtrin(double s, double p, double pr) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double ompr = 1.0 - pr
        double result, bound
        int status = 10
        char *argnames[5]

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "s"
    argnames[3] = "pr"
    argnames[4] = "ompr"

    result, status, bound = cdfnbn_which3(p, q, s, pr, ompr)
    return get_result("nbdtrin", argnames, result, status, bound, 1)


cdef inline double ncfdtri(double dfn, double dfd, double nc, double p) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[5]

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[1] = "dfn"
    argnames[2] = "dfd"
    argnames[3] = "nc"

    result, status, bound = cdffnc_which2(p, q, dfn, dfd, nc)
    return get_result("ncfdtri", argnames, result, status, bound, 1)
