#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "numpy/arrayobject.h"
#include "ccallback.h"
#include "src/lsoda.h"
#include <math.h>

#define PyArray_MAX(a,b) (((a)>(b))?(a):(b))

#ifdef HAVE_BLAS_ILP64
#define F_INT npy_int64
#define F_INT_NPY NPY_INT64
#else
#define F_INT int
#define F_INT_NPY NPY_INT
#endif


#define PYERR(errobj,message) {\
    PyErr_SetString(errobj,message); \
    goto fail; \
}
#define PYERR2(errobj,message) { \
    PyErr_Print(); \
    PyErr_SetString(errobj, message); \
    goto fail; \
}

static PyObject *odepack_error;

static char doc_odeint[] =
    "[y,{infodict,}istate] = odeint(fun, y0, t, args=(), Dfun=None, "
    "col_deriv=0, ml=, mu=, full_output=0, rtol=, atol=, tcrit=, h0=0.0, "
    "hmax=0.0, hmin=0.0, ixpr=0.0, mxstep=0.0, mxhnil=0, mxordn=0, "
    "mxords=0)\n  yprime = fun(y,t,...)";

/*
 *  Copy a contiguous matrix at `c` to a Fortran-ordered matrix at `f`.
 *  `ldf` is the leading dimension of the Fortran array at `f`.
 *  `nrows` and `ncols` are the number of rows and columns of the matrix, resp.
 *  If `transposed` is 0, c[i, j] is *(c + ncols*i + j).
 *  If `transposed` is nonzero, c[i, j] is *(c + i + nrows*j)  (i.e. `c` is
 *  stored in F-contiguous order).
 */
static void
copy_array_to_fortran(double *f, int ldf, int nrows, int ncols,
                      double *c, int transposed)
{
    int i, j;
    int row_stride, col_stride;

    /* The strides count multiples of sizeof(double), not bytes. */
    if (transposed) {
        row_stride = 1;
        col_stride = nrows;
    }
    else {
        row_stride = ncols;
        col_stride = 1;
    }
    for (i = 0; i < nrows; ++i) {
        for (j = 0; j < ncols; ++j) {
            double value;
            /* value = c[i,j] */
            value = *(c + row_stride*i + col_stride*j);
            /* f[i,j] = value */
            *(f + ldf*j + i) = value;
        }
    }
}


// Callback Infrastructure

typedef struct {
    PyObject *ode_function;
    PyObject *jac_function;
    PyObject *func_args;
    int jac_transpose;
    int jac_type;
    int tfirst;               // Controls argument order: 0=(y,t,...), 1=(t,y,...) but why?
} odepack_callback_t;

// Thread-local storage for callbacks
static SCIPY_TLS odepack_callback_t* current_odepack_callback = NULL;


/*
 * Modern ODE function thunk - replaces the global variable approach
 * Called from C/Fortran code, interfaces with Python function
 */
static void
ode_function_thunk(int *n, double *t, double *y, double *ydot)
{
    if (!current_odepack_callback || !current_odepack_callback->ode_function) { *n = -1; return; }
    npy_intp dims[1] = {*n};

    PyObject *capsule = PyCapsule_New((void*)y, NULL, NULL);
    if (!capsule) { *n = -1; return; }

    PyObject* py_y = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, (void*)y);
    if (!py_y) {
        Py_DECREF(capsule);
        *n = -1;
        return;
    }

    if (PyArray_SetBaseObject((PyArrayObject*)py_y, capsule) < 0) {
        Py_DECREF(py_y);
        *n = -1;
        return;
    }

    PyObject *args_tuple;
    if (current_odepack_callback->tfirst == 0) {
        // Order: (y, t, *extra_args)
        if (current_odepack_callback->func_args) {
            Py_ssize_t extra_size = PyTuple_Size(current_odepack_callback->func_args);
            args_tuple = PyTuple_New(2 + extra_size);
            if (!args_tuple) {
                *n = -1;
                return;
            }
            PyTuple_SetItem(args_tuple, 0, py_y);
            PyTuple_SetItem(args_tuple, 1, PyFloat_FromDouble(*t));
            // Copy extra arguments
            for (Py_ssize_t i = 0; i < extra_size; i++) {
                PyObject *item = PyTuple_GetItem(current_odepack_callback->func_args, i);
                Py_INCREF(item);
                PyTuple_SetItem(args_tuple, 2 + i, item);
            }
        } else {
            // No extra arguments
            args_tuple = PyTuple_New(2);
            if (!args_tuple) {
                *n = -1;
                return;
            }
            PyTuple_SetItem(args_tuple, 0, py_y);
            PyTuple_SetItem(args_tuple, 1, PyFloat_FromDouble(*t));
        }
    } else {
        // Alternative order: (t, y, *extra_args)
        if (current_odepack_callback->func_args) {
            Py_ssize_t extra_size = PyTuple_Size(current_odepack_callback->func_args);
            args_tuple = PyTuple_New(2 + extra_size);
            if (!args_tuple) {
                *n = -1;
                return;
            }
            PyTuple_SetItem(args_tuple, 0, PyFloat_FromDouble(*t));
            PyTuple_SetItem(args_tuple, 1, py_y);
            // Copy extra arguments
            for (Py_ssize_t i = 0; i < extra_size; i++) {
                PyObject *item = PyTuple_GetItem(current_odepack_callback->func_args, i);
                Py_INCREF(item);
                PyTuple_SetItem(args_tuple, 2 + i, item);
            }
        } else {
            // No extra arguments
            args_tuple = PyTuple_New(2);
            if (!args_tuple) {
                *n = -1;
                return;
            }
            PyTuple_SetItem(args_tuple, 0, PyFloat_FromDouble(*t));
            PyTuple_SetItem(args_tuple, 1, py_y);
        }
    }

    // Call Python function
    PyObject *result = PyObject_CallObject(current_odepack_callback->ode_function, args_tuple);
    Py_DECREF(args_tuple);
    if (!result) { *n = -1; return; }

    // Process result and validate dimensions
    PyArrayObject *result_array = (PyArrayObject*)PyArray_FROM_OTF(result, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (!result_array) {
        Py_DECREF(result);
        *n = -1;
        return;
    }

    if (PyArray_NDIM(result_array) > 1) {
        PyErr_Format(PyExc_RuntimeError,
                "The array returned by func must be one-dimensional, but got ndim=%d.",
                PyArray_NDIM(result_array));
        Py_DECREF(result_array);
        Py_DECREF(result);
        *n = -1;
        return;
    }

    if (PyArray_Size((PyObject *)result_array) != *n) {
        PyErr_Format(PyExc_RuntimeError,
            "The size of the array returned by func (%ld) does not match "
            "the size of y0 (%d).",
            PyArray_Size((PyObject *)result_array), *n);
        Py_DECREF(result_array);
        Py_DECREF(result);
        *n = -1;
        return;
    }

    // Copy result to output array - manual loop instead of memcpy
    double *result_data = (double*)PyArray_DATA(result_array);
    for (int i = 0; i < *n; i++) {
        ydot[i] = result_data[i];
    }

    Py_DECREF(result_array);
    Py_DECREF(result);
}


static void
ode_jacobian_thunk(int *n, double *t, double *y, int *ml, int *mu, double *pd, int *nrowpd)
{
    if (!current_odepack_callback || !current_odepack_callback->jac_function) {
        *n = -1;
        return;
    }

    npy_intp dims[1] = {*n};

    // Create a PyCapsule to manage the data lifetime
    PyObject *capsule = PyCapsule_New((void*)y, NULL, NULL);
    if (!capsule) {
        *n = -1;
        return;
    }

    PyObject* py_y = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, (void*)y);
    if (!py_y) {
        Py_DECREF(capsule);
        *n = -1;
        return;
    }

    if (PyArray_SetBaseObject((PyArrayObject*)py_y, capsule) < 0) {
        Py_DECREF(py_y);
        *n = -1;
        return;
    }

    PyObject *jac_args_tuple;
    if (current_odepack_callback->tfirst == 0) {
        // Standard ODEPACK order: (y, t, *extra_args)
        if (current_odepack_callback->func_args) {
            // There are extra arguments
            Py_ssize_t extra_size = PyTuple_Size(current_odepack_callback->func_args);
            jac_args_tuple = PyTuple_New(2 + extra_size);
            if (!jac_args_tuple) {
                *n = -1;
                return;
            }
            PyTuple_SetItem(jac_args_tuple, 0, py_y);
            PyTuple_SetItem(jac_args_tuple, 1, PyFloat_FromDouble(*t));
            // Copy extra arguments
            for (Py_ssize_t i = 0; i < extra_size; i++) {
                PyObject *item = PyTuple_GetItem(current_odepack_callback->func_args, i);
                Py_INCREF(item);
                PyTuple_SetItem(jac_args_tuple, 2 + i, item);
            }
        } else {
            // No extra arguments
            jac_args_tuple = PyTuple_New(2);
            if (!jac_args_tuple) {
                *n = -1;
                return;
            }
            PyTuple_SetItem(jac_args_tuple, 0, py_y);
            PyTuple_SetItem(jac_args_tuple, 1, PyFloat_FromDouble(*t));
        }
    } else {
        // Alternative order: (t, y, *extra_args)
        if (current_odepack_callback->func_args) {
            // There are extra arguments
            Py_ssize_t extra_size = PyTuple_Size(current_odepack_callback->func_args);
            jac_args_tuple = PyTuple_New(2 + extra_size);
            if (!jac_args_tuple) {
                *n = -1;
                return;
            }
            PyTuple_SetItem(jac_args_tuple, 0, PyFloat_FromDouble(*t));
            PyTuple_SetItem(jac_args_tuple, 1, py_y);
            // Copy extra arguments
            for (Py_ssize_t i = 0; i < extra_size; i++) {
                PyObject *item = PyTuple_GetItem(current_odepack_callback->func_args, i);
                Py_INCREF(item);
                PyTuple_SetItem(jac_args_tuple, 2 + i, item);
            }
        } else {
            // No extra arguments
            jac_args_tuple = PyTuple_New(2);
            if (!jac_args_tuple) {
                *n = -1;
                return;
            }
            PyTuple_SetItem(jac_args_tuple, 0, PyFloat_FromDouble(*t));
            PyTuple_SetItem(jac_args_tuple, 1, py_y);
        }
    }

    // Call Python Jacobian function
    PyObject *result = PyObject_CallObject(current_odepack_callback->jac_function, jac_args_tuple);
    Py_DECREF(jac_args_tuple);
    if (!result) {
        *n = -1;
        return;
    }

    PyArrayObject *result_array = (PyArrayObject*)PyArray_FROM_OTF(result, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (!result_array) {
        Py_DECREF(result);
        *n = -1;
        return;
    }

    // Calculate expected dimensions based on Jacobian type
    npy_intp expected_nrows, expected_ncols;
    expected_ncols = *n;
    if (current_odepack_callback->jac_type == 4) {
        // Banded Jacobian
        expected_nrows = *ml + *mu + 1;
    } else {
        // Full Jacobian
        expected_nrows = *n;
    }

    // Handle transpose case
    if (!current_odepack_callback->jac_transpose) {
        npy_intp tmp = expected_nrows;
        expected_nrows = expected_ncols;
        expected_ncols = tmp;
    }

    // Validate result dimensions
    npy_intp ndim = PyArray_NDIM(result_array);
    if (ndim > 2) {
        PyErr_Format(PyExc_RuntimeError,
            "The Jacobian array must be two dimensional, but got ndim=%ld.", ndim);
        Py_DECREF(result_array);
        Py_DECREF(result);
        *n = -1;
        return;
    }

    npy_intp *dims_result = PyArray_DIMS(result_array);
    int dim_error = 0;
    if (ndim == 0) {
        if ((expected_nrows != 1) || (expected_ncols != 1)) {
            dim_error = 1;
        }
    } else if (ndim == 1) {
        if ((expected_nrows != 1) || (dims_result[0] != expected_ncols)) {
            dim_error = 1;
        }
    } else if (ndim == 2) {
        if ((dims_result[0] != expected_nrows) || (dims_result[1] != expected_ncols)) {
            dim_error = 1;
        }
    }

    if (dim_error) {
        const char *jac_type_str = (current_odepack_callback->jac_type == 4) ? "banded " : "";
        PyErr_Format(PyExc_RuntimeError,
            "Expected a %sJacobian array with shape (%ld, %ld), but got (%ld, %ld)",
            jac_type_str, expected_nrows, expected_ncols,
            (ndim >= 1) ? dims_result[0] : 1,
            (ndim >= 2) ? dims_result[1] : ((ndim == 1) ? dims_result[0] : 1));
        Py_DECREF(result_array);
        Py_DECREF(result);
        *n = -1;
        return;
    }

    // Copy result to output array with proper layout handling
    if ((current_odepack_callback->jac_type == 1) && !current_odepack_callback->jac_transpose) {
        // Full Jacobian, no transpose needed, use manual loop
        double *src_data = (double*)PyArray_DATA(result_array);
        for (ssize_t i = 0; i < (*n) * (*nrowpd); i++) {
            pd[i] = src_data[i];
        }
    } else {
        // Need to copy with proper Fortran layout
        npy_intp m = (current_odepack_callback->jac_type == 4) ? (*ml + *mu + 1) : *n;
        copy_array_to_fortran(pd, *nrowpd, m, *n,
                             (double*)PyArray_DATA(result_array),
                             !current_odepack_callback->jac_transpose);
    }

    Py_DECREF(result_array);
    Py_DECREF(result);
}

// ================================================================================


static int
setup_odepack_callback(odepack_callback_t *callback,
                       PyObject *ode_func, PyObject *jac_func, PyObject *extra_args,
                       int jac_transpose, int jac_type, int tfirst)
{
    // Initialize callback structure
    callback->ode_function = ode_func;
    callback->jac_function = jac_func;
    callback->jac_transpose = jac_transpose;
    callback->jac_type = jac_type;
    callback->tfirst = tfirst;
    callback->func_args = NULL;

    // Increment reference counts for Python objects
    if (ode_func) {
        Py_INCREF(ode_func);
    }
    if (jac_func && jac_func != Py_None) {
        Py_INCREF(jac_func);
    }

    /*
     *
     * Both ode_function and jac_function will be called with:
     *   - tfirst=0: f = func(y, t, *extra_args)
     *   - tfirst=1: f = func(t, y, *extra_args)
     *
     */
    if (extra_args && PyTuple_Check(extra_args)) {
        Py_INCREF(extra_args);
        callback->func_args = extra_args;
    } else {
        // No extra arguments - leave as NULL
        callback->func_args = NULL;
    }

    return 0;
}

// Clean up callback
static void
cleanup_odepack_callback(odepack_callback_t *callback)
{
    if (callback->func_args) {
        Py_DECREF(callback->func_args);
        callback->func_args = NULL;
    }

    if (callback->ode_function) {
        Py_DECREF(callback->ode_function);
        callback->ode_function = NULL;
    }

    if (callback->jac_function && callback->jac_function != Py_None) {
        Py_DECREF(callback->jac_function);
        callback->jac_function = NULL;
    }
}

// Activate callback for use in thunks - sets thread-local storage
static void
activate_odepack_callback(odepack_callback_t *callback) { current_odepack_callback = callback; }


static void
deactivate_odepack_callback(void) { current_odepack_callback = NULL; }

// ================================================================================


int
setup_extra_inputs(PyArrayObject **ap_rtol, PyObject *o_rtol,
                   PyArrayObject **ap_atol, PyObject *o_atol,
                   PyArrayObject **ap_tcrit, PyObject *o_tcrit,
                   long *numcrit, int neq)
{
    int itol = 0;
    double tol = 1.49012e-8;
    npy_intp one = 1;

    /* Setup tolerances */
    if (o_rtol == NULL) {
        *ap_rtol = (PyArrayObject *) PyArray_SimpleNew(1, &one, NPY_DOUBLE);
        if (*ap_rtol == NULL) {
            PYERR2(odepack_error,"Error constructing relative tolerance.");
        }
        *(double *) PyArray_DATA(*ap_rtol) = tol;                /* Default */
    }
    else {
        *ap_rtol = (PyArrayObject *) PyArray_ContiguousFromObject(o_rtol,
                                                            NPY_DOUBLE, 0, 1);
        if (*ap_rtol == NULL) {
            PYERR2(odepack_error,"Error converting relative tolerance.");
        }
        /* XXX Fix the following. */
        if (PyArray_NDIM(*ap_rtol) == 0); /* rtol is scalar */
        else if (PyArray_DIMS(*ap_rtol)[0] == neq) {
            itol |= 2;      /* Set rtol array flag */
        }
        else {
            PYERR(odepack_error, "Tolerances must be an array of the same length as the\n     number of equations or a scalar.");
        }
    }

    if (o_atol == NULL) {
        *ap_atol = (PyArrayObject *) PyArray_SimpleNew(1, &one, NPY_DOUBLE);
        if (*ap_atol == NULL) {
            PYERR2(odepack_error,"Error constructing absolute tolerance");
        }
        *(double *)PyArray_DATA(*ap_atol) = tol;
    }
    else {
        *ap_atol = (PyArrayObject *) PyArray_ContiguousFromObject(o_atol, NPY_DOUBLE, 0, 1);
        if (*ap_atol == NULL) {
            PYERR2(odepack_error,"Error converting absolute tolerance.");
        }
        /* XXX Fix the following. */
        if (PyArray_NDIM(*ap_atol) == 0); /* atol is scalar */
        else if (PyArray_DIMS(*ap_atol)[0] == neq) {
            itol |= 1;        /* Set atol array flag */
        }
        else {
            PYERR(odepack_error,"Tolerances must be an array of the same length as the\n     number of equations or a scalar.");
        }
    }
    itol++;             /* increment to get correct value */

    /* Setup t-critical */
    if (o_tcrit != NULL) {
        *ap_tcrit = (PyArrayObject *) PyArray_ContiguousFromObject(o_tcrit, NPY_DOUBLE, 0, 1);
        if (*ap_tcrit == NULL) {
            PYERR2(odepack_error,"Error constructing critical times.");
        }
        *numcrit = PyArray_Size((PyObject *) (*ap_tcrit));
    }
    return itol;

fail:       /* Needed for use of PYERR */
    return -1;
}


int
compute_lrw_liw(int *lrw, int *liw, int neq, int jt, int ml, int mu,
                int mxordn, int mxords)
{
    int lrn, lrs, nyh, lmat;

    if (jt == 1 || jt == 2) {
        lmat = neq*neq + 2;
    }
    else if (jt == 4 || jt == 5) {
        lmat = (2*ml + mu + 1)*neq + 2;
    }
    else {
        PYERR(odepack_error,"Incorrect value for jt.");
    }

    if (mxordn < 0) {
        PYERR(odepack_error,"Incorrect value for mxordn.");
    }
    if (mxords < 0) {
        PYERR(odepack_error,"Incorrect value for mxords.");
    }
    nyh = neq;

    lrn = 20 + nyh*(mxordn+1) + 3*neq;
    lrs = 20 + nyh*(mxords+1) + 3*neq + lmat;

    *lrw = PyArray_MAX(lrn,lrs);
    *liw = 20 + neq;
    return 0;

fail:
    return -1;
}



static PyObject *
odepack_odeint(PyObject *dummy, PyObject *args, PyObject *kwdict)
{
    PyObject *fcn, *y0, *p_tout, *o_rtol = NULL, *o_atol = NULL;
    PyArrayObject *ap_y = NULL, *ap_yout = NULL;
    PyArrayObject *ap_rtol = NULL, *ap_atol = NULL;
    PyArrayObject *ap_tout = NULL;
    PyObject *extra_args = NULL;
    PyObject *Dfun = Py_None;
    int neq, itol = 1, itask = 1, istate = 1, iopt = 0, lrw, *iwork, liw, jt = 4;
    double *y, t, *tout, *rtol, *atol, *rwork;
    double h0 = 0.0, hmax = 0.0, hmin = 0.0;
    long ixpr = 0, mxstep = 0, mxhnil = 0, mxordn = 12, mxords = 5, ml = -1, mu = -1;
    long tfirst;
    PyObject *o_tcrit = NULL;
    PyArrayObject *ap_tcrit = NULL;
    PyArrayObject *ap_hu = NULL, *ap_tcur = NULL, *ap_tolsf = NULL, *ap_tsw = NULL;
    PyArrayObject *ap_nst = NULL, *ap_nfe = NULL, *ap_nje = NULL, *ap_nqu = NULL;
    PyArrayObject *ap_mused = NULL;
    long imxer = 0, lenrw = 0, leniw = 0, col_deriv = 0;
    npy_intp out_sz = 0, dims[2];
    long k, ntimes, crit_ind = 0;
    long allocated = 0, full_output = 0, numcrit = 0;
    long t0count;
    double *yout, *yout_ptr, *tout_ptr, *tcrit = NULL;
    double *wa;
    odepack_callback_t callback;  /* Stack-allocated callback structure */
    static char *kwlist[] = {"fun", "y0", "t", "args", "Dfun", "col_deriv",
                             "ml", "mu", "full_output", "rtol", "atol", "tcrit",
                             "h0", "hmax", "hmin", "ixpr", "mxstep", "mxhnil",
                             "mxordn", "mxords", "tfirst", NULL};


    if (!PyArg_ParseTupleAndKeywords(args, kwdict, "OOO|OOllllOOOdddllllll", kwlist,
                                     &fcn, &y0, &p_tout, &extra_args, &Dfun,
                                     &col_deriv, &ml, &mu, &full_output, &o_rtol, &o_atol,
                                     &o_tcrit, &h0, &hmax, &hmin, &ixpr, &mxstep, &mxhnil,
                                     &mxordn, &mxords, &tfirst)) {
        return NULL;
    }

    if (o_tcrit == Py_None) {
        o_tcrit = NULL;
    }
    if (o_rtol == Py_None) {
        o_rtol = NULL;
    }
    if (o_atol == Py_None) {
        o_atol = NULL;
    }

    /* Set up jt, ml, and mu */
    if (Dfun == Py_None) {
        /* set jt for internally generated */
        jt++;
    }
    if (ml < 0 && mu < 0) {
        /* neither ml nor mu given, mark jt for full jacobian */
        jt -= 3;
    }
    if (ml < 0) {
        /* if one but not both are given */
        ml = 0;
    }
    if (mu < 0) {
        mu = 0;
    }

    /* Setup modern callback infrastructure to replace global variables */
    memset(&callback, 0, sizeof(odepack_callback_t));

    if (extra_args == NULL) {
        if ((extra_args = PyTuple_New(0)) == NULL) {
            goto fail;
        }
    }
    else {
        Py_INCREF(extra_args);   /* We decrement on exit. */
    }
    if (!PyTuple_Check(extra_args)) {
        PYERR(odepack_error, "Extra arguments must be in a tuple.");
    }
    if (!PyCallable_Check(fcn) || (Dfun != Py_None && !PyCallable_Check(Dfun))) {
        PYERR(odepack_error, "The function and its Jacobian must be callable functions.");
    }

    /* Initialize modern callback structure */
    if (setup_odepack_callback(&callback, fcn, Dfun, extra_args,
                               !(col_deriv), jt, tfirst) < 0) {
        PYERR(odepack_error, "Failed to setup callback infrastructure.");
    }

    /* Initial input vector */
    ap_y = (PyArrayObject *) PyArray_ContiguousFromObject(y0, NPY_DOUBLE, 0, 0);
    if (ap_y == NULL) {
        goto fail;
    }
    if (PyArray_NDIM(ap_y) > 1) {
        PyErr_SetString(PyExc_ValueError, "Initial condition y0 must be one-dimensional.");
        goto fail;
    }
    y = (double *) PyArray_DATA(ap_y);
    neq = PyArray_Size((PyObject *) ap_y);
    dims[1] = neq;

    /* Set of output times for integration */
    ap_tout = (PyArrayObject *) PyArray_ContiguousFromObject(p_tout, NPY_DOUBLE, 0, 0);
    if (ap_tout == NULL) {
        goto fail;
    }
    if (PyArray_NDIM(ap_tout) > 1) {
        PyErr_SetString(PyExc_ValueError, "Output times t must be one-dimensional.");
        goto fail;
    }
    tout = (double *) PyArray_DATA(ap_tout);
    ntimes = PyArray_Size((PyObject *)ap_tout);
    dims[0] = ntimes;

    t0count = 0;
    if (ntimes > 0) {
        /* Copy tout[0] to t, and count how many times it occurs. */
        t = tout[0];
        t0count = 1;
        while ((t0count < ntimes) && (tout[t0count] == t)) {
            ++t0count;
        }
    }

    /* Set up array to hold the output evaluations*/
    ap_yout= (PyArrayObject *) PyArray_SimpleNew(2,dims,NPY_DOUBLE);
    if (ap_yout== NULL) {
        goto fail;
    }
    yout = (double *) PyArray_DATA(ap_yout);

    /* Copy initial vector into first row(s) of output */
    yout_ptr = yout;
    for (k = 0; k < t0count; ++k) {
        for (int i = 0; i < neq; i++) {
            yout_ptr[i] = y[i];
        }
        yout_ptr += neq;
    }

    itol = setup_extra_inputs(&ap_rtol, o_rtol, &ap_atol, o_atol, &ap_tcrit,
                              o_tcrit, &numcrit, neq);
    if (itol < 0 ) {
        goto fail;  /* Something didn't work */
    }
    rtol = (double *) PyArray_DATA(ap_rtol);
    atol = (double *) PyArray_DATA(ap_atol);

    /* Find size of working arrays*/
    if (compute_lrw_liw(&lrw, &liw, neq, jt, ml, mu, mxordn, mxords) < 0) {
        goto fail;
    }

    if ((wa = (double *)calloc(lrw*sizeof(double) + liw*sizeof(int), 1))==NULL) {
        PyErr_NoMemory();
        goto fail;
    }
    allocated = 1;
    rwork = wa;
    iwork = (int *)(wa + lrw);

    iwork[0] = ml;
    iwork[1] = mu;

    if (h0 != 0.0 || hmax != 0.0 || hmin != 0.0 || ixpr != 0 || mxstep != 0 ||
            mxhnil != 0 || mxordn != 0 || mxords != 0) {
        rwork[4] = h0;
        rwork[5] = hmax;
        rwork[6] = hmin;
        iwork[4] = ixpr;
        iwork[5] = mxstep;
        iwork[6] = mxhnil;
        iwork[7] = mxordn;
        iwork[8] = mxords;
        iopt = 1;
    }
    istate = 1;
    k = t0count;

    /* If full output make some useful output arrays */
    if (full_output) {
        out_sz = ntimes-1;
        ap_hu = (PyArrayObject *) PyArray_SimpleNew(1, &out_sz, NPY_DOUBLE);
        ap_tcur = (PyArrayObject *) PyArray_SimpleNew(1, &out_sz, NPY_DOUBLE);
        ap_tolsf = (PyArrayObject *) PyArray_SimpleNew(1, &out_sz, NPY_DOUBLE);
        ap_tsw = (PyArrayObject *) PyArray_SimpleNew(1, &out_sz, NPY_DOUBLE);
        ap_nst = (PyArrayObject *) PyArray_SimpleNew(1, &out_sz, NPY_INT);
        ap_nfe = (PyArrayObject *) PyArray_SimpleNew(1, &out_sz, NPY_INT);
        ap_nje = (PyArrayObject *) PyArray_SimpleNew(1, &out_sz, NPY_INT);
        ap_nqu = (PyArrayObject *) PyArray_SimpleNew(1, &out_sz, NPY_INT);
        ap_mused = (PyArrayObject *) PyArray_SimpleNew(1, &out_sz, NPY_INT);
        if (ap_hu == NULL || ap_tcur == NULL || ap_tolsf == NULL ||
                ap_tsw == NULL || ap_nst == NULL || ap_nfe == NULL ||
                ap_nje == NULL || ap_nqu == NULL || ap_mused == NULL) {
            goto fail;
        }
    }

    if (o_tcrit != NULL) {
        /* There are critical points */
        itask = 4;
        tcrit = (double *)(PyArray_DATA(ap_tcrit));
        rwork[0] = *tcrit;
    }

    /* Activate modern callback infrastructure for thread-safe operation */
    activate_odepack_callback(&callback);

    while (k < ntimes && istate > 0) {    /* loop over desired times */

        tout_ptr = tout + k;
        /* Use tcrit if relevant */
        if (itask == 4) {
            if (!tcrit) {
                PYERR(odepack_error, "Internal error - tcrit must be defined!");
            }
            if (*tout_ptr > *(tcrit + crit_ind)) {
                crit_ind++;
                rwork[0] = *(tcrit+crit_ind);
            }
        }
        if (crit_ind >= numcrit) {
            itask = 1;  /* No more critical values */
        }

        lsoda(ode_function_thunk, neq, y, &t, tout_ptr, itol, rtol, atol, &itask,
              &istate, &iopt, rwork, lrw, iwork, liw,
              ode_jacobian_thunk, (const int)jt);
        if (full_output) {
            *((double *)PyArray_DATA(ap_hu) + (k-1)) = rwork[10];
            *((double *)PyArray_DATA(ap_tcur) + (k-1)) = rwork[12];
            *((double *)PyArray_DATA(ap_tolsf) + (k-1)) = rwork[13];
            *((double *)PyArray_DATA(ap_tsw) + (k-1)) = rwork[14];
            *((int *)PyArray_DATA(ap_nst) + (k-1)) = iwork[10];
            *((int *)PyArray_DATA(ap_nfe) + (k-1)) = iwork[11];
            *((int *)PyArray_DATA(ap_nje) + (k-1)) = iwork[12];
            *((int *)PyArray_DATA(ap_nqu) + (k-1)) = iwork[13];
            if (istate == -5 || istate == -4) {
                imxer = iwork[15];
            }
            else {
                imxer = -1;
            }
            lenrw = iwork[16];
            leniw = iwork[17];
            *((int *)PyArray_DATA(ap_mused) + (k-1)) = iwork[18];
        }
        if (PyErr_Occurred()) {
            goto fail;
        }
        /* copy integration result to output - manual loop instead of memcpy */
        for (int i = 0; i < neq; i++) {
            yout_ptr[i] = y[i];
        }
        yout_ptr += neq;
        k++;
    }

    /* Deactivate callback infrastructure and clean up */
    deactivate_odepack_callback();
    cleanup_odepack_callback(&callback);

    Py_DECREF(extra_args);
    Py_DECREF(ap_atol);
    Py_DECREF(ap_rtol);
    Py_XDECREF(ap_tcrit);
    Py_DECREF(ap_y);
    Py_DECREF(ap_tout);
    free(wa);

    /* Do Full output */
    if (full_output) {
        return Py_BuildValue("N{s:N,s:N,s:N,s:N,s:N,s:N,s:N,s:N,s:l,s:l,s:l,s:N}l",
                    PyArray_Return(ap_yout),
                    "hu", PyArray_Return(ap_hu),
                    "tcur", PyArray_Return(ap_tcur),
                    "tolsf", PyArray_Return(ap_tolsf),
                    "tsw", PyArray_Return(ap_tsw),
                    "nst", PyArray_Return(ap_nst),
                    "nfe", PyArray_Return(ap_nfe),
                    "nje", PyArray_Return(ap_nje),
                    "nqu", PyArray_Return(ap_nqu),
                    "imxer", imxer,
                    "lenrw", lenrw,
                    "leniw", leniw,
                    "mused", PyArray_Return(ap_mused),
                    (long)istate);
    }
    else {
        return Py_BuildValue("Nl", PyArray_Return(ap_yout), (long)istate);
    }

fail:
    /* Clean up modern callback infrastructure */
    deactivate_odepack_callback();
    cleanup_odepack_callback(&callback);

    Py_XDECREF(extra_args);
    Py_XDECREF(ap_y);
    Py_XDECREF(ap_rtol);
    Py_XDECREF(ap_atol);
    Py_XDECREF(ap_tcrit);
    Py_XDECREF(ap_tout);
    Py_XDECREF(ap_yout);
    if (allocated) {
        free(wa);
    }
    if (full_output) {
        Py_XDECREF(ap_hu);
        Py_XDECREF(ap_tcur);
        Py_XDECREF(ap_tolsf);
        Py_XDECREF(ap_tsw);
        Py_XDECREF(ap_nst);
        Py_XDECREF(ap_nfe);
        Py_XDECREF(ap_nje);
        Py_XDECREF(ap_nqu);
        Py_XDECREF(ap_mused);
    }
    return NULL;
}


static struct PyMethodDef odepacklib_module_methods[] = {
    {"odeint", (PyCFunction) odepack_odeint, METH_VARARGS|METH_KEYWORDS, doc_odeint},
    {NULL, NULL, 0, NULL}
};


static int
odepacklib_module_exec(PyObject *module)
{
    if (_import_array() < 0) { return -1; }

    odepack_error = PyErr_NewException("_odepack.error", NULL, NULL);
    if (odepack_error == NULL) { return -1; }

    if (PyModule_AddObject(module, "error", odepack_error) < 0) {
        Py_DECREF(odepack_error);
        return -1;
    }

    return 0;
}


static struct PyModuleDef_Slot odepacklib_module_slots[] = {
    {Py_mod_exec, odepacklib_module_exec},
#if PY_VERSION_HEX >= 0x030c00f0  // Python 3.12+
    // signal that this module can be imported in isolated subinterpreters
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#endif
#if PY_VERSION_HEX >= 0x030d00f0  // Python 3.13+
    // signal that this module supports running without an active GIL
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL},
};


static struct
PyModuleDef moduledef = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "_odepack",
    .m_size = 0,
    .m_methods = odepacklib_module_methods,
    .m_slots = odepacklib_module_slots,
};


PyMODINIT_FUNC
PyInit__odepack(void)
{
    return PyModuleDef_Init(&moduledef);
}
