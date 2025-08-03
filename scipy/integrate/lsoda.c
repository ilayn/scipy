#include "lsoda.h"
#include <math.h>

/**
 * @brief Computes the weighted max-norm of a vector.
 *
 * This function computes the weighted max-norm of the vector of length n contained in the array v,
 * with weights contained in the array w of length n.
 *
 * The result is:
 *   vmnorm = max_{i=1,...,n} (abs(v[i]) * w[i])
 *
 * @param n Length of the vector.
 * @param v Input vector of length n.
 * @param w Weights vector of length n.
 * @return The weighted max-norm.
 */
static double
vmnorm(const int n, double* restrict v, double* restrict w)
{
    double vm = 0.0;
    for (int i = 0; i < n; i++)
    {
        vm = fmax(vm, fabs(v[i]) * w[i]);
    }
    return vm;
}



/**
 * @brief Computes and processes the matrix P = I - h*el[0]*J, where J is an approximation to the Jacobian.
 *
 * prja is called by stoda to compute and process the matrix P.
 * - J is computed by the user-supplied routine jac if miter == 0 or 3, or by finite differencing if miter == 1 or 4.
 * - J, scaled by -h*el[0], is stored in wm.
 * - The norm of J (matrix norm consistent with the weighted max-norm on vectors given by vmnorm) is computed, and J is overwritten by P.
 * - P is then subjected to LU decomposition in preparation for later solution of linear systems with P as coefficient matrix.
 *   This is done by dgetrf if miter == 0 or 1, and by dgbtrf if miter == 3 or 4.
 *
 * In addition to variables described previously, communication with prja uses the following:
 * @param y      Array containing predicted values on entry.
 * @param ftem   Work array of length n (acor in stoda).
 * @param savf   Array containing f evaluated at predicted y.
 * @param wm     Real work space for matrices. On output, it contains the LU decomposition of P.
 *               Storage of matrix elements starts at wm[2].
 *               wm also contains the following matrix-related data:
 *               - wm[0] = sqrt(uround), used in numerical Jacobian increments.
 * @param iwm    Integer work space containing pivot information, starting at iwm[20].
 *               iwm also contains the band parameters ml = iwm[0] and mu = iwm[1] if miter == 4 or 5.
 * @param el0    el[0] (input).
 * @param pdnorm Norm of Jacobian matrix (output).
 * @param ierpj  Output error flag: 0 if no trouble, >0 if P matrix found to be singular.
 * @param jcur   Output flag: 1 to indicate that the Jacobian matrix (or approximation) is now current.
 *
 * This routine also uses the common variables el0, h, tn, uround, miter, n, nfe, and nje.
 */
static void
prja(int* neq, double* y, double* yh, int nyh, double* ewt, double* ftem, double* savf, double* wm, int* iwm, f, jac)
{

}


/**
 * @brief Manages the solution of the linear system arising from a chord iteration.
 *
 * This routine is called if miter != 0.
 * - If miter is 0 or 1, it calls dgetrs to solve the system.
 * - If miter == 2, it updates the coefficient h*el0 in the diagonal matrix, and then computes the solution.
 * - If miter is 3 or 4, it calls dgbtrs.
 *
 * Communication with solsy uses the following variables:
 * @param wm   Real work space containing the inverse diagonal matrix if miter == 3,
 *             and the LU decomposition of the matrix otherwise.
 *             Storage of matrix elements starts at wm[2].
 *             wm also contains the following matrix-related data:
 *             - wm[0] = sqrt(uround) (not used here)
 *             - wm[1] = hl0, the previous value of h*el0, used if miter == 3
 * @param iwm  Integer work space containing pivot information, starting at iwm[20],
 *             if miter is 1, 2, 4, or 5. iwm also contains band parameters
 *             ml = iwm[0] and mu = iwm[1] if miter is 4 or 5.
 * @param x    The right-hand side vector on input, and the solution vector on output, of length n.
 * @param tem  Vector of work space of length n, not used in this version.
 * @param iersl Output flag (in common). iersl = 0 if no trouble occurred,
 *              iersl = 1 if a singular matrix arose with miter == 3.
 *
 * This routine also uses the common variables el0, h, miter, and n.
 */
static void
solsy(double* wm, int* iwm, double* x, double* tem)
{
    int int1 = 1, ierr = 0;
    iersl = 0;

    if ((miter == 0) || (miter == 1))
    {
        dgetrs_("N", &n, &int1, &wm[2], &n, &iwm[20], x, &n, &ierr);
        return;
    } else if (miter == 2) {
        phl0 = wm[1];
        hl0 = h * el0;
        if (hl0 != phl0)
        {
            r = hl0 / phl0;
            for (int i = 0; i < n; i++)
            {
                double di = 1.0 - r * (1.0 - 1.0 / wm[i + 2]);
                if (dabs(di) == 0.0)
                {
                    iersl = 1;
                    return;
                }
                wm[i + 2] = 1.0 / di;
            }
            // 320
        }
        for (int i = 0; i < n; i++)
        {
            x[i] = wm[i + 2] * x[i];
        }
        // 340
    } else if ((miter == 3) || (miter == 4)) {
        ml = iwm[0];
        mu = iwm[1];
        meband = 2*ml + mu + 1;
        dgbtrs_("N", &n, &ml, &mu, &int1, &wm[2], &meband, &iwm[20], x, &n, &ierr);
    }

    return;
}


