#include "_arpack_s_double.h"

typedef int ARPACK_compare_rfunc(const double, const double);

static int sortr_LM(const double, const double);
static int sortr_SM(const double, const double);
static int sortr_LA(const double, const double);
static int sortr_SA(const double, const double);

static const double unfl = 2.2250738585072014e-308;
static const double ovfl = 1.0 / 2.2250738585072014e-308;
static const double ulp = 2.220446049250313e-16;


void
dsaitr(struct ARPACK_dsaupd_variables *V, double* resid, double* rnorm,
       double* v, int ldv, double* h, int ldh, int* ipntr, double* workd)
{
    int i, infol, iter, ipj, irj, ivj, jj, n, tmp_int;
    double smlnum = unfl * ( V->n / ulp);
    double xtemp[2] = { 0.0 };
    const double sq2o2 = sqrt(2.0) / 2.0;

    char *MTYPE = "G", *TRANS = "T", *NORM = "1";
    int int1 = 1, int0 = 0;
    double dbl1 = 1.0, dbl0 = 0.0, temp1, tmp_dbl, tst1;

    n = V->n;  // n is constant, this is just for typing convenience
    ipj = 0;
    irj = ipj + n;
    ivj = irj + n;

    if (V->ido == ido_FIRST)
    {
         /*-----------------------------*
         | Initial call to this routine |
         *-----------------------------*/
        V->info = 0;
        V->saitr_step3 = 0;
        V->saitr_step4 = 0;
        V->saitr_orth1 = 0;
        V->saitr_orth2 = 0;
        V->saitr_restart = 0;
        V->saitr_j = V->nev;
    }

     /*------------------------------------------------*
     | When in reverse communication mode one of:      |
     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
     | will be .true. when ....                        |
     | STEP3: return from computing OP*v_{j}.          |
     | STEP4: return from computing B-norm of OP*v_{j} |
     | ORTH1: return from computing B-norm of r_{j+1}  |
     | ORTH2: return from computing B-norm of          |
     |        correction to the residual vector.       |
     | RSTART: return from OP computations needed by   |
     |         dgetv0.                                 |
     *------------------------------------------------*/
    if (V->saitr_step3) { goto LINE50; }
    if (V->saitr_step4) { goto LINE60; }
    if (V->saitr_orth1) { goto LINE70; }
    if (V->saitr_orth2) { goto LINE90; }
    if (V->saitr_restart) { goto LINE30; }

     /*----------------------------*
     | Else this is the first step |
     *----------------------------*/

     /*-------------------------------------------------------------*
     |                                                              |
     |        A R N O L D I     I T E R A T I O N     L O O P       |
     |                                                              |
     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
     *-------------------------------------------------------------*/

LINE1000:


         /*--------------------------------------------------------*
         | Check for exact zero. Equivalent to determining whether |
         | a j-step Arnoldi factorization is present.              |
         *--------------------------------------------------------*/

        if (*rnorm > 0.0) { goto LINE40; }

         /*--------------------------------------------------*
         | Invariant subspace found, generate a new starting |
         | vector which is orthogonal to the current Arnoldi |
         | basis and continue the iteration.                 |
         *--------------------------------------------------*/

         /*--------------------------------------------*
         | ITRY is the loop variable that controls the |
         | maximum amount of times that a restart is   |
         | attempted. NRSTRT is used by stat.h         |
         *--------------------------------------------*/


}


void
dsapps(int n, int* kev, int np, double* shift, double* v, int ldv, double* h, int ldh,
       double* resid, double* q, int ldq, double* workd)
{
    int i, iend, istart, itop, j, jj, kplusp, tmp_int, int1 = 1;
    double a1, a2, a3, a4, big, c, f, g, r, s, dbl0 = 0.0, dbl1 = 1.0, dblm1 = -1.0;

    itop = 0;
    kplusp = *kev + np;

     /*-------------------------------------------*
     | Initialize Q to the identity to accumulate |
     | the rotations and reflections              |
     *-------------------------------------------*/
    for (i = 0; i < kplusp; i++)
    {
        for (j = 0; j < kplusp; j++)
        {
            q[j + ldq*i] = 0.0;
            if (i == j) { q[j + ldq*i] = 1.0; }
        }
    }

     /*---------------------------------------------*
     | Quick return if there are no shifts to apply |
     *---------------------------------------------*/
    if (np == 0) { return; }

    for (jj = 0; jj < np; jj++)
    {
        istart = itop;
         /*---------------------------------------------------------*
         | Check for splitting and deflation. Currently we consider |
         | an off-diagonal element h(i+1,1) negligible if           |
         |         h(i+1,1) .le. epsmch*( |h(i,2)| + |h(i+1,2)| )   |
         | for i=1:KEV+NP-1.                                        |
         | If above condition tests true then we set h(i+1,1) = 0.  |
         | Note that h(1:KEV+NP,1) are assumed to be non negative.  |
         *---------------------------------------------------------*/

        do
        {
             /*-----------------------------------------------*
             | The following loop exits early if we encounter |
             | a negligible off diagonal element.             |
             *-----------------------------------------------*/

            for (i = istart; i < kplusp - 1; i++)
            {
                big = fabs(h[i + ldh]) + fabs(h[i + 1 + ldh]);
                if (h[i + 1] <= ulp*big)
                {
                    h[i + 1] = 0.0;
                    iend = i;
                    break;
                }
            }
            // No break
            if (i == kplusp - 1) { iend = i; }

            if (istart < iend)
            {
                 /*-------------------------------------------------------*
                 | Construct the plane rotation G'(istart,istart+1,theta) |
                 | that attempts to drive h(istart+1,1) to zero.          |
                 *-------------------------------------------------------*/
                f = h[istart + ldh] - shift[jj];
                g = h[istart + 1];
                dlartg_(&f, &g, &c, &s, &r);

                 /*------------------------------------------------------*
                 | Apply rotation to the left and right of H;            |
                 | H <- G' * H * G,  where G = G(istart,istart+1,theta). |
                 | This will create a "bulge".                           |
                 *------------------------------------------------------*/
                a1 = c*h[istart + ldh]     + s*h[istart+1];
                a2 = c*h[istart + 1]       + s*h[istart+1 + ldh];
                a4 = c*h[istart + 1 + ldh] - s*h[istart+1];
                a3 = c*h[istart + 1]       - s*h[istart + ldh];
                h[istart + ldh]     = c*a1 + s*a2;
                h[istart + 1 + ldh] = c*a4 - s*a3;
                h[istart + 1]       = c*a3 + s*a4;

                 /*---------------------------------------------------*
                 | Accumulate the rotation in the matrix Q;  Q <- Q*G |
                 *---------------------------------------------------*/
                tmp_int = (istart + jj > kplusp - 1 ? kplusp - 1 : istart + jj);
                for (j = 0; j < tmp_int; j++)
                {
                    a1                    =   c*q[j + istart*ldq] + s*q[j + (istart+1)*ldq];
                    q[j + (istart+1)*ldq] = - s*q[j + istart*ldq] + c*q[j + (istart+1)*ldq];
                    q[j + istart*ldq]     = a1;
                }

                 /*---------------------------------------------*
                 | The following loop chases the bulge created. |
                 | Note that the previous rotation may also be  |
                 | done within the following loop. But it is    |
                 | kept separate to make the distinction among  |
                 | the bulge chasing sweeps and the first plane |
                 | rotation designed to drive h(istart+1,1) to  |
                 | zero.                                        |
                 *---------------------------------------------*/

                for (i = istart + 1; i < iend - 1; i++)
                {
                     /*---------------------------------------------*
                     | Construct the plane rotation G'(i,i+1,theta) |
                     | that zeros the i-th bulge that was created   |
                     | by G(i-1,i,theta). g represents the bulge.   |
                     *---------------------------------------------*/
                    f = h[i];
                    g = s*h[i+1];

                     /*---------------------------------*
                     | Final update with G(i-1,i,theta) |
                     *---------------------------------*/
                    h[i + 1] = c*h[i + 1];
                    dlartg_(&f, &g, &c, &s, &r);

                     /*------------------------------------------*
                     | The following ensures that h(1:iend-1,1), |
                     | the first iend-2 off diagonal of elements |
                     | H, remain non negative.                   |
                     *------------------------------------------*/
                    if (r < 0)
                    {
                        r = -r;
                        c = -c;
                        s = -s;
                    }

                     /*-------------------------------------------*
                     | Apply rotation to the left and right of H; |
                     | H <- G * H * G',  where G = G(i,i+1,theta) |
                     *-------------------------------------------*/

                    h[i] = r;

                    a1 = c*h[i + ldh]     + s*h[i + 1];
                    a2 = c*h[i + 1]       + s*h[i + 1 + ldh];
                    a3 = c*h[i + 1]       - s*h[i + ldh];
                    a4 = c*h[i + 1 + ldh] - s*h[i + 1];

                    h[i + ldh]     = c*a1 + s*a2;
                    h[i + 1 + ldh] = c*a4 - s*a3;
                    h[i + 1]       = c*a3 + s*a4;

                     /*---------------------------------------------------*
                     | Accumulate the rotation in the matrix Q;  Q <- Q*G |
                     *---------------------------------------------------*/
                    tmp_int = (i + jj > kplusp - 1 ? kplusp - 1 : i + jj);
                    for (j = 0; j < tmp_int; j++)
                    {
                        a1               =   c*q[j + i*ldq] + s*q[j + (i+1)*ldq];
                        q[j + (i+1)*ldq] = - s*q[j + i*ldq] + c*q[j + (i+1)*ldq];
                        q[j + i*ldq  ]   = a1;
                    }
                    // 50
                }
                // 70
            }

             /*-------------------------*
             | Update the block pointer |
             *-------------------------*/

            istart = iend + 1;

             /*-----------------------------------------*
             | Make sure that h(iend,1) is non-negative |
             | If not then set h(iend,1) <-- -h(iend,1) |
             | and negate the last column of Q.         |
             | We have effectively carried out a        |
             | similarity on transformation H           |
             *-----------------------------------------*/

            if (h[iend] < 0.0)
            {
                h[iend] = -h[iend];
                dscal_(&kplusp, &dblm1, &q[ldq*iend], &int1);
            }

             /*-------------------------------------------------------*
             | Apply the same shift to the next block if there is any |
             *-------------------------------------------------------*/
        } while (iend < kplusp - 1);

         /*----------------------------------------------------*
         | Check if we can increase the the start of the block |
         *----------------------------------------------------*/
        for (i = itop; i < kplusp - 1; i++)
        {
            if (h[i+1] > 0.0) { break; }
            itop += 1;
        }
         /*----------------------------------*
         | Finished applying the jj-th shift |
         *----------------------------------*/
    }
    // 90

     /*-----------------------------------------*
     | All shifts have been applied. Check for  |
     | more possible deflation that might occur |
     | after the last shift is applied.         |
     *-----------------------------------------*/
    for (i = itop; i < kplusp - 1; i++)
    {
        big = fabs(h[i + ldh]) + fabs(h[i+1 + ldh]);
        if (h[i+1] <= ulp*big)
        {
            h[i+1] = 0.0;
        }
    }
    // 100

     /*------------------------------------------------*
     | Compute the (kev+1)-st column of (V*Q) and      |
     | temporarily store the result in WORKD(N+1:2*N). |
     | This is not necessary if h(kev+1,1) = 0.        |
     *------------------------------------------------*/
    if (h[*kev] > 0.0)
    {
        dgemv_("N", &n, &kplusp, &dbl1, v, ldv, &q[ldq*(*kev)], &int1, &dbl0, &workd[n], &int1);
    }

     /*------------------------------------------------------*
     | Compute column 1 to kev of (V*Q) in backward order    |
     | taking advantage that Q is an upper triangular matrix |
     | with lower bandwidth np.                              |
     | Place results in v(:,kplusp-kev:kplusp) temporarily.  |
     *------------------------------------------------------*/

    for (i = 0; i < *kev; i++)
    {
        tmp_int = kplusp - i;
        dgemv_("N", &n, &tmp_int, &dbl1, v, ldv, &q[ldq*(*kev-i)], &int1, &dbl0, workd, &int1);
        dcopy_(&n, workd, &int1, &v[ldv*(kplusp-i)], &int1);
    }
    // 130

     /*------------------------------------------------*
     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
     *------------------------------------------------*/
    for (i = 0; i < *kev; i++)
    {
        dcopy_(&n, &v[ldv*(np+i)], &int1, &v[ldv*i], &int1);
    }
    // 140

    if (h[*kev] > 0.0)
    {
        dcopy_(&n, &workd[n], &int1, &v[ldv*(*kev)], &int1);
    }

     /*------------------------------------*
     | Update the residual vector:         |
     |    r <- sigmak*r + betak*v(:,kev+1) |
     | where                               |
     |    sigmak = (e_{kev+p}'*Q)*e_{kev}  |
     |    betak = e_{kev+1}'*H*e_{kev}     |
     *------------------------------------*/

    dscal_(&n, &q[kplusp-1, *kev-1], resid, &int1);
    if (h[*kev] > 0.0)
    {
        daxpy_(&n, &h[*kev], &v[ldv*(*kev)], &int1, resid, &int1);
    }

    return;
}


void
dsgets(struct ARPACK_dsaupd_variables *V, int* kev, int* np, double* ritz, double* bounds, double* shifts)
{
    int kevd2, tmp1, tmp2, int1 = 1;
    if (V->which == which_BE)
    {
         /*----------------------------------------------------*
         | Both ends of the spectrum are requested.            |
         | Sort the eigenvalues into algebraically increasing  |
         | order first then swap high end of the spectrum next |
         | to low end in appropriate locations.                |
         | NOTE: when np < floor(kev/2) be careful not to swap |
         | overlapping locations.                              |
         *----------------------------------------------------*/
        dsortr(which_LA, 1, *kev + *np, ritz, bounds);
        kevd2 = *kev / 2;
        if (*kev > 1)
        {
            tmp1 = (kevd2 > *np ? *np : kevd2);
            tmp2 = (kevd2 > *np ? kevd2 : *np);
            dswap_(&tmp1, ritz, &int1, &ritz[tmp2], &int1);
            dswap_(&tmp1, bounds, &int1, &bounds[tmp2], &int1);
        }
    } else {
         /*---------------------------------------------------*
         | LM, SM, LA, SA case.                               |
         | Sort the eigenvalues of H into the desired order   |
         | and apply the resulting order to BOUNDS.           |
         | The eigenvalues are sorted so that the wanted part |
         | are always in the last KEV locations.              |
         *---------------------------------------------------*/
        dsortr(V->which, 1, *kev + *np, ritz, bounds);
    }

    if ((V->ishift == 1) && (*np > 0))
    {
         /*------------------------------------------------------*
         | Sort the unwanted Ritz values used as shifts so that  |
         | the ones with largest Ritz estimates are first.       |
         | This will tend to minimize the effects of the         |
         | forward instability of the iteration when the shifts  |
         | are applied in subroutine dsapps.                     |
         *------------------------------------------------------*/
        dsortr(which_SM, 1, *np, bounds, ritz);
        dcopy_(*np, ritz, &int1, shifts, &int1);
    }
}


void
dsortr(const enum ARPACK_which w, const int apply, const int n, double* x1, double* x2)
{
    int i, igap, j;
    double temp;
    ARPACK_compare_rfunc *f;

    switch (w)
    {
        case which_LM:
            f = sortr_LM;
            break;
        case which_SM:
            f = sortr_SM;
            break;
        case which_LA:
            f = sortr_LA;
            break;
        case which_SA:
            f = sortr_SA;
            break;
        default:
            f = sortr_LM;
            break;
    }

    igap = n / 2;

    while (igap != 0)
    {
        j = 0;
        for (i = igap; i < n; i++)
        {
            while (f(x1[j], x2[j]))
            {
                if (j < 0) { break; }
                temp = x1[j];
                x1[j] = x1[j+igap];
                x1[j+igap] = temp;

                if (apply)
                {
                    temp = x2[j];
                    x2[j] = x2[j+igap];
                    x2[j+igap] = temp;
                }
                j -= igap;
            }
            j = i - igap + 1;
        }
        igap = igap / 2;
    }
}


int
sortr_LM(const double x1, const double x2)
{
    return (fabs(x1) > fabs(x2));
}


int
sortr_SM(const double x1, const double x2)
{
    return (fabs(x1) < fabs(x2));
}


int
sortr_LA(const double x1, const double x2)
{
    return (x1 > x2);
}


int
sortr_SA(const double x1, const double x2)
{
    return (x1 < x2);
}
