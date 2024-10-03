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
dsaitr(struct ARPACK_arnoldi_update_vars_d *V, double* resid, double* rnorm,
       double* v, int ldv, double* h, int ldh, int* ipntr, double* workd)
{
    int i, infol, ipj, irj, ivj, jj, n, tmp_int;
    double smlnum = unfl * ( V->n / ulp);
    double xtemp[2] = { 0.0 };
    const double sq2o2 = sqrt(2.0) / 2.0;

    char *MTYPE = "G", *TRANS = "T", *NORM = "1";
    int int1 = 1, int0 = 0;
    double dbl1 = 1.0, dbl0 = 0.0, dblm1 = -1.0, temp1, tmp_dbl, tst1;

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
        V->aitr_step3 = 0;
        V->aitr_step4 = 0;
        V->aitr_orth1 = 0;
        V->aitr_orth2 = 0;
        V->aitr_restart = 0;
        V->aitr_j = V->nev;
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
    if (V->aitr_step3) { goto LINE50; }
    if (V->aitr_step4) { goto LINE60; }
    if (V->aitr_orth1) { goto LINE70; }
    if (V->aitr_orth2) { goto LINE90; }
    if (V->aitr_restart) { goto LINE30; }

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

    V->getv0_itry = 1;
LINE20:
    V->aitr_restart = 1;
    V->ido = ido_FIRST;
LINE30:
    dgetv0(V, 0, n, V->aitr_j, v, ldv, resid, rnorm, ipntr, workd, &V->aitr_ierr);
    if (V->ido != ido_DONE) { return; }
    if (V->aitr_ierr < 0)
    {
        V->getv0_itry += 1;
        if (V->getv0_itry <= 3) { goto LINE20; }
         /*-----------------------------------------------*
         | Give up after several restart attempts.        |
         | Set INFO to the size of the invariant subspace |
         | which spans OP and exit.                       |
         *-----------------------------------------------*/
        V->info = V->aitr_j;
        V->ido = ido_DONE;
        return;
    }

LINE40:
     /*--------------------------------------------------------*
     | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
     | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
     | when reciprocating a small RNORM, test against lower    |
     | machine bound.                                          |
     *--------------------------------------------------------*/
    dcopy_(&n, resid, &int1, &v[ldv*(V->aitr_j)], &int1);
    if (*rnorm >= unfl)
    {
        temp1 = 1.0 / *rnorm;
        dscal_(&n, &temp1, &v[ldv*(V->aitr_j)], &int1);
        dscal_(&n, &temp1, &workd[ipj], &int1);
    } else {
        dlascl_(MTYPE, &i, &i, rnorm, &dbl1, &n, &int1, &v[ldv*(V->aitr_j)], &n, &infol);
        dlascl_(MTYPE, &i, &i, rnorm, &dbl1, &n, &int1, &workd[ipj], &n, &infol);
    }

     /*-----------------------------------------------------*
     | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
     | Note that this is not quite yet r_{j}. See STEP 4    |
     *-----------------------------------------------------*/
    V->aitr_step3 = 1;
    dcopy_(&n, &v[ldv*(V->aitr_j)], &int1, &workd[ivj], &int1);
    ipntr[0] = ivj;
    ipntr[1] = irj;
    ipntr[2] = ipj;
    V->ido = ido_OPX;

     /*----------------------------------*
     | Exit in order to compute OP*v_{j} |
     *----------------------------------*/
    return;

LINE50:
     /*---------------------------------*
     | Back from reverse communication; |
     | WORKD(IRJ:IRJ+N-1) := OP*v_{j}   |
     *---------------------------------*/
    V->aitr_step3 = 0;

     /*-----------------------------------------*
     | Put another copy of OP*v_{j} into RESID. |
     *-----------------------------------------*/
    dcopy_(&n, &workd[irj], &int1, resid, &int1);

     /*------------------------------------------*
     | STEP 4:  Finish extending the symmetric   |
     |          Arnoldi to length j. If MODE = 2 |
     |          then B*OP = B*inv(B)*A = A and   |
     |          we don't need to compute B*OP.   |
     | NOTE: If MODE = 2 WORKD(IVJ:IVJ+N-1) is   |
     | assumed to have A*v_{j}.                  |
     *------------------------------------------*/
    if (V->mode == 2) { goto LINE65; }
    if (V->bmat)
    {
        ipntr[0] = irj;
        ipntr[1] = ipj;
        V->ido = ido_BX;
         /*------------------------------------*
         | Exit in order to compute B*OP*v_{j} |
         *------------------------------------*/
        return;
    } else {
        dcopy_(&n, resid, &int1, &workd[ipj], &int1);
    }

LINE60:

     /*---------------------------------*
     | Back from reverse communication; |
     | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j} |
     | if step4 = .true.                |
     *---------------------------------*/
    V->aitr_step4 = 0;

     /*------------------------------------*
     | The following is needed for STEP 5. |
     | Compute the B-norm of OP*v_{j}.     |
     *------------------------------------*/

LINE65:

    if (V->mode == 2)
    {
         /*---------------------------------*
         | Note that the B-norm of OP*v_{j} |
         | is the inv(B)-norm of A*v_{j}.   |
         *---------------------------------*/
        V->aitr_wnorm = ddot_(&n, resid, &int1, &workd[ivj], &int1);
        V->aitr_wnorm = sqrt(fabs(V->aitr_wnorm));
    } else if (V->bmat) {
        V->aitr_wnorm = ddot_(&n, resid, &int1, &workd[ipj], &int1);
        V->aitr_wnorm = sqrt(fabs(V->aitr_wnorm));
    } else {
        V->aitr_wnorm = dnrm2_(&n, resid, &int1);
    }

     /*----------------------------------------*
     | Compute the j-th residual corresponding |
     | to the j step factorization.            |
     | Use Classical Gram Schmidt and compute: |
     | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
     | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
     *----------------------------------------*/

     /*-----------------------------------------*
     | Compute the j Fourier coefficients w_{j} |
     | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
     *-----------------------------------------*/

    if (V->mode != 2)
    {
        dgemv_("T", &n, &V->aitr_j, &dbl1, v, &ldv, &workd[ipj], &int1, &dbl0, &workd[irj], &int1);
    } else {
        dgemv_("T", &n, &V->aitr_j, &dbl1, v, &ldv, &workd[ivj], &int1, &dbl0, &workd[irj], &int1);
    }

     /*-------------------------------------*
     | Orthgonalize r_{j} against V_{j}.    |
     | RESID contains OP*v_{j}. See STEP 3. |
     *-------------------------------------*/
    dgemv("N", &n, &V->aitr_j, &dblm1, v, &ldv, &workd[irj], &int1, &dbl0, resid, &int1);

     /*-------------------------------------*
     | Extend H to have j rows and columns. |
     *-------------------------------------*/
    h[V->aitr_j + ldh] = workd[irj + V->aitr_j + 1];

    if ((V->aitr_j == 0) || (V->aitr_restart))
    {
        h[V->aitr_j] = 0.0;
    } else {
        h[V->aitr_j] = *rnorm;
    }

    V->aitr_orth1 = 1;
    V->aitr_iter = 0;

    if (V->bmat)
    {
        dcopy_(&n, resid, &int1, &workd[irj], &int1);
        ipntr[0] = irj;
        ipntr[1] = ipj;
        V->ido = ido_BX;
         /*---------------------------------*
         | Exit in order to compute B*r_{j} |
         *---------------------------------*/
        return;
    } else {
        dcopy_(&n, resid, &int1, &workd[ipj], &int1);
    }

LINE70:

     /*--------------------------------------------------*
     | Back from reverse communication if ORTH1 = .true. |
     | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
     *--------------------------------------------------*/

    V->aitr_orth1 = 0;

     /*-----------------------------*
     | Compute the B-norm of r_{j}. |
     *-----------------------------*/
    if (V->bmat)
    {
        *rnorm = ddot_(&n, resid, &int1, &workd[ipj], &int1);
        *rnorm = sqrt(fabs(*rnorm));
    } else {
        *rnorm = dnrm2_(&n, resid, &int1);
    }

     /*----------------------------------------------------------*
     | STEP 5: Re-orthogonalization / Iterative refinement phase |
     | Maximum NITER_ITREF tries.                                |
     |                                                           |
     |          s      = V_{j}^T * B * r_{j}                     |
     |          r_{j}  = r_{j} - V_{j}*s                         |
     |          alphaj = alphaj + s_{j}                          |
     |                                                           |
     | The stopping criteria used for iterative refinement is    |
     | discussed in Parlett's book SEP, page 107 and in Gragg &  |
     | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
     | Determine if we need to correct the residual. The goal is |
     | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
     *----------------------------------------------------------*/

    if (*rnorm > sq2o2) { goto LINE100; }
    V->aitr_iter = 0;

     /*--------------------------------------------------*
     | Enter the Iterative refinement phase. If further  |
     | refinement is necessary, loop back here. The loop |
     | variable is ITER. Perform a step of Classical     |
     | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
     *--------------------------------------------------*/
LINE80:

     /*---------------------------------------------------*
     | Compute V_{j}^T * B * r_{j}.                       |
     | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
     *---------------------------------------------------*/

    dgemv_("T", &n, &V->aitr_j, &dbl1, v, &ldv, &workd[ipj], &int1, &dbl0, &workd[irj], &int1);

     /*--------------------------------------------*
     | Compute the correction to the residual:     |
     | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1). |
     | The correction to H is v(:,1:J)*H(1:J,1:J)  |
     | + v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j.         |
     *--------------------------------------------*/

    dgemv_("N", &n, &V->aitr_j, &dblm1, v, &ldv, &workd[irj], &int1, &dbl1, resid, &int1);

    if ((V->aitr_j == 0) || (V->aitr_restart))
    {
        h[V->aitr_j] = 0.0;
    }
    h[V->aitr_j + ldh] = workd[irj + V->aitr_j + 1];
    V->aitr_orth2 = 1;

    if (V->bmat)
    {
        dcopy_(&n, resid, &int1, &workd[irj], &int1);
        ipntr[0] = irj;
        ipntr[1] = ipj;
        V->ido = ido_BX;
         /*----------------------------------*
         | Exit in order to compute B*r_{j}. |
         | r_{j} is the corrected residual.  |
         *----------------------------------*/
        return;
    } else {
        dcopy_(&n, resid, &int1, &workd[ipj], &int1);
    }

LINE90:
     /*--------------------------------------------------*
     | Back from reverse communication if ORTH2 = .true. |
     *--------------------------------------------------*/

     /*----------------------------------------------------*
     | Compute the B-norm of the corrected residual r_{j}. |
     *----------------------------------------------------*/
    if (V->bmat)
    {
        V->aitr_rnorm1 = ddot_(&n, resid, &int1, &workd[ipj], &int1);
        V->aitr_rnorm1 = sqrt(fabs(V->aitr_rnorm1));
    } else {
        V->aitr_rnorm1 = dnrm2_(&n, resid, &int1);
    }

     /*----------------------------------------*
     | Determine if we need to perform another |
     | step of re-orthogonalization.           |
     *----------------------------------------*/
    if (V->aitr_rnorm1 > sq2o2)
    {
         /*--------------------------------------*
         | No need for further refinement.       |
         *--------------------------------------*/
        *rnorm = V->aitr_rnorm1;

    } else {
         /*------------------------------------------*
         | Another step of iterative refinement step |
         | is required.                              |
         *------------------------------------------*/
        *rnorm = V->aitr_rnorm1;
        V->aitr_iter += 1;
        if (V->aitr_iter < 2) { goto LINE80; }

         /*------------------------------------------------*
         | Otherwise RESID is numerically in the span of V |
         *------------------------------------------------*/
        for (jj = 0; jj < n; jj++)
        {
            resid[jj] = 0.0;
        }
        // 95
        *rnorm = 0.0;
    }

     /*---------------------------------------------*
     | Branch here directly if iterative refinement |
     | wasn't necessary or after at most NITER_REF  |
     | steps of iterative refinement.               |
     *---------------------------------------------*/
LINE100:

    V->aitr_restart = 0;
    V->aitr_orth2 = 0;

     /*---------------------------------------------------------*
     | Make sure the last off-diagonal element is non negative  |
     | If not perform a similarity transformation on H(1:j,1:j) |
     | and scale v(:,j) by -1.                                  |
     *---------------------------------------------------------*/
    if (h[V->aitr_j] < 0.0)
    {
        h[V->aitr_j] = -h[V->aitr_j];
        if (V->aitr_j < V->nev + V->np)
        {
            dscal_(&n, &dblm1, &v[V->aitr_j + 1], &int1);
        } else {
            dscal_(&n, &dblm1, resid, &int1);
        }
    }

     /*-----------------------------------*
     | STEP 6: Update  j = j+1;  Continue |
     *-----------------------------------*/
    V->aitr_j += 1;
    if (V->aitr_j > V->nev + V->np)
    {
        V->ido = ido_DONE;
        return;
    }

     /*-------------------------------------------------------*
     | Loop back to extend the factorization by another step. |
     *-------------------------------------------------------*/

    goto LINE1000;

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
dsgets(struct ARPACK_arnoldi_update_vars_d *V, int* kev, int* np, double* ritz, double* bounds, double* shifts)
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
