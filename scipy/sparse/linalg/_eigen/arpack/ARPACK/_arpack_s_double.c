#include "_arpack_s_double.h"

typedef int ARPACK_compare_rfunc(const double, const double);

static int sortr_LM(const double, const double);
static int sortr_SM(const double, const double);
static int sortr_LA(const double, const double);
static int sortr_SA(const double, const double);


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
        dcopy_(np, ritz, &int1, shifts, &int1);
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
