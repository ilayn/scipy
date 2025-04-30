#include <math.h>




void fpbspl(double* t, int n, int k, double x, int l, double* h)
{
    // subroutine fpbspl evaluates the (k+1) non-zero b-splines of degree k
    // at t(l) <= x < t(l+1) using the stable recurrence relation of de boor and cox.
    // Travis Oliphant  2007
    //    changed so that weighting of 0 is used when knots with multiplicity are present.
    //    Also, notice that l+k <= n and 1 <= l+1-k or else the routine will be
    //    accessing memory outside t.  Thus it is imperative that k <= l <= n-k
    //    but this is not checked.
    int li, lj;
    double f;
    double hh[19] = { 0.0 };

    h[0] = 1.0;
    for (int j = 0; j < k; j++) {
        for (int i = 0; i <= j; i++)
        {
            hh[i] = h[i];
        }
        // 10
        h[0] = 0.0;
        for (int i = 0; i <= j; i++)
        {
            li = l + i + 1;
            lj  = li - j - 1;
            if (t[li] != t[lj])
            {
                f = hh[i] / (t[li] - t[lj]);
                h[i] = h[i] + f * (t[li] - x);
                h[i + 1] = f * (x - t[lj]);
            }
            else {
                h[i + 1] = 0.0;
            }
        }
    }
    // 20
}


void fpdeno(int maxtr, int* up, int* left, int* right, int nbind, int* merk)
{

}




void fpdisc(double* t, int n, int k2, double* b, int nest)
{

}


void fpgivs(double piv, double* ww, double* c, double* s)
{
    // fpgivs calculates the parameters of a givens transformation.
    // Copied from dlartg function of LAPACK.
    const double safmin = 2.2250738585072014e-308;
    const double safmax = 1.7976931348623157e+308;
    const double rtmin = sqrt(safmin);
    const double rtmax = sqrt(safmax/2);
    double f1 = fabs(piv);
    double g1 = fabs(*ww);

    if(*ww == 0.0)
    {
       *c = 1.0;
       *s = 0.0;
       *ww = piv;
    }
    else if (piv == 0.0)
    {
       *c = 0.0;
       *s = copysign(1.0, *ww);
       *ww = g1;
    }
    else if (((f1 > rtmin) && (f1 < rtmax)) && ((g1 > rtmin) && (g1 < rtmax )))
    {
       double d = hypot(piv, *ww);
       double orig_ww = *ww;
       *c = f1 / d;
       *ww = copysign(d, piv);
       *s = orig_ww / *ww;
    }
    else
    {
       double u = fmin(safmax, fmax(safmin, fmax(f1, g1)));
       double fs = piv / u;
       double gs = *ww / u;
       double d = hypot(fs, gs);
       *c = fabs(fs) / d;
       *ww = copysign(d, piv);
       *s = gs / *ww;
       *ww = (*ww)*u;
    }
    return;
}


void fporde(double* x, double* y, const int m, const int kx, const int ky,
            double* tx, const int nx, double* ty, const int ny, int* nummer,
            int* indx, const int nreg)
{
    // subroutine fporde sorts the data points (x(i),y(i)),i=1,2,...,m
    // according to the panel tx(l)<=x<tx(l+1),ty(k)<=y<ty(k+1), they belong
    // to. for each panel a stack is constructed  containing the numbers
    // of data points lying inside; index(j),j=1,2,...,nreg points to the
    // first data point in the jth panel while nummer(i),i=1,2,...,m gives
    // the number of the next data point in the panel.

}


void fprota(double c, double s, double *a, double *b)
{
    // fprota applies a givens rotation to the pair (a,b).
    double temp = c * (*a) - s * (*b);
    *b = c * (*b) + s * (*a);
    *a = temp;
    return;
}
