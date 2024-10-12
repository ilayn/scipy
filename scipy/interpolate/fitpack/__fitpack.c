void
fpbspl(const double* t, const int k, const double x, const int l, double* h)
{
    double f, hh[19];
    int i, j, li, lj;

    h[0] = 1.0;
    for (j = 0; j < k; j++)
    {
        for (i = 0; i <= j; i++) { hh[i] = h[i]; }
        h[0] = 0.0;
        for (i = 0; i <= j; i++)
        {
            li = l + i;
            lj = li - j;
            if (t[li] == t[lj])
            {
                h[i+1] = 0.0;
            } else {
                f = hh[i] / (t[li] - t[lj]);
                h[i] = h[i] + f*(t[li] - x);
                h[i+1] = f * (x - t[lj]);
            }
        }
    }
    return;
}