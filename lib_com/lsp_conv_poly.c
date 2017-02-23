/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*-------------------------------------------------------------------*
 * Local constants
 *-------------------------------------------------------------------*/

/* The conversion modes. */
#define DOWNCONV   0
#define UPCONV     1
/* The cap of the inverse power spectrum. */
#define MAXPOWERSPECT   1e-5f
#define N50             GRID50_POINTS
#define N40             GRID40_POINTS


/*-------------------------------------------------------------------*
 * Local functions
 *-------------------------------------------------------------------*/

static void powerspect( const float x[], const short N, float  R[], float  S[], float G[], const short mode );

static void spectautocorr( const float x[], const short N, const float G[], float  r[] );

static void zeros2poly( float x[], float R[], float S[] );

static void polydecomp( float A[], float P[], float Q[] );

static void cheb2poly( float P[] );


/*---------------------------------------------------------------------*
 *  lsp_convert_poly()
 *
 *  Converts the LP filter estimated at 16.0 kHz sampling rate down
 *  12.8 kHz frequency scale or alternatively from 12.8 kHz up to
 *  16.0 kHz.  The former is called down conversation and latter up
 *  conversion.  The resulting LP filter is characterized with its
 *  line spectrum pairs.  The original Lp filter can be either in
 *  its immittance, used for the AMR-WB IO mode, or line spectrum
 *  pair representation.
 *
 *  The conversion is based the autocorrelation computed from the
 *  power spectrum of the LP filter that is truncated or extrapolated
 *  to the desired frequency scale.
 *---------------------------------------------------------------------*/

short lsp_convert_poly(
    float  w[],           /* i/o: LSP or ISP parameters          */
    const short L_frame,        /* i  : flag for up or down conversion */
    const short Opt_AMRWB       /* i  : flag for the AMR-WB IO mode    */
)
{
    float epsP[M+1];
    float G[GRID50_POINTS], r[M+1], A[M+1], A_old[M+1], R[NC+1], S[NC+1];
    short i;
    short flag;

    /*---------------------------------------------------------------------*
     * Because AMR-WB IO uses immittance spectrum frequency representation
     * instead of line spectrum frequency representation, the input
     * parameters do not give the zeros of the polynomials R(x) and S(x).
     * Hence R(x) and S(x) are formed via the polynomial A(z) of the linear
     * prediction filter.
     *---------------------------------------------------------------------*/

    if( Opt_AMRWB )
    {
        isp2a( w, A_old, M );
        polydecomp( A_old, R, S );
    }

    /*---------------------------------------------------------------------*
     * Form the polynomials R(x) and S(x) from their zeros that are the
     * line spectrum pairs of A(z).  The polynomial coefficients can be
     * scaled for convenience, because scaling will not affect the
     * resulting LP coefficients.  Scaling by 128 gives the correct offset
     * to the power spectrum for n = 16.
     *---------------------------------------------------------------------*/

    else
    {
        zeros2poly( w, R, S );
        for (i = 0; i <= NC; i++)
        {
            R[i] *= 128.0f;
            S[i] *= 128.0f;
        }
        lsp2a_stab( w, A_old, M );
    }

    /*---------------------------------------------------------------------*
     * Conversion from 16.0 kHz down to 12.8 kHz.  The power spectrum
     * needs to be computed only up to 6.4 kHz, because the upper band
     * is omitted.
     *---------------------------------------------------------------------*/

    if( L_frame == L_FRAME )
    {
        powerspect( grid50, N50, R, S, G, DOWNCONV );
        spectautocorr( grid40, N40, G, r );
    }

    /*---------------------------------------------------------------------*
     * Conversion from 12.8 kHz up to 16.0 kHz.
     * Compute the power spectrum of the LP filter, extrapolate the
     * power spectrum from 6.4 kHz to 8.0 kHz, and compute auto-
     * correlation on this power spectrum.
     *---------------------------------------------------------------------*/

    else
    {
        powerspect( grid40, N40, R, S, G, UPCONV );
        for (i = N40; i < N50; i++)
        {
            G[i] = G[N40-1];
        }

        spectautocorr( grid50, N50, G, r );
    }

    /*---------------------------------------------------------------------*
     * Compute the linear prediction coefficients from the autocorrelation
     * and convert to line spectrum pairs.
     *---------------------------------------------------------------------*/

    flag = lev_dur( A, r, M, epsP );
    a2lsp_stab( A, w, stable_LSP );

    return(flag);
}


/*---------------------------------------------------------------------*
 *  powerspect()
 *
 *  Computes the power spectrum G(w) = 1/|A(w)|^2 at N points on
 *  the real axis x = cos w by utilizing the line spectrum frequency
 *  decomposition
 *
 *     A(z) = (P(z) + Q(z))/2,
 *
 *  where assuming A(z) of an even degree n,
 *
 *     P(z) = [A(z) + z^(n+1) A(1/z)]/(1/z + 1),
 *     Q(z) = [A(z) - z^(n+1) A(1/z)]/(1/z - 1).
 *
 *  The zeros of these polynomials give the line spectrum frequencies
 *  of A(z).  It can be shown that for an even n,
 *
 *     |A(x)|^2 = 2 (1 + x) R(x)^2 + 2 (1 - x) S(x)^2,
 *
 *  where x = cos w, and R(x) and S(x) are the direct polynomials
 *  resulting from the Chebyshev series representation of P(z)
 *  and Q(z).
 *
 *  This routine assumes the grid X = 1, x[0], x[1], .., x[m-1],
 *  -, ..., -x[1], -x[0], -1 such that x[i] = cos((i+1)*pi/N) for
 *  evaluating the power spectrum.  Only m = (N-1)/2 - 1 grid points
 *  need to be stored, because cos(0) and cos(pi/2) are trivial,
 *  and the points above pi/2 are obtained readily using the symmetry
 *  of cosine.
 *
 *  The power spectrum can be scaled as a*G[], where a is chosen
 *  for convenience. This is because the scaling has no impact on
 *  the LP coefficients to be determined based on the power spectrum.
 *---------------------------------------------------------------------*/

static void powerspect(
    const float x[],        /* i: Grid points x[0:N-1]             */
    const short N,          /* i: Number of grid points            */
    float  R[],             /* i: Coefficients of R(x) in R[0:NC]  */
    float  S[],             /* i: Coefficients of S(x) in S[0:NC]  */
    float  G[],             /* o: Power spectrum G[0:N]            */
    const short mode        /* i: Flag for up or down conversion   */
)
{
    float  re, ro, se, so;
    float  s0, s1, r0, r1, x0, x1, x2;
    Word16 i, j;
    Word16 iuni, imid;

    /* init buffer */
    for(i=0; i<N; i++)
    {
        G[i] = 0;
    }
    /*---------------------------------------------------------------------*
     * Down conversion yields iuni unique grid points that do not have
     * symmetric counterparts above x = cos(pi/2) = 0.
     * Set the mid point of the frequency grid.
     *---------------------------------------------------------------------*/

    if (mode == DOWNCONV)
    {
        iuni = (GRID50_POINTS - 1)/5 - 1;
        imid = (GRID50_POINTS - 1)/2;
    }

    /*---------------------------------------------------------------------*
    * Power spectrum at x = cos(pi) = -1 that is not needed in down
    * conversion. Set the mid point of the frequency grid.
    *---------------------------------------------------------------------*/

    else
    {
        iuni = 0;
        imid = (GRID40_POINTS - 1)/2;

        s0 = S[0];
        for (j = 1; j <= NC; j++)
        {
            s0 = S[j] - s0;
        }
        G[N-1] = 4.0f*s0*s0;
    }

    /*---------------------------------------------------------------------*
     * Power spectrum at x = cos(0) = 1.
     *---------------------------------------------------------------------*/

    r0 = R[0];
    for (j = 1; j <= NC; j++)
    {
        r0 += R[j];
    }
    r0 = max(r0, 0.000000953674316f);
    G[0] = 4.0f*r0*r0;

    /*---------------------------------------------------------------------*
     * Power spectrum at x = cos(pi/2) = 0.
     *---------------------------------------------------------------------*/

    r0 = R[NC];
    s0 = S[NC];
    r0 = r0*r0;
    s0 = s0*s0;

    G[imid] = 2.0f*(r0 + s0);

    /*---------------------------------------------------------------------*
     * Power spectrum at unique points that do not have symmetric
     * counterparts at x > cos(pi/2) = 0.
     *---------------------------------------------------------------------*/

    for (i = 1; i <= iuni; i++)
    {
        x0 = x[i-1];
        r0 = R[0];
        s0 = S[0];

        for (j = 1; j <= NC; j++)
        {
            r0 = x0*r0 + R[j];
            s0 = x0*s0 + S[j];
        }
        r0 = (1.0f + x0)*r0*r0;
        s0 = (1.0f - x0)*s0*s0;

        G[i] = 2.0f*(r0 + s0);
    }

    /*---------------------------------------------------------------------*
     * Power spectrum at points other than x = -1, 0, and 1 and unique
     * points is computed using the anti-symmetry of the grid relative
     * to the midpoint x = 0 in order to reduce looping.
     *---------------------------------------------------------------------*/

    for ( ; i < imid; i++)
    {
        x0 = x[i-1];
        x2 = x0*x0;

        re = R[0];
        se = S[0];
        ro = R[1];
        so = S[1];

        for (j = 2; j < NC; j+=2)
        {
            re = x2*re + R[j];
            ro = x2*ro + R[j+1];
            se = x2*se + S[j];
            so = x2*so + S[j+1];
        }

        re  = x2*re + R[j];
        ro *= x0;
        se  = x2*se + S[j];
        so *= x0;

        r0 = re + ro;
        s0 = se + so;
        r1 = re - ro;
        s1 = se - so;

        x1 = 1.0f + x0;
        x2 = 1.0f - x0;

        r0 = x1*r0*r0;
        s0 = x2*s0*s0;
        r1 = x2*r1*r1;
        s1 = x1*s1*s1;

        G[i]     = 2.0f*(r0 + s0);
        G[N-i-1] = 2.0f*(r1 + s1);
    }

    /*---------------------------------------------------------------------*
     * Compute the power spectrum 1/|A(x)|^2 from |A(x)|^2 with logic
     * to prevent division by small number and upper bound the spectrum.
     * This upper bound is implicit in fixed point arithmetic, but is
     * needed explicitly in floating point implementations to avoid numeric
     * problems.
     *---------------------------------------------------------------------*/

    for(i=0; i<N; i++)
    {
        if (G[i] < MAXPOWERSPECT)
        {
            G[i] = MAXPOWERSPECT;
        }
        G[i] = 1.0f/G[i];
    }
    return;
}

/*---------------------------------------------------------------------*
 *  spectautocorr()
 *
 *  Computes the autocorrelation r[j] for j = 0, 1, ..., M from
 *  the power spectrum P(w) by using the rectangle rule to
 *  approximate the integral
 *
 *             1     pi
 *     r[j] = ---    I  P(w) cos(j*w) dw.
 *            2*pi -pi
 *
 *  It is sufficient to evaluate the integrand only from w = 0 to
 *  w = pi due to the symmetry P(-w) = P(w).  We can further
 *  employ the relation
 *
 *      cos(j*(pi - w)) = (-1)^j cos(j*w)
 *
 *  to use symmetries relative to w = pi/2 for w in (0, pi/2).
 *
 *  When applying the rectangle rule, it is convenient to evaluate
 *  w = 0, w = pi/2, and w = pi separately.  By using a frequency
 *  grid of N points, we can express the rectangle rule as
 *
 *     r[j] = G[0] + 2*a*G[(N-1)/2] + b*G[N-1]
 *
 *                      M
 *                 + 2 sum (G[i] - G[N-i-1]) cos(j*x[i])
 *                     i=1
 *
 *  where G[i] is the power spectrum at x[i] = i*pi/(N-1) and
 *  M = (N-1)/2 - 1 is the number of the grid points in the
 *  interval(0, pi/2).  The number of grid points N is assumed odd.
 *
 *  The coefficients
 *
 *     b = (-1)^j
 *     a = (1 + (-1)^(j+1))(-1)^floor(j/2)
 *
 *  follow from the properties of the cosine function.  The
 *  computation further uses the recursion
 *
 *     cos(j*w) = 2*cos(w)*cos((j-1)*w) - cos((j-2)*w)
 *
 *  Note that the autocorrelation can be scaled for convenience,
 *  because this scaling has no impact on the LP coefficients to be
 *  calculated from the autocorrelation. The expression of r[j] thus
 *  omits the division by N.
 *
 *  See the powerspect function on the definition of the grid.
 *
 *  References
 *  J. Makhoul, "Spectral linear prediction: properties and
 *  applications," IEEE Trans. on Acoustics, Speech and Signal
 *  Processing, Vol. 23, No. 3, pp.283-296, June 1975
 *---------------------------------------------------------------------*/

static void spectautocorr(
    const float x[],        /* i: Grid points x[0:N-1]             */
    const short N,          /* i: Number of grid points            */
    const float G[],        /* i: Power spectrum G[0:N-1]          */
    float r[]         /* o: Autocorrelation r[0:M]           */
)
{
    float  c[M+1];          /* c[j] = cos(j*w) */
    float  gp, gn, c2;
    Word16 i, j;
    Word16 imid;

    /*---------------------------------------------------------------------*
     * The midpoint of the cosine table x of N entries. Only the entries
     * x[0] = cos(d), x[1] = cos(2*d), ..., x[imid-1] = cos((imid-1)*d)
     * need to be stored due to trivial cos(0), cos(pi/2), cos(pi), and
     * symmetry relative to pi/2. Here d = pi/(N - 1).
     *---------------------------------------------------------------------*/

    imid = (N - 1)/2;

    /*---------------------------------------------------------------------*
     * Autocorrelation r[j] at zero lag j = 0 for the upper half of the
     * unit circle, but excluding the points x = cos(0) and x = cos(pi).
     *---------------------------------------------------------------------*/

    r[0] = G[1];
    for (i = 2; i < N-1; i++)
    {
        r[0] += G[i];
    }

    /*---------------------------------------------------------------------*
     * Initialize the autocorrelation r[j] at lags greater than zero
     * by adding the midpoint x = cos(pi/2) = 0.
     *---------------------------------------------------------------------*/

    r[1] = 0.0f;
    r[2] = -G[imid];

    for (j = 3; j < M; j+=2)
    {
        r[j]   =  0.0f;
        r[j+1] = -r[j-1];
    }

    /*---------------------------------------------------------------------*
     * Autocorrelation r[j] at lags j = 1, 2, ..., M.  The computation
     * employes the relation cos(j*(pi - w)) = (-1)^j cos(j*w) and
     * cos(j*w) = 2*cos(w)*cos((j-1)*w) - cos((j-2)*w) for obtaining
     * the cosine c[j] = cos(j*w).
     *---------------------------------------------------------------------*/

    c[0] = 1.0f;

    for (i = 1; i < imid; i++)
    {
        gp = G[i] + G[N-i-1];
        gn = G[i] - G[N-i-1];

        c[1] = x[i-1];
        c2 = 2.0f*c[1];

        r[1] += gn*c[1];

        for (j = 2; j < M; j+=2)
        {
            c[j]    = c2*c[j-1] - c[j-2];
            r[j]   += gp*c[j];
            c[j+1]  = c2*c[j]   - c[j-1];
            r[j+1] += gn*c[j+1];
        }

        c[j] = c2*c[j-1] - c[j-2];
        r[j] += gp*c[j];
    }

    /*---------------------------------------------------------------------*
     * Add the endpoints x = cos(0) = 1 and x = cos(pi) = -1 as
     * well as the lower half of the unit circle.
     *---------------------------------------------------------------------*/

    gp = G[0] + G[N-1];
    gn = G[0] - G[N-1];

    for (j = 0; j < M; j+=2)
    {
        r[j]   = 2.0f*r[j]   + gp;
        r[j+1] = 2.0f*r[j+1] + gn;
    }
    r[j] = 2.0f*r[j] + gp;

    return;
}

/*---------------------------------------------------------------------*
 *  zeros2poly()
 *
 *  Computes the coefficients of the polynomials
 *
 *      R(x) =   prod  (x - x[i]),
 *             i = 0,2,4,...
 *
 *      S(x) =   prod  (x - x[i]),
 *             i = 1,3,5,...
 *
 *  given their zeros x[i] for i = 0, 1, ..., n-1. The routine
 *  assumes n = 1 or even n greater than or equal to 4.
 *
 *  The polynomial coefficients are returned in R[0:n/2-1] and
 *  S[0:n/2-1]. The leading coefficients are in R[0] and S[0].
 *
 *  In this implementation, n is set to a compile time constant.
 *---------------------------------------------------------------------*/

static void zeros2poly(
    float x[],           /* i:  Zeros of R(x) and S(x)      */
    float R[],           /* o:  Coefficients of R(x)        */
    float S[]            /* o:  Coefficients of S(x)        */
)
{
    float   xr, xs;
    Word16  i, j;

    R[0] = 1.0f;
    R[1] = -*x++;
    S[0] = 1.0f;
    S[1] = -*x++;

    for (i = 2; i <= NC; i++)
    {
        xr = -*x++;
        xs = -*x++;

        R[i] = xr*R[i-1];
        S[i] = xs*S[i-1];

        for (j = i-1; j > 0; j--)
        {
            R[j] += xr*R[j-1];
            S[j] += xs*S[j-1];
        }
    }

    return;
}

/*---------------------------------------------------------------------*
 *  polydecomp()
 *
 *  Computes the coefficients of the symmetric and antisymmetric
 *  polynomials P(z) and Q(z) that define the line spectrum pair
 *  decomposition of a given polynomial A(z) of order n. For even n,
 *
 *      P(z) = [A(z) + z^(n+1) A(1/z)]/(1/z + 1),
 *      Q(z) = [A(z) - z^(n+1) A(1/z)]/(1/z - 1),
 *
 *  These polynomials are then expressed in their direct form,
 *  respectively, R(x) and S(x), on the real axis x = cos w using
 *  explicit Chebyshev polynomials of the first kind.
 *
 *  The coefficients of the polynomials R(x) and S(x) are returned
 *  in R[0:n/2] and S[0:n/2] for the given linear prediction
 *  coefficients A[0:n/2]. Note that R(x) and S(x) are formed in
 *  place such that P(z) is stored in the same array than R(x),
 *  and Q(z) is stored in the same array than S(x).
 *
 *  The routines assumes n = 16.
 *---------------------------------------------------------------------*/

static void polydecomp(
    float A[],            /* i: linear prediction coefficients */
    float R[],            /* o: coefficients of R(x)           */
    float S[]             /* o: coefficients of S(x)           */
)
{
    float *P = &R[0], *Q = &S[0];
    Word16 i, j;

    P[0] = A[0];
    Q[0] = A[0];
    for (i = 1, j = M; i <= NC; i++, j--)
    {
        P[i] = A[i] + A[j] - P[i-1];
        Q[i] = A[i] - A[j] + Q[i-1];
    }

    cheb2poly(P);
    cheb2poly(Q);

    return;
}

/*---------------------------------------------------------------------*
 *  cheb2poly()
 *
 *  Computes the coefficients of the explicit Chebyshev polynomial
 *  P(x) = P[0]*x^n + P[1]*x^(n-1) + ... + P[n] given the coefficients
 *  of the series
 *
 *     C(x) = C[0]*T_n(x) + C[1]*T_n-1(x) + ... + C[n]*T_0(x),
 *
 *  where T_n(x) is the nth Chebyshev polynomial of the first kind.
 *  This implementation assumes C[0] = 1. Only the value n = 8 is
 *  supported.
 *
 *  The conversion from C(x) to P(x) is done in place such that the
 *  coefficients of C(x) are given in P[0:8] and those of P(x) are
 *  returned in the same array.
 *---------------------------------------------------------------------*/

static void cheb2poly(
    float P[]       /* i/o: The coefficients of C(x) and P(x) */
)
{
    float c1, c2, c3, c4, c5, c6, c7, c8;

    c1 = P[1];
    c2 = P[2];
    c3 = P[3];
    c4 = P[4];
    c5 = P[5];
    c6 = P[6];
    c7 = P[7];
    c8 = P[8];

    P[0] = 128.0f;
    P[1] =  64.0f*c1;
    P[2] =  32.0f*c2 - 256.0f;
    P[3] =  16.0f*c3 - 112.0f*c1;
    P[4] = 160.0f    -  48.0f*c2 + 8.0f*c4;
    P[5] =  56.0f*c1 -  20.0f*c3 + 4.0f*c5;
    P[6] =  18.0f*c2 -  32.0f    - 8.0f*c4 + 2.0f*c6;
    P[7] =   5.0f*c3 -   7.0f*c1 - 3.0f*c5 +      c7;
    P[8] =   1.0f     -       c2 +      c4 -      c6 + 0.5f*c8;

    return;
}
