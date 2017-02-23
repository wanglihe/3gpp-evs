/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include "basop_proto_func.h"

/*-------------------------------------------------------------------*
 * Local functions
 *-------------------------------------------------------------------*/

static float chebps2(const float x, const float *f, const short n);
static float LPC_chebyshev (float x, float f[], int n);
static void get_isppol( const float *isp, float f[], const short n );

/*---------------------------------------------------------------------*
 * a2isp()
 *
 * Compute the ISPs from the LPC coefficients a[] using Chebyshev
 * polynomials. The found ISPs are in the cosine domain with values
 * in the range from 1 down to -1.
 * The table grid[] contains the points (in the cosine domain) at
 * which the polynomials are evaluated.
 *
 * The ISPs are the roots of the two polynomials F1(z) and F2(z)
 * defined as
 *               F1(z) = [A(z) + z^-M A(z^-1)]
 *  and          F2(z) = [A(z) - z^-M A(z^-1)]/ (1-z^-2)
 *
 * For an even order M=2N, F1(z) has M/2 conjugate roots on the unit circle
 * and F2(z) has M/2-1 conjugate roots on the unit circle in addition to two
 * roots at 0 and pi.
 *
 * For a 16th order LP analysis (M=16), F1(z) and F2(z) can be written as
 *
 *   F1(z) = (1 + a[16]) *  PRODUCT  (1 - 2 cos(w_i) z^-1 + z^-2 )
 *                         i=0,2,...,14
 *
 *   F2(z) = (1 - a[16]) *  PRODUCT  (1 - 2 cos(w_i) z^-1 + z^-2 )
 *                         i=1,3,...,13
 *
 * The ISPs are frequencies w_i, i=0...M-2 plus the last predictor
 * coefficient a[M].
 *---------------------------------------------------------------------*/

void a2isp(
    const float *a,        /* i:   LP filter coefficients     */
    float *isp,      /* o:   Immittance spectral pairs  */
    const float *old_isp   /* i:   ISP vector from past frame */
)
{
    short j, i, nf, ip, order;
    float xlow,ylow,xhigh,yhigh,xmid,ymid,xint;
    float f1[NC+1], f2[NC];
    float *coef;

    /*-------------------------------------------------------------*
     * find the sum and diff polynomials F1(z) and F2(z)
     *     F1(z) = [A(z) + z^M A(z^-1)]
     *     F2(z) = [A(z) - z^M A(z^-1)]/(1-z^-2)
     *
     * for (i=0; i<NC; i++)
     * {
     *     f1[i] = a[i] + a[M-i];
     *     f2[i] = a[i] - a[M-i];
     * }
     * f1[NC] = 2.0*a[NC];
     *
     * for (i=2; i<NC; i++)            Divide by (1-z^-2)
     *     f2[i] += f2[i-2];
     *-------------------------------------------------------------*/

    for (i=0; i<NC; i++)
    {
        f1[i] = a[i] + a[M-i];
        f2[i] = a[i] - a[M-i];
    }

    f1[NC] = 2.0f*a[NC];
    for (i=2; i<NC; i++)        /* divide by (1 - z^-2) */
    {
        f2[i] += f2[i-2];
    }

    /*-----------------------------------------------------------------*
     * Find the ISPs (roots of F1(z) and F2(z) ) using the
     * Chebyshev polynomial evaluation.
     * The roots of F1(z) and F2(z) are alternatively searched.
     * We start by finding the first root of F1(z) then we switch
     * to F2(z) then back to F1(z) and so on until all roots are found.
     *
     *  - Evaluate Chebyshev pol. at grid points and check for sign change.
     *  - If sign change track the root by subdividing the interval
     *    4 times and ckecking sign change.
     *-----------------------------------------------------------------*/

    nf = 0;      /* number of found frequencies */
    ip = 0;      /* flag to first polynomial   */

    coef = f1;   /* start with F1(z) */
    order = NC;

    xlow = grid100[0];
    ylow = chebps2( xlow, coef, order );

    j = 0;

    while ( (nf < M-1) && (j < GRID100_POINTS) )
    {
        j++;
        xhigh = xlow;
        yhigh = ylow;
        xlow = grid100[j];
        ylow = chebps2( xlow, coef, order );
        if ( ylow*yhigh <= 0.0f )  /* if sign changes new root exists */
        {
            j--;

            /* divide the interval of sign change by NO_ITER */
            for (i = 0; i < NO_ITER; i++)
            {
                xmid = (float)(0.5f * (xlow + xhigh));
                ymid = chebps2(xmid,coef,order);

                if (ylow*ymid <= 0.0f)
                {
                    yhigh = ymid;
                    xhigh = xmid;
                }
                else
                {
                    ylow = ymid;
                    xlow = xmid;
                }
            }

            /*--------------------------------------------------------*
             * Linear interpolation
             * xint = xlow - ylow*(xhigh-xlow)/(yhigh-ylow)
             *--------------------------------------------------------*/

            xint = xlow - ylow*(xhigh-xlow)/(yhigh-ylow);
            isp[nf] = xint;    /* new root */
            nf++;

            ip = 1 - ip;               /* flag to other polynomial    */
            coef = ip ? f2 : f1;        /* pointer to other polynomial */

            order = ip ? (NC-1) : NC;  /* order of other polynomial   */

            xlow = xint;
            ylow = chebps2(xlow,coef,order);

        }
    }

    isp[M-1] = a[M];

    /*-----------------------------------------------------------------*
     * Check if m-1 roots found, if not use the ISPs from previous frame
     *-----------------------------------------------------------------*/

    if( (nf < M-1) || (a[M] > 1.0f) || (a[M] < -1.0f) )
    {
        for( i=0; i<M; i++ )
        {
            isp[i] = old_isp[i];
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * isp2a()
 *
 * Convert ISPs to predictor coefficients a[]
 *-------------------------------------------------------------------*/

void isp2a(
    const float *isp,  /* i  : ISP vector (in the cosine domain) */
    float *a,    /* o  : LP filter coefficients            */
    const short m      /* i  : order of LP analysis              */
)
{
    float f1[NC16k+1], f2[NC16k];
    short i, j;
    short nc;

    nc = m/2;

    /*-----------------------------------------------------------------*
     *  Find the polynomials F1(z) and F2(z)                           *
     *-----------------------------------------------------------------*/

    get_isppol(&isp[0], f1, nc);
    get_isppol(&isp[1], f2, (short)(nc-1));

    /*-----------------------------------------------------------------*
     *  Multiply F2(z) by (1 - z^-2)                                   *
     *-----------------------------------------------------------------*/

    for (i = nc-1; i > 1; i--)
    {
        f2[i] -= f2[i-2];
    }

    /*-----------------------------------------------------------------*
     *  Scale F1(z) by (1+isp[m-1])  and  F2(z) by (1-isp[m-1])        *
     *-----------------------------------------------------------------*/

    for (i = 0; i < nc; i++)
    {
        f1[i] *= (1.0f + isp[m-1]);
        f2[i] *= (1.0f - isp[m-1]);
    }

    /*-----------------------------------------------------------------*
     *  A(z) = (F1(z)+F2(z))/2                                         *
     *  F1(z) is symmetric and F2(z) is asymmetric                     *
     *-----------------------------------------------------------------*/

    a[0] = 1.0f;
    for (i=1, j=m-1; i<nc; i++, j--)
    {
        a[i] = (float)(0.5f * (f1[i] + f2[i]));
        a[j] = (float)(0.5f * (f1[i] - f2[i]));
    }
    a[nc] = (float)(0.5f * f1[nc] * (1.0f + isp[m-1]));
    a[m] = isp[m-1];

    return;
}

/*-------------------------------------------------------------------*
 * a2lsp()
 *
 * Convert predictor coefficients a[] to LSPs
 *-------------------------------------------------------------------*/

short a2lsp (
    float *freq,        /* o  : LSP vector                  */
    const float *a_in,        /* i  : predictor coefficients      */
    const short order         /* i  : order of LP analysis        */
)
{
    short i, iswitch, offset, STEPindex;
    short lspnumber, root, notlast, order_by_2;
    float temp, temp2;
    float q[20], prev[2];
    float frequency, LastFreq, STEP;
    const float *a;
    float space_min;

    a = &(a_in[1]);
    order_by_2 = order / 2;
    LastFreq = 0;

    /* calculate q[z] and p[z] , they are all stored in q */
    offset = order_by_2;
    q[0] = (float)(a[0] + a[order - 1] - 1.0);
    q[offset] = (float)(a[0] - a[order - 1] + 1.0);
    for (i = 1; i < order_by_2; i++)
    {
        q[i] = a[i] + a[order - 1 - i] - q[i - 1];
        q[i + offset] = a[i] - a[order - 1 - i] + q[i - 1 + offset];
    }

    q[order_by_2 - 1] = q[order_by_2 - 1] / 2;
    q[order_by_2 - 1 + offset] = q[order_by_2 - 1 + offset] / 2;

    prev[0] = 9e9f;
    prev[1] = 9e9f;
    lspnumber = 0;
    notlast = 1;
    iswitch = 0;
    frequency = 0;

    while ( notlast )
    {
        root = 1;
        offset = iswitch * order_by_2;
        STEPindex = 0;        /* Start with low resolution grid */
        STEP = STEPS[STEPindex];
        while (root)
        {
            temp = (float)cos (frequency * 6.2832);
            if (order >= 4)
            {
                temp2 = LPC_chebyshev (temp, q + offset, order_by_2);
            }
            else
            {
                temp2 = temp + q[0 + offset];
            }

            if ((temp2 * prev[iswitch]) <= 0.0 || frequency >= 0.5)
            {
                if (STEPindex == STEPSNUM - 1)
                {
                    if (fabs (temp2) < fabs (prev[iswitch]))
                    {
                        freq[lspnumber] = frequency;
                    }
                    else
                    {
                        freq[lspnumber] = frequency - STEP;
                    }
                    if ((prev[iswitch]) < 0.0)
                    {
                        prev[iswitch] = 9e9f;
                    }
                    else
                    {
                        prev[iswitch] = -9e9f;
                    }
                    root = 0;
                    frequency = LastFreq;
                    STEPindex = 0;
                }
                else
                {
                    if (STEPindex == 0)
                    {
                        LastFreq = frequency;
                    }
                    frequency -= STEPS[++STEPindex];    /* Go back one grid step */
                    STEP = STEPS[STEPindex];
                }
            }
            else
            {
                prev[iswitch] = temp2;
                frequency += STEP;
            }
        } /* while(root) */

        lspnumber++;
        if (lspnumber > order - 1)
        {
            notlast = 0;
        }
        iswitch = 1 - iswitch;
    } /* while (notlast) */

    /* stability check */
    space_min = 1;
    for( i=1; i < order; i++ )
    {
        space_min = ((freq[i] - freq[i-1]) < space_min)?(freq[i] - freq[i-1]):space_min;
    }

    if( space_min <= 0 )
    {
        return 0;
    }

    return 1;
}

/*-------------------------------------------------------------------*
 * lsp2a()
 *
 * Convert LSPs to predictor coefficients a[]
 *-------------------------------------------------------------------*/

void lsp2a (
    float *pc_in,       /* i/o: predictor coefficients  */
    float *freq,        /* i/o: LSP coefficients        */
    const short order         /* i  : order of LP analysis    */
)
{
    float p[LPC_SHB_ORDER], q[LPC_SHB_ORDER];
    float a[LPC_SHB_ORDER], a1[LPC_SHB_ORDER], a2[LPC_SHB_ORDER];
    float b[LPC_SHB_ORDER], b1[LPC_SHB_ORDER], b2[LPC_SHB_ORDER];

    float xx;
    int i, k;
    int lspflag;
    float *pc;
    int order_by_2;

    lspflag = 1;
    pc = &(pc_in[1]);

    order_by_2 = order / 2;

    /* check input for ill-conditioned cases */
    if ((freq[0] <= 0.0) || (freq[order - 1] >= 0.5))
    {
        lspflag = 0;

        if (freq[0] <= 0)
        {
            freq[0] = 0.022f;
        }

        if (freq[order - 1] >= 0.5)
        {
            freq[order - 1] = 0.499f;
        }
    }

    if (!lspflag)
    {
        xx = (freq[order - 1] - freq[0]) * recip_order[order];
        for (i = 1; i < order; i++)
        {
            freq[i] = freq[i - 1] + xx;
        }
    }

    for (i = 0; i <= order_by_2; i++)
    {
        a[i] = 0.;
        a1[i] = 0.;
        a2[i] = 0.;
        b[i] = 0.;
        b1[i] = 0.;
        b2[i] = 0.;
    }

    for (i = 0; i < order_by_2; i++)
    {
        p[i] = (float)cos (6.2832 * freq[2 * i]);
        q[i] = (float)cos (6.2832 * freq[2 * i + 1]);
    }

    a[0] = 0.25f;
    b[0] = 0.25f;

    for (i = 0; i < order_by_2; i++)
    {
        a[i + 1] = a[i] - 2 * p[i] * a1[i] + a2[i];
        b[i + 1] = b[i] - 2 * q[i] * b1[i] + b2[i];
        a2[i] = a1[i];
        a1[i] = a[i];
        b2[i] = b1[i];
        b1[i] = b[i];
    }
    a[0] = 0.25f;
    b[0] = -0.25f;

    for (i = 0; i < order_by_2; i++)
    {
        a[i + 1] = a[i] - 2 * p[i] * a1[i] + a2[i];
        b[i + 1] = b[i] - 2 * q[i] * b1[i] + b2[i];
        a2[i] = a1[i];
        a1[i] = a[i];
        b2[i] = b1[i];
        b1[i] = b[i];
    }

    pc[0] = 2 * (a[order_by_2] + b[order_by_2]);

    for (k = 2; k <= order; k++)
    {
        a[0] = 0.0f;
        b[0] = 0.0f;

        for (i = 0; i < order_by_2; i++)
        {
            a[i + 1] = a[i] - 2 * p[i] * a1[i] + a2[i];
            b[i + 1] = b[i] - 2 * q[i] * b1[i] + b2[i];
            a2[i] = a1[i];
            a1[i] = a[i];
            b2[i] = b1[i];
            b1[i] = b[i];
        }
        pc[k - 1] = 2 * (a[order_by_2] + b[order_by_2]);
    }

    return;
}

/*-----------------------------------------------------------*
 * get_lsppol()
 *
 * Find the polynomial F1(z) or F2(z) from the LSPs.
 * This is performed by expanding the product polynomials:
 *
 * F1(z) =   product   ( 1 - 2 LSF_i z^-1 + z^-2 )
 *         i=0,2,4,..,n-2
 * F2(z) =   product   ( 1 - 2 LSF_i z^-1 + z^-2 )
 *         i=1,3,5,..,n-1
 *
 * where LSP_i are the LSPs in the cosine domain.
 *-----------------------------------------------------------*/

static void get_lsppol(
    const float lsp[],    /* i  : line spectral freq. (cosine domain)  */
    float f[],            /* o  : the coefficients of F1 or F2         */
    const short n,        /* i  : no of coefficients (m/2)             */
    short flag            /* i  : 1 ---> F1(z) ; 2 ---> F2(z)          */
)
{
    float b;
    const float *plsp;
    int i,j;

    plsp = lsp + flag - 1;

    f[0] = 1.0f;
    b = -2.0f * *plsp;
    f[1] = b;

    for (i = 2; i <= n; i++)
    {
        plsp += 2;
        b = -2.0f * *plsp;
        f[i] = b*f[i-1] + 2.0f*f[i-2];

        for (j = i-1; j > 1; j--)
        {
            f[j] += b*f[j-1] + f[j-2];
        }

        f[1] += b;
    }

    return;
}

/*---------------------------------------------------------------------*
 * a2lsp_stab()
 *
 * Compute the LSPs from the LPC coefficients a[] using Chebyshev
 * polynomials. The found LSPs are in the cosine domain with values
 * in the range from 1 down to -1.
 * The table grid[] contains the points (in the cosine domain) at
 * which the polynomials are evaluated.
 *
 * The ISPs are the roots of the two polynomials F1(z) and F2(z)
 * defined as
 *               F1(z) = [A(z) + z^-M A(z^-1)]/ (1+z^-1)
 *  and          F2(z) = [A(z) - z^-M A(z^-1)]/ (1-z^-1)
 *---------------------------------------------------------------------*/

void a2lsp_stab(
    const float *a,        /* i:   LP filter coefficients     */
    float *lsp,      /* o:   LSP vector                 */
    const float *old_lsp   /* i:   LSP vector from past frame */
)
{
    short j, i, nf, ip;
    float xlow, ylow, xhigh, yhigh, xmid, ymid, xint;
    float *pf1, *pf2;
    const float *pa1, *pa2;
    float f1[NC+1], f2[NC+1];

    /*-------------------------------------------------------------*
     * find the sum and diff polynomials F1(z) and F2(z)           *
     *      F1(z) = [A(z) + z^11 A(z^-1)]/(1+z^-1)                 *
     *      F2(z) = [A(z) - z^11 A(z^-1)]/(1-z^-1)                 *
     *-------------------------------------------------------------*/

    pf1 = f1;                               /* Equivalent code using indices   */
    pf2 = f2;
    *pf1++ = 1.0f;                          /* f1[0] = 1.0;                    */
    *pf2++ = 1.0f;                          /* f2[0] = 1.0;                    */
    pa1 = a + 1;
    pa2 = a + M;

    for( i=0; i<=NC-1; i++ )                /* for (i=1, j=M; i<=NC; i++, j--) */
    {
        *pf1 = *pa1 + *pa2 - *(pf1-1);    /* f1[i] = a[i]+a[j]-f1[i-1];   */
        *pf2 = *pa1++ - *pa2-- + *(pf2-1);/* f2[i] = a[i]-a[j]+f2[i-1];   */
        pf1++;
        pf2++;
    }

    /*---------------------------------------------------------------------*
     * Find the LSPs (roots of F1(z) and F2(z) ) using the                 *
     * Chebyshev polynomial evaluation.                                    *
     * The roots of F1(z) and F2(z) are alternatively searched.            *
     * We start by finding the first root of F1(z) then we switch          *
     * to F2(z) then back to F1(z) and so on until all roots are found.    *
     *                                                                     *
     *  - Evaluate Chebyshev pol. at grid points and check for sign change.*
     *  - If sign change track the root by subdividing the interval        *
     *    4 times and ckecking sign change.                                *
     *---------------------------------------------------------------------*/

    nf=0;      /* number of found frequencies */
    ip=0;      /* flag to first polynomial   */

    pf1 = f1;  /* start with F1(z) */

    xlow = grid100[0];
    ylow = chebps2( xlow, pf1, NC );

    j = 0;
    while( (nf < M) && (j < GRID100_POINTS) )
    {
        j++;
        xhigh = xlow;
        yhigh = ylow;
        xlow = grid100[j];
        ylow = chebps2( xlow, pf1, NC );

        if( ylow*yhigh <= 0.0 )  /* if sign change new root exists */
        {
            j--;
            /* divide the interval of sign change by NO_ITER */
            for (i = 0; i < NO_ITER; i++)
            {
                xmid = 0.5f * ( xlow + xhigh );
                ymid = chebps2( xmid, pf1, NC );
                if( ylow*ymid <= 0.0 )
                {
                    yhigh = ymid;
                    xhigh = xmid;
                }
                else
                {
                    ylow = ymid;
                    xlow = xmid;
                }
            }

            /* linear interpolation for evaluating the root */
            xint = xlow - ylow*(xhigh-xlow)/(yhigh-ylow);
            lsp[nf] = xint;    /* new root */
            nf++;
            ip = 1 - ip;         /* flag to other polynomial    */
            pf1 = ip ? f2 : f1;  /* pointer to other polynomial */
            xlow = xint;
            ylow = chebps2( xlow, pf1, NC );
        }
    }

    /* Check if M roots found */
    /* if not use the LSPs from previous frame */
    if( nf < M )
    {

        for( i=0; i<M; i++ )
        {
            lsp[i] = old_lsp[i];
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * lsp2a_stab()
 *
 * Convert LSPs to predictor coefficients A[]
 *-------------------------------------------------------------------*/

void lsp2a_stab(
    const float *lsp,  /* i  : LSF vector (in the cosine domain) */
    float *a,          /* o  : LP filter coefficients            */
    const short m      /* i  : order of LP analysis              */
)
{
    float f1[NC+1], f2[NC+1];
    short i, k,nc;
    float *pf1, *pf2, *pf1_1, *pf2_1, *pa1, *pa2;

    nc=m/2;

    /*-----------------------------------------------------*
     *  Find the polynomials F1(z) and F2(z)               *
     *-----------------------------------------------------*/

    get_lsppol(lsp,f1,nc,1);
    get_lsppol(lsp,f2,nc,2);

    /*-----------------------------------------------------*
     *  Multiply F1(z) by (1+z^-1) and F2(z) by (1-z^-1)   *
     *-----------------------------------------------------*/

    pf1 = f1 + nc;
    pf1_1 = pf1 - 1;
    pf2 = f2 + nc;                      /* Version using indices            */
    pf2_1 = pf2 - 1;
    k = nc-1;
    for (i = 0; i <= k; i++)            /* for (i = NC; i > 0; i--)         */
    {
        *pf1-- += *pf1_1--;             /*   f1[i] += f1[i-1];              */
        *pf2-- -= *pf2_1--;             /*   f2[i] -= f2[i-1];              */
    }

    /*-----------------------------------------------------*
     *  A(z) = (F1(z)+F2(z))/2                             *
     *  F1(z) is symmetric and F2(z) is antisymmetric      *
     *-----------------------------------------------------*/

    pa1 = a;
    *pa1++ = 1.0;                       /* a[0] = 1.0;                      */
    pa2 = a + m;
    pf1 = f1 + 1;
    pf2 = f2 + 1;
    for (i = 0; i <= k; i++)            /* for (i=1, j=M; i<=NC; i++, j--)  */
    {
        *pa1++ = 0.5f*(*pf1 + *pf2);    /*   a[i] = 0.5*(f1[i] + f2[i]);    */
        *pa2-- = 0.5f*(*pf1++ - *pf2++);/*   a[j] = 0.5*(f1[i] - f2[i]);    */
    }

    return;
}

/*---------------------------------------------------------------------------
 * reorder_lsf()
 *
 * To make sure that the LSFs are properly ordered and to keep a certain
 * minimum distance between consecutive LSFs.
 *--------------------------------------------------------------------------*/

void reorder_lsf(
    float *lsf,            /* i/o: LSF vector */
    const float min_dist,  /* i  : minimum required distance */
    const short n,         /* i  : LPC order                 */
    const float fs         /* i  : sampling frequency        */
)
{
    short i;
    float lsf_min;
    float lsf_max;

    /*-----------------------------------------------------------------*
     * Verify the LSF ordering and minimum GAP
     *-----------------------------------------------------------------*/

    lsf_min = min_dist;
    for (i = 0; i < n; i++)
    {
        if (lsf[i] < lsf_min)
        {
            lsf[i] = lsf_min;
        }
        lsf_min = lsf[i] + min_dist;
    }

    /*------------------------------------------------------------------------------------------*
     * Reverify the LSF ordering and minimum GAP in the reverse order (security)
     *------------------------------------------------------------------------------------------*/

    lsf_max = fs/2.0f - min_dist;

    if( lsf[n-1] > lsf_max )        /* If danger of unstable filter in case of resonance in HF */
    {
        for (i = n-1; i >=0; i--)   /* Reverify the minimum LSF gap in the reverse sens */
        {
            if (lsf[i] > lsf_max)
            {
                lsf[i] = lsf_max;
            }

            lsf_max = lsf[i] - min_dist;
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * space_lsfs()
 *
 *-------------------------------------------------------------------*/

void space_lsfs (
    float *lsfs,        /* i/o: Line spectral frequencies   */
    const short order         /* i  : order of LP analysis        */
)
{
    float delta;
    short i, flag = 1;

    while (flag == 1)
    {
        flag = 0;
        for (i = 0; i <= order; i++)
        {
            delta = (float)( i == 0 ? lsfs[0] : (i == order ? 0.5f - lsfs[i - 1] : (lsfs[i] - lsfs[i -1])));
            if (delta < SPC)
            {
                flag = 1;
                delta -= SPC_plus;
                if (i == order)
                {
                    lsfs[i - 1] += delta;
                }
                else
                {
                    if (i == 0)
                    {
                        lsfs[i] -= delta;
                    }
                    else
                    {
                        delta *= 0.5f;
                        lsfs[i - 1] += delta;
                        lsfs[i] -= delta;
                    }
                }
            }
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * lsp_weights()
 *
 *-------------------------------------------------------------------*/

void lsp_weights (
    const float *lsps,        /* i  : Line spectral pairs         */
    float *weight,      /* o  : weights                     */
    const short order         /* i  : order of LP analysis        */
)
{
    short i;
    float delta1, delta2;

    for (i = 0; i < order; i++)
    {
        delta1 = (float)((i == 0) ? lsps[i] : lsps[i] - lsps[i - 1]);
        delta2 = (float)((i == order - 1) ? 0.5f - lsps[i] : lsps[i + 1] - lsps[i]);

        weight[i] = 250 * root_a_over_b(ALPHA_SQ, delta1 * delta2);
    }

    if (order != LPC_SHB_ORDER_WB)
    {
        weight[3] *= 1.1f;
        weight[4] *= 1.1f;
    }

    return;
}


/*-----------------------------------------------------------------------*
 * isf2isp()
 *
 * Transformation of ISF to ISP
 *
 * ISP are immitance spectral pairs in cosine domain (-1 to 1).
 * ISF are immitance spectral pairs in frequency domain (0 to fs/2).
 *-----------------------------------------------------------------------*/

void isf2isp(
    const float isf[],  /* i  : isf[m] normalized (range: 0<=val<=fs/2) */
    float isp[],  /* o  : isp[m] (range: -1<=val<1)               */
    const short m,      /* i  : LPC order                               */
    const float fs      /* i  : sampling frequency                      */
)
{
    short i;

    for(i=0; i<m-1; i++)
    {
        isp[i] = (float)cos(isf[i] * EVS_PI / (fs/2));
    }
    isp[m-1] = (float)cos(isf[m-1] * EVS_PI / (fs/2) * 2.0);

    return;
}


/*-----------------------------------------------------------------------*
 * isp2isf()
 *
 * Transformation of ISP to ISF
 *
 * ISP are immitance spectral pair in cosine domain (-1 to 1).
 * ISF are immitance spectral pair in frequency domain (0 to fs/2).
 *-----------------------------------------------------------------------*/

void isp2isf(
    const float isp[],  /* i  : isp[m] (range: -1<=val<1)               */
    float isf[],  /* o  : isf[m] normalized (range: 0<=val<=fs/2) */
    const short m,      /* i  : LPC order                               */
    const float fs      /* i  : sampling frequency                      */
)
{
    short i;

    /* convert ISPs to frequency domain 0..6400 */
    for(i=0; i<m-1; i++)
    {
        isf[i] = (float)(acos(isp[i]) * ((fs/2) / EVS_PI));
    }
    isf[m-1] = (float)(acos(isp[m-1]) * ((fs/2) / EVS_PI) * 0.5f);

    return;
}

/*-------------------------------------------------------------------*
 * a2rc()
 *
 * Convert from LPC to reflection coeff
 *-------------------------------------------------------------------*/

void a2rc(
    const float *a,           /* i  : LPC coefficients            */
    float *refl,        /* o  : Reflection co-efficients    */
    const short lpcorder      /* i  : LPC order                   */
)
{
    float f[M];
    short m, j, n;
    float km, denom, x;

    for (m = 0; m < lpcorder; m++)
    {
        f[m] = -a[m];
    }

    /* Initialization */
    for (m = lpcorder - 1; m >= 0; m--)
    {
        km = f[m];
        if (km <= -1.0 || km >= 1.0)
        {
            for ( j = 0; j < lpcorder; j++ )
            {
                refl[j] = 0.f;
            }

            return;
        }

        refl[m] = -km;
        denom = 1.0f / (1.0f - km * km);

        for (j = 0; j < m / 2; j++)
        {
            n = m - 1 - j;
            x = denom * f[j] + km * denom * f[n];
            f[n] = denom * f[n] + km * denom * f[j];
            f[j] = x;
        }

        if (m & 1)
        {
            f[j] = denom * f[j] + km * denom * f[j];
        }
    }

    return;
}


/*----------------------------------------------------------------------------------*
 * vq_dec_lvq()
 *
 * Multi-stage VQ decoder for LSF parameters, last stage LVQ
 *----------------------------------------------------------------------------------*/

short vq_dec_lvq (
    short sf_flag,        /* i  : safety net flag                           */
    float x[],            /* o  : Decoded vector                            */
    short indices[],      /* i  : Indices                                   */
    short stages,         /* i  : Number of stages                          */
    short N,              /* i  : Vector dimension                          */
    short mode,           /* (i): mode_lvq, or mode_lvq_p */
    short no_bits,        /* (i): no. bits for lattice */
    unsigned int *p_offset_scale1,
    unsigned int *p_offset_scale2,
    unsigned int *p_offset_scale1_p,
    unsigned int *p_offset_scale2_p,
    short *p_no_scales,
    short *p_no_scales_p
)
{
    float  x_lvq[16];
    short i;
    short ber_flag;
    /* clear vector */
    set_f(x, 0, N);

    /*-----------------------------------------------*
     * add contribution of each stage
     *-----------------------------------------------*/

    if (sf_flag == 1)
    {
        for(i=0; i<stages-1; i++)
        {
            v_add(x, &Quantizers[CB[mode]+i][indices[i]*N], x, N);
        }

        ber_flag =
            deindex_lvq(&indices[stages-1], x_lvq, mode, sf_flag, no_bits, p_offset_scale1, p_offset_scale2, p_no_scales);
    }
    else
    {
        for(i=0; i<stages-1; i++)
        {
            v_add(x, &Quantizers_p[CB_p[mode]+i][indices[i]*N], x, N);
        }

        ber_flag =
            deindex_lvq(&indices[stages-1], x_lvq, mode, sf_flag, no_bits, p_offset_scale1_p, p_offset_scale2_p, p_no_scales_p);
    }

    v_add(x, x_lvq, x, N);

    return ber_flag;
}


/*----------------------------------------------------------------------------------*
 * lsf_allocate()
 *
 * Calculate number of stages and levels for each stage based on the allowed bit allocation
 *----------------------------------------------------------------------------------*/

void lsf_allocate(
    const short nBits,       /* i  : Number of bits to use for quantization     */
    const short framemode,   /* i  : LSF quantizer mode                         */
    const short framemode_p, /* i  : ISF quantizer mode predmode (mode_lvq_p)   */
    short *stages0,    /* o  : Number of stages for safety-net quantizer  */
    short *stages1,    /* o  : Number of stages for predictive quantizer  */
    short levels0[],   /* o  : Number of vectors for each stage for SFNET */
    short levels1[],   /* o  : Number of vectors for each stage for pred  */
    short bits0[],     /* o  : Number of bits for each stage safety net   */
    short bits1[]      /* o  : Number of bits for each stage pred         */
)
{
    short i;
    short cumleft;
    short bits_lvq,  n_stages,nbits0;

    /* VOICED@16kHz */
    if( framemode == 14 )
    {
        return;
    }

    /*---------------------------------------------------*
     * Calculate bit allocation for safety-net quantizer
     *---------------------------------------------------*/

    cumleft = BitsVQ[framemode];
    bits_lvq = nBits-cumleft;
    nbits0 = CBbits[framemode];
    if (nbits0 > -1)
    {
        if (nbits0>0) /* */
        {
            n_stages = 2;
            levels0[0] = CBsizes[nbits0];
            bits0[0] = nbits0;
            bits0[1] = cumleft-nbits0;

            if ( bits0[1] == 0 )
            {
                n_stages--;
            }
            else
            {
                levels0[1] = CBsizes[cumleft-nbits0];
            }
        }
        else /* no bits for VQ stage */
        {
            n_stages = 0;
        }

        *stages0 = n_stages;
        if(bits_lvq > 0)
        {
            bits0[n_stages] = bits_lvq;
            levels0[n_stages] = bits_lvq; /* this is number of bits, not levels */
            *stages0 = n_stages+1;
        }
    }
    else
    {
        *stages0 = 0;
    }

    /*---------------------------------------------------*
     * Calculate bit allocation for predictive quantizer
     *---------------------------------------------------*/

    if ( framemode_p > -1 )
    {
        cumleft = BitsVQ_p[framemode_p];
        bits_lvq = nBits - cumleft;
        nbits0 = CBbits_p[framemode_p];

        if (nbits0 > -1)
        {
            if ( nbits0 > 0 )
            {
                if ( framemode_p == 7 )
                {
                    /* for UNVOICED_WB only */
                    n_stages = 3;
                    for( i=0; i<n_stages; i++ )
                    {
                        levels1[i] = CBsizes[nbits0];
                        bits1[i] = nbits0;
                    }
                    bits1[n_stages] = bits_lvq;
                    levels1[n_stages] = bits_lvq;
                    *stages1 = n_stages+1;
                }
                else
                {
                    n_stages = 1;
                    levels1[0] = CBsizes[nbits0];
                    bits1[0] = nbits0;
                    nbits0 = cumleft-nbits0;
                    if (nbits0>0)
                    {
                        levels1[1] = CBsizes[nbits0];
                        bits1[1] = nbits0;
                        n_stages = 2;
                    }

                    levels1[n_stages] = bits_lvq; /* this is number of bits, not levels */
                    bits1[n_stages] = bits_lvq;
                    *stages1 = n_stages+1;
                }
            }
            else
            {
                *stages1 = 1;
                bits1[0] = bits_lvq;
                levels1[0] = bits_lvq;
            }
        }
        else
        {
            fprintf(stderr, "lsf_allocate(): invalid number of bits in used predictive mode\n");
            exit(-1);
        }
    }

    return;
}

/*----------------------------------------------------------------------------------*
 * find_pred_mode()
 *
 *
 *----------------------------------------------------------------------------------*/

short find_pred_mode(
    const short coder_type,
    const short bwidth,
    const float int_fs,
    short * p_mode_lvq,
    short * p_mode_lvq_p,
    short core_brate
)
{
    short idx, predmode;

    idx = bwidth;
    if (idx>1)
    {
        idx = 1;
    }
    if ((int_fs == INT_FS_16k))
    {
        idx = 2;
    }
    else
    {
        if ((core_brate>=GENERIC_MA_LIMIT)&&(coder_type==GENERIC)
                &&(idx==1)
           )
        {
            idx = 3;
        }
    }

    predmode = predmode_tab[idx][coder_type];

    if (idx<=2)
    {
        *p_mode_lvq = NO_CODING_MODES*idx + coder_type;
        if (predmode>0)
        {
            *p_mode_lvq_p = *p_mode_lvq;
        }
        else
        {
            *p_mode_lvq_p = -1;
        }
    }
    else  /* WB 12.8 with MA pred in GENERIC*/
    {
        *p_mode_lvq = NO_CODING_MODES + coder_type;
        if (coder_type == GENERIC)
        {
            *p_mode_lvq_p = 18;
        }
        else
        {
            if (predmode>0)
            {
                *p_mode_lvq_p = *p_mode_lvq;
            }
            else
            {
                *p_mode_lvq_p = -1;
            }
        }
    }


    if (predmode==-1)
    {
        fprintf (stderr, "find_pred_mode(): incorrect coder_type specification: %d\n", coder_type);
        exit(-1);
    }

    return predmode;
}


/*---------------------------------------------------------------------------
 * reorder_isf()
 *
 * To make sure that the ISFs are properly ordered and to keep a certain
 * minimum distance between consecutive ISFs.
 *--------------------------------------------------------------------------*/

void reorder_isf(
    float *isf,      /* i/o: ISF vector */
    const float min_dist,  /* i  : minimum required distance */
    const short n,         /* i  : LPC order                 */
    const float fs         /* i  : sampling frequency        */
)
{
    short i;
    float isf_min;
    float isf_max;

    /*-----------------------------------------------------------------*
     * Verify the ISF ordering and minimum GAP
     *-----------------------------------------------------------------*/

    isf_min = min_dist;
    for (i = 0; i < n-1; i++)
    {
        if (isf[i] < isf_min)
        {
            isf[i] = isf_min;
        }

        isf_min = isf[i] + min_dist;
    }

    /*------------------------------------------------------------------------------------------*
     * Reverify the ISF ordering and minimum GAP in the reverse order (security)
     *------------------------------------------------------------------------------------------*/

    isf_max = fs/2.0f - min_dist;

    if( isf[n-2] > isf_max )        /* If danger of unstable filter in case of resonance in HF */
    {
        for (i = n-2; i >=0; i--)   /* Reverify the minimum ISF gap in the reverse sens */
        {
            if (isf[i] > isf_max)
            {
                isf[i] = isf_max;
            }

            isf_max = isf[i] - min_dist;
        }
    }

    return;
}

/*----------------------------------------------------------------------------------*
 * lsf_stab()
 *
 * Check LSF stability (distance between old LSFs and current LSFs)
 *----------------------------------------------------------------------------------*/

float lsf_stab(               /* o  : LP filter stability */
    const float *lsf,         /* i  : LSF vector          */
    const float *lsfold,      /* i  : old LSF vector      */
    const short Opt_AMR_WB,   /* i  : flag indicating AMR-WB IO mode */
    const short L_frame       /* i  : frame length */
)
{
    short i, m;
    float scale, stab_fac, tmp;

    tmp = 0.0f;
    if ( Opt_AMR_WB )
    {
        m = M-1;
    }
    else
    {
        m = M;
    }

    for (i=0; i<m; i++)
    {
        tmp += (lsf[i] - lsfold[i]) * (lsf[i] - lsfold[i]);
    }

    scale = (float)L_FRAME / (float)L_frame;
    scale *= scale;

    stab_fac = (float)(1.25f - (scale * tmp/400000.0f));

    if (stab_fac > 1.0f)
    {
        stab_fac = 1.0f;
    }

    if (stab_fac < 0.0f)
    {
        stab_fac = 0.0f;
    }
    return stab_fac;
}


/*---------------------------------------------------------------------
 * get_isppol()
 *
 * Find the polynomial F1(z) or F2(z) from the ISPs.
 * This is performed by expanding the product polynomials:
 *
 * F1(z) =   PRODUCT   ( 1 - 2 isp_i z^-1 + z^-2 )
 *         i=0,2,...,14
 * F2(z) =   PRODUCT   ( 1 - 2 isp_i z^-1 + z^-2 )
 *         i=1,3,...,13
 *
 * where isp_i are the ISPs in the cosine domain.
 *---------------------------------------------------------------------*/

static void get_isppol(
    const float *isp,  /* i  : Immitance spectral pairs (cosine domaine) */
    float f[],   /* o  : the coefficients of F1 or F2              */
    const short n      /* i  : nb. of coefficients (m/2)                 */
)
{
    float b;
    short i,j;

    f[0] = 1.0f;
    b = (float)(-2.0f * *isp);
    f[1] = b;
    for (i = 2; i <= n; i++)
    {
        isp += 2;
        b = (float)(-2.0f * *isp);
        f[i] = (float)(b * f[i-1] + 2.0f * f[i-2]);

        for (j = i-1; j > 1; j--)
        {
            f[j] += b*f[j-1] + f[j-2];
        }
        f[1] += b;
    }

    return;
}


/*---------------------------------------------------------------------
 * chebps2()
 *
 * Evaluates the Chebyshev polynomial series
 *
 * The polynomial order is
 *     n = m/2   (m is the prediction order)
 * The polynomial is given by
 *     C(x) = f(0)T_n(x) + f(1)T_n-1(x) + ... +f(n-1)T_1(x) + f(n)/2
 *---------------------------------------------------------------------*/

static float chebps2(  /* o:   the value of the polynomial C(x)    */
    const float x,     /* i:   value of evaluation; x=cos(freq)    */
    const float *f,    /* i:   coefficients of sum or diff polyn.  */
    const short n      /* i:   order of polynomial                 */
)
{
    float b1, b2, b0, x2;
    short i;


    x2 = (float)(2.0f * x);
    b2 = f[0];

    b1 = x2*b2 + f[1];

    for (i=2; i<n; i++)
    {
        b0 = x2*b1 - b2 + f[i];
        b2 = b1;
        b1 = b0;
    }

    return (float)(x * b1 - b2 + 0.5f * f[n]);
}


/*---------------------------------------------------------------------
 * LPC_chebyshev()
 *
 * Evaluate a series expansion in Chebyshev polynomials
 *
 *  The polynomial order is
 *     n = m/2   (m is the prediction order)
 *  The polynomial is given by
 *    C(x) = T_n(x) + f(1)T_n-1(x) + ... +f(n-1)T_1(x) + f(n)/2
 *---------------------------------------------------------------------*/

static float LPC_chebyshev (float x, float f[], int n)
{
    float b1, b2, b0, x2, val;
    int i;

    x2 = 2.0f * x;
    b2 = 1.0f;
    b1 = x2 + f[0];

    for (i = 2; i < n; i++)
    {
        b0 = x2 * b1 - b2 + f[i - 1];

        /* was previously casting x2 into float to have this
        equation evaluated in float to be same
        same as EVRC-B only code which has 2.0 in equation
        instead of a float x2. This was causing non-bit-exactness
        in a very very very rare corner case.
        Doesnt really matter, but just to be picky! */
        b2 = b1;
        b1 = b0;
    }

    val = (x * b1 - b2 + f[i - 1]);

    return val;
}


/*-------------------------------------------------------------------*
 * lsp2isp()
 *
 * Convert LSPs to ISPs via predictor coefficients A[]
 *-------------------------------------------------------------------*/

void lsp2isp(
    const float *lsp,           /* i  : LSP vector                        */
    float *isp,           /* o  : ISP filter coefficients           */
    float *stable_isp,    /* i/o: ISP filter coefficients           */
    const short m               /* i  : order of LP analysis              */
)
{
    float a[M+1];

    /* LSP --> A */
    lsp2a_stab( lsp, a, m );

    /* A --> ISP */
    a2isp( a, isp, stable_isp );

    /* Update to latest stable ISP */
    mvr2r( isp, stable_isp, M );
}

/*-------------------------------------------------------------------*
 * isp2lsp()
 *
 * Convert ISPs to LSPs via predictor coefficients A[]
 *-------------------------------------------------------------------*/

void isp2lsp(
    const float *isp,           /* i  : LSP vector                        */
    float *lsp,           /* o  : ISP filter coefficients           */
    float *stable_lsp,    /* i/o: stable LSP filter coefficients    */
    const short m               /* i  : order of LP analysis              */
)
{
    float a[M+1];

    /* ISP --> A */
    isp2a( isp, a, m );

    /* A --> LSP */
    a2lsp_stab( a, lsp, stable_lsp );

    /* Update to latest stable LSP */
    mvr2r( lsp, stable_lsp, M );
}


/*-------------------------------------------------------------------*
 * lsf2isf()
 *
 * Convert LSPs to ISPs
 *-------------------------------------------------------------------*/

void lsf2isf(
    const float *lsf,           /* i  : LSF vector                        */
    float *isf,           /* o  : ISF vector                        */
    float *stable_isp,    /* i/o: stable ISP filter coefficients    */
    const short m,              /* i  : order of LP analysis              */
    const float int_fs
)
{
    float tmp_lsp[M];
    float tmp_isp[M];

    /* LSF --> LSP */
    lsf2lsp( lsf, tmp_lsp, m, int_fs );

    /* LSP --> ISP */
    lsp2isp( tmp_lsp, tmp_isp, stable_isp, m );

    /* ISP --> ISF */
    isp2isf( tmp_isp, isf, m, int_fs );

    return;
}

/*-------------------------------------------------------------------*
 * isf2lsf()
 *
 * Convert ISFs to LSFs
 *-------------------------------------------------------------------*/

void isf2lsf(
    const float *isf,           /* i  : ISF vector                        */
    float *lsf,           /* o  : LSF vector                        */
    float *stable_lsp,    /* i/o: stable LSP filter coefficients    */
    const short m,              /* i  : order of LP analysis              */
    const float int_fs
)
{
    float tmp_isp[M];
    float tmp_lsp[M];

    /* ISF --> ISP */
    isf2isp( isf, tmp_isp, m, int_fs );

    /* ISP --> LSP */
    isp2lsp( tmp_isp, tmp_lsp, stable_lsp, m );

    /* LSP --> LSF */
    lsp2lsf( tmp_lsp, lsf, m, int_fs );

    return;
}

/*-----------------------------------------------------------------------*
 * lsp2lsf()
 *
 * Transformation of LSPs to LSFs
 *
 * LSP are line spectral pair in cosine domain (-1 to 1).
 * LSF are line spectral frequencies (0 to fs/2).
 *-----------------------------------------------------------------------*/

void lsp2lsf(
    const float lsp[],  /* i  : isp[m] (range: -1<=val<1)               */
    float lsf[],  /* o  : isf[m] normalized (range: 0<=val<=fs/2) */
    const short m,      /* i  : LPC order                               */
    const float fs      /* i  : sampling frequency                      */
)
{
    short i;

    /* convert LSPs to LSFs */
    for( i=0; i<m; i++ )
    {
        lsf[i] = (float)(acos(lsp[i]) * ((fs/2) / EVS_PI));
    }

    return;
}

/*-----------------------------------------------------------------------*
* lsf2lsp()
 *
 * Transformation of LSFs to LSPs
 *
 * LSP are line spectral pairs in cosine domain (-1 to 1).
 * LSF are line spectral frequencies (0 to fs/2).
 *-----------------------------------------------------------------------*/

void lsf2lsp(
    const float lsf[],  /* i  : isf[m] normalized (range: 0<=val<=fs/2) */
    float lsp[],        /* o  : isp[m] (range: -1<=val<1)               */
    const short m,      /* i  : LPC order                               */
    const float fs      /* i  : sampling frequency                      */
)
{
    short i;

    /* convert LSFs to LSPs */
    for( i=0; i<m; i++ )
    {
        lsp[i] = (float)cos(lsf[i] * EVS_PI / (fs/2));
    }

    return;
}


/*---------------------------------------------------------------------------
  * tcvq_Dec()
  *
  *
  *--------------------------------------------------------------------------*/

static void tcvq_Dec(short *ind, float *d_out, short safety_net)
{
    short i;
    short index[9];
    short stage, state, branch[N_STAGE], codeword[N_STAGE];
    short fins, iwd;
    float pred[N_DIM];
    float D[N_STAGE_VQ][N_DIM];
    const float (*TCVQ_CB_SUB1)[128][2], (*TCVQ_CB_SUB2)[64][2], (*TCVQ_CB_SUB3)[32][2];
    const float (*IntraCoeff)[2][2];

    mvs2s(ind, index, 9);

    if(safety_net)
    {
        TCVQ_CB_SUB1 = SN_TCVQ_CB_SUB1;
        TCVQ_CB_SUB2 = SN_TCVQ_CB_SUB2;
        TCVQ_CB_SUB3 = SN_TCVQ_CB_SUB3;
        IntraCoeff   = SN_IntraCoeff;
    }
    else
    {
        TCVQ_CB_SUB1 = AR_TCVQ_CB_SUB1;
        TCVQ_CB_SUB2 = AR_TCVQ_CB_SUB2;
        TCVQ_CB_SUB3 = AR_TCVQ_CB_SUB3;
        IntraCoeff   = AR_IntraCoeff;
    }

    /* Decode Final state */
    fins = index[0] & 15;

    /* Decode Branch info */
    branch[0] = index[1] >> 4;
    branch[1] = index[2] >> 4;

    for(stage = 2; stage < N_STAGE_VQ-4; stage++)
    {
        branch[stage] = index[stage+1] >> 3;
    }

    branch[4] = fins & 0x1;
    branch[5] = (fins >> 1) & 0x1;
    branch[6] = (fins >> 2) & 0x1;
    branch[7] = (fins >> 3) & 0x1;

    /* Decode Codeword info */
    for(stage = 0; stage < 2; stage++)
    {
        codeword[stage] = (index[stage+1] & 15) << 3;
    }

    for(stage = 2; stage < N_STAGE_VQ-4; stage++)
    {
        codeword[stage] = (index[stage+1] & 7) << 3;
    }

    for(stage = N_STAGE_VQ-4; stage < N_STAGE_VQ; stage++)
    {
        codeword[stage] = (index[stage+1] & 3) << 3;
    }

    state = (fins >> 2) << 2;

    /* stage #1 */
    iwd    = NTRANS2[branch[0]+2][state] + codeword[0];
    D[0][0] = TCVQ_CB_SUB1[0][iwd][0];
    D[0][1] = TCVQ_CB_SUB1[0][iwd][1];
    state  = NTRANS2[branch[0]][state];

    /* stage #2 */
    pred[0] = IntraCoeff[0][0][0] * D[0][0] + IntraCoeff[0][0][1]*D[0][1];
    pred[1]  = IntraCoeff[0][1][0] * D[0][0] + IntraCoeff[0][1][1]*D[0][1];

    iwd    = NTRANS2[branch[1]+2][state] + codeword[1];
    D[1][0] = TCVQ_CB_SUB1[1][iwd][0] + pred[0];
    D[1][1] = TCVQ_CB_SUB1[1][iwd][1] + pred[1];
    state  = NTRANS2[branch[1]][state];

    /* stage #3 - #4 */
    for(stage = 2; stage < N_STAGE_VQ-4; stage++)
    {
        pred[0]     = IntraCoeff[stage-1][0][0] * D[stage-1][0] + IntraCoeff[stage-1][0][1]*D[stage-1][1];
        pred[1]      = IntraCoeff[stage-1][1][0] * D[stage-1][0] + IntraCoeff[stage-1][1][1]*D[stage-1][1];

        iwd      = NTRANS2[branch[stage]+2][state] + codeword[stage];
        D[stage][0] = TCVQ_CB_SUB2[stage-2][iwd][0] + pred[0];
        D[stage][1] = TCVQ_CB_SUB2[stage-2][iwd][1] + pred[1];
        state    = NTRANS2[branch[stage]][state];
    }

    /* stage #5 - #8 */
    for(stage = N_STAGE_VQ-4; stage < N_STAGE_VQ; stage++)
    {
        pred[0]    = IntraCoeff[stage-1][0][0] * D[stage-1][0] + IntraCoeff[stage-1][0][1]*D[stage-1][1];
        pred[1]    = IntraCoeff[stage-1][1][0] * D[stage-1][0] + IntraCoeff[stage-1][1][1]*D[stage-1][1];

        iwd      = NTRANS2[branch[stage]+2][state] + codeword[stage];
        D[stage][0] = TCVQ_CB_SUB3[stage-4][iwd][0] + pred[0];
        D[stage][1] = TCVQ_CB_SUB3[stage-4][iwd][1] + pred[1];
        state    = NTRANS2[branch[stage]][state];
    }

    for(stage = 0; stage < N_STAGE_VQ; stage++)
    {
        for(i = 0; i < N_DIM; i++)
        {
            d_out[(N_DIM*stage) + i] = D[stage][i];
        }
    }
    return;
}

/*---------------------------------------------------------------------------
  * qlsf_ARSN_tcvq_Dec_16k()
  *
  * Predictive BC-TCQ encoder for LSF quantization
  *--------------------------------------------------------------------------*/

short qlsf_ARSN_tcvq_Dec_16k (
    float *y,           /* o  : Quantized LSF vector    */
    short *indice,      /* i  : Indices                 */
    const short nBits         /* i  : number of bits          */
)
{
    short i;
    float error_svq_q[M];
    short safety_net;

    /* Select Mode */
    safety_net = indice[0];


    if(safety_net == 1)
    {
        tcvq_Dec(&indice[1], y, safety_net);

        if(nBits>30)
        {
            for(i = 0; i < 8; i++)
            {
                error_svq_q[i] = AR_SVQ_CB1[indice[10]][i];
                error_svq_q[i+8] = AR_SVQ_CB2[indice[11]][i];
            }

            for(i = 0; i < M; i++)
            {
                y[i] = y[i] + error_svq_q[i] * scale_ARSN[i];
            }
        }
    }
    else
    {
        tcvq_Dec(&indice[1], y, safety_net);

        if(nBits>30)
        {
            for(i = 0; i < 8; i++)
            {
                error_svq_q[i] = AR_SVQ_CB1[indice[10]][i];
                error_svq_q[i+8] = AR_SVQ_CB2[indice[11]][i];
            }

            for(i = 0; i < M; i++)
            {
                y[i] = y[i] + error_svq_q[i];
            }
        }
    }

    return safety_net;
}


/*-------------------------------------------------------------------*
 * lsf_syn_mem_backup()
 *
 * back-up synthesis filter memory and LSF qunatizer memories (used in SC-VBR)
 *-------------------------------------------------------------------*/

void lsf_syn_mem_backup(
    Encoder_State *st,                        /* i: state structure                                       */
    LPD_state* LPDmem,                 /* i: LPD state memory structure                            */
    float *btilt_code,                /* i: tilt code                                             */
    float *bgc_threshold,             /* i:                                                       */
    float *clip_var_bck,              /* o:                                                       */
    short *next_force_sf_bck,           /* 0:                                                       */

    float *lsp_new,                   /* i: LSP vector to quantize                                */
    float *lsp_mid,                   /* i: mid-frame LSP vector                                  */
    float *clip_var,                  /* o: pitch clipping state var                              */
    float *mem_AR,                    /* o: quantizer memory for AR model                         */
    float *mem_MA,                    /* o: quantizer memory for AR model                         */
    float *lsp_new_bck,               /* o: LSP vector to quantize- backup                        */
    float *lsp_mid_bck,               /* o: mid-frame LSP vector - backup                         */
    short *mCb1,                      /* o: counter for stationary frame after a transition frame */
    float *Bin_E,                     /* o: FFT Bin energy 128 *2 sets                            */
    float *Bin_E_old,                 /* o: FFT Bin energy 128 sets                               */
    float *mem_syn_bck,               /* o: synthesis filter memory                               */
    float *mem_w0_bck,                /* o: memory of the weighting filter                        */
    float *streaklimit,
    short *pstreaklen
)
{
    short i;

    *clip_var = st->clip_var[0];

    for(i=0; i<M; i++)
    {
        mem_AR[i] = st->mem_AR[i];
        mem_MA[i] = st->mem_MA[i];
        lsp_new_bck[i] = lsp_new[i];
        lsp_mid_bck[i] = lsp_mid[i];
    }

    *mCb1 = st->mCb1;
    *streaklimit = st->streaklimit;
    *pstreaklen = st->pstreaklen;

    for(i=0; i<L_FFT; i++)
    {
        Bin_E[i]=st->Bin_E[i];
    }

    for(i=0; i<(L_FFT/2); i++)
    {
        Bin_E_old[i]=st->Bin_E_old[i];
    }

    /* back-up memories */
    for(i=0; i<M; i++)
    {
        mem_syn_bck[i] = st->LPDmem.mem_syn[i];
    }

    *mem_w0_bck = st->LPDmem.mem_w0;

    *btilt_code = LPDmem->tilt_code;
    *bgc_threshold = LPDmem->gc_threshold;
    mvr2r( st->clip_var, clip_var_bck, 6 );
    *next_force_sf_bck = st->next_force_safety_net;


    return;
}


/*-------------------------------------------------------------------*
 * lsf_syn_mem_restore()
 *
 * restore synthesis filter memory and LSF quantizer memories
 *-------------------------------------------------------------------*/

void lsf_syn_mem_restore(
    Encoder_State *st,                        /* o: state structure                                        */
    LPD_state* LPDmem,                /* o: LPD_state vewctor                                      */
    float btilt_code,                 /* i:                                                        */
    float gc_threshold,               /* i:                                                        */
    float *clip_var_bck,              /* i:                                                        */
    short next_force_sf_bck,          /* i:                                                        */

    float *lsp_new,                   /* o: LSP vector to quantize                                 */
    float *lsp_mid,                   /* o: mid-frame LSP vector                                   */
    float clip_var,                   /* i: pitch clipping state var                               */
    float *mem_AR,                    /* i: quantizer memory for AR model                          */
    float *mem_MA,                    /* i: quantizer memory for AR model                          */
    float *lsp_new_bck,               /* i: LSP vector to quantize- backup                         */
    float *lsp_mid_bck,               /* i: mid-frame LSP vector - backup                          */
    short mCb1,                       /* i: counter for stationary frame after a transition frame  */
    float *Bin_E,                     /* i: FFT Bin energy 128 *2 sets                             */
    float *Bin_E_old,                 /* i: FFT Bin energy 128 sets                                */
    float *mem_syn_bck,               /* i: synthesis filter memory                                */
    float mem_w0_bck,                 /* i: memory of the weighting filter                         */
    float streaklimit,
    short pstreaklen
)
{
    short i;

    /* restore lsf memories */
    st->clip_var[0] = clip_var;

    for(i=0; i<M; i++)
    {
        st->mem_AR[i] = mem_AR[i];
        st->mem_MA[i] = mem_MA[i];
        lsp_new[i] = lsp_new_bck[i];
        lsp_mid[i] = lsp_mid_bck[i];
    }

    st->mCb1 = mCb1;
    st->streaklimit = streaklimit;
    st->pstreaklen = pstreaklen;

    for(i=0; i<L_FFT; i++)
    {
        st->Bin_E[i] = Bin_E[i];
    }

    for(i=0; i<(L_FFT/2); i++)
    {
        st->Bin_E_old[i]=Bin_E_old[i];
    }

    /* restoring memories */
    st->LPDmem.mem_w0 = mem_w0_bck;

    for(i=0; i<M; i++)
    {
        st->LPDmem.mem_syn[i] = mem_syn_bck[i];
    }

    LPDmem->tilt_code = btilt_code;
    LPDmem->gc_threshold = gc_threshold;
    mvr2r( clip_var_bck, st->clip_var, 6 );
    st->next_force_safety_net = next_force_sf_bck;


    return;
}

/*--------------------------------------------------------------------------*
 * lsf_update_memory()
 *
 *
 *--------------------------------------------------------------------------*/

void lsf_update_memory(
    int narrowband,           /* i  : narrowband flag                 */
    const float qlsf[],       /* i  : quantized lsf coefficients      */
    float old_mem_MA[],       /* i  : MA memory                       */
    float mem_MA[]            /* o  : updated MA memory               */
)
{
    int i;

    for (i=0; i<M; ++i)
    {
        mem_MA[i] = qlsf[i] - lsf_means[narrowband][i] - MU_MA * old_mem_MA[i];
    }

    return;
}

/*--------------------------------------------------------------------------*
 * tcxlpc_get_cdk()
 *
 *
 *--------------------------------------------------------------------------*/

int tcxlpc_get_cdk(         /* Returns: codebook index */
    int coder_type            /* (I) GC/VC indicator   */
)
{
    int cdk;

    cdk = (coder_type == VOICED);

    return cdk;
}

/*--------------------------------------------------------------------------*
 * msvq_dec()
 *
 *
 *--------------------------------------------------------------------------*/


void msvq_dec
(
    const float *const *cb,   /* i  : Codebook (indexed cb[*stages][levels][p])            */
    const int dims[],         /* i  : Dimension of each codebook stage (NULL: full dim.)   */
    const int offs[],         /* i  : Starting dimension of each codebook stage (NULL: 0)  */
    int stages,               /* i  : Number of stages                                     */
    int N,                    /* i  : Vector dimension                                     */
    int maxN,                 /* i  : Codebook dimension                                   */
    const int Idx[],          /* i  : Indices                                              */
    float *uq,                /* o  : quantized vector                                     */
    Word16 *uq_ind            /* o  : quantized vector (fixed point)                       */
)
{
    int i;
    int n, maxn, start;
    Word16 j;
    set_zero( uq, N );

    IF( uq_ind )
    {
        FOR (i=0; i<N; ++i)
        {
            uq_ind[i] = 0;
            move16();
        }
    }
    for( i=0; i<stages; i++ )
    {
        if (dims)
        {
            n = dims[i];
            maxn = n;
        }
        else
        {
            n = N;
            maxn = maxN;
        }
        if (offs)
        {
            start = offs[i];
        }
        else
        {
            start = 0;
        }

        v_add( uq+start, cb[i]+Idx[i]*maxn, uq+start, n );

        IF( uq_ind )
        {
            FOR (j=0; j<n; ++j)
            {
                move16();
                uq_ind[start+j] = add(uq_ind[start+j], (Word16)(cb[i][Idx[i]*maxn+j]*2.0f*1.28f));
            }
        }
    }

    return;
}



/*--------------------------------------------------------------------------*
 * lsf_update_memory()
 *
 *
 *--------------------------------------------------------------------------*/

static void spec2isf(
    float/*double*/	spec_r[],   /* input spectrum real part (only left half + one zero)*/
    float/*double*/	spec_i[],   /* input spectrum imag part (only left half+right halt with zeros)*/
    int/*short*/   speclen,    /* length of spectrum (only left half)*/
    float /*double*/	lsf[],       /* locations of LSFs (buffer must be sufficiently long) */ /*15Q16*/
    const float /*double*/	old_lsf[]       /* locations of LSFs (buffer must be sufficiently long) */ /*15Q16*/
)
{

    /*spec_r[] needs a 0 in the end!*/
    float s;
    int specix,lsfix, i;

    specix = lsfix = 0;
    s = spec_r[specix++];

    while ((specix < speclen) && lsfix <= 15)
    {

        /*check for next zero crossing*/
        for (; s*spec_r[specix] >= 0; specix++);

        lsf[lsfix++] = (specix-1 + spec_r[specix-1]/(spec_r[specix-1]-spec_r[specix]))*(12800/256);

        /*check for the next zero crossing*/
        for (; s*spec_i[specix] >= 0; specix++);

        lsf[lsfix++] = (specix-1 + spec_i[specix-1]/(spec_i[specix-1]-spec_i[specix]))*(12800/256);

        spec_r[speclen] = s;
        spec_i[speclen] = s;

        s =-s;
    }

    if (lsfix < 16)
    {
        for(i=0; i<16; i++)
        {
            lsf[i] = old_lsf[i];
        }
    }

    return;
}

/*--------------------------------------------------------------------------*
 * a2isf()
 *
 *
 *--------------------------------------------------------------------------*/

#define SCALE1_HALF  1018.59161376953f

typedef struct
{
    float re;
    float im;
} Pfloat;

void a2isf(
    float *a,
    float *isf,
    const float *old_isf,
    short lpcOrder
)
{
    float RealFFT[128];
    float ImagFFT[128];
    float RealOut[130];
    float ImagOut[130];
    float *ptrReal;
    float *ptrImag;
    int n, i, j;
    const Pfloat *ptwiddle;
    Pfloat *pwn17, pwn[128], *pwn15, tmpw15;
    int N = 256;
    float s[4];
    float L_tmp, L_tmp1;
    float lpc[19], fftTmpRe[16], fftTmpIm[16];
    Pfloat twid[128];
    float c;

    set_zero(fftTmpRe, 16);
    set_zero(fftTmpIm, 16);

    /* half length FFT */
    /*c = [sum(a) ((-1).^(1:length(a)))*a];*/

    L_tmp = 0;
    for(j=0; j<=lpcOrder; j++)
    {
        L_tmp += a[j];
    }

    L_tmp1 = 0;
    for(j=0; j<lpcOrder/2; j++)
    {
        L_tmp1 -= a[2*j];
        L_tmp1 += a[2*j+1];
    }
    L_tmp1 -= a[2*j];



    /*s = [1 -2*(c(1)-c(2))/(c(1)+c(2)) 1]';*/
    s[0]  = 1;
    if((L_tmp+L_tmp1) != 0)
    {
        s[1]  = -2*((L_tmp-L_tmp1)/(L_tmp+L_tmp1));
    }
    else
    {
        s[1] = 1;
    }
    s[2]  = 1;

    lpc[0] = a[0] * s[0];
    L_tmp =  a[1] * s[0];
    lpc[1] = L_tmp + (a[1-1] * s[1]);


    for (n = 2; n < 17; n++)
    {
        L_tmp = a[n] * s[0];
        L_tmp += (a[n - 1] * s[1]);
        lpc[n] = L_tmp + (a[n - 2] * s[2]);
    }
    lpc[18] = a[16] * s[0];
    L_tmp =   a[15] * s[0];
    lpc[17] = L_tmp + (a[16] * s[1]);


    ptrReal = RealFFT;
    ptrImag = ImagFFT;

    for(j=0; j<9; j++)
    {
        fftTmpRe[j] = lpc[2*j];
        fftTmpIm[j] = lpc[2*j+1];
    }
    fftTmpRe[j] = lpc[2*j];
    fftTmpIm[j] = 0;
    j++;

    for(; j<16; j++)
    {
        fftTmpRe[j] = 0;
        fftTmpIm[j] = 0;
    }

    DoRTFTn(fftTmpRe, fftTmpIm, 16);

    for(j=0; j<16; j++)
    {
        ptrReal[j*8] = fftTmpRe[j];
        ptrImag[j*8] = fftTmpIm[j];
    }

    ptrReal++;
    ptrImag++;

    for(i=1; i<8; i++)
    {
        /*X=x(i:8:M/8) .* exp(-j*2*pi*i*(0:M/8-1)/M);*/
        ptwiddle = (const Pfloat *)w_a[i-1];

        fftTmpRe[0] = lpc[0];
        fftTmpIm[0] = lpc[1];

        for(j=1; j<9; j++)
        {
            fftTmpRe[j] = (lpc[2*j] * ptwiddle->re) - (lpc[2*j+1] * ptwiddle->im);
            fftTmpIm[j] = (lpc[2*j+1] * ptwiddle->re) + (lpc[2*j] * ptwiddle->im);
            ptwiddle++;
        }

        fftTmpRe[j] = lpc[2*j]*ptwiddle->re;
        fftTmpIm[j] = lpc[2*j]*ptwiddle->im;
        ptwiddle++;
        j++;
        for(; j<16; j++)
        {
            fftTmpRe[j] = 0;
            fftTmpIm[j] = 0;
            ptwiddle++;
        }

        DoRTFTn(fftTmpRe, fftTmpIm, 16);

        for(j=0; j<16; j++)
        {
            ptrReal[j*8] = fftTmpRe[j];
            ptrImag[j*8] = fftTmpIm[j];
        }

        ptrReal++;
        ptrImag++;

    }

    c = EVS_PI / ( 2.0f * (float)128 );

    for ( i = 1 ; i < 128 ; i++ )
    {
        twid[i-1].im = (float)sin( c * ( 2.0f * (float)i ) );
        twid[i-1].re = (float)cos( c * ( 2.0f * (float)i ) );
    }
    ptwiddle = twid;

    /* pre-twiddle */
    for ( i = 1 ; i < 128 ; i++ )
    {
        pwn[i-1].im = (float)sin( c * ( 18.0f * (float)i ) );
        pwn[i-1].re = (float)cos( c * ( 18.0f * (float)i ) );
    }

    pwn17 = pwn;
    pwn15 = &tmpw15;

    /*Y(1) = real(X(1)) + imag(X(1));*/
    RealOut[0] = (RealFFT[0] + ImagFFT[0]);
    ImagOut[0]  = 0;

    /*Y(N/2+1) = 0.5*(X(1) + conj(X(1))).*exp(pi*i*128*(18)/N) - i*0.5*(X(1) - conj(X(1))).*exp(pi*i*128*(16)/N);*/
    RealOut[128] = 0;
    ImagOut[128] = (RealFFT[0] + RealFFT[0]) - (ImagFFT[0] + ImagFFT[0]);

    /*Y(2:N/2) = (0.5*(X(2:N/2) + conj(X(N/2:-1:2))) - i*0.5*(X(2:N/2) - conj(X(N/2:-1:2))).*exp(-pi*i*r*(2)/N)).*exp(pi*i*r*(18)/N);*/
    for(i=1; i<N/2; i++)
    {
        float ReAr = (RealFFT[i] + RealFFT[N/2-i]);
        float ReBr = (RealFFT[N/2-i] - RealFFT[i]);
        float ImAr = (ImagFFT[i] - ImagFFT[N/2-i]);
        float ImBr =-(ImagFFT[i] + ImagFFT[N/2-i]);

        tmpw15.re = (ptwiddle->re * pwn17->re) + (ptwiddle->im * pwn17->im);
        tmpw15.im = (ptwiddle->re * pwn17->im) - (ptwiddle->im * pwn17->re);

        /*RealOut[i] = mac_r(L_msu(L_msu(L_mult(ReAr, pwn17->re),ImAr, pwn17->im), ReBr, pwn15->v.im), ImBr, pwn15->re); move16();*/
        RealOut[i] = (ReAr * pwn17->re) - (ImAr * pwn17->im) - ((ReBr * pwn15->im) + (ImBr * pwn15->re));
        ImagOut[i] = (ReAr * pwn17->im) + (ImAr * pwn17->re) +  (ReBr * pwn15->re) - (ImBr * pwn15->im);

        ptwiddle++;
        pwn17++;
    }

    spec2isf(RealOut, ImagOut, 128, isf, old_isf);

    isf[lpcOrder - 1] = (Float32)(acos(a[lpcOrder]) * SCALE1_HALF);

    return;
}

