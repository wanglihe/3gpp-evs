/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <assert.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*-----------------------------------------------------------------*
 * Local constants
 *-----------------------------------------------------------------*/

#define MAX_LEN_LP   960
#define SCALE1_LPC   2037.1832275390625f /* LSP to LSF conversion factor */


/*---------------------------------------------------------------------*
 * autocorr()
 *
 * Compute autocorrelations of input signal
 *---------------------------------------------------------------------*/

void autocorr(
    const float *x,        /* i  : input signal               */
    float *r,        /* o  : autocorrelations vector    */
    const short m,         /* i  : order of LP filter         */
    const short len,       /* i  : window size                */
    const float *wind,     /* i  : window                     */
    const short rev_flag,  /* i  : flag to reverse window     */
    const short sym_flag,  /* i  : symmetric window flag      */
    const short no_thr     /* i  : flag to avoid thresholding */
)
{
    float t[MAX_LEN_LP];
    float s;
    short i, j;

    /* Windowing of signal */
    if (rev_flag == 1)
    {
        /* time reversed window */
        for (i = 0; i < len; i++)
        {
            t[i] = x[i] * wind[len-i-1];
        }
    }
    else if( sym_flag == 1 )
    {
        /* symmetric window of even length */
        for( i=0; i<len/2; i++ )
        {
            t[i] = x[i] * wind[i];
        }

        for( ; i<len; i++ )
        {
            t[i] = x[i] * wind[len-1-i];
        }
    }
    else  /* assymetric window */
    {
        for (i = 0; i < len; i++)
        {
            t[i] = x[i] * wind[i];
        }
    }

    /* Compute r[1] to r[m] */
    for (i = 0; i <= m; i++)
    {
        s = t[0] * t[i];
        for (j = 1; j < len-i; j++)
        {
            s += t[j] * t[i+j];
        }
        r[i] = s;
    }

    if ( r[0] < 100.0f && no_thr == 0 )
    {
        r[0] = 100.0f;
    }

    return;
}

/*---------------------------------------------------------------------*
 * lev_dur()
 *
 * Wiener-Levinson-Durbin algorithm to compute LP parameters from the autocorrelations
 * of input signal
 *---------------------------------------------------------------------*/

short lev_dur(          /* o:   energy of prediction error   */
    float       *a,     /* o:   LP coefficients (a[0] = 1.0) */
    const float *r,     /* i:   vector of autocorrelations   */
    const short m,      /* i:   order of LP filter           */
    float       epsP[]  /* o:   prediction error energy      */
)
{
    short i, j, l;
    float buf[TCXLTP_LTP_ORDER];
    float *rc;    /* reflection coefficients  0,...,m-1 */
    float s, at, err;
    short flag=0;

    rc = &buf[0];
    rc[0] = (-r[1])/r[0];
    a[0] = 1.0f;
    a[1] = rc[0];
    err = r[0] + r[1]*rc[0];
    if ( epsP != NULL)
    {
        epsP[0] = r[0];
        epsP[1] = err;
    }

    for ( i = 2; i <= m; i++ )
    {
        s = 0.0f;
        for ( j = 0; j < i; j++ )
        {
            s += r[i-j] * a[j];
        }

        rc[i-1]= (-s) / err;

        if (fabs(rc[i-1]) > 0.99945f)
        {
            flag=1;/* Test for unstable filter. If unstable keep old A(z) */
        }

        for ( j = 1; j <= i/2; j++ )
        {
            l = i-j;
            at = a[j] + rc[i-1] * a[l];
            a[l] += rc[i-1] * a[j];
            a[j] = at;
        }

        a[i] = rc[i-1];

        err += rc[i-1] * s;
        if ( err <= 0.0f )
        {
            err = 0.01f;
        }

        if ( epsP != NULL)
        {
            epsP[i] = err;
        }
    }

    return (flag);
}


/*---------------------------------------------------------------------*
 * E_LPC_int_lpc_tcx()
 *
 *
 *---------------------------------------------------------------------*/

void E_LPC_int_lpc_tcx(
    const float lsf_old[],    /* input : LSFs from past frame            */
    const float lsf_new[],    /* input : LSFs from present frame         */
    float a[]                 /* output: interpolated LP coefficients    */
)
{
    int i;
    float lsf[M];

    for (i = 0; i < M; i++)
    {
        lsf[i] = lsf_old[i]*0.125f + lsf_new[i]*0.875f;
    }

    lsp2a_stab(lsf, a, M);

    return;
}

/*---------------------------------------------------------------------*
 * lsp_reorder()
 *
 *
 *---------------------------------------------------------------------*/

static void lsp_reorder(
    float *lsp,                /* (I/O): lsp vector (acos() domain) */
    float min_dist,            /* (I): minimum required distance    */
    int lpcorder               /* (I): LPC order                    */
)
{
    short i;
    float lsp_min, lsp_max;

    /* Verify the LSF ordering and minimum GAP */
    lsp_min = min_dist;

    for (i=0; i<lpcorder; ++i)
    {
        if (lsp[i] < lsp_min)
        {
            lsp[i] = lsp_min;
        }
        lsp_min = lsp[i] + min_dist;
    }

    /* Reverify the LSF ordering and minimum GAP in the reverse order (security) */
    lsp_max = EVS_PI - min_dist;

    /* If danger of unstable filter in case of resonance in HF */
    if (lsp[lpcorder-1] > lsp_max)
    {
        /* Reverify the minimum LSF gap in the reverse sense */
        for (i = lpcorder-1; i>=0; --i)
        {
            if (lsp[i] > lsp_max)
            {
                lsp[i] = lsp_max;
            }
            lsp_max = lsp[i] - min_dist;
        }
    }

    return;
}


/*---------------------------------------------------------------------*
 * E_LPC_lsp_unweight()
 *
 * Approximate unweighting
 *---------------------------------------------------------------------*/

int E_LPC_lsp_unweight(
    float lsp_w[],             /* (I): weighted lsp             */
    float lsp_uw[],            /* (O): unweighted lsp           */
    float lsf_uw[],            /* (O): unweighted lsf           */
    float inv_gamma            /* (I): inverse weighting factor */
)
{
    float lsp_w_orig[M], lsp_w_diff[M], mean, step;
    const lsp_unw_triplet *unw_coeffs = NULL;
    short i;

    /* Table selection */
    if ((float)fabs(inv_gamma - 1.0f / 0.94f) < 0.0001f)
    {
        unw_coeffs = p16_gamma0_94to1;
    }
    else if ((float)fabs(inv_gamma - 1.0f / 0.92f) < 0.0001f)
    {
        unw_coeffs = p16_gamma0_92to1;
    }
    else
    {
        assert(0);
    }

    step = EVS_PI/(float)(M+1);
    mean = step;

    /* Apply acos() and get mean removed version */
    for (i=0; i<M; ++i)
    {
        lsp_w_orig[i] = (float)acos(lsp_w[i]);
        lsp_w_diff[i] = lsp_w_orig[i] - mean;
        mean += step;
    }

    /* Approximate unweighting by 3-tap FIR */

    lsp_uw[0] = lsp_w_orig[0] + unw_coeffs[0][1] * lsp_w_diff[0] + unw_coeffs[0][2] * lsp_w_diff[1];
    for (i=1; i<M-1; ++i)
    {
        lsp_uw[i] = lsp_w_orig[i] + unw_coeffs[i][0] * lsp_w_diff[i-1] + unw_coeffs[i][1] * lsp_w_diff[i] + unw_coeffs[i][2] * lsp_w_diff[i+1];
    }

    lsp_uw[M-1] = lsp_w_orig[M-1] + unw_coeffs[M-1][0] * lsp_w_diff[M-2] + unw_coeffs[M-1][1] * lsp_w_diff[M-1];

    /* Reorder */
    lsp_reorder( lsp_uw, 50.0f/SCALE1_LPC, M );

    /* Convert to lsf, apply cos() */
    for( i=0; i<M; ++i )
    {
        lsf_uw[i] = lsp_uw[i] * SCALE1_LPC;
        lsp_uw[i] = (float)cos(lsp_uw[i]);
    }

    return 0;
}

