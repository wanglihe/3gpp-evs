/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "basop_proto_func.h"
#include "control.h"
#include "basop_util.h"

#define NC_MAX 8

static Word16 E_LPC_f_lsp_pol_get(const Word16 lsp[], Word32 f[],  const Word16 n, const Word16 past_Ovf, const Word16 isMODE1);


/*
 * E_LPC_f_lsp_a_conversion
 *
 * Parameters:
 *    lsp            I: Line spectral pairs          Q15
 *    a              O: Predictor coefficients (order = m)  Qx (The Q factor of the output to be deduced from a(0))
 *    m              I: order of LP filter
 *
 * Function:
 *    Convert ISPs to predictor coefficients a[]
 *
 * Returns:
 *    void
 */
void basop_E_LPC_f_lsp_a_conversion(const Word16 *lsp, Word16 *a,  const Word16 m)
{
    Word16 i, j, k;
    Word32 f1[NC_MAX+1], f2[NC_MAX+1];
    Word16 nc;
    Word32 t0;
    Word16 Ovf, Ovf2;


    /*-----------------------------------------------------*
     *  Find the polynomials F1(z) and F2(z)               *
     *-----------------------------------------------------*/

    nc = shr(m, 1);

    assert(m == 16 || m == 10);

    Ovf = 0;
    move16();
    Ovf = E_LPC_f_lsp_pol_get(&lsp[0], f1, nc, Ovf, 1);
    Ovf2 = E_LPC_f_lsp_pol_get(&lsp[1], f2, nc, Ovf, 1);
    IF(sub(Ovf2,Ovf) !=0)
    {
        /* to ensure similar scaling for f1 and f2 in case
          an overflow would be detected only in f2,
          but this case never happen on my dtb */
        E_LPC_f_lsp_pol_get(&lsp[0], f1, nc, s_max(Ovf2,Ovf), 1);
    }
    /*-----------------------------------------------------*
     *  Multiply F1(z) by (1+z^-1) and F2(z) by (1-z^-1)   *
     *-----------------------------------------------------*/
    /*modification*/
    k = sub(nc,1);
    FOR (i = 0; i <= k; i++)
    {
        f1[nc-i] = L_add(f1[nc-i],f1[nc-i-1]);
        move32();
        f2[nc-i] = L_sub(f2[nc-i],f2[nc-i-1]);
        move32();
    }

    /*-----------------------------------------------------*
     *  A(z) = (F1(z)+F2(z))/2                             *
     *  F1(z) is symmetric and F2(z) is antisymmetric      *
     *-----------------------------------------------------*/

    t0 = L_deposit_l(0);
    FOR (i = 1; i <= nc; i++)
    {
        t0 = L_max( t0, L_abs(L_add(f1[i], f2[i])) );
        t0 = L_max( t0, L_abs(L_sub(f1[i], f2[i])) );
    }
    k = s_min( norm_l(t0), 6 );
    a[0] = shl( 256, k );
    move16();
    test();
    IF( Ovf || Ovf2)
    {
        a[0] = shl( 256, sub(k,Ovf) );
        move16();
    }
    j = m;
    FOR (i = 1;  i <= nc; i++)
    {
        /* a[i] = 0.5*(f1[i] + f2[i]) */
        t0 = L_add(f1[i],f2[i]);
        t0 = L_shl(t0, k);
        a[i] = round_fx(t0);             /* from Q23 to Qx and * 0.5 */

        /* a[j] = 0.5*(f1[i] - f2[i]) */
        t0 = L_sub(f1[i],f2[i]);
        t0 = L_shl(t0, k);
        a[j] = round_fx(t0);             /* from Q23 to Qx and * 0.5 */
        j--;
    }

    return;
}


/*---------------------------------------------------------------------------
* procedure  reorder_lsf()
*
* To make sure that the  lsfs are properly ordered and to keep a certain
* minimum distance between consecutive lsfs.
*--------------------------------------------------------------------------*/
void basop_reorder_lsf(
    Word16 *lsf,      /* i/o: LSFs in the frequency domain (0..0.5)   Q(x2.56)*/
    const Word16 min_dist,  /* i  : minimum required distance               x2.56*/
    const Word16 n,         /* i  : LPC order                               */
    const Word32 fs         /* i  : sampling frequency                      */
)
{
    Word16 i, lsf_min, n_m_1;
    Word16 lsf_max;

    lsf_min = min_dist;
    move16();

    /*-----------------------------------------------------------------------*
     * Verify the LSF ordering and minimum GAP
     *-----------------------------------------------------------------------*/

    FOR (i = 0; i < n; i++)
    {
        if (sub(lsf[i], lsf_min) < 0)
        {
            lsf[i] = lsf_min;
            move16();
        }
        lsf_min = add(lsf[i], min_dist);
    }

    /*-----------------------------------------------------------------------*
     * Reverify the LSF ordering and minimum GAP in the reverse order (security)
     *-----------------------------------------------------------------------*/
    lsf_max = round_fx(L_sub(L_shr(L_mult0(extract_l(L_shr(fs,1)), 1311),9-16), L_deposit_h(min_dist))); /* Q0 + Q9 , 1311 is 2.56 in Q9 */
    n_m_1 = sub(n,1);
    IF (sub(lsf[n_m_1], lsf_max) > 0)    /* If danger of unstable filter in case of resonance in HF */
    {
        FOR (i = n_m_1; i >= 0; i--) /* Reverify the minimum LSF gap in the reverse direction */
        {
            if (sub(lsf[i], lsf_max) > 0)
            {
                lsf[i] = lsf_max;
                move16();
            }
            lsf_max = sub(lsf[i], min_dist);
        }
    }
}


/*
 * E_LPC_f_lsp_pol_get
 *
 * Parameters:
 *    lsp/isp        I: Line spectral pairs (cosine domaine)      Q15
 *    f              O: the coefficients of F1 or F2        Q23
 *    n              I: no of coefficients (m/2)
 *                      == NC for F1(z); == NC-1 for F2(z)
 *    fact           I: scaling factor
 *
 *-----------------------------------------------------------*
 * procedure E_LPC_f_lsp_pol_get:                            *
 *           ~~~~~~~~~~~                                     *
 *   Find the polynomial F1(z) or F2(z) from the LSPs.       *
 * This is performed by expanding the product polynomials:   *
 *                                                           *
 * F1(z) =   product   ( 1 - 2 LSF_i z^-1 + z^-2 )           *
 *         i=0,2,4,6,8                                       *
 * F2(z) =   product   ( 1 - 2 LSF_i z^-1 + z^-2 )           *
 *         i=1,3,5,7,9                                       *
 *                                                           *
 * where LSP_i are the LSPs in the cosine domain.            *
 *                                                           *
 *-----------------------------------------------------------*
 *   R.A.Salami    October 1990                              *
 *-----------------------------------------------------------*
 */
static
Word16 E_LPC_f_lsp_pol_get(const Word16 lsp[], Word32 f[],  const Word16 n, const Word16 past_Ovf, const Word16 isMODE1)
{
    /* All computation in Q23 */
    const Word16 *plsp;
    Word16 i, j;
    Word16 b;
    Word32 b32;
    Word16 Ovf = 0;
    Word16 Q_out;
    Word16 m2;


    Q_out = 31-23;
    move16();
    Ovf = past_Ovf;
    move16();

    test();
    if(past_Ovf && isMODE1) /* Currently this feature is implemented only in MODE1 */
    {
        /* In some NB cases, overflow where detectected
            in f1 or f2 polynomial computation when it
            happen we reduce the precision of the computing
            to limit the risk of saturation*/
        Q_out = add(Q_out, past_Ovf);
    }
    Overflow = 0;
    move16();
    plsp = lsp;
    f[0] = L_shl(1, sub(31, Q_out));
    move32();
    /*b = -2.0f * *plsp;*/
    b = *plsp;
    move16();
    m2 = shl(-2, sub(15, Q_out));
    f[1] = L_mult(b, m2);
    move32();

    FOR (i = 2; i <= n; i++)
    {
        plsp += 2;
        /*b = 2.0f * *plsp;*/
        move16();
        b = *plsp;
        b32 = L_mult(b, m2);

        /*f[i] = -b*f[i-1] + 2.0f*f[i-2];*/
        move32();
        f[i] = L_shl(L_sub(f[i-2], Mpy_32_16(f[i-1], b)),1);

        FOR (j = i-1; j > 1; j--)
        {
            /*f[j] += b*f[j-1] + f[j-2];*/
            move32();
            f[j] = L_add(f[j], L_sub(f[j-2], L_shl(Mpy_32_16(f[j-1], b),1)));
        }
        move32();
        f[1] = L_add(f[1], b32);
    }


    test();
    IF (Overflow>0 && isMODE1)
    {
        /* If an overflow is detected, redo the computation with 1 bit less */
        Ovf = add(Ovf,1);
        Ovf = E_LPC_f_lsp_pol_get(lsp, f, n ,Ovf, isMODE1);
    }
    return Ovf;
}
