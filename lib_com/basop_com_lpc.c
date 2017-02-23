/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/


#include <math.h>
#include <memory.h>
#include <assert.h>
#include <stdio.h>
#include "typedef.h"
#include "basop_proto_func.h"
#include "cnst.h"

#include "basop_util.h"
#include "stl.h"

#define UNROLL_CHEBYSHEV_INNER_LOOP
#define NC_MAX 8
#define GUESS_TBL_SZ 256

#define Madd_32_16(accu, x, y) L_add(accu, Mpy_32_16(x, y))
#define Msub_32_16(accu, x, y) L_sub(accu, Mpy_32_16(x, y))


/*
 * weight_a
 *
 * Parameters:
 *    a              I: LP filter coefficients			Q12
 *    ap             O: weighted LP filter coefficients Q12
 *    gamma          I: weighting factor				Q15
 *
 * Function:
 *    Weighting of LP filter coefficients, ap[i] = a[i] * (gamma^i).
 *
 * Returns:
 *    void
 */
void basop_weight_a(const Word16 *a, Word16 *ap, const Word16 gamma)
{
    Word16 i, fac;
    Word32 Amax;
    Word16 shift;



    fac = gamma;
    Amax = L_mult( 16384, a[0] );
    FOR (i = 1; i < M; i++)
    {
        Amax = L_max( Amax, L_abs( L_mult0( fac, a[i] ) ) );
        fac = mult_r( fac, gamma );
    }
    Amax = L_max( Amax, L_abs( L_mult0( fac, a[M] ) ) );
    shift = norm_l( Amax );
    fac = gamma;
    ap[0] = shl( a[0], sub(shift,1) );
    move16();
    FOR (i = 1; i < M; i++)
    {
        ap[i] = round_fx(L_shl(L_mult0(a[i], fac),shift));
        move16();
        fac = mult_r( fac, gamma );
    }
    ap[M] = round_fx(L_shl(L_mult0(a[M], fac),shift));
    move16();


    return;
}

/*
 * weight_a_inv
 *
 * Parameters:
 *    a              I: LP filter coefficients			Q12
 *    ap             O: weighted LP filter coefficients Q12
 *    inv_gamma      I: inverse weighting factor				Q14
 *
 * Function:
 *    Weighting of LP filter coefficients, ap[i] = a[i] * (inv_gamma^i).
 *
 * Returns:
 *    void
 */
void basop_weight_a_inv(const Word16 *a, Word16 *ap, const Word16 inv_gamma)
{
    Word16 i;
    static const Word16 inv_gamma_tab_12k8[16] = { 17809, 19357, 21041, 22870, 24859, 27020, 29370, 31924,   /* Q14 */
                                                   17350, 18859, 20499, 22281, 24219, 26325, 28614, 31102
                                                 }; /* Q13 */
    static const Word16 inv_gamma_tab_16k[16]  = { 17430, 18542, 19726, 20985, 22324, 23749, 25265, 26878,   /* Q14 */
                                                   14297, 15209, 16180, 17213, 18312, 19480, 20724, 22047
                                                 }; /* Q13 */
    const Word16 *inv_gamma_tab;
    Word32 L_tmp;
    Word32 Amax;
    Word16 shift;


    IF (inv_gamma == 16384)
    {
        FOR (i = 0; i <= M; i++)
        {
            ap[i] = a[i];
            move16();
        }
        return;
    }

    assert( inv_gamma==GAMMA1_INV || inv_gamma==GAMMA16k_INV );

    inv_gamma_tab = inv_gamma_tab_12k8;
    move16();
    if (sub(inv_gamma,GAMMA16k_INV) == 0)
    {
        inv_gamma_tab = inv_gamma_tab_16k;
        move16();
    }

    Amax = L_mult( 16384, a[0] );
    FOR (i = 1; i < 9; i++)
    {
        Amax = L_max( Amax, L_abs( L_mult( a[i], inv_gamma_tab[i-1] ) ) );
    }
    FOR (i = 9; i < 17; i++)
    {
        Amax = L_max( Amax, L_abs( L_shl( L_mult( a[i], inv_gamma_tab[i-1] ), 1 ) ) );
    }
    shift = norm_l( Amax );
    ap[0] = shl( a[0], sub(shift,1) );
    move16();
    FOR (i = 1; i < 9; i++)
    {
        L_tmp = L_mult( a[i], inv_gamma_tab[i-1] );
        ap[i] = round_fx( L_shl( L_tmp, shift ) );
        move16();
    }
    shift = add(shift,1);
    FOR (i = 9; i < 17; i++)
    {
        L_tmp = L_mult( a[i], inv_gamma_tab[i-1] );
        ap[i] = round_fx( L_shl( L_tmp, shift ) );
        move16();
    }


    return;
}

/*
 * basop_E_LPC_a_add_tilt
 *
 * Parameters:
 *    a              I: LP filter coefficients (m+1 coeffs)
 *    ap             O: modified LP filter coefficients (m+2 coeffs)
 *    gamma          I: tilt factor
 *
 * Function:
 *    Modified LP filter by adding 1st order pre-premphasis, Ap(z)=A(z).(1-gamma.z^(-1))
 *
 * Returns:
 *    void
 */
void basop_E_LPC_a_add_tilt(const Word16 *a, Word16 *ap, Word16 gamma)
{
    Word16 i;
    Word32 Amax, Atmp[M+2];
    Word16 shift;




    Amax = L_mult( 16384, a[0] );
    FOR (i = 1; i <= M; i++)
    {
        Atmp[i] = L_sub( L_mult(16384, a[i]), L_mult0(gamma, a[i-1]) );
        move32();
        Amax = L_max( Amax, L_abs( Atmp[i] ) );
    }
    Atmp[M+1] = L_negate( L_mult0(gamma, a[M]) );
    move32();
    Amax = L_max( Amax, L_abs( Atmp[M+1] ) );
    shift = norm_l( Amax );
    ap[0] = shl( a[0], sub(shift,1) );
    move16();
    FOR (i = 1; i <= M; i++)
    {
        ap[i] = round_fx( L_shl( Atmp[i], shift ) );
        move16();
    }
    ap[M+1] = round_fx( L_shl( Atmp[M+1], shift ) );
    move16();

    return;
}



static Word16 xsf_to_xsp(Word16 xsf)
{
    /* xsp = cos(xsf * 3.1415/6400); */
    return getCosWord16R2(xsf);
}

/*
 * lsf2lsp
 *
 * Parameters:
 *    lsf            I: lsf[m] normalized (range: 0 <= val <= 0.5)	x2.56
 *    lsp            O: lsp[m] (range: -1 <= val < 1)				Q15
 *
 * Function:
 *    Transformation lsf to lsp
 *
 *    LSF are line spectral pair in frequency domain (0 to 6400).
 *    LSP are line spectral pair in cosine domain (-1 to 1).
 *
 * Returns:
 *    void
 */
void basop_lsf2lsp(const Word16 lsf[], Word16 lsp[])
{
    Word16 i;



    /* convert ISFs to the cosine domain */
    FOR (i = 0; i < M; i++)
    {
        *lsp++ = xsf_to_xsp(*lsf++);
        move16();
    }


    return;
}
