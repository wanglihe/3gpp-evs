/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "options.h"
#include "cnst.h"
#include "basop_proto_func.h"
#include "stl.h"
#include "prot.h"
#include "rom_com.h"

/* compare two positive normalized 16 bit mantissa/exponent values */
/* return value: positive if first value greater, negative if second value greater, zero if equal */
static Word16 compMantExp16Unorm(Word16 m1, Word16 e1, Word16 m2, Word16 e2)
{
    Word16 tmp;

    assert((m1 >= 0x4000) && (m2 >= 0x4000)); /* comparisons below work only for normalized mantissas */

    tmp = sub(e1, e2);
    if (tmp == 0) tmp = sub(m1, m2);

    return tmp;
}

void basop_lpc2mdct(Word16 *lpcCoeffs, Word16 lpcOrder,
                    Word16 *mdct_gains, Word16 *mdct_gains_exp,
                    Word16 *mdct_inv_gains, Word16 *mdct_inv_gains_exp)
{
    Word32 RealData[FDNS_NPTS];
    Word32 ImagData[FDNS_NPTS];
    Word16 i, j, k, step, scale, s, tmp16;
    Word16 g, g_e, ig, ig_e;
    Word32 tmp32;
    const PWord16 *ptwiddle;



    /* short-cut, to avoid calling of BASOP_getTables() */
    ptwiddle = SineTable512_fx;
    step = 8;

    /*ODFT*/
    assert(lpcOrder < FDNS_NPTS);

    /* pre-twiddle */
    FOR (i=0; i<=lpcOrder; i++)
    {
        RealData[i] = L_mult(lpcCoeffs[i], ptwiddle->v.re);
        move32();
        ImagData[i] = L_negate(L_mult(lpcCoeffs[i], ptwiddle->v.im));
        move32();
        ptwiddle += step;
    }

    /* zero padding */
    FOR ( ; i<FDNS_NPTS; i++)
    {
        RealData[i] = 0;
        move32();
        ImagData[i] = 0;
        move32();
    }

    /* half length FFT */
    scale = add(norm_s(lpcCoeffs[0]),1);
    move16();
    BASOP_cfft( RealData, ImagData, 1, &scale );    /* sizeOfFft == FDNS_NPTS == 8 */


    /*Get amplitude*/
    j = FDNS_NPTS - 1;
    k = 0;
    move16();

    FOR (i=0; i<FDNS_NPTS/2; i++)
    {
        s = sub(norm_l(L_max(L_abs(RealData[i]), L_abs(ImagData[i]))), 1);

        tmp16 = extract_h(L_shl(RealData[i], s));
        tmp32 = L_mult(tmp16, tmp16);

        tmp16 = extract_h(L_shl(ImagData[i], s));
        tmp16 = mac_r(tmp32, tmp16, tmp16);

        s = shl(sub(scale, s), 1);

        if (tmp16 == 0)
        {
            s = -16;
            move16();
        }
        if (tmp16 == 0)
        {
            tmp16 = 1;
            move16();
        }

        BASOP_Util_Sqrt_InvSqrt_MantExp(tmp16, s, &g, &g_e, &ig, &ig_e);

        if (mdct_gains != 0)
        {
            mdct_gains[k] = g;
            move16();
        }

        if (mdct_gains_exp != 0)
        {
            mdct_gains_exp[k] = g_e;
            move16();
        }

        if (mdct_inv_gains != 0)
        {
            mdct_inv_gains[k] = ig;
            move16();
        }

        if (mdct_inv_gains_exp != 0)
        {
            mdct_inv_gains_exp[k] = ig_e;
            move16();
        }

        k = add(k, 1);


        s = sub(norm_l(L_max(L_abs(RealData[j]), L_abs(ImagData[j]))), 1);

        tmp16 = extract_h(L_shl(RealData[j], s));
        tmp32 = L_mult(tmp16, tmp16);

        tmp16 = extract_h(L_shl(ImagData[j], s));
        tmp16 = mac_r(tmp32, tmp16, tmp16);

        s = shl(sub(scale, s), 1);

        if (tmp16 == 0)
        {
            s = -16;
            move16();
        }
        if (tmp16 == 0)
        {
            tmp16 = 1;
            move16();
        }

        BASOP_Util_Sqrt_InvSqrt_MantExp(tmp16, s, &g, &g_e, &ig, &ig_e);

        if (mdct_gains != 0)
        {
            mdct_gains[k] = g;
            move16();
        }

        if (mdct_gains_exp != 0)
        {
            mdct_gains_exp[k] = g_e;
            move16();
        }

        if (mdct_inv_gains != 0)
        {
            mdct_inv_gains[k] = ig;
            move16();
        }

        if (mdct_inv_gains_exp != 0)
        {
            mdct_inv_gains_exp[k] = ig_e;
            move16();
        }

        j = sub(j, 1);
        k = add(k, 1);
    }

}


void basop_mdct_noiseShaping_interp(Word32 x[], Word16 lg, Word16 gains[], Word16 gains_exp[])
{
    Word16 i, j, jp, jn, k, l;
    Word16 g, pg, ng, e, tmp;


    assert(lg % FDNS_NPTS == 0);
    k = shr(lg, 6); /* FDNS_NPTS = 64 */

    IF (gains)
    {
        /* Linear interpolation */
        IF (sub(k, 4) == 0)
        {
            jp = 0;
            move16();
            j = 0;
            move16();
            jn = 1;
            move16();

            FOR (i = 0; i < lg; i += 4)
            {
                pg = gains[jp];
                move16();
                g  = gains[j];
                move16();
                ng = gains[jn];
                move16();

                /* common exponent for pg and g */
                tmp = sub(gains_exp[j], gains_exp[jp]);
                if (tmp > 0) pg = shr(pg, tmp);
                if (tmp < 0) g = shl(g, tmp);
                e = s_max(gains_exp[j], gains_exp[jp]);

                tmp = mac_r(L_mult(pg, FL2WORD16(0.375f)), g, FL2WORD16(0.625f));
                x[i] = L_shl(Mpy_32_16(x[i], tmp), e);
                move32();

                tmp = mac_r(L_mult(pg, FL2WORD16(0.125f)), g, FL2WORD16(0.875f));
                x[i+1] = L_shl(Mpy_32_16(x[i+1], tmp), e);
                move32();

                /* common exponent for g and ng */
                g = gains[j];
                move16();
                tmp = sub(gains_exp[j], gains_exp[jn]);
                if (tmp > 0) ng = shr(ng, tmp);
                if (tmp < 0) g = shl(g, tmp);
                e = s_max(gains_exp[j], gains_exp[jn]);

                tmp = mac_r(L_mult(g, FL2WORD16(0.875f)), ng, FL2WORD16(0.125f));
                x[i+2] = L_shl(Mpy_32_16(x[i+2], tmp), e);
                move32();

                tmp = mac_r(L_mult(g, FL2WORD16(0.625f)), ng, FL2WORD16(0.375f));
                x[i+3] = L_shl(Mpy_32_16(x[i+3], tmp), e);
                move32();

                jp = j;
                move16();
                j = jn;
                move16();
                jn = s_min(add(jn, 1), FDNS_NPTS-1);
            }
        }
        ELSE IF (sub(k, 5) == 0)
        {
            jp = 0;
            move16();
            j = 0;
            move16();
            jn = 1;
            move16();

            FOR (i = 0; i < lg; i += 5)
            {
                pg = gains[jp];
                move16();
                g  = gains[j];
                move16();
                ng = gains[jn];
                move16();

                /* common exponent for pg and g */
                tmp = sub(gains_exp[j], gains_exp[jp]);
                if (tmp > 0) pg = shr(pg, tmp);
                if (tmp < 0) g = shl(g, tmp);
                e = s_max(gains_exp[j], gains_exp[jp]);

                tmp = mac_r(L_mult(pg, FL2WORD16(0.40f)), g, FL2WORD16(0.60f));
                x[i]   = L_shl(Mpy_32_16(x[i], tmp), e);
                move32();

                tmp = mac_r(L_mult(pg, FL2WORD16(0.20f)), g, FL2WORD16(0.80f));
                x[i+1] = L_shl(Mpy_32_16(x[i+1], tmp), e);
                move32();


                x[i+2] = L_shl(Mpy_32_16(x[i+2], gains[j]), gains_exp[j]);
                move32();

                /* common exponent for g and ng */
                g = gains[j];
                move16();
                tmp = sub(gains_exp[j], gains_exp[jn]);
                if (tmp > 0) ng = shr(ng, tmp);
                if (tmp < 0) g = shl(g, tmp);
                e = s_max(gains_exp[j], gains_exp[jn]);

                tmp = mac_r(L_mult(g, FL2WORD16(0.80f)), ng, FL2WORD16(0.20f));
                x[i+3] = L_shl(Mpy_32_16(x[i+3], tmp), e);
                move32();

                tmp = mac_r(L_mult(g, FL2WORD16(0.60f)), ng, FL2WORD16(0.40f));
                x[i+4] = L_shl(Mpy_32_16(x[i+4], tmp), e);
                move32();

                jp = j;
                move16();
                j = jn;
                move16();
                jn = s_min(add(jn, 1), FDNS_NPTS-1);
            }
        }
        ELSE   /* no interpolation */
        {
            FOR (i = 0; i < FDNS_NPTS; i++)
            {
                FOR (l = 0; l < k; l++)
                {
                    *x = L_shl(Mpy_32_16(*x, *gains), *gains_exp);
                    move32();
                    x++;
                }

                gains++;
                gains_exp++;
            }
        }
    }

}


void basop_PsychAdaptLowFreqDeemph(Word32 x[],
                                   const Word16 lpcGains[], const Word16 lpcGains_e[],
                                   Word16 lf_deemph_factors[]
                                  )
{
    Word16 i;
    Word16 max, max_e, fac, min, min_e, tmp, tmp_e;
    Word32 L_tmp;



    assert(lpcGains[0] >= 0x4000);

    max = lpcGains[0];
    move16();
    max_e = lpcGains_e[0];
    move16();
    min = lpcGains[0];
    move16();
    min_e = lpcGains_e[0];
    move16();

    /* find minimum (min) and maximum (max) of LPC gains in low frequencies */
    FOR (i = 1; i < 9; i++)
    {
        IF (compMantExp16Unorm(lpcGains[i], lpcGains_e[i], min, min_e) < 0)
        {
            min = lpcGains[i];
            move16();
            min_e = lpcGains_e[i];
            move16();
        }

        IF (compMantExp16Unorm(lpcGains[i], lpcGains_e[i], max, max_e) > 0)
        {
            max = lpcGains[i];
            move16();
            max_e = lpcGains_e[i];
            move16();
        }
    }

    min_e = add(min_e, 5); /* min *= 32.0f; */

    test();
    IF ((compMantExp16Unorm(max, max_e, min, min_e) < 0) && (min > 0))
    {
        /* fac = tmp = (float)pow(max / min, 0.0078125f); */
        tmp_e = min_e;
        move16();
        tmp = Inv16(min, &tmp_e);
        L_tmp = L_shl(L_mult(tmp, max), add(tmp_e, max_e)); /* Q31 */
        L_tmp = BASOP_Util_Log2(L_tmp); /* Q25 */
        L_tmp = L_shr(L_tmp, 7); /* 0.0078125f = 1.f/(1<<7) */
        L_tmp = BASOP_Util_InvLog2(L_tmp); /* Q31 */
        tmp = round_fx(L_tmp); /* Q15 */
        fac = tmp; /* Q15 */                                                        move16();

        /* gradual lowering of lowest 32 bins; DC is lowered by (max/tmp)^1/4 */
        FOR (i = 31; i >= 0; i--)
        {
            x[i] = Mpy_32_16(x[i], fac);
            move32();
            if (lf_deemph_factors != NULL)
            {
                lf_deemph_factors[i] = mult_r(lf_deemph_factors[i], fac);
                move16();
            }
            fac = mult_r(fac, tmp);
        }
    }

}



