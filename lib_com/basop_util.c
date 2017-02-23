/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/


#include <stdio.h>
#include <assert.h>
#include "options.h"
#include "basop_util.h"
#include "rom_com.h"
#include "basop_settings.h"
#include "basop_mpy.h"
#include "control.h"
#include "cnst.h"

extern const Word32 SqrtTable[32];
extern const Word16 SqrtDiffTable[32];

extern const Word32 ISqrtTable[32];
extern const Word16 ISqrtDiffTable[32];

extern const Word32 InvTable[32];
extern const Word16 InvDiffTable[32];


Word32 BASOP_Util_Log2(Word32 x)
{
    Word32  exp;
    Word16  exp_e;
    Word16  nIn;
    Word16  accuSqr;
    Word32  accuRes;

    assert(x >= 0);

    if (x == 0)
    {
        return ((Word32)MIN_32);
    }

    /* normalize input, calculate integer part */
    exp_e = norm_l(x);
    x = L_shl(x,exp_e);
    exp = L_deposit_l(exp_e);

    /* calculate (1-normalized_input) */
    nIn = extract_h(L_sub(MAX_32,x));

    /* approximate ln() for fractional part (nIn *c0 + nIn^2*c1 + nIn^3*c2 + ... + nIn^8 *c7) */

    /* iteration 1, no need for accumulation */
    accuRes = L_mult(nIn,ldCoeff[0]);             /* nIn^i * coeff[0] */
    accuSqr = mult(nIn,nIn);                      /* nIn^2, nIn^3 .... */

    /* iteration 2 */
    accuRes = L_mac(accuRes,accuSqr,ldCoeff[1]);  /* nIn^i * coeff[1] */
    accuSqr = mult(accuSqr,nIn);                  /* nIn^2, nIn^3 .... */

    /* iteration 3 */
    accuRes = L_mac(accuRes,accuSqr,ldCoeff[2]);  /* nIn^i * coeff[2] */
    accuSqr = mult(accuSqr,nIn);                  /* nIn^2, nIn^3 .... */

    /* iteration 4 */
    accuRes = L_mac(accuRes,accuSqr,ldCoeff[3]);  /* nIn^i * coeff[3] */
    accuSqr = mult(accuSqr,nIn);                  /* nIn^2, nIn^3 .... */

    /* iteration 5 */
    accuRes = L_mac(accuRes,accuSqr,ldCoeff[4]);  /* nIn^i * coeff[4] */
    accuSqr = mult(accuSqr,nIn);                  /* nIn^2, nIn^3 .... */

    /* iteration 6 */
    accuRes = L_mac(accuRes,accuSqr,ldCoeff[5]);  /* nIn^i * coeff[5] */
    accuSqr = mult(accuSqr,nIn);                  /* nIn^2, nIn^3 .... */

    /* iteration 7, no need to calculate accuSqr any more */
    accuRes = L_mac(accuRes,accuSqr,ldCoeff[6]);  /* nIn^i * coeff[6] */

    /* ld(fractional part) = ln(fractional part)/ln(2), 1/ln(2) = (1 + 0.44269504) */
    accuRes = L_mac0(L_shr(accuRes,1),extract_h(accuRes),14506);

    accuRes = L_shr(accuRes,LD_DATA_SCALE-1);     /* fractional part/LD_DATA_SCALE */
    exp = L_shl(exp,(31-LD_DATA_SCALE));          /* integer part/LD_DATA_SCALE */
    accuRes = L_sub(accuRes,exp);                 /* result = integer part + fractional part */

    return (accuRes);
}


Word32 BASOP_Util_InvLog2(Word32 x)
{
    Word16  frac;
    Word16  exp;
    Word32  retVal;
    UWord32 index3;
    UWord32 index2;
    UWord32 index1;
    UWord32 lookup3f;
    UWord32 lookup12;
    UWord32 lookup;

    if ( x < FL2WORD32(-31.0/64.0) )
    {
        return 0;
    }
    test();
    if ( (L_sub(x,FL2WORD32(31.0/64.0)) >= 0) || (x == 0) )
    {
        return 0x7FFFFFFF;
    }

    frac   = extract_l(L_and(x,0x3FF));

    index3 = L_and(L_shr(x,10),0x1F);
    index2 = L_and(L_shr(x,15),0x1F);
    index1 = L_and(L_shr(x,20),0x1F);

    exp = extract_l(L_shr(x,25));
    if ( x > 0 )
    {
        exp = sub(31,exp);
    }
    if ( x < 0 )
    {
        exp = negate(exp);
    }

    lookup3f = L_add(exp2x_tab_long[index3],L_shr(Mpy_32_16(0x0016302F,frac),1));
    lookup12 = Mpy_32_32(exp2_tab_long[index1],exp2w_tab_long[index2]);
    lookup   = Mpy_32_32(lookup12, lookup3f);

    retVal = L_shr(lookup,sub(exp,3));

    return retVal;
}


Word16 BASOP_Util_Add_MantExp                    /*!< Exponent of result        */
(Word16   a_m,       /*!< Mantissa of 1st operand a */
 Word16   a_e,       /*!< Exponent of 1st operand a */
 Word16   b_m,       /*!< Mantissa of 2nd operand b */
 Word16   b_e,       /*!< Exponent of 2nd operand b */
 Word16  *ptrSum_m)  /*!< Mantissa of result */
{
    Word32 L_lm, L_hm;
    Word16 shift;

    /* Compare exponents: the difference is limited to +/- 15
       The Word16 mantissa of the operand with higher exponent is moved into the low
     part of a Word32 and shifted left by the exponent difference. Then, the
     unshifted mantissa of the operand with the lower exponent is added to the lower
     16 bits. The addition result is normalized and the upper Word16 of the result represents
     the mantissa to return. The returned exponent takes into account all shift operations
     including the final 16-bit extraction.
     Note: The resulting mantissa may be inaccurate in the case, where the mantissa of the operand
           with higher exponent is not really left-aligned, while the mantissa of the operand with
    	   lower exponent is so. If in such a case, the difference in exponents is more than 15,
    	   an inaccuracy is introduced.
    	   Example:
    	   A: a_e = 20, a_m = 0x0001
    	   B: b_e =  0, b_m = 0x4000
             correct:      A+B=1*2^20+1*2^14=0x0010.0000+0x0000.4000=0x0010.4000=0x4100*2^6
    	   previously:   A+B=1*2^20+1*2^14=0x0001+0x0000=0x0001*2^20
    	   this version: A+B=1*2^20+1*2^14=0x0000.8000+0x0000.4000=0x6000*2^6
    */

    shift = sub(a_e, b_e);
    if (shift >= 0)
        shift = s_min( 15,shift);

    if (shift < 0)
        shift = s_max(-15,shift);
    a_e = s_max(a_e,b_e);
    L_hm = L_deposit_l(a_m);       /* mantissa belonging to higher exponent */
    L_lm = L_deposit_l(a_m);       /* mantissa belonging to lower exponent */
    if (shift >= 0)
        L_lm = L_deposit_l(b_m);
    if (shift < 0)
        L_hm = L_deposit_l(b_m);

    if (shift > 0)
        shift = negate(shift);

    L_hm = L_shr(L_hm,shift);    /* shift left due to negative shift parameter */
    a_e = add(a_e, shift);
    L_hm = L_add(L_hm, L_lm);
    shift = norm_l(L_hm);
    L_hm = L_shl(L_hm,shift);
    *ptrSum_m = extract_h(L_hm);
    move16();

    a_e = sub(a_e,shift);
    if (L_hm)
        a_e = add(a_e,16);

    return (a_e);
}


/* local function for Sqrt16 */
static Word16 Sqrt16_common(Word16 m,
                            Word16 e)
{
    Word16 index, frac;

    assert((m >= 0x4000) || (m == 0));

    /* get table index (upper 6 bits minus 32) */
    /* index = (m >> 9) - 32; */
    index = mac_r(-32768 - (32 << 16), m, 1 << 6);

    /* get fractional part for interpolation (lower 9 bits) */
    frac = s_and(m, 0x1FF); /* Q9 */

    /* interpolate */
    if (m != 0)
    {
        BASOP_SATURATE_WARNING_OFF;
        m = mac_r(SqrtTable[index], SqrtDiffTable[index], frac);
        BASOP_SATURATE_WARNING_ON;
    }

    /* handle odd exponents */
    if (s_and(e, 1) != 0) m = mult_r(m, 0x5a82);

    return m;
}


Word16 Sqrt16(                  /*!< output mantissa */
    Word16 mantissa,  /*!< input mantissa */
    Word16 *exponent  /*!< pointer to exponent */
)
{
    Word16 preShift, e;

    assert(mantissa >= 0);

    /* normalize */
    preShift = norm_s(mantissa);

    e = sub(*exponent, preShift);
    mantissa = shl(mantissa, preShift);

    /* calc mantissa */
    mantissa = Sqrt16_common(mantissa, e);

    /* e = (e + 1) >> 1 */
    *exponent = mult_r(e, 1 << 14);
    move16();

    return mantissa;
}

Word16 Inv16(                  /*!< output mantissa */
    Word16 mantissa,  /*!< input mantissa */
    Word16 *exponent  /*!< pointer to exponent */
)
{
    Word16 index, frac;
    Word16 preShift;
    Word16 m, e;

    assert(mantissa != 0);

    /* absolute */
    BASOP_SATURATE_WARNING_OFF;
    m = abs_s(mantissa);
    BASOP_SATURATE_WARNING_ON;

    /* normalize */
    preShift = norm_s(m);

    e = sub(*exponent, preShift);
    m = shl(m, preShift);

    /* get table index (upper 6 bits minus 32) */
    /* index = (m >> 9) - 32; */
    index = mac_r(-32768 - (32 << 16), m, 1 << 6);

    /* get fractional part for interpolation (lower 9 bits) */
    frac = shl(s_and(m, 0x1FF), 1); /* Q10 */

    /* interpolate */
    m = msu_r(InvTable[index], InvDiffTable[index], frac);

    /* restore sign */
    if (mantissa < 0) m = negate(m);

    /* e = 1 - e */
    *exponent = sub(1, e);
    move16();

    return m;
}


void BASOP_Util_Sqrt_InvSqrt_MantExp (Word16 mantissa,      /*!< mantissa */
                                      Word16 exponent,      /*!< expoinent */
                                      Word16 *sqrt_mant,    /*!< Pointer to sqrt mantissa */
                                      Word16 *sqrt_exp,     /*!< Pointer to sqrt exponent */
                                      Word16 *isqrt_mant,   /*!< Pointer to 1/sqrt mantissa */
                                      Word16 *isqrt_exp     /*!< Pointer to 1/sqrt exponent */
                                     )
{
    Word16 index, frac;
    Word16 preShift;
    Word16 m, mi, e_odd;

    assert(mantissa > 0);

    /* normalize */
    preShift = norm_s(mantissa);

    exponent = sub(exponent, preShift);
    mantissa = shl(mantissa, preShift);

    /* get table index (upper 6 bits minus 32) */
    /* index = (m >> 9) - 32; */
    index = mac_r(-32768 - (32 << 16), mantissa, 1 << 6);

    /* get fractional part for interpolation (lower 9 bits) */
    frac = s_and(mantissa, 0x1FF); /* Q9 */

    /* interpolate */
    BASOP_SATURATE_WARNING_OFF;
    m = mac_r(SqrtTable[index], SqrtDiffTable[index], frac);
    mi = msu_r(ISqrtTable[index], ISqrtDiffTable[index], frac);
    BASOP_SATURATE_WARNING_ON;

    /* handle even/odd exponents */
    e_odd = s_and(exponent, 1);
    if (e_odd != 0) m = mult_r(m, 0x5a82);
    if (e_odd == 0) mi = mult_r(mi, 0x5a82);

    /* e = (e + 1) >> 1 */
    *sqrt_exp = mult_r(exponent, 1 << 14);
    move16();

    /* e = (2 - e) >> 1 */
    *isqrt_exp = msu_r(1L << 15, exponent, 1 << 14);
    move16();

    /* Write result */
    *sqrt_mant = m;
    move16();
    *isqrt_mant = mi;
    move16();

}

/********************************************************************/
/*!
  \brief   Calculates the scalefactor needed to normalize input array

    The scalefactor needed to normalize the Word32 input array is returned <br>
    If the input array contains only '0', a scalefactor 0 is returned <br>
    Scaling factor is determined wrt a normalized target x: 1073741824 <= x <= 2147483647 for positive x <br>
    and   -2147483648 <= x <= -1073741824 for negative x
*/

Word16 getScaleFactor32(                 /* o: measured headroom in range [0..31], 0 if all x[i] == 0 */
    const Word32 *x,      /* i: array containing 32-bit data */
    const Word16 len_x)   /* i: length of the array to scan  */
{
    Word16 i, i_min, i_max;
    Word32 x_min, x_max;

    x_max = L_add(0, 0);
    x_min = L_add(0, 0);
    FOR (i = 0; i < len_x; i++)
    {
        if (x[i] >= 0)
            x_max = L_max(x_max,x[i]);
        if (x[i] < 0)
            x_min = L_min(x_min,x[i]);
    }

    i_max = 0x20;
    move16();
    i_min = 0x20;
    move16();

    if (x_max != 0)
        i_max = norm_l(x_max);

    if (x_min != 0)
        i_min = norm_l(x_min);

    i = s_and(s_min(i_max, i_min),0x1F);

    return i;
}


Word16 BASOP_Util_Divide1616_Scale(Word16 x, Word16 y, Word16 *s)
{
    Word16 z;
    Word16 sx;
    Word16 sy;
    Word16 sign;

    /* assert (x >= (Word16)0); */
    assert (y != (Word16)0);

    sign = 0;
    move16();

    IF ( x < 0 )
    {
        x = negate(x);
        sign = s_xor(sign,1);
    }

    IF ( y < 0 )
    {
        y = negate(y);
        sign= s_xor(sign,1);
    }

    IF ( x == (Word16)0 )
    {
        move16();
        *s = 0;

        return ((Word16)0);
    }

    sx = norm_s(x);
    x  = shl(x,sx);
    x  = shr(x,1);
    move16();
    *s = sub(1,sx);

    sy = norm_s(y);
    y  = shl(y,sy);
    move16();
    *s = add(*s,sy);

    z = div_s(x,y);

    if ( sign != 0 )
    {
        z = negate(z);
    }

    return z;
}


void set_val_Word16(Word16 X[], const Word16 val, Word16 n)
{
    Word16 i;


    FOR (i = 0; i < n; i++)
    {
        X[i] = val;
        move16();
    }


    return;
}

void set_val_Word32(Word32 X[], const Word32 val, Word16 n)
{
    Word16 i;


    FOR (i = 0; i < n; i++)
    {
        X[i] = val;
        move32();
    }


    return;
}

Word16 mult0 (  Word16 x, Word16 y)
{
    return extract_l(L_mult0(x,y));
}


#define SINETAB SineTable512_fx
#define LD 9

/*
 * Calculates coarse lookup values for sine/cosine and residual angle.
 * \param x angle in radians with exponent = 2 or as radix 2 with exponent = 0.
 * \param scale shall always be 2
 * \param sine pointer to where the sine lookup value is stored into
 * \param cosine pointer to where the cosine lookup value is stored into
 * \param flag_radix2 flag indicating radix 2 angle if non-zero.
 */
static Word16 fixp_sin_cos_residual_16(Word16 x, const Word16 scale, Word16 *sine, Word16 *cosine, Word8 flag_radix2)
{
    Word16 residual;
    Word16 s;
    Word16 ssign;
    Word16 csign;
    Word16 tmp, cl = 0, sl = 0;
    const Word16 shift = 15-LD-1-scale;

    if (flag_radix2 == 0)
    {
        x = mult_r(x, FL2WORD16(1.0/EVS_PI));
    }
    s = shr(x, shift);

    residual = s_and(x, (1<<shift)-1);
    /* We assume "2+scale" is a constant */
    residual = shl(residual,2+scale);
    residual = mult_r(residual,FL2WORD16(EVS_PI/4.0));

    /* Sine sign symmetry */
    ssign = s_and(s, (1<<LD)<<1);

    /* Cosine sign symmetry */
    csign = s_and(add(s, (1<<LD)), (1<<LD)<<1);

    /* Modulo EVS_PI */
    s = s_and(s, (2<<LD)-1);

    /* EVS_PI/2 symmetry */
    s = s_min(s, sub(2<<LD, s));

    {
        tmp = s_min(sub(1<<LD, s), s);
        s = sub(tmp,s);

        if ( ! s)
        {
            move16();
            sl = SINETAB[tmp].v.im;
        }
        if (! s)
        {
            move16();
            cl = SINETAB[tmp].v.re;
        }
        if (s)
        {
            move16();
            sl = SINETAB[tmp].v.re;
        }
        if (s)
        {
            move16();
            cl = SINETAB[tmp].v.im;
        }

        if (ssign)
        {
            sl = negate(sl);
        }
        if (csign)
        {
            cl = negate(cl);
        }

        move16();
        move16();
        *sine   = sl;
        *cosine = cl;
    }

    return residual;
}


Word16 getCosWord16R2(Word16 theta)
{
    Word16 result, residual, sine, cosine;

    residual = fixp_sin_cos_residual_16(theta, 1, &sine, &cosine, 1);
    /* This negation prevents the subsequent addition from overflow */
    /* The negation cannot overflow, sine is in range [0x0..0x7FFF] */
    BASOP_SATURATE_WARNING_OFF
    sine = negate(sine);
    result = msu_r(L_mult(sine, residual), cosine, -32768);
    BASOP_SATURATE_WARNING_ON

    return result;
}


Word16 idiv1616U(Word16 x, Word16 y)
{
    Word16 s;


    /* make y > x */
    s = add(sub(norm_s(y), norm_s(x)), 1);
    s = s_max(s, 0);

    BASOP_SATURATE_WARNING_OFF
    y = shl(y, s);
    BASOP_SATURATE_WARNING_ON

    /* divide and shift */
    y = shr(div_s(x, y), sub(15, s));


    return y;
}

Word32 BASOP_util_Pow2(
    const Word32 exp_m, const Word16 exp_e,
    Word16 *result_e
)
{
    static const Word16 pow2Coeff[8] =
    {
        FL2WORD16(0.693147180559945309417232121458177),    /* ln(2)^1 /1! */
        FL2WORD16(0.240226506959100712333551263163332),    /* ln(2)^2 /2! */
        FL2WORD16(0.0555041086648215799531422637686218),   /* ln(2)^3 /3! */
        FL2WORD16(0.00961812910762847716197907157365887),  /* ln(2)^4 /4! */
        FL2WORD16(0.00133335581464284434234122219879962),  /* ln(2)^5 /5! */
        FL2WORD16(1.54035303933816099544370973327423e-4),  /* ln(2)^6 /6! */
        FL2WORD16(1.52527338040598402800254390120096e-5),  /* ln(2)^7 /7! */
        FL2WORD16(1.32154867901443094884037582282884e-6)   /* ln(2)^8 /8! */
    };

    Word32 frac_part = 0, tmp_frac, result_m;
    Word16 int_part = 0;

    IF (exp_e > 0)
    {
        /* "+ 1" compensates L_shr(,1) of the polynomial evaluation at the loop end. */

        int_part = add(1,extract_l(L_shr(exp_m, sub(31, exp_e))));
        frac_part = L_lshl(exp_m, exp_e);
        frac_part = L_and(0x7FFFFFFF, frac_part);
    }
    if (exp_e <= 0)
        frac_part = L_shl(exp_m, exp_e);
    if (exp_e <= 0)
    {
        int_part = 1;
        move16();
    }

    /* Best accuracy is around 0, so try to get there with the fractional part. */
    IF( (tmp_frac = L_sub(frac_part,FL2WORD32(0.5))) >= 0)
    {
        int_part = add(int_part, 1);
        frac_part = L_sub(tmp_frac,FL2WORD32(0.5));
    }
    ELSE IF( (tmp_frac = L_add(frac_part,FL2WORD32(0.5))) < 0)
    {
        int_part = sub(int_part, 1);
        frac_part = L_add(tmp_frac,FL2WORD32(0.5));
    }

    /* Evaluate taylor polynomial which approximates 2^x */
    {
        Word32 p;
        Word16 i;


        /* First taylor series coefficient a_0 = 1.0, scaled by 0.5 due to L_shr(,1). */
        result_m = L_add(FL2WORD32(1.0/2.0),L_shr(Mpy_32_16(frac_part, pow2Coeff[0]), 1));
        p = Mpy_32_32(frac_part, frac_part);
        FOR (i = 1; i < 7; i++)
        {
            /* next taylor series term: a_i * x^i, x=0 */
            result_m = L_add(result_m, L_shr(Mpy_32_16(p, pow2Coeff[i]), 1));
            p = Mpy_32_32(p, frac_part);
        }
        result_m = L_add(result_m, L_shr(Mpy_32_16(p, pow2Coeff[i]), 1));
    }
    *result_e = int_part;
    move16();

    return result_m;
}

Word16 BASOP_Util_Divide3216_Scale(     /* o: result of division x/y, not normalized  */
    Word32 x,                             /* i: numerator, signed                       */
    Word16 y,                             /* i: denominator, signed                     */
    Word16 *s)                            /* o: scaling, 0, if x==0                     */
{
    Word16 z;
    Word16 sx;
    Word16 sy;
    Word16 sign;

    /*assert (x > (Word32)0);
    assert (y >= (Word16)0);*/

    /* check, if numerator equals zero, return zero then */
    IF ( x == (Word32)0 )
    {
        *s = 0;
        move16();

        return ((Word16)0);
    }

    sign = s_xor(extract_h(x),y);  /* just to exor the sign bits */
    BASOP_SATURATE_WARNING_OFF
    x = L_abs(x);
    y = abs_s(y);
    BASOP_SATURATE_WARNING_ON
    sx = sub(norm_l(x),1);
    x  = L_shl(x,sx);
    sy = norm_s(y);
    y  = shl(y,sy);
    *s = sub(sy,sx);
    move16();

    z = div_s(round_fx(x),y);

    if ( sign < 0 )                /* if sign bits differ, negate the result */
    {
        z = negate(z);
    }

    return z;
}


static const Word16 table_pow2[32] =
{
    16384, 16743, 17109, 17484, 17867, 18258, 18658, 19066, 19484, 19911,
    20347, 20792, 21247, 21713, 22188, 22674, 23170, 23678, 24196, 24726,
    25268, 25821, 26386, 26964, 27554, 28158, 28774, 29405, 30048, 30706,
    31379, 32066
};
/* table of table_pow2[i+1] - table_pow2[i] */
static const Word16 table_pow2_diff_x32[32] =
{
    11488, 11712, 12000, 12256, 12512, 12800, 13056, 13376, 13664, 13952,
    14240, 14560, 14912, 15200, 15552, 15872, 16256, 16576, 16960, 17344,
    17696, 18080, 18496, 18880, 19328, 19712, 20192, 20576, 21056, 21536,
    21984, 22432
};

Word32 Pow2(                              /* (o) Q0  : result       (range: 0<=val<=0x7fffffff) */
    Word16 exponant,                      /* (i) Q0  : Integer part.      (range: 0<=val<=30)   */
    Word16 fraction                       /* (i) Q15 : Fractional part.   (range: 0.0<=val<1.0) */
)
{
    Word16 exp, i, a;
    Word32 L_x;

    i = mac_r(-32768, fraction, 32);         /* Extract b10-b16 of fraction */
    a = s_and(fraction, 0x3ff);              /* Extract  b0-b9  of fraction */

    L_x = L_deposit_h(table_pow2[i]);           /* table[i] << 16   */
    L_x = L_mac(L_x, table_pow2_diff_x32[i], a);/* L_x -= diff*a*2  */

    exp = sub(30, exponant);

    L_x = L_shr_r(L_x, exp);

    return L_x;
}

/*************************************************************************
 *
 *   FUNCTION:   Log2_norm()
 *
 *   PURPOSE:   Computes log2(L_x, exp),  where   L_x is positive and
 *              normalized, and exp is the normalisation exponent
 *              If L_x is negative or zero, the result is 0.
 *
 *   DESCRIPTION:
 *        The function Log2(L_x) is approximated by a table and linear
 *        interpolation. The following steps are used to compute Log2(L_x)
 *
 *           1- exponent = 30-norm_exponent
 *           2- i = bit25-b31 of L_x;  32<=i<=63  (because of normalization).
 *           3- a = bit10-b24
 *           4- i -=32
 *           5- fraction = table[i]<<16 - (table[i] - table[i+1]) * a * 2
 *
 *************************************************************************/

static const Word32 L_table[32] =
{
    -32768L,   95322112L,  187793408L,  277577728L,
    364871680L,  449740800L,  532381696L,  612859904L,
    691306496L,  767787008L,  842432512L,  915308544L,
    986546176L, 1056210944L, 1124302848L, 1190887424L,
    1256095744L, 1319993344L, 1382580224L, 1443921920L,
    1504083968L, 1563131904L, 1621000192L, 1677885440L,
    1733722112L, 1788510208L, 1842380800L, 1895399424L,
    1947435008L, 1998618624L, 2049015808L, 2098626560L
};

static const Word16 table_diff[32] =
{
    1455, 1411,  1370, 1332,  1295, 1261,  1228, 1197,
    1167, 1139,  1112, 1087,  1063, 1039,  1016,  995,
    975,  955,   936,  918,   901,  883,   868,  852,
    836,  822,   809,  794,   781,  769,   757,  744
};

Word16 Log2_norm_lc (   /* (o) : Fractional part of Log2. (range: 0<=val<1)  */
    Word32 L_x          /* (i) : input value (normalized)                    */
)
{
    Word16 i, a;
    Word16 y;


    L_x = L_shr (L_x, 9);
    a = extract_l (L_x);                      /* Extract b10-b24 of fraction */
    a = lshr(a, 1);

    i = mac_r(L_x, -32*2-1, 16384);           /* Extract b25-b31 minus 32 */

    y = mac_r(L_table[i], table_diff[i], a);/* table[i] << 16 - diff*a*2 */


    return y;
}


Word32 BASOP_Util_fPow(
    Word32 base_m, Word16 base_e,
    Word32 exp_m, Word16 exp_e,
    Word16 *result_e
)
{

    Word16 ans_lg2_e, base_lg2_e;
    Word32 base_lg2_m, ans_lg2_m, result_m;
    Word16 shift;

    test();
    IF ((base_m == 0) && (exp_m != 0))
    {
        *result_e = 0;
        move16();
        return 0;
    }
    /* Calc log2 of base */
    shift      = norm_l(base_m);
    base_m     = L_shl(base_m, shift);
    base_e     = sub(base_e, shift);
    base_lg2_m = BASOP_Util_Log2(base_m);

    /* shift: max left shift such that neither base_e or base_lg2_m saturate. */
    shift = sub(s_min(norm_s(base_e), WORD16_BITS-1-LD_DATA_SCALE), 1);
    /* Compensate shift into exponent of result. */
    base_lg2_e = sub(WORD16_BITS-1, shift);
    base_lg2_m = L_add(L_shr(base_lg2_m, sub(WORD16_BITS-1-LD_DATA_SCALE, shift)), L_deposit_h(shl(base_e, shift)));

    /* Prepare exp */
    shift = norm_l(exp_m);
    exp_m = L_shl(exp_m, shift);
    exp_e = sub(exp_e, shift);

    /* Calc base pow exp */
    ans_lg2_m = Mpy_32_32(base_lg2_m, exp_m);
    ans_lg2_e = add(exp_e, base_lg2_e);

    /* Calc antilog */
    result_m = BASOP_util_Pow2(ans_lg2_m, ans_lg2_e, result_e);

    return result_m;
}

Word32 BASOP_Util_Add_Mant32Exp                  /*!< o: normalized result mantissa */
(Word32   a_m,       /*!< i: Mantissa of 1st operand a  */
 Word16   a_e,       /*!< i: Exponent of 1st operand a  */
 Word32   b_m,       /*!< i: Mantissa of 2nd operand b  */
 Word16   b_e,       /*!< i: Exponent of 2nd operand b  */
 Word16  *ptr_e)     /*!< o: exponent of result         */
{
    Word32 L_tmp;
    Word16 shift;

    /* Compare exponents: the difference is limited to +/- 30
       The Word32 mantissa of the operand with lower exponent is shifted right by the exponent difference.
       Then, the unshifted mantissa of the operand with the higher exponent is added. The addition result
       is normalized and the result represents the mantissa to return. The returned exponent takes into
       account all shift operations.
    */

    if (!a_m)
        a_e = add(b_e,0);

    if (!b_m)
        b_e = add(a_e,0);

    shift = sub(a_e, b_e);
    shift = s_max(-31,shift);
    shift = s_min(31, shift);
    if (shift < 0)
    {
        /* exponent of b is greater than exponent of a, shr a_m */
        a_m = L_shl(a_m,shift);
    }
    if (shift > 0)
    {
        /* exponent of a is greater than exponent of b */
        b_m = L_shr(b_m,shift);
    }
    a_e = add(s_max(a_e,b_e),1);
    L_tmp = L_add(L_shr(a_m,1),L_shr(b_m,1));
    shift = norm_l(L_tmp);
    if (shift)
        L_tmp = L_shl(L_tmp,shift);
    if (L_tmp == 0)
        a_e = add(0,0);
    if (L_tmp != 0)
        a_e = sub(a_e,shift);
    *ptr_e = a_e;

    return (L_tmp);
}

static const Word32 L_table_isqrt[48] =
{
    2147418112L,  2083389440L,  2024669184L,  1970667520L,
    1920794624L,  1874460672L,  1831403520L,  1791098880L,
    1753415680L,  1717960704L,  1684602880L,  1653145600L,
    1623326720L,  1595080704L,  1568276480L,  1542782976L,
    1518469120L,  1495334912L,  1473183744L,  1451950080L,
    1431633920L,  1412169728L,  1393491968L,  1375469568L,
    1358168064L,  1341521920L,  1325465600L,  1309933568L,
    1294991360L,  1280507904L,  1266548736L,  1252982784L,
    1239875584L,  1227161600L,  1214775296L,  1202847744L,
    1191182336L,  1179910144L,  1168965632L,  1158283264L,
    1147863040L,  1137770496L,  1127940096L,  1118306304L,
    1108934656L,  1099825152L,  1090912256L,  1082261504L
};
/* table of table_isqrt[i] - table_isqrt[i+1] */
static const Word16 table_isqrt_diff[48] =
{
    977,   896,   824,   761,   707,   657,   615,   575,
    541,   509,   480,   455,   431,   409,   389,   371,
    353,   338,   324,   310,   297,   285,   275,   264,
    254,   245,   237,   228,   221,   213,   207,   200,
    194,   189,   182,   178,   172,   167,   163,   159,
    154,   150,   147,   143,   139,   136,   132,   130
};
static const Word16 shift[] = {9,10};
Word32 Isqrt_lc1(
    Word32 frac,  /* (i)   Q31: normalized value (1.0 < frac <= 0.5) */
    Word16 * exp  /* (i/o)    : exponent (value = frac x 2^exponent) */
)
{
    Word16 i, a;
    Word32 L_tmp;

    IF (frac <= (Word32) 0)
    {
        *exp = 0;
        move16();
        return 0x7fffffff; /*0x7fffffff*/
    }

    /* If exponant odd -> shift right by 10 (otherwise 9) */
    L_tmp = L_shr(frac, shift[s_and(*exp, 1)]);

    /* 1) -16384 to shift left and change sign                 */
    /* 2) 32768 to Add 1 to Exponent like it was divided by 2  */
    /* 3) We let the mac_r add another 0.5 because it imitates */
    /*    the behavior of shr on negative number that should   */
    /*    not be rounded towards negative infinity.            */
    /* It replaces:                                            */
    /*    *exp = negate(shr(sub(*exp, 1), 1));   move16();     */
    *exp = mac_r(32768, *exp, -16384);
    move16();

    a = extract_l(L_tmp);                           /* Extract b10-b24 */
    a = lshr(a, 1);

    i = mac_r(L_tmp, -16*2-1, 16384);               /* Extract b25-b31 minus 16 */

    L_tmp = L_msu(L_table_isqrt[i], table_isqrt_diff[i], a);/* table[i] << 16 - diff*a*2 */

    return L_tmp;
}


