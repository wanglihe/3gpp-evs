/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#ifndef __BASOP_UTIL_H__
#define __BASOP_UTIL_H__

#include "basop_settings.h"
#include "typedef.h"
#include "basop32.h"
#include "basop_mpy.h"


#define LD_DATA_SCALE     (6)

#define LD_DATA_SHIFT_I5  (7)

#define modDiv2(x)        sub(x,shl(shr(x,1),1))
#define modDiv8(x)        L_sub(x,L_shl(L_shr(x,3),3))

#ifndef CHEAP_NORM_SIZE
#define CHEAP_NORM_SIZE 161
#endif

static __inline Word16 limitScale16( Word16 s)
{
    /* It is assumed, that s is calculated just before, therefore we can switch upon sign */
    if (s >= 0)
        s = s_min(s,WORD16_BITS-1);
    if (s < 0)
        s = s_max(s,1-WORD16_BITS);
    return (s);
}

static __inline Word16 limitScale32( Word16 s)
{
    /* It is assumed, that s is calculated just before, therefore we can switch upon sign */
    if (s >= 0)
        s = s_min(s, WORD32_BITS-1);
    if (s < 0)
        s = s_max(s, 1-WORD32_BITS);
    return (s);
}


/*!**********************************************************************
   \brief   Add two values given by mantissa and exponent.

   Mantissas are in 16-bit-fractional format with values between 0 and 1. <br>
   The base for exponents is 2.  Example:  \f$  a = a\_m * 2^{a\_e}  \f$<br>

************************************************************************/
Word16 BASOP_Util_Add_MantExp                    /*!< Exponent of result        */
(Word16   a_m,       /*!< Mantissa of 1st operand a */
 Word16   a_e,       /*!< Exponent of 1st operand a */
 Word16   b_m,       /*!< Mantissa of 2nd operand b */
 Word16   b_e,       /*!< Exponent of 2nd operand b */
 Word16  *ptrSum_m); /*!< Mantissa of result */



/************************************************************************/
/*!
  \brief   Calculate the squareroot of a number given by mantissa and exponent

  Mantissa is in 16/32-bit-fractional format with values between 0 and 1. <br>
  For *norm versions mantissa has to be between 0.5 and 1. <br>
  The base for the exponent is 2.  Example:  \f$  a = a\_m * 2^{a\_e}  \f$<br>
  The exponent is addressed via pointers and will be overwritten with the result.
*/
Word16 Sqrt16(                  /*!< output mantissa */
    Word16 mantissa,  /*!< input mantissa */
    Word16 *exponent  /*!< pointer to exponent */
);

/*****************************************************************************/
/*!
  \brief   Calculate the inverse of a number given by mantissa and exponent

  Mantissa is in 16-bit-fractional format with values between 0 and 1. <br>
  The base for the exponent is 2.  Example:  \f$  a = a\_m * 2^{a\_e}  \f$<br>
  The operand is addressed via pointers and will be overwritten with the result.

  The function uses a table lookup and a newton iteration.
*/
Word16 Inv16(                  /*!< output mantissa */
    Word16 mantissa,  /*!< input mantissa */
    Word16 *exponent  /*!< pointer to exponent */
);

/******************************************************************************/
/*!
  \brief   Calculate the squareroot and inverse of squareroot of a number given by mantissa and exponent

  Mantissa is in 16-bit-fractional format with values between 0 and 1. <br>
  The base for the exponent is 2.  Example:  \f$  a = a\_m * 2^{a\_e}  \f$<br>
*/
void BASOP_Util_Sqrt_InvSqrt_MantExp (Word16 mantissa,      /*!< mantissa */
                                      Word16 exponent,      /*!< expoinent */
                                      Word16 *sqrt_mant,    /*!< Pointer to sqrt mantissa */
                                      Word16 *sqrt_exp,     /*!< Pointer to sqrt exponent */
                                      Word16 *isqrt_mant,   /*!< Pointer to 1/sqrt mantissa */
                                      Word16 *isqrt_exp     /*!< Pointer to 1/sqrt exponent */
                                     );



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
    const Word16 len_x);  /* i: length of the array to scan  */

/************************************************************************/
/*!
  \brief 	Binary logarithm with 7 iterations

  \param   x

  \return log2(x)/64
 */
/************************************************************************/
Word32 BASOP_Util_Log2(Word32 x);



/****************************************************************************/
/*!
  \brief   Does fractional division of Word16 arg1 by Word16 arg2


  \return fractional Q15 Word16 z = arg1(Q15)/arg2(Q15)  with scaling s
*/
Word16 BASOP_Util_Divide1616_Scale( Word16 x,   /*!< i  : Numerator*/
                                    Word16 y,   /*!< i  : Denominator*/
                                    Word16 *s); /*!< o  : Additional scalefactor difference*/

/****************************************************************************/
/*!
  \brief   Does fractional integer division of Word32 arg1 by Word16 arg2


  \return fractional Word16 integer z = arg1(32bits)/arg2(16bits) , z not normalized
*/
Word16 BASOP_Util_Divide3216_Scale(     Word32 x,   /*!< i  : Numerator  */
                                        Word16 y,   /*!< i  : Denominator*/
                                        Word16 *s); /*!< o  : Additional scalefactor difference*/


/************************************************************************/
/*!
 * \brief 	Binary logarithm with 5 iterations
 *
 * \param[i] val
 *
 * \return basop_log2(val)/128
 */
/************************************************************************/
Word32 BASOP_Util_log2_i5(Word32 val);

/************************************************************************/
/*!
  \brief 	Binary power

  Date: 06-JULY-2012 Arthur Tritthart, IIS Fraunhofer Erlangen

  Version with 3 table lookup and 1 linear interpolations

  Algorithm: compute power of 2, argument x is in Q7.25 format
             result = 2^(x/64)
             We split exponent (x/64) into 5 components:
             integer part:      represented by b31..b25  (exp)
             fractional part 1: represented by b24..b20  (lookup1)
             fractional part 2: represented by b19..b15  (lookup2)
             fractional part 3: represented by b14..b10  (lookup3)
             fractional part 4: represented by b09..b00  (frac)
             => result = (lookup1*lookup2*(lookup3+C1*frac)<<3)>>exp

  Due to the fact, that all lookup values contain a factor 0.5
  the result has to be shifted by 3 to the right also.
  Table exp2_tab_long contains the log2 for 0 to 1.0 in steps
  of 1/32, table exp2w_tab_long the log2 for 0 to 1/32 in steps
  of 1/1024, table exp2x_tab_long the log2 for 0 to 1/1024 in
  steps of 1/32768. Since the 2-logarithm of very very small
  negative value is rather linear, we can use interpolation.

  Limitations:

  For x <= 0, the result is fractional positive
  For x > 0, the result is integer in range 1...7FFF.FFFF
  For x < -31/64, we have to clear the result
  For x = 0, the result is ~1.0 (0x7FFF.FFFF)
  For x >= 31/64, the result is 0x7FFF.FFFF

  \param  x

  \return pow(2,(x/64))
 */
/************************************************************************/
Word32 BASOP_Util_InvLog2(Word32 x);


/****************************************************************************/
/*!
  \brief Sets Array Word16 arg1 to value Word16 arg2 for Word16 arg3 elements
*/
void set_val_Word16(    Word16 X[],         /*!< Address of array               */
                        const Word16 val,   /*!< Value to copy into array       */
                        Word16 n);          /*!< Number of elements to process  */


/****************************************************************************/
/*!
  \brief Sets Array Word32 arg1 to value Word32 arg2 for Word16 arg3 elements
*/
void set_val_Word32(    Word32 X[],         /*!< Address of array               */
                        const Word32 val,   /*!< Value to copy into array       */
                        Word16 n);          /*!< Number of elements to process  */

/****************************************************************************/
/*!
  \brief    Does a multiplication of Word16 input values

  \return   z16 = x16 * y16
*/
Word16 mult0 (  Word16 x,   /*!< i  : Multiplier    */
                Word16 y);  /*!< i  : Multiplicand  */

/**
 * \brief calculate cosine of angle. Tuned for ISF domain.
 * \param theta Angle normalized to radix 2, theta = (angle in radians)*2.0/pi
 * \return result with exponent 0.
 */
Word16 getCosWord16R2(Word16 theta);

/****************************************************************************/
/*!
  \brief    16/16->16 unsigned integer division

  x and y have to be positive, x has to be < 16384

  \return 16/16->16 integer
 */

Word16 idiv1616U(Word16 x, Word16 y);


/**
 * \brief return 2 ^ (exp * 2^exp_e)
 * \param exp_m mantissa of the exponent to 2.0f
 * \param exp_e exponent of the exponent to 2.0f
 * \param result_e pointer to a INT where the exponent of the result will be stored into
 * \return mantissa of the result
 */
Word32 BASOP_util_Pow2(
    const Word32 exp_m, const Word16 exp_e,
    Word16 *result_e
);


Word32 Isqrt_lc1(
    Word32 frac,  /* (i)   Q31: normalized value (1.0 < frac <= 0.5) */
    Word16 * exp  /* (i/o)    : exponent (value = frac x 2^exponent) */
);

/*****************************************************************************/
/*!

   \brief Calculates pow(2,x)
  ___________________________________________________________________________
 |                                                                           |
 |   Function Name : Pow2()                                                  |
 |                                                                           |
 |     L_x = pow(2.0, exponant.fraction)         (exponent = interger part)  |
 |         = pow(2.0, 0.fraction) << exponent                                |
 |---------------------------------------------------------------------------|
 |  Algorithm:                                                               |
 |                                                                           |
 |   The function Pow2(L_x) is approximated by a table and linear            |
 |   interpolation.                                                          |
 |                                                                           |
 |   1- i = bit10-b15 of fraction,   0 <= i <= 31                            |
 |   2- a = bit0-b9   of fraction                                            |
 |   3- L_x = table[i]<<16 - (table[i] - table[i+1]) * a * 2                 |
 |   4- L_x = L_x >> (30-exponant)     (with rounding)                       |
 |___________________________________________________________________________|
*/
Word32 Pow2(                              /*!< (o) Q0  : result       (range: 0<=val<=0x7fffffff) */
    Word16 exponant,                      /*!< (i) Q0  : Integer part.      (range: 0<=val<=30)   */
    Word16 fraction                       /*!< (i) Q15 : Fractionnal part.  (range: 0.0<=val<1.0) */
);

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

Word16 Log2_norm_lc ( /* (o) : Fractional part of Log2. (range: 0<=val<1)  */
    Word32 L_x        /* (i) : input value (normalized)                    */
);

/*************************************************************************
 *
 *   FUNCTION:   BASOP_Util_fPow()
 */
/**
 * \brief BASOP_Util_fPow
 *
 *   PURPOSE:   Computes pow(base_m, base_e, exp_m, exp_e), where base_m and base_e
 *              specify the base, and exp_m and exp_e specify the exponent.
 *              The result is returned in a mantissa and exponent representation.
 *
 *   DESCRIPTION:
 *        The function BASOP_Util_fPow(L_x) calculates the power function by
 *        calculating 2 ^ (log2(base)*exp)
 *
 * \param base_m mantissa of base
 * \param base_e exponent of base
 * \param exp_m mantissa of exponent
 * \param exp_e exponent of exponent
 * \param result_e pointer to exponent of result
 * \return Word32 mantissa of result
 *
 *************************************************************************/

Word32 BASOP_Util_fPow(               /* (o) : mantissa of result                                              */
    Word32 base_m, Word16 base_e, /* (i) : input value for base (mantissa and exponent)                    */
    Word32 exp_m, Word16 exp_e,   /* (i) : input value for exponent (mantissa and exponent)                */
    Word16 *result_e              /* (o) : output pointer to exponent of result                            */
);

/*!**********************************************************************
   \brief   Add two values given by mantissa and exponent.

   Mantissas are in 32-bit-fractional format with values between 0 and 1. <br>
   The base for exponents is 2.  Example:  \f$  a = a\_m * 2^{a\_e}  \f$<br>

************************************************************************/
Word32 BASOP_Util_Add_Mant32Exp                  /*!< o: normalized result mantissa */
(Word32   a_m,       /*!< i: Mantissa of 1st operand a  */
 Word16   a_e,       /*!< i: Exponent of 1st operand a  */
 Word32   b_m,       /*!< i: Mantissa of 2nd operand b  */
 Word16   b_e,       /*!< i: Exponent of 2nd operand b  */
 Word16  *ptr_e);    /*!< o: exponent of result         */

/****************************************************************************/
/*!
  \brief   Accumulates multiplications

	Accumulates the elementwise multiplications of Word32 Array X with Word16 Array Y
	pointed to by arg1 and arg2 with specified headroom. Length of to be multiplied arrays is arg3,
    headroom with has to be taken into account is specified in arg4

  \return Word32 result of accumulated multiplications over Word32 array arg1 and Word16 array arg2 and Word16 pointer
          to exponent correction factor which needs to be added to the exponent of the result vector
*/
Word32 dotWord32_16_guards(const Word32 * X, const Word16 * Y, Word16 n, Word16 hr, Word16 * shift);


#endif /* __BASOP_UTIL_H__ */
