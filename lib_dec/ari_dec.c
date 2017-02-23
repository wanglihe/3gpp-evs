/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdio.h>
#include "cnst.h"
#include "prot.h"
#include "stat_com.h"
#include "assert.h"
#include "basop_util.h"



/*---------------------------------------------------------------
 * Ari decode 14 bits routines
  -------------------------------------------------------------*/

/**
 * \brief 	Start ArCo decoding
 *
 * \param[i/o] st
 * \param[o] s
 *
 * \return bit consumption
 */
void ari_start_decoding_14bits(Decoder_State *st,Tastat *s)
{
    long	val;

    val = get_next_indice(st, cbitsnew);

    s->low = 0;
    s->high = ari_q4new;
    s->vobf=val;
}

/**
 * \brief 	Start ArCo decoding
 *
 * \param[o] prm
 * \param[i] bp
 * \param[o] s
 *
 * \return bit consumption
 */
long ari_start_decoding_14bits_prm(const int *ptr,long bp,Tastat *s)
{
    long	val;
    long	i;
    const int *p;



    val = 0;

    p=ptr+bp;
    for (i=0; i<cbitsnew; i++)
    {
        val = (val<<1) | *(p+i);
    }

    s->low = 0;
    s->high = ari_q4new;
    s->vobf=val;

    bp=bp+i;


    return bp;
}



/**
 * \brief Only for 17 symbols with new extended Tables:
 * based on tri's optimization
 * based on extended Tables which need less branches for coding
 *
 * \param[i/o] st
 * \param[o]   res
 * \param[i/o] s
 * \param[i]   cum_freq
 *
 * \return bit consumption
 */
void ari_decode_14bits_s17_ext(Decoder_State *st,
                               int *res,
                               Tastat *s,
                               const unsigned short *cum_freq)
{
    unsigned long symbol;
    unsigned long low, high, range, value;
    unsigned long cum;
    unsigned short const  *p;




    /* read s->low,high,vobf sequentially */
    low = s->low;
    high = s->high;
    value = s->vobf;

    range = high-low+1;                                             /* keep: tmp=low-1 */
    cum =((((int) (value-low+1))<<stat_bitsnew)-((int) 1));


    p = cum_freq;

    /* Note: For each indirect addressing p[i], we assume a tmp pointer init followed by a costfree reading the value */
    /* If the value multiplied by range is greater than cum, the pointer p is set to the tmp pointer                  */
    /*    tmp_p = p+8; if (tmp_p[0]*range>cum) p = tmp_p;                                                             */

    /*  */
    if(p[8]*range>cum)
    {
        p+= 8;
    }
    /*  */
    if(p[4]*range>cum)
    {
        p+= 4;
    }
    /*  */
    if(p[2]*range>cum)
    {
        p+= 2;
    }
    /*  */
    if(p[1]*range>cum)
    {
        p+= 1;
        if((p==cum_freq+15) && (p[1]*range>cum))
        {
            p+=1;
        }
    }

    symbol = p-cum_freq;

    high  = low + mul_sbc_14bits(range,cum_freq[symbol]) - 1;
    low  += mul_sbc_14bits(range,cum_freq[symbol+1]);

    for (;;)
    {
        if ( high>=ari_q2new )
        {
            if ( low>=ari_q2new )
            {
                value -= ari_q2new;
                low -= ari_q2new;
                high -= ari_q2new;
            }
            else
            {
                if ( low>=ari_q1new && high<ari_q3new )
                {
                    value -= ari_q1new;
                    low -= ari_q1new;
                    high -= ari_q1new;
                }
                else
                {
                    break;
                }
            }
        }
        low  += low;                                                  /*  */
        high += high+1;                                                /* ------- */

        value = (value<<1) | get_next_indice_1(st);                     /*  */
    }

    s->low = low;
    s->high = high;
    s->vobf = value;

    *res=symbol;
}

/**
 * \brief Only for 27 symbols with new extended Tables:
 * based on tri's optimization
 * based on extended Tables which need less branches for coding
 * copied from ari_decode_14bits_s17_ext, with changes marked
 *
 * \param[i/o] st
 * \param[o]   res
 * \param[i/o] s
 * \param[i]   cum_freq
 */
void ari_decode_14bits_s27_ext(Decoder_State *st,
                               int *res,
                               Tastat *s,
                               const unsigned short *cum_freq)
{
    unsigned long symbol;
    unsigned long low, high, range, value;
    unsigned long cum;
    unsigned long il, ih, im;




    /* read s->low,high,vobf sequentially */
    low = s->low;
    high = s->high;
    value = s->vobf;

    range = high-low+1;                                             /* keep: tmp=low-1 */
    cum =((((int) (value-low+1))<<stat_bitsnew)-((int) 1));


    /* Note: For each indirect addressing p[i], we assume a tmp pointer init followed by a costfree reading the value */
    /* If the value multiplied by range is greater than cum, the pointer p is set to the tmp pointer                  */
    /*    tmp_p = p+8; if (tmp_p[0]*range>cum) p = tmp_p;                                                             */

    /* begin change when compared with ari_decode_14bits_s17_ext,
       starting with line: if (p[8] * range > cum) { */
    il = 0;
    ih = 27;

    /* do a five step binary search, using the interval [il, ih) */
    im = 13; /* (il + ih) >> 1 */
    if (cum_freq[im] * range > cum)
    {
        il = im;
    }
    else
    {
        ih = im;
    }

    im = (il + ih) >> 1;
    if (cum_freq[im] * range > cum)
    {
        il = im;
    }
    else
    {
        ih = im;
    }

    im = (il + ih) >> 1;
    if (cum_freq[im] * range > cum)
    {
        il = im;
    }
    else
    {
        ih = im;
    }

    im = (il + ih) >> 1;
    if (cum_freq[im] * range > cum)
    {
        il = im;
    }
    else
    {
        ih = im;
    }

    if (ih - il > 1)   /* if the interval has more than one symbol */
    {
        /* here, only ih == il + 2 is possible, which means two symbols in the interval */
        im = il + 1; /* (il + ih) >> 1 */
        if (cum_freq[im] * range > cum)
        {
            il = im;
        }
    }

    symbol = il;
    /* end change when compared with ari_decode_14bits_s17_ext,
       ending with line: symbol = p - cum_freq; */

    high  = low + mul_sbc_14bits(range,cum_freq[symbol]) - 1;
    low  += mul_sbc_14bits(range,cum_freq[symbol+1]);

    for (;;)
    {
        if ( high>=ari_q2new )
        {
            if ( low>=ari_q2new )
            {
                value -= ari_q2new;
                low -= ari_q2new;
                high -= ari_q2new;
            }
            else
            {
                if ( low>=ari_q1new && high<ari_q3new )
                {
                    value -= ari_q1new;
                    low -= ari_q1new;
                    high -= ari_q1new;
                }
                else
                {
                    break;
                }
            }
        }
        low  += low;                                                  /*  */
        high += high+1;                                                /* ------- */

        value = (value<<1) | get_next_indice_1(st);                     /*  */
    }

    s->low = low;
    s->high = high;
    s->vobf = value;

    *res=symbol;
}

/**
 * \brief Only for decoding one bit with uniform probability:
 * based on tri's optimization
 * copied from ari_decode_14bits_s17_ext, with changes marked
 * the equivalent cum_freq table used is {16384, 8192, 0}
 *
 * \param[i/o] st
 * \param[o]   res
 * \param[i/o] s
 */
void ari_decode_14bits_bit_ext(Decoder_State *st,
                               int *res,
                               Tastat *s
                              )
{
    unsigned long symbol;
    unsigned long low, high, range, value;
    unsigned long cum;




    /* read s->low,high,vobf sequentially */
    low = s->low;
    high = s->high;
    value = s->vobf;

    range = high-low+1;                                             /* keep: tmp=low-1 */
    cum =((((int) (value-low+1))<<stat_bitsnew)-((int) 1));


    /* Note: For each indirect addressing p[i], we assume a tmp pointer init followed by a costfree reading the value */
    /* If the value multiplied by range is greater than cum, the pointer p is set to the tmp pointer                  */
    /*    tmp_p = p+8; if (tmp_p[0]*range>cum) p = tmp_p;                                                             */

    /* begin change when compared with ari_decode_14bits_s17_ext,
       starting with line: if (p[8] * range > cum) { */
    symbol = 0;

    if ((range << 13) > cum)
    {
        symbol = 1;
    }

    if (symbol == 0)
    {
        /* high is unchanged */
        low = low + (range >> 1);
    }
    else
    {
        high = low + (range >> 1) - 1;
        /* low is unchanged */
    }
    /* end change when compared with ari_decode_14bits_s17_ext,
       ending with line: low += mul_sbc_14bits(range, cum_freq[symbol + 1]); */


    for (;;)
    {
        if ( high>=ari_q2new )
        {
            if ( low>=ari_q2new )
            {
                value -= ari_q2new;
                low -= ari_q2new;
                high -= ari_q2new;
            }
            else
            {
                if ( low>=ari_q1new && high<ari_q3new )
                {
                    value -= ari_q1new;
                    low -= ari_q1new;
                    high -= ari_q1new;
                }
                else
                {
                    break;
                }
            }
        }
        low  += low;                                                  /*  */
        high += high+1;                                                /* ------- */

        value = (value<<1) | get_next_indice_1(st);                     /*  */
    }

    s->low = low;
    s->high = high;
    s->vobf = value;

    *res=symbol;
}

/*------------------------------------------------------------------------
 * Function: ari_decode_14bits_pow
 *
 * Decode a symbol which follows the exponential distribution. That is,
 * symbols are in the following intervals
 *
 * p(x = 0) = 1 - exp(- 0.5 * base * 2)
 * p(x = q>0) = exp(- (q-0.5)*base* 2) - exp(- (q+0.5)*base*2 )
 *
 *-------------------------------------------------------------------------*/
long ari_decode_14bits_pow(const int *ptr, long bp, long bits, int *res, Tastat *s,	unsigned base)
{
    unsigned long symbol;
    unsigned long low, high, range, value;
    unsigned long cum;
    Word16 pows[12];	/* "base" to the power of 2, 4, 8,... 2^12 */
    Word16 lowlim, highlim=0, testval;
    int k;



    low = s->low;
    high = s->high + 1;
    value = s->vobf;
    lowlim = 0;
    symbol = 0;

    range = high - low;

    /* the value read from bit-stream */
    assert(value >= low);
    cum = ((((int) (value - low)) << stat_bitsnew) + (1 << stat_bitsnew) - 1);

    /* search for the interval where "cum" fits */
    if (((base >> 1) * range) > cum)	/* below pow-1 */
    {
        pows[0] = testval = base;
        /* increase exponent until it is smaller than "cum" */
        for (k=1; k<12; k++)
        {
            highlim = testval;
            pows[k] = mult_r(pows[k-1], pows[k-1]);
            move16();
            testval = mult_r(pows[k], base);
            if (((testval >> 1) * range) <= cum) /* found! big range is [lowlim,testval], (now narrow it down) */
            {
                lowlim = testval;
                k--;
                symbol = 1<<k;
                break;
            }
        }
        assert(k<12);	/* maximum 2^10-1*/
        /* narrow the range down */
        for (k--; k>0; k--)
        {
            testval = mult_r(highlim, pows[k]);
            if (((testval >> 1) * range) <= cum)
            {
                lowlim = testval;
                symbol -= 1<<(k-1);
            }
            else
            {
                highlim = testval;
            }


        }
        highlim >>= 1;
        lowlim >>= 1;
    }
    else     /* trivial case, above pow-1, that is, first symbol */
    {
        symbol = 0;
        lowlim = base >> 1;
        highlim = 16384;
    }


    high = low + mul_sbc_14bits(range, highlim);

    low += mul_sbc_14bits(range, lowlim);

    /*ptr init for ptr*/
    for (; bp < bits;)
    {
        if (high > ari_q2new)
        {
            if (low >= ari_q2new)
            {
                value -= ari_q2new;
                low   -= ari_q2new;
                high  -= ari_q2new;
            }
            else
            {
                if (low >= ari_q1new && high <= ari_q3new)
                {
                    value -= ari_q1new;
                    low   -= ari_q1new;
                    high  -= ari_q1new;
                }
                else
                {
                    break;
                }
            }
        }
        low  += low;
        high += high;

        value = (value << 1) | ptr[bp++];
    }

    if( !(bp != bits || ! ((s->low == (int) low) && (s->high == (int) high) && (s->vobf == (int) value)) ) )
    {
        /* This should not happen except of bit errors. */
        s->high = s->low = 0;
        *res = 0;
        return -1;
    }

    s->low = low;
    s->high = high - 1;
    s->vobf = value;

    *res = symbol;


    return bp;
}



/*------------------------------------------------------------------------
 * Function: ari_decode_14bits_sign
 *
 * Decode a sign with equal probabilities.
 *-------------------------------------------------------------------------*/
long ari_decode_14bits_sign(const int *ptr, long bp, long bits, int *res, Tastat *s)
{
    unsigned long symbol;
    unsigned long low, high, range, value;
    unsigned long cum;



    low = s->low;
    high = s->high + 1;
    value = s->vobf;

    range = high - low;

    if (bp < bits)
    {
        assert(value >= low);
        cum = ((((int) (value - low)) << stat_bitsnew) + (1 << stat_bitsnew) - 1);
        if (8192 * range > cum)
        {
            symbol = 2;
            high = low + (range >> 1);
        }
        else
        {
            symbol = 1;
            low += range >> 1;
        }

        /*ptr init for ptr*/
        for (; bp < bits;)
        {
            if (high > ari_q2new)
            {
                if (low >= ari_q2new)
                {
                    value -= ari_q2new;
                    low   -= ari_q2new;
                    high  -= ari_q2new;
                }
                else
                {
                    if (low >= ari_q1new && high <= ari_q3new)
                    {
                        value -= ari_q1new;
                        low   -= ari_q1new;
                        high  -= ari_q1new;
                    }
                    else
                    {
                        break;
                    }
                }
            }
            low  += low;
            high += high;

            value = (value << 1) | ptr[bp++];
        }
    }
    else
    {
        cum = value-low;
        range >>= 1;
        if (range > cum)
        {
            symbol = 2;
            high = low + range;
        }
        else
        {
            symbol = 1;
            low += range;
        }
    }

    s->low = low;
    s->high = high - 1;
    s->vobf = value;

    *res = symbol;


    return bp;
}

