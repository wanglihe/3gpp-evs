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




/**
 * \brief 	Copy state
 *
 * \param[i] source
 * \param[o] dest
 *
 * \return none
 */
void ari_copy_states(Tastat *source, Tastat *dest)
{
    dest->low  = source->low;
    dest->high = source->high;
    dest->vobf = source->vobf;
}

/*---------------------------------------------------------------
  Ari encoder 14 bits routines
  -------------------------------------------------------------*/

/**
 * \brief 	Start ArCo encoding
 *
 * \param[o] s
 *
 * \return none
 */
void ari_start_encoding_14bits(Tastat *s)
{
    /* : addressing is made with walking pointer s */
    s->low  = 0;
    s->high = ari_q4new;
    s->vobf = 0;
}

/**
 * \brief 	Finish ArCo encoding
 *
 * \param[o] ptr
 * \param[i] bp
 * \param[i] s
 *
 * \return bit consumption
 */
long ari_done_encoding_14bits(int *ptr,long bp,Tastat *s)
{
    long	low;
    long	bits_to_follow;




    /*  not needed, s points to s->low */
    low = s->low;
    bits_to_follow = s->vobf+1;

    if ( low < ari_q1new )
    {
        ptr[bp++] = 0; /*send a zero*/
        for(; bits_to_follow>0; bits_to_follow--)
        {
            ptr[bp++] = 1; /*send a one*/
        }
    }
    else
    {
        ptr[bp++] = 1;  /*send a one*/
        for(; bits_to_follow>0; bits_to_follow--)
        {
            ptr[bp++] = 0;  /*send a zero*/
        }
    }

    /*It is done so no need to save values-> no counting*/
    /*s->low = low;
    s->vobf = bits_to_follow;*/


    return bp;
}


/**
 * \brief encode function for extended proba tables: less branches needed for coding
 *
 * \param[o]   ptr
 * \param[i]   bp
 * \param[i/o] s
 * \param[i]   symbol
 * \param[i]   cum_freq
 *
 * \return bit consumption
 */
long ari_encode_14bits_ext(int *ptr,long bp,Tastat *s,long symbol,const unsigned short *cum_freq)
{
    long low, high, range;
    long  bits_to_follow;



    /*for all operation using bit_ptr=&ptr[bp]   */
    /* for reading s->high,low,vobf sequentially */
    high=s->high;
    low =s->low;
    range = high-low+1;

    high  = low + mul_sbc_14bits(range,cum_freq[symbol]) - 1;
    low  += mul_sbc_14bits(range,cum_freq[symbol+1]);

    bits_to_follow = s->vobf;

    for (;;)
    {
        if ( high<ari_q2new )
        {
            ptr[bp++] = 0; /*send a zero*/
            for(; bits_to_follow>0; bits_to_follow--)
            {
                ptr[bp++] = 1; /*send a one*/
            }
        }
        else
        {
            if ( low>=ari_q2new )
            {
                ptr[bp++] = 1; /*send a one*/
                for(; bits_to_follow>0; bits_to_follow--)
                {
                    ptr[bp++] = 0; /*send a zero*/
                }
                low -= ari_q2new;
                high -= ari_q2new; /* Subtract offset to top.  */
            }
            else
            {
                /* Output an opposite bit   */
                if ( low>=ari_q1new && high<ari_q3new )   /* Output an opposite bit   */
                {
                    /* later if in middle half. */
                    bits_to_follow += 1;
                    low -= ari_q1new;    /* Subtract offset to middle*/
                    high -= ari_q1new;
                }
                else
                {
                    break;  /* Otherwise exit loop. */
                }
            }
        }
        low  += low;
        high += high+1;                    /* Scale up code range.     */
    }

    s->low  = low;
    s->high = high;
    s->vobf = bits_to_follow;


    return bp;
}

/*------------------------------------------------------------------------
 * Function: ari_encode_14bits_range
 *
 * Encode an cumulative frequency interval.
 *-------------------------------------------------------------------------*/

long ari_encode_14bits_range(int *ptr, long bp, long bits, Tastat *s, unsigned short cum_freq_low, unsigned short cum_freq_high)
{
    long low, high, range;
    long bits_to_follow;



    /*  not needed, s points to s->low */
    high = s->high;
    high++;
    low = s->low;
    range = high - low;

    high = low + mul_sbc_14bits(range, cum_freq_high);
    low += mul_sbc_14bits(range, cum_freq_low);

    bits_to_follow = s->vobf;

    /* while there are more than 16 bits left */
    for (; bp + 16 + bits_to_follow - bits < 0;)
    {
        if (high <= ari_q2new)
        {
            ptr[bp++] = 0; /*send a zero*/
            for (; bits_to_follow > 0; bits_to_follow--)
            {
                ptr[bp++] = 1; /*send a one*/
            }
        }
        else if (low >= ari_q2new)
        {
            /* to reach this branch */
            ptr[bp++] = 1; /*send a one*/
            for (; bits_to_follow > 0; bits_to_follow--)
            {
                ptr[bp++] = 0; /*send a zero*/
            }
            low -= ari_q2new;
            high -= ari_q2new; /* Subtract offset to top.  */
        }
        else if (low >= ari_q1new && high <= ari_q3new)
        {
            /* to reach this branch */
            /* Output an opposite bit   */
            /* later if in middle half. */
            bits_to_follow += 1;
            low -= ari_q1new; /* Subtract offset to middle*/
            high -= ari_q1new;
        }
        else
        {
            /* to reach this branch */
            break; /* Otherwise exit loop.     */
        }

        low += low;
        high += high; /* Scale up code range.     */
    }
    /* if there are <= 16 bits left */
    if (bp + 16 + bits_to_follow - bits >= 0)
    {
        /* No need to do anyhing, but let's keep a place for a breakpoint */
        s->vobf = -1;
    }

    s->low = low;
    s->high = high - 1;
    s->vobf = bits_to_follow;

    return bp;
}

/*------------------------------------------------------------------------
 * Function: ari_encode_14bits_sign
 *
 * Encode a sign with equal probabilities.
 *-------------------------------------------------------------------------*/
long ari_encode_14bits_sign(int *ptr, long bp, long bits, Tastat *s, long sign)
{
    long low, high, range;
    long bits_to_follow;



    /*  not needed, s points to s->low */
    high = s->high;
    high++;
    low = s->low;
    range = high - low;

    if (sign)
    {
        high = low + (range >> 1);
    }
    else
    {
        low += range >> 1;
    }

    bits_to_follow = s->vobf;

    /* while there are more than 16 bits left */
    for (; bp + 16 + bits_to_follow - bits < 0;)
    {
        if (high <= ari_q2new)
        {
            ptr[bp++] = 0; /*send a zero*/
            for (; bits_to_follow > 0; bits_to_follow--)
            {
                ptr[bp++] = 1; /*send a one*/
            }
        }
        else if (low >= ari_q2new)
        {
            /* to reach this branch */
            ptr[bp++] = 1; /*send a one*/
            for (; bits_to_follow > 0; bits_to_follow--)
            {
                ptr[bp++] = 0; /*send a zero*/
            }
            low -= ari_q2new;
            high -= ari_q2new; /* Subtract offset to top.  */
        }
        else if (low >= ari_q1new && high <= ari_q3new)
        {
            /* to reach this branch */
            /* Output an opposite bit   */
            /* later if in middle half. */
            bits_to_follow += 1;
            low -= ari_q1new; /* Subtract offset to middle*/
            high -= ari_q1new;
        }
        else
        {
            /* to reach this branch */
            break; /* Otherwise exit loop.     */
        }

        low += low;
        high += high; /* Scale up code range.     */
    }

    s->low = low;
    s->high = high - 1;
    s->vobf = bits_to_follow;

    return bp;
}

/*------------------------------------------------------------------------
 * Function: ari_done_cbr_encoding_14bits
 *
 * Finish up encoding in CBR mode.
 *-------------------------------------------------------------------------*/
long ari_done_cbr_encoding_14bits(int *ptr, long bp, long bits, Tastat *s)
{

    long high;
    long bits_to_follow;
    int k;


    while (bits - bp - 16 - s->vobf > 0)
    {
        bp = ari_encode_14bits_sign(ptr, bp, bits, s, 0);
    }

    /*  not needed, s points to s->low */
    high = s->high;
    bits_to_follow = s->vobf;

    if (bits_to_follow)
    {
        /* If in upper half, then output a one, bits_to_follow zeros, and the remaining bits, except the first one */
        if (high < 0x8000)
        {
            ptr[bp++] = 0; /*send a zero*/
            for (; bits_to_follow > 0; bits_to_follow--)
            {
                ptr[bp++] = 1; /*send a one*/
            }
        }
        else
        {
            ptr[bp++] = 1; /*send a one*/
            for (; bits_to_follow > 0; bits_to_follow--)
            {
                ptr[bp++] = 0; /*send a zero*/
            }
        }
        /* write remaining bits */
        for (k = 0x4000; k > 0; k >>= 1)
        {
            if (k & high)
            {
                ptr[bp++] = 1; /*send a one*/
            }
            else
            {
                ptr[bp++] = 0; /*send a zero*/
            }
            if (bp >= bits)
                break;
        }

    }
    else
    {
        /* no carry-bits, just write all bits */
        for (k = 0x8000; k > 0; k >>= 1)
        {
            if (k & high)
            {
                ptr[bp++] = 1; /*send a one*/
            }
            else
            {
                ptr[bp++] = 0; /*send a zero*/
            }
            if (bp >= bits)
                break;
        }

    }


    return bp;
}

