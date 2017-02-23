/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * Local functions
 *-------------------------------------------------------------------*/

static void rc_enc_shift(Encoder_State *st);
static void rc_enc_write(Encoder_State *st, short byte, short bits);


/*-------------------------------------------------------------------*
 * rc_enc_init()
 *
 *  Initalize range coder
 *-------------------------------------------------------------------*/

void rc_enc_init(
    Encoder_State *st,          /* i/o: Encoder state       */
    short tot_bits              /* i  : Total bit budget    */
)
{
    st->rc_low = 0;
    st->rc_range = 0xffffffff;
    st->rc_cache = -1;
    st->rc_carry = 0;
    st->rc_carry_count = 0;
    st->rc_num_bits = 0;
    st->rc_tot_bits = tot_bits;
    st->rc_offset = 0;

    return;
}

/*-------------------------------------------------------------------*
 * rc_encode()
 *
 *  Encode symbol with range coder
 *-------------------------------------------------------------------*/

void rc_encode(
    Encoder_State *st,              /* i/o: Encoder state                       */
    unsigned int cum_freq,          /* i  : Cumulative frequency up to symbol   */
    unsigned int sym_freq,          /* i  : Symbol probability                  */
    unsigned int tot                /* i  : Total cumulative frequency          */
)
{
    unsigned int r, tmp;
    unsigned int inv_tot;
    short exp;

    inv_tot = UL_inverse(tot, &exp);
    tmp = UMult_32_32(st->rc_range, inv_tot);
    r = tmp >> (exp - 32);
    tmp = r*cum_freq;

    st->rc_low = st->rc_low + tmp;
    if (st->rc_low < tmp)
    {
        st->rc_carry = 1;
    }

    st->rc_range = r*sym_freq;

    while (st->rc_range < 1<<24)
    {
        st->rc_range = st->rc_range << 8;
        st->rc_num_bits += 8;
        rc_enc_shift(st);
    }

    return;
}

/*-------------------------------------------------------------------*
 * rc_enc_finish()
 *
 *  Finalize range coder
 *-------------------------------------------------------------------*/

void rc_enc_finish(
    Encoder_State *st           /* i/o: Encoder state       */
)
{
    unsigned int val, mask, high;
    short bits, over1, over2;

    bits = norm_ul(st->rc_range) + 1;  /* aligned to BASOP */

    mask = 0xffffffff >> bits;

    val = st->rc_low + mask;
    high = st->rc_low + st->rc_range;

    over1 = val < st->rc_low;
    over2 = high < st->rc_low;

    val = val & ~mask;

    if ( !(over1 ^ over2) )
    {
        if ((val + mask) >= high)
        {
            bits++;
            mask >>= 1;
            val = (st->rc_low + mask) & ~mask;
        }

        if (val < st->rc_low)
        {
            st->rc_carry = 1;
        }
    }

    st->rc_low = val;

    if (bits > st->rc_tot_bits - st->rc_num_bits)
    {
        bits = st->rc_tot_bits - st->rc_num_bits;

    }

    st->rc_num_bits += bits;
    while (bits > 0)
    {
        rc_enc_shift(st);
        bits -= 8;
    }

    bits += 8;

    if ( st->rc_carry_count > 0 )
    {
        rc_enc_write(st, st->rc_cache + st->rc_carry, 8);

        while (st->rc_carry_count > 1)
        {
            rc_enc_write(st, (st->rc_carry + 0xff), 8);
            st->rc_carry_count--;
        }
        rc_enc_write(st, (st->rc_carry + 0xff) & ((1<<bits) - 1), bits);
    }
    else
    {
        rc_enc_write(st, (st->rc_cache + st->rc_carry)>>(8-bits), bits);
    }

    bits = st->rc_num_bits;
    while (bits < st->rc_tot_bits-16)
    {
        rc_enc_write(st, 0, 16);
        bits += 16;
    }

    bits = st->rc_tot_bits - bits;
    if (bits > 0)
    {
        rc_enc_write(st, 0, bits);
    }

    return;
}

/*-------------------------------------------------------------------*
 * rc_enc_shift()
 *
 * Shift a byte out to bitstream
 *-------------------------------------------------------------------*/

static void rc_enc_shift(
    Encoder_State *st           /* i/o: Encoder state       */
)
{
    if ( st->rc_low < (0xffUL << 24) || st->rc_carry )
    {
        if (st->rc_cache >= 0)
        {
            rc_enc_write(st, st->rc_cache + st->rc_carry, 8);
        }

        while (st->rc_carry_count > 0)
        {
            rc_enc_write(st, (st->rc_carry + 0xff) & 255, 8);
            st->rc_carry_count--;
        }

        st->rc_cache = st->rc_low >> 24;
        st->rc_carry = 0;
    }
    else
    {
        st->rc_carry_count++;
    }

    st->rc_low = st->rc_low << 8;

    return;
}

/*-------------------------------------------------------------------*
 * rc_enc_bits()
 *
 *
 *-------------------------------------------------------------------*/

void rc_enc_bits(
    Encoder_State *st,          /* i/o: Encoder state       */
    unsigned int value,         /* i  : Value to encode     */
    short bits                  /* i  : Number of bits used */
)
{
    if ( rc_get_bits2(st->rc_num_bits, st->rc_range) + bits <= st->rc_tot_bits)
    {
        st->rc_num_bits += bits;
        if ( bits > 16 )
        {
            push_indice(st, IND_RC_END - st->rc_offset, value >> 16, bits - 16);
            st->rc_offset++;

            push_indice(st, IND_RC_END - st->rc_offset, value & ((1 << 16) - 1), 16);
            st->rc_offset++;
        }
        else
        {
            push_indice(st, IND_RC_END - st->rc_offset++, value, bits);
        }
    }
    else
    {
    }

    return;
}

/*-------------------------------------------------------------------*
 * rc_enc_uniform()
 *
 * Encode with uniform distribution
 *-------------------------------------------------------------------*/

void rc_enc_uniform(
    Encoder_State *st,          /* i/o: Encoder state       */
    unsigned int value,         /* i  : Value to encode     */
    unsigned int tot            /* i  : Maximum value       */
)
{
    short n;

    n = 32 - norm_ul(tot-1);   /* aligned to BASOP */

    if (n <= 8)
    {
        rc_encode(st, value, 1, tot);
    }
    else
    {
        n -= 8;
        rc_encode(st, value >> n, 1, (tot >> n) + 1);
        rc_enc_bits(st, value & ((1<<n)-1), n);
    }

    return;
}

/*-------------------------------------------------------------------*
 * rc_enc_write()
 *
 *  Write a byte to bitstream
 *-------------------------------------------------------------------*/

static void rc_enc_write(
    Encoder_State *st,        /* i/o: Encoder state         */
    short byte,               /* i  : Byte to write         */
    short bits                /* i  : Number of bits        */
)
{
    push_indice(st, IND_RC_START, byte, bits);

    return;
}

