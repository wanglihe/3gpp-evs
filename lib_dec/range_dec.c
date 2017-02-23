/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


static short rc_dec_read( Decoder_State *st );

/*-------------------------------------------------------------------*
 * rc_dec_init()
 *
 * Initialize range coder
 *-------------------------------------------------------------------*/

void rc_dec_init(
    Decoder_State *st,          /* i/o: Decoder State       */
    short tot_bits              /* i  : Total bit budget    */
)
{
    short i;

    st->rc_low = 0;
    st->rc_range = 0xffffffff;
    st->rc_num_bits = 0;
    st->rc_offset = tot_bits + st->next_bit_pos;
    st->rc_end = st->rc_offset;

    for (i = 0; i < 4; i++)
    {
        st->rc_low = (st->rc_low << 8) + rc_dec_read(st);
    }

    return;
}

/*-------------------------------------------------------------------*
 * rc_decode()
 *
 *  Decode symbol
 *-------------------------------------------------------------------*/

unsigned int rc_decode(         /* o  : Decoded cumulative frequency    */
    Decoder_State *st,          /* i/o: Decoder State                   */
    unsigned int tot            /* i  : Total cumulative frequency      */
)
{
    unsigned int inv, val;
    short exp;

    inv = UL_inverse(tot, &exp);
    st->rc_help = UMult_32_32(st->rc_range, inv);
    st->rc_help = st->rc_help >> (exp - 32);

    /* safety check in case of bit errors */
    val = st->rc_low/st->rc_help;
    if (val > tot)
    {
        st->BER_detect = 1;
        return 0;
    }
    return val;
}

/*-------------------------------------------------------------------*
 * rc_dec_update()
 *
 *  Update range coder
 *-------------------------------------------------------------------*/

void rc_dec_update(
    Decoder_State *st,              /* i/o: Decoder State           */
    unsigned int cum_freq,          /* i  : Cumulative frequency    */
    unsigned int sym_freq           /* i  : Symbol frequency        */
)
{
    st->rc_low = st->rc_low - cum_freq*st->rc_help;
    st->rc_range = st->rc_help*sym_freq;

    while (st->rc_range < (1<<24))
    {
        st->rc_num_bits += 8;
        st->rc_low = (st->rc_low << 8) + rc_dec_read(st);
        st->rc_range <<= 8;
    }

    return;
}

/*-------------------------------------------------------------------*
 * rc_dec_bits()
 *
 *  Encode bits
 *-------------------------------------------------------------------*/

unsigned int rc_dec_bits(       /* i  : Decoded value   */
    Decoder_State *st,          /* i/o: Decoder State   */
    short bits                  /* i  : Number of bits  */
)
{
    unsigned int value;

    st->rc_num_bits += bits;

    if (bits > 16)
    {
        st->rc_offset -= bits - 16;
        value = get_indice(st, st->rc_offset, bits - 16) << 16;
        st->rc_offset -= 16;
        value |= get_indice(st, st->rc_offset, 16);
    }
    else
    {
        st->rc_offset -= bits;
        value = get_indice(st, st->rc_offset, bits);

    }

    return value;
}

/*-------------------------------------------------------------------*
 * rc_dec_uniform()
 *
 * Encode with uniform distribution
 *-------------------------------------------------------------------*/

unsigned int rc_dec_uniform(    /* i  : Decoded value   */
    Decoder_State *st,          /* i/o: Decoder State   */
    unsigned int tot            /* i  : Maximum value   */
)
{
    unsigned int value;
    short n;
    n = 32 - norm_ul(tot - 1);   /* aligned to BASOP */
    if (n <= 8)
    {
        value = rc_decode(st, tot);
        rc_dec_update(st, value, 1);
    }
    else
    {
        n -= 8;
        value = rc_decode(st, (tot >> n) + 1);
        rc_dec_update(st, value, 1);
        value <<= n;
        value |= rc_dec_bits(st, n);
    }

    return value;
}

/*-------------------------------------------------------------------*
 * rc_dec_finish()
 *
 *  Finalize range decoder
 *-------------------------------------------------------------------*/

void rc_dec_finish(Decoder_State *st)
{
    st->next_bit_pos = st->rc_end;
}


/*-------------------------------------------------------------------*
 * rc_dec_read()
 *
 *  Read a byte from bit stream
 *-------------------------------------------------------------------*/

static short rc_dec_read(Decoder_State *st)
{
    return (short)get_next_indice(st, 8);
}

