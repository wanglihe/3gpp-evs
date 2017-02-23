/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"


/*----------------------------------------------------------------------------------*
 * dec_acelp_2t32()
 *
 * 12 bits algebraic codebook decoder.
 * 2 track x 32 positions per track = 64 samples.
 *
 * 12 bits --> 2 pulses in a frame of 64 samples.
 *
 * All pulses can have two (2) possible amplitudes: +1 or -1.
 * Each pulse can have 32 possible positions.
 *
 * See cod2t32.c for more details of the algebraic code.
 *----------------------------------------------------------------------------------*/

void dec_acelp_2t32(
    Decoder_State *st,      /* i/o: decoder state structure   */
    float code[]    /* o:   algebraic (fixed) codebook excitation */
)
{
    short index, i0, i1;

    index = (short) get_next_indice( st, 12 );

    set_f( code, 0.0f, L_SUBFR );

    /*-----------------------------------------------------------------*
     * decode the positions and signs of pulses and build the codeword
     *-----------------------------------------------------------------*/

    i0 = ((index>>6) & (NB_POS_FCB_2T-1)) * NB_TRACK_FCB_2T;
    i1 = ((index & (NB_POS_FCB_2T-1)) * NB_TRACK_FCB_2T) + 1;
    code[i0] = -1.0f;
    if ((index & 0x800) == 0)
    {
        code[i0] = 1.0f;
    }

    code[i1] = -1.0f;
    if ((index & 0x20) == 0)
    {
        code[i1] = 1.0f;
    }

    return;
}


/*----------------------------------------------------------------------------------*
 * dec_acelp_1t64()
 *
 * 7 bits algebraic codebook.
 * 1 track x 64 positions per track = 64 samples.
 *
 * The pulse can have 64 possible positions and two (2) possible amplitudes: +1 or -1.
 *----------------------------------------------------------------------------------*/

void dec_acelp_1t64(
    Decoder_State *st,       /* i/o: decoder state structure   */
    float code[]    /* o:   algebraic (fixed) codebook excitation */
)
{
    short pos, sgn;

    /*-----------------------------------------------------------------*
     * decode the positions and signs of pulses and build the codeword
     *-----------------------------------------------------------------*/

    pos = (short)get_next_indice( st, 7 );

    sgn = -1;
    if( pos >= L_SUBFR )
    {
        pos -= L_SUBFR;
        sgn = 1;
    }

    set_f( code, 0.0f, L_SUBFR );
    code[pos] = sgn;

    return;
}
