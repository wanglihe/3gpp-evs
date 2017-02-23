/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*-----------------------------------------------------------------*
 * Transform domain contribution decoding
 *-----------------------------------------------------------------*/

void transf_cdbk_dec(
    Decoder_State *st,            /* i/o: decoder state structure   */
    const long  core_brate,     /* i  : core bitrate                                    */
    const short coder_type,     /* i  : coding type                                     */
    const short harm_flag_acelp,/* i  : harmonic flag for higher rates ACELP            */
    const short i_subfr,        /* i  : subframe index                                  */
    const short tc_subfr,       /* i  : TC subframe index                               */
    const float Es_pred,        /* i  : predicited scaled innovation energy             */
    const float gain_code,      /* i  : innovative excitation gain                      */
    float *mem_preemp,    /* i/o: dequantizer preemhasis memory                   */
    float *gain_preQ,     /* o  : prequantizer excitation gain                    */
    float *norm_gain_preQ,/* o  : normalized prequantizer excitation gain         */
    float code_preQ[],    /* o  : prequantizer excitation                         */
    short *unbits         /* o  : number of AVQ unused bits                       */
)
{
    short i, index, nBits, Nsv;
    int x_norm[L_SUBFR];
    float Ecode;
    short nq[L_SUBFR/WIDTH_BAND];

    /*--------------------------------------------------------------*
     * Set bit-allocation
     *--------------------------------------------------------------*/

    Nsv = 8;
    nBits = AVQ_bits_16kHz_tbl[BIT_ALLOC_IDX_16KHZ(core_brate, coder_type, i_subfr, TC_SUBFR2IDX_16KHZ(tc_subfr))];

    /* increase # of AVQ allocated bits by unused bits from the previous subframe */
    nBits += (*unbits);

    /*--------------------------------------------------------------*
     * Dequantize prequantizer excitation gain
     *--------------------------------------------------------------*/

    index = (short)get_next_indice( st, G_AVQ_BITS );

    if( coder_type == INACTIVE )
    {
        if( core_brate == ACELP_64k )
        {
            *gain_preQ = usdequant( index, G_AVQ_MIN_INACT_64k, G_AVQ_DELTA_INACT_64k );
        }
        else if( core_brate == ACELP_48k )
        {
            *gain_preQ = usdequant( index, G_AVQ_MIN_INACT_48k, G_AVQ_DELTA_INACT_48k );
        }
        else
        {
            *gain_preQ = usdequant( index, G_AVQ_MIN_INACT, G_AVQ_DELTA_INACT );
        }

        *gain_preQ *= gain_code;
    }
    else
    {
        if( core_brate > ACELP_24k40 && core_brate <= ACELP_32k )
        {
            *gain_preQ = gain_dequant( index, 0.1f*G_AVQ_MIN, G_AVQ_MAX, G_AVQ_BITS );
        }
        else
        {
            *gain_preQ = gain_dequant( index, G_AVQ_MIN, G_AVQ_MAX, G_AVQ_BITS );
        }
        if( Es_pred < 0  )
        {
            *gain_preQ *= 0.25f * fabs(Es_pred);
        }
        else
        {
            *gain_preQ *= Es_pred;
        }
    }

    /*--------------------------------------------------------------*
     * Encode and multiplex subvectors into bit-stream
     *--------------------------------------------------------------*/

    AVQ_demuxdec( st, x_norm, &nBits, Nsv, nq );

    /* save # of AVQ unused bits for next subframe */
    *unbits = nBits;

    /*--------------------------------------------------------------*
     * iDCT transform
     *--------------------------------------------------------------*/

    set_f( code_preQ, 0.0f, L_SUBFR );
    for( i=0; i<Nsv*WIDTH_BAND; i++ )
    {
        code_preQ[i] = (float) (x_norm[i]);
    }

    if( coder_type == INACTIVE || core_brate > ACELP_32k || harm_flag_acelp )
    {
        edct2( L_SUBFR, 1, code_preQ, code_preQ, ip_edct2_64, w_edct2_64 );
    }

    /*--------------------------------------------------------------*
     * Preemphasise
     *--------------------------------------------------------------*/

    /* in extreme cases at subframe boundaries, lower the preemphasis memory to avoid a saturation */
    if( (nq[7] != 0) && (st->last_nq_preQ - nq[0] > 7) )
    {
        *mem_preemp /= 16;
    }

    st->last_nq_preQ = nq[7];

    preemph( code_preQ, FAC_PRE_AVQ, L_SUBFR, mem_preemp );

    /*--------------------------------------------------------------*
     * Compute normalized prequantizer excitation gain for FEC
     *--------------------------------------------------------------*/

    Ecode = (sum2_f( code_preQ, L_SUBFR ) + 0.01f) / L_SUBFR;

    /* somewhat attenuate pre-quantizer normalized gain for FEC */
    *norm_gain_preQ = 0.8f * (*gain_preQ) * (float)sqrt( Ecode );

    st->use_acelp_preq = 1;

    return;

}

/*-----------------------------------------------------------*
 * gain_dequant()
 *
 * Returns decoded gain quantized between the specified
 * range using the specified number of levels.
 *-----------------------------------------------------------*/

float gain_dequant(             /* o:   decoded gain                  */
    short index,          /* i:   quantization index            */
    const float min,            /* i:   value of lower limit          */
    const float max,            /* i:   value of upper limit          */
    const short bits            /* i:   number of bits to dequantize  */
)
{
    float gain,  c_min, c_mult;
    short levels;

    levels = 1<<bits;

    c_min = (float)log10(min);
    c_mult = (float) ((levels-1)/(log10(max)-c_min));

    gain = (float)pow( 10.0, (((float)index)/c_mult) + c_min );

    return( gain );

}
