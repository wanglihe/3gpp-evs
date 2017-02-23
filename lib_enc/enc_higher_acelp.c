/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*---------------------------------------------------------------------*
 * Local functions
 *---------------------------------------------------------------------*/

static void find_cn( const float xn[], const float Ap[], const float *p_Aq, float cn[] );

/*-----------------------------------------------------------------*
 * Transform domain contribution encoding
 *-----------------------------------------------------------------*/

void transf_cdbk_enc(
    Encoder_State *st,            /* i/o: encoder state structure      */
    const long  core_brate,     /* i  : core bitrate                                    */
    const short extl,           /* i  : extension layer                                 */
    const short coder_type,     /* i  : coding type                                     */
    const short harm_flag_acelp,/* i  : harmonic flag for higher rates ACELP            */
    const short i_subfr,        /* i  : subframe index                                  */
    const short tc_subfr,       /* i  : TC subframe index                               */
    float cn[],           /* i/o: target vector in residual domain                */
    float exc[],          /* i/o: pointer to excitation signal frame              */
    const float *p_Aq,          /* i  : 12k8 Lp coefficient                             */
    const float Ap[],           /* i  : weighted LP filter coefficients                 */
    const float h1[],           /* i  : weighted filter input response                  */
    float xn[],           /* i/o: target vector                                   */
    float xn2[],          /* i/o: target vector for innovation search             */
    float y1[],           /* i/o: zero-memory filtered adaptive excitation        */
    const float y2[],           /* i  : zero-memory filtered innovative excitation      */
    const float Es_pred,        /* i  : predicited scaled innovation energy             */
    float *gain_pit,      /* i/o: adaptive excitation gain                        */
    const float gain_code,      /* i  : innovative excitation gain                      */
    float g_corr[],       /* o  : ACELP correlation values                        */
    const short clip_gain,      /* i  : adaptive gain clipping flag                     */
    float *mem_deemp,     /* i/o: prequantizer deemhasis memory                   */
    float *mem_preemp,    /* i/o: prequantizer preemhasis memory                  */
    float *gain_preQ,     /* o  : prequantizer excitation gain                    */
    float code_preQ[],    /* o  : prequantizer excitation                         */
    short *unbits         /* o  : number of AVQ unused bits                       */
)
{
    short i, index, nBits, Nsv;
    float x_in[L_SUBFR], x_tran[L_SUBFR], temp;
    int x_norm[L_SUBFR+L_SUBFR/WIDTH_BAND];
    float corr, ener;
    short nq[L_SUBFR/WIDTH_BAND];

    /*--------------------------------------------------------------*
     * Set bit-allocation
     *--------------------------------------------------------------*/

    Nsv = 8;
    nBits = AVQ_bits_16kHz_tbl[BIT_ALLOC_IDX_16KHZ(core_brate, coder_type, i_subfr, TC_SUBFR2IDX_16KHZ(tc_subfr))];

    /* increase # of AVQ allocated bits by unused bits from the previous subframe */
    nBits += (*unbits);

    /*--------------------------------------------------------------*
     * Compute/Update target
     * For inactive frame, find target in residual domain
     * Deemphasis
     *--------------------------------------------------------------*/

    if ( coder_type == INACTIVE )
    {
        for( i=0; i<L_SUBFR; i++ )
        {
            x_tran[i] = xn[i] - *gain_pit * y1[i] - gain_code * y2[i];
        }

        find_cn( x_tran, Ap, p_Aq, x_in );
    }
    else
    {
        updt_tar( cn, x_in, &exc[i_subfr], *gain_pit, L_SUBFR );
    }

    deemph( x_in, FAC_PRE_AVQ, L_SUBFR, mem_deemp );

    /*--------------------------------------------------------------*
     * DCT-II
     *--------------------------------------------------------------*/

    if( coder_type != INACTIVE && core_brate > ACELP_24k40 && core_brate <= ACELP_32k && !harm_flag_acelp )
    {
        mvr2r( x_in, x_tran, L_SUBFR );
    }
    else
    {
        edct2( L_SUBFR, -1, x_in, x_tran, ip_edct2_64, w_edct2_64 );
    }

    /*--------------------------------------------------------------*
     * Split algebraic vector quantizer based on RE8 lattice
     *--------------------------------------------------------------*/

    AVQ_cod( x_tran, x_norm, nBits, Nsv );

    /*--------------------------------------------------------------*
     * Find prequantizer excitation gain
     * Quantize the gain
     *--------------------------------------------------------------*/

    corr = 0;
    ener = 1e-6f;

    for (i = 0; i < Nsv*8; i++)
    {
        corr += x_tran[i]*(float)x_norm[i];
        ener += (float)x_norm[i]*(float)x_norm[i];
    }

    *gain_preQ = corr/ener;

    if ( coder_type == INACTIVE )
    {
        *gain_preQ /= gain_code;

        if( core_brate == ACELP_64k )
        {
            index = usquant( *gain_preQ, gain_preQ, G_AVQ_MIN_INACT_64k, G_AVQ_DELTA_INACT_64k, (1 << G_AVQ_BITS) );
        }
        else if( core_brate == ACELP_48k )
        {
            index = usquant( *gain_preQ, gain_preQ, G_AVQ_MIN_INACT_48k, G_AVQ_DELTA_INACT_48k, (1 << G_AVQ_BITS) );
        }
        else
        {
            index = usquant( *gain_preQ, gain_preQ, G_AVQ_MIN_INACT, G_AVQ_DELTA_INACT, (1 << G_AVQ_BITS) );
        }

        *gain_preQ *= gain_code;
    }
    else
    {
        if( Es_pred < 0  )
        {
            temp = 0.25f * fabs(Es_pred);
        }
        else
        {
            temp = Es_pred;
        }
        *gain_preQ /= temp;

        if( core_brate > ACELP_24k40 && core_brate <= ACELP_32k )
        {
            index = gain_quant( gain_preQ, 0.1f*G_AVQ_MIN, G_AVQ_MAX, G_AVQ_BITS );
        }
        else
        {
            index = gain_quant( gain_preQ, G_AVQ_MIN, G_AVQ_MAX, G_AVQ_BITS );
        }
        *gain_preQ *= temp;
    }

    push_indice( st, IND_AVQ_GAIN, index, G_AVQ_BITS );

    /*--------------------------------------------------------------*
     * Encode and multiplex subvectors into bit-stream
     *--------------------------------------------------------------*/

    AVQ_encmux( st, -1, x_norm, &nBits, Nsv, nq );

    /* save # of AVQ unused bits for next subframe */
    *unbits = nBits;

    /* at the last subframe, write AVQ unused bits */
    if( i_subfr == 4*L_SUBFR && extl != SWB_BWE_HIGHRATE && extl != FB_BWE_HIGHRATE )
    {
        while( *unbits > 0 )
        {
            i = min(*unbits, 16);
            push_indice( st, IND_UNUSED, 0, i );
            *unbits -= i;
        }
    }

    /*--------------------------------------------------------------*
     * DCT transform
     *--------------------------------------------------------------*/

    for( i=0; i<Nsv*WIDTH_BAND; i++ )
    {
        x_tran[i] = (float) (x_norm[i]);
    }

    set_f( x_tran+Nsv*WIDTH_BAND, 0.0f, L_SUBFR-WIDTH_BAND*Nsv );

    if( coder_type != INACTIVE && core_brate > ACELP_24k40 && core_brate <= ACELP_32k && !harm_flag_acelp )
    {
        mvr2r( x_tran, code_preQ, L_SUBFR );
    }
    else
    {
        edct2( L_SUBFR, 1, x_tran, code_preQ, ip_edct2_64, w_edct2_64 );
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
     * For inactive segments
     * - Zero-memory filtered pre-filter excitation
     * - Update of targets and gain_pit
     * For active segments
     * - Update xn[L_subfr-1] for updating the memory of the weighting filter
     *--------------------------------------------------------------*/

    if ( coder_type == INACTIVE )
    {
        temp = code_preQ[0] * h1[L_SUBFR-1];

        for( i=1; i<L_SUBFR; i++ )
        {
            temp += code_preQ[i] * h1[L_SUBFR-1-i];
        }

        xn[L_SUBFR-1] -= *gain_preQ * temp;
    }
    else
    {
        conv( code_preQ, h1, x_tran, L_SUBFR );

        updt_tar( cn, cn, code_preQ, *gain_preQ, L_SUBFR );
        updt_tar( xn, xn, x_tran, *gain_preQ, L_SUBFR );

        *gain_pit = corr_xy1( xn, y1, g_corr, L_SUBFR, 0);

        /* clip gain if necessary to avoid problems at decoder */
        if( clip_gain == 1 && *gain_pit > 0.95f )
        {
            *gain_pit = 0.95f;
        }

        updt_tar( xn, xn2, y1, *gain_pit, L_SUBFR );
    }

    st->use_acelp_preq = 1;

    return;

}

/*-------------------------------------------------------------------*
 * Find target in residual domain - cn[]
 *-------------------------------------------------------------------*/

static void find_cn(
    const float xn[],           /* i  : target signal                                   */
    const float Ap[],           /* i  : weighted LP filter coefficients                 */
    const float *p_Aq,          /* i  : 12k8 LP coefficients                            */
    float cn[]            /* o  : target signal in residual domain                */
)
{
    float tmp, tmp_fl[L_SUBFR+M];

    set_f( tmp_fl, 0, M );
    mvr2r( xn, tmp_fl+M, L_SUBFR );
    tmp = 0.0f;

    preemph( tmp_fl+M, PREEMPH_FAC_16k, L_SUBFR, &tmp );
    syn_filt( Ap, M, tmp_fl+M, tmp_fl+M, L_SUBFR, tmp_fl, 0 );
    residu( p_Aq, M, tmp_fl+M, cn, L_SUBFR );

    return;
}

/*---------------------------------------------------------------*
 * gain_quant()
 *
 * Quantization of gains between the specified range
 * using the specified number of levels.
 *---------------------------------------------------------------*/

short gain_quant(           /* o:   quantization index            */
    float *gain,      /* i/o: quantized gain                */
    const float min,        /* i:   value of lower limit          */
    const float max,        /* i:   value of upper limit          */
    const short bits        /* i:   number of bits to quantize    */
)
{
    short index, levels;
    float tmp, c_min, c_mult;

    levels = 1<<bits;

    if( *gain < FLT_MIN )
    {
        *gain = FLT_MIN;
    }

    c_min = (float)log10(min);
    c_mult = (float) ((levels-1)/(log10(max)-c_min));

    tmp = c_mult * ((float)log10(*gain) - c_min);
    index = (short)(tmp + 0.5f);

    if( index < 0 )
    {
        index = 0;
    }
    if( index > levels-1 )
    {
        index = levels-1;
    }

    *gain = (float)pow( 10.0, (((float)index)/c_mult) + c_min );

    return(index);
}
