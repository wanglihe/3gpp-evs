/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*---------------------------------------------------------------------*
 * Es_pred_dec()
 *
 * Decoding of scaled predicted innovation energy to be used in all subframes
 *---------------------------------------------------------------------*/

void Es_pred_dec(
    float *Es_pred,     /* o  : predicited scaled innovation energy       */
    const int   enr_idx,      /* i  : indice                                    */
    const short nb_bits,      /* i  : number of bits                             */
    const short no_ltp        /* i  : no LTP flag                                */
)
{
    if( !no_ltp )
    {
        switch ( nb_bits )
        {
        case 5:
            *Es_pred = Es_pred_qua_5b[enr_idx];
            break;
        case 4:
            *Es_pred = Es_pred_qua_4b[enr_idx];
            break;
        case 3:
            *Es_pred = Es_pred_qua_3b[enr_idx];
            break;
        default:
            *Es_pred = Es_pred_qua_5b[enr_idx];
            break;
        }
    }
    else
    {
        *Es_pred = Es_pred_qua_4b_no_ltp[enr_idx];
    }

    return;
}


/*--------------------------------------------------------------------------*
 * lp_gain_updt()
 *
 * Update of LP pitch and code gains (FEC)
 *-------------------------------------------------------------------------*/

void lp_gain_updt(
    const short i_subfr,          /* i  :  subframe number            */
    const float gain_pit,         /* i  : Decoded gain pitch          */
    const float norm_gain_code,   /* i  : Normalised gain code        */
    float *lp_gainp,              /* i/o: LP-filtered pitch gain(FEC) */
    float *lp_gainc,              /* i/o: LP-filtered code gain (FEC) */
    const short L_frame           /* i  : length of the frame         */
)
{
    if( L_frame == L_FRAME )
    {
        if( i_subfr == 0 )
        {
            *lp_gainp = 0.1f * gain_pit;
            *lp_gainc = 0.1f * norm_gain_code;
        }
        else if( i_subfr == L_SUBFR )
        {
            *lp_gainp += 0.2f * gain_pit;
            *lp_gainc += 0.2f * norm_gain_code;
        }
        else if( i_subfr == 2*L_SUBFR )
        {
            *lp_gainp += 0.3f * gain_pit;
            *lp_gainc += 0.3f * norm_gain_code;
        }
        else  /* i_subfr == 3*L_SUBFR */
        {
            *lp_gainp += 0.4f * gain_pit;
            *lp_gainc += 0.4f * norm_gain_code;
        }
    }
    else
    {
        if( i_subfr == 0 )
        {
            *lp_gainp = (1.0f/15.0f) * gain_pit;
            *lp_gainc = (1.0f/15.0f) * norm_gain_code;
        }
        else if( i_subfr == L_SUBFR )
        {
            *lp_gainp += (2.0f/15.0f) * gain_pit;
            *lp_gainc += (2.0f/15.0f) * norm_gain_code;
        }
        else if( i_subfr == 2*L_SUBFR )
        {
            *lp_gainp += (3.0f/15.0f) * gain_pit;
            *lp_gainc += (3.0f/15.0f) * norm_gain_code;
        }
        else if( i_subfr == 3*L_SUBFR )
        {
            *lp_gainp += (4.0f/15.0f) * gain_pit;
            *lp_gainc += (4.0f/15.0f) * norm_gain_code;
        }
        else  /* i_subfr == 4*L_SUBFR */
        {
            *lp_gainp += (5.0f/15.0f) * gain_pit;
            *lp_gainc += (5.0f/15.0f) * norm_gain_code;
        }
    }

    return;
}


/*---------------------------------------------------------------------*
 * gain_dec_tc()
 *
 * Decoding of pitch and codebook gains and updating long term energies
 *---------------------------------------------------------------------*/

void gain_dec_tc(
    Decoder_State *st,            /* i/o: decoder state structure             */
    const long  core_brate,     /* i  : core bitrate                        */
    const short L_frame,        /* i  : length of the frame                 */
    const short i_subfr,        /* i  : subframe number                     */
    const short tc_subfr,       /* i  :  TC subframe index                  */
    const float Es_pred,        /* i  :  predicted scaled innov. energy     */
    const float *code,          /* i  : algebraic code excitation           */
    float *gain_pit,      /* o  : pitch gain                          */
    float *gain_code,     /* o  : Quantized codeebook gain            */
    float *gain_inov,     /* o  : unscaled innovation gain            */
    float *norm_gain_code /* o  : norm. gain of the codebook excit.   */
)
{
    short index, nBits;
    float Ecode, gcode0;
    float Ei;

    *gain_pit = 0;

    /*----------------------------------------------------------------*
     * find number of bits for gain dequantization
     *----------------------------------------------------------------*/

    if( L_frame == L_FRAME )
    {
        nBits = gain_bits_tbl[BIT_ALLOC_IDX(core_brate, TRANSITION, i_subfr, TC_SUBFR2IDX(tc_subfr))];
    }
    else  /* L_frame == L_FRAME16k */
    {
        nBits = gain_bits_16kHz_tbl[BIT_ALLOC_IDX_16KHZ(core_brate, TRANSITION, i_subfr, TC_SUBFR2IDX_16KHZ(tc_subfr))];
    }

    /*-----------------------------------------------------------------*
     * calculate the predicted gain code
     *-----------------------------------------------------------------*/

    Ecode = (dotp( code, code, L_SUBFR) + 0.01f) / L_SUBFR;
    *gain_inov = 1.0f / (float)sqrt( Ecode );
    Ei = 10 * (float)log10( Ecode );
    gcode0 = (float) pow(10, 0.05 * (Es_pred - Ei));

    /*------------------------------------------------------------------------------------------*
     * Select the gain quantization table and dequantize the gain
     *------------------------------------------------------------------------------------------*/

    index = (short)get_next_indice( st, nBits );

    if( nBits > 3 )
    {
        *gain_code = gain_dequant( index, G_CODE_MIN, G_CODE_MAX, nBits );
    }
    else /* nBits == 3 */
    {
        *gain_code = tbl_gain_code_tc[index];
    }

    /*-----------------------------------------------------------------*
     * decode normalized codebook gain
     *-----------------------------------------------------------------*/

    *gain_code *= gcode0;
    *norm_gain_code = *gain_code / *gain_inov;

    return;
}

/*---------------------------------------------------------------------*
 * gain_dec_amr_wb()
 *
 * Decoding of pitch and fixed codebook gains (used also in AMR-WB IO mode)
 *---------------------------------------------------------------------*/

void gain_dec_amr_wb(
    Decoder_State *st,              /* i/o: decoder state structure               */
    const long  core_brate,       /* i  : core bitrate                          */
    float *gain_pit,        /* o  : Quantized pitch gain                  */
    float *gain_code,       /* o  : Quantized codeebook gain              */
    float *past_qua_en,     /* i/o: gain quantization memory (4 words)    */
    float *gain_inov,       /* o  : unscaled innovation gain              */
    const float *code,            /* i  : algebraic code excitation             */
    float *norm_gain_code   /* o  : norm. gain of the codebook excitation */
)
{
    short i, index;
    short nbits;
    float gcode0, qua_en;
    const float *t_qua_gain;

    *gain_inov = 1.0f/ (float)sqrt( ( dotp( code, code, L_SUBFR ) + 0.01f ) / L_SUBFR );

    /*-----------------------------------------------------------------*
     * Select the gain quantization table
     *-----------------------------------------------------------------*/

    if( core_brate < ACELP_12k65)
    {
        nbits = 6;
        t_qua_gain = t_qua_gain6b;
    }
    else
    {
        nbits = 7;
        t_qua_gain = t_qua_gain7b;
    }

    /*-----------------------------------------------------------------*
     * predicted code gain
     *-----------------------------------------------------------------*/

    /* start with predicting code energy in dB */
    gcode0 = MEAN_ENER;
    for (i=0; i<GAIN_PRED_ORDER; i++)
    {
        gcode0 += pred_gain[i] * past_qua_en[i];
    }
    gcode0 += (float)(20.0 * log10( *gain_inov ) );

    /* convert from energy in dB to gain */
    gcode0 = (float)pow( 10.0, gcode0/20.0);

    /*-----------------------------------------------------------------*
     * Decode pitch gain
     *-----------------------------------------------------------------*/

    index = (short)get_next_indice( st, nbits );
    *gain_pit = t_qua_gain[index * 2];

    /*-----------------------------------------------------------------*
     * Decode code gain
     *-----------------------------------------------------------------*/

    qua_en = t_qua_gain[index * 2 + 1];
    *gain_code = qua_en * gcode0;

    /*-----------------------------------------------------------------*
     * update table of past quantized energies
     *-----------------------------------------------------------------*/

    for (i=GAIN_PRED_ORDER-1; i>0; i--)
    {
        past_qua_en[i] = past_qua_en[i-1];
    }
    past_qua_en[0] = (float)(20.0*log10(qua_en));

    /*-----------------------------------------------------------------*
     * Normalized code gain
     *-----------------------------------------------------------------*/

    *norm_gain_code = *gain_code / *gain_inov;

    return;
}

/*--------------------------------------------------------------------------*
 * gain_dec_mless()
 *
 * Decoding of pitch and codebook gains without updating long term energies
 *-------------------------------------------------------------------------*/

void gain_dec_mless(
    Decoder_State *st,              /* i/o: decoder state structure               */
    const long  core_brate,       /* i  : core bitrate                          */
    const short L_frame,          /* i  : length of the frame                   */
    const short coder_type,       /* i  : coding type                           */
    const short i_subfr,          /* i  : subframe number                       */
    const short tc_subfr,         /* i  : TC subframe index                     */
    const float *code,            /* i  : algebraic code excitation             */
    const float Es_pred,          /* i  : predicted scaled innov. energy        */
    float *gain_pit,        /* o  : Quantized pitch gain                  */
    float *gain_code,       /* o  : Quantized codeebook gain              */
    float *gain_inov,       /* o  : unscaled innovation gain              */
    float *norm_gain_code   /* o  : norm. gain of the codebook excitation */
)
{
    short index, nBits;
    float gcode0, Ei, Ecode;
    const float *qua_table;

    /*-----------------------------------------------------------------*
     * decode pitch gain
     *-----------------------------------------------------------------*/

    if( L_frame == L_FRAME )
    {
        nBits = gain_bits_tbl[BIT_ALLOC_IDX(core_brate, coder_type, i_subfr, TC_SUBFR2IDX(tc_subfr))];
    }
    else  /* L_frame == L_FRAME16k */
    {
        nBits = gain_bits_16kHz_tbl[BIT_ALLOC_IDX_16KHZ(core_brate, coder_type, i_subfr, TC_SUBFR2IDX_16KHZ(tc_subfr))];
    }

    if( (tc_subfr == 3*L_SUBFR && i_subfr == 3*L_SUBFR && L_frame == L_FRAME) ||
            (tc_subfr == 4*L_SUBFR && i_subfr == 4*L_SUBFR && L_frame == L_FRAME16k) )
    {
        /* decode pitch gain */
        index = (short)get_next_indice( st, nBits >> 1 );
        Ei = (G_PITCH_MAX_TC192 - G_PITCH_MIN_TC192) / ((1 << (nBits >> 1)) - 1);  /* set quantization step */
        *gain_pit = usdequant( index, G_PITCH_MIN_TC192, Ei );

        /* calculate the predicted gain code */
        Ecode = (dotp( code, code, L_SUBFR) + 0.01f) / L_SUBFR;
        *gain_inov = 1.0f / (float)sqrt( Ecode );
        Ei = 10 * (float)log10( Ecode );
        gcode0 = (float) pow(10, 0.05 * (Es_pred - Ei));

        /* decode normalized codebook gain */
        index = (short)get_next_indice( st, (nBits+1)>>1 );
        *gain_code = gain_dequant( index, G_CODE_MIN_TC192, G_CODE_MAX_TC192, (nBits+1)>>1 );
        *gain_code *= gcode0;
    }
    else
    {
        switch ( nBits )
        {
        case 7:
        {
            qua_table = gain_qua_mless_7b;
            break;
        }
        case 6:
        {
            qua_table = gain_qua_mless_6b;
            break;
        }
        case 5:
        {
            qua_table = gain_qua_mless_5b;
            break;
        }
        default:
        {
            qua_table = gain_qua_mless_6b;
            break;
        }
        }

        if( coder_type == INACTIVE && nBits == 6  )
        {
            nBits--;
        }

        index = (short)get_next_indice( st, nBits );

        *gain_pit = qua_table[index * 2];
        Ecode = (dotp( code, code, L_SUBFR) + 0.01f) / L_SUBFR;
        *gain_inov = 1.0f / (float)sqrt( Ecode );
        Ei = 10 * (float)log10( Ecode );

        /*-----------------------------------------------------------------*
         * calculate the predicted gain code
         *-----------------------------------------------------------------*/

        gcode0 = (float) pow(10, 0.05 * (Es_pred - Ei));

        /*-----------------------------------------------------------------*
         * decode normalized codebook gain
         *-----------------------------------------------------------------*/

        *gain_code = qua_table[index * 2 + 1] * gcode0;
    }

    *norm_gain_code = *gain_code / *gain_inov;

    return;
}

/*--------------------------------------------------------------------------*
 * gain_dec_lbr()
 *
 * Decoding of pitch and codebook gains in ACELP at 6.6 and 7.5 kbps
 *-------------------------------------------------------------------------*/

void gain_dec_lbr(
    Decoder_State *st,               /* i/o: decoder state structure                          */
    const long  core_brate,        /* i  : core bitrate                                     */
    const short coder_type,        /* i  : coding type                                      */
    const short i_subfr,           /* i  : subframe index                                   */
    const float *code,             /* i  : algebraic excitation                             */
    float *gain_pit,         /* o  : quantized pitch gain                             */
    float *gain_code,        /* o  : quantized codebook gain                          */
    float *gain_inov,        /* o  : gain of the innovation (used for normalization)  */
    float *norm_gain_code,   /* o  : norm. gain of the codebook excitation            */
    float gains_mem[]        /* i/o: pitch gain and code gain from previous subframes */
)
{
    short index, nBits, n_pred, ctype;
    float gcode0, aux[10], Ecode;
    const float *b, *cdbk = 0;

    Ecode = ( dotp( code, code, L_SUBFR ) + 0.01f ) / L_SUBFR;
    *gain_inov = 1.0f / (float)sqrt(Ecode);

    /*-----------------------------------------------------------------*
     * select the codebook, size and number of bits
     * set the gains searching range
     *-----------------------------------------------------------------*/

    nBits = gain_bits_tbl[BIT_ALLOC_IDX(core_brate, coder_type, i_subfr, 0)];
    ctype = 2*(coder_type - 1);

    /*-----------------------------------------------------------------*
     * calculate prediction of gcode
     * search for the best codeword
     *-----------------------------------------------------------------*/

    if (i_subfr == 0)
    {
        b = b_1sfr;
        n_pred = 2;

        switch ( nBits )
        {
        case 8:
        {
            cdbk = gp_gamma_1sfr_8b;
            break;
        }
        case 7:
        {
            cdbk = gp_gamma_1sfr_7b;
            break;
        }
        case 6:
        {
            cdbk = gp_gamma_1sfr_6b;
            break;
        }
        }

        /* calculate predicted gain */
        aux[0] = 1.0f;
        aux[1] = ctype;
        gcode0 = (float)pow(10, dotp(b, aux, n_pred) - 0.5f * (float)log10(Ecode));

        /* retrieve the codebook index and calculate both gains */
        index = (short)get_next_indice( st, nBits );
        *gain_pit = cdbk[index * 2];
        *gain_code = cdbk[index * 2 + 1] * gcode0;
        gains_mem[0] = *gain_code;
        gains_mem[3] = *gain_pit;
    }
    else if (i_subfr == L_SUBFR)
    {
        b = b_2sfr;
        n_pred = 4;

        switch ( nBits )
        {
        case 7:
        {
            cdbk = gp_gamma_2sfr_7b;
            break;
        }
        case 6:
        {
            cdbk = gp_gamma_2sfr_6b;
            break;
        }
        }

        /* calculate predicted gain */
        aux[0] = 1.0f;
        aux[1] = ctype;
        aux[2] = (float)log10(gains_mem[0]);
        aux[3] = gains_mem[3];
        gcode0 = (float)pow(10, dotp(b, aux, n_pred));

        /* retrieve the codebook index and calculate both gains */
        index = (short)get_next_indice( st, nBits );
        *gain_pit = cdbk[index * 2];
        *gain_code = cdbk[index * 2 + 1] * gcode0;
        gains_mem[1] = *gain_code;
        gains_mem[4] = *gain_pit;
    }
    else if (i_subfr == 2*L_SUBFR)
    {
        b = b_3sfr;
        n_pred = 6;

        cdbk = gp_gamma_3sfr_6b;

        /* calculate predicted gain */
        aux[0] = 1.0f;
        aux[1] = ctype;
        aux[2] = (float)log10(gains_mem[0]);
        aux[3] = (float)log10(gains_mem[1]);
        aux[4] = gains_mem[3];
        aux[5] = gains_mem[4];
        gcode0 = (float)pow(10, dotp(b, aux, n_pred));

        /* retrieve the codebook index and calculate both gains */
        index = (short)get_next_indice( st, nBits );
        *gain_pit = cdbk[index * 2];
        *gain_code = cdbk[index * 2 + 1] * gcode0;

        gains_mem[2] = *gain_code;
        gains_mem[5] = *gain_pit;
    }
    else if (i_subfr == 3*L_SUBFR)
    {
        b = b_4sfr;
        n_pred = 8;

        cdbk = gp_gamma_4sfr_6b;

        /* calculate predicted gain */
        aux[0] = 1.0f;
        aux[1] = ctype;
        aux[2] = (float)log10(gains_mem[0]);
        aux[3] = (float)log10(gains_mem[1]);
        aux[4] = (float)log10(gains_mem[2]);
        aux[5] = gains_mem[3];
        aux[6] = gains_mem[4];
        aux[7] = gains_mem[5];
        gcode0 = (float)pow(10, dotp(b, aux, n_pred));

        /* retrieve the codebook index and calculate both gains */
        index = (short)get_next_indice( st, nBits );
        *gain_pit = cdbk[index * 2];
        *gain_code = cdbk[index * 2 + 1] * gcode0;
    }

    *norm_gain_code = *gain_code / *gain_inov;

    return;
}

/*--------------------------------------------------------------------------*
 * gain_dec_SQ()
 *
 * Decoding of pitch and codebook gains using scalar quantizers
 *-------------------------------------------------------------------------*/

void gain_dec_SQ(
    Decoder_State *st,              /* i/o: decoder state structure               */
    const long  core_brate,       /* i  : core bitrate                          */
    const short coder_type,       /* i  : coding type                           */
    const short i_subfr,          /* i  : subframe number                       */
    const short tc_subfr,         /* i  : TC subframe index                     */
    const float *code,            /* i  : algebraic code excitation             */
    const float Es_pred,          /* i  : predicted scaled innov. energy        */
    float *gain_pit,        /* o  : Quantized pitch gain                  */
    float *gain_code,       /* o  : Quantized codeebook gain              */
    float *gain_inov,       /* o  : unscaled innovation gain              */
    float *norm_gain_code   /* o  : norm. gain of the codebook excitation */
)
{
    short index, nBits;
    float gcode0, Ei, Ecode;
    short tmp16;
    /*-----------------------------------------------------------------*
     * get number of bits
     *-----------------------------------------------------------------*/

    nBits = gain_bits_16kHz_tbl[ BIT_ALLOC_IDX_16KHZ(core_brate, coder_type, i_subfr, TC_SUBFR2IDX_16KHZ(tc_subfr)) ];

    /*-----------------------------------------------------------------*
     * decode pitch gain
     *-----------------------------------------------------------------*/

    index = (short)get_next_indice( st, nBits>>1 );
    /*Ei = (G_PITCH_MAX - G_PITCH_MIN) / ((1 << (nBits>>1)) - 1);*/  /* set quantization step */
    tmp16 = div_s(1,((1 << (nBits>>1)) - 1) );      /* Q15*/
    Ei = (float)mult_r((short)(G_PITCH_MAX*8192.0f+0.5f),tmp16)/8192.0f;

    *gain_pit = usdequant( index, G_PITCH_MIN, Ei );

    /*-----------------------------------------------------------------*
     * calculate the predicted gain code
     *-----------------------------------------------------------------*/

    Ecode = (dotp( code, code, L_SUBFR) + 0.01f) / L_SUBFR;
    *gain_inov = 1.0f / (float)sqrt( Ecode );
    Ei = 10 * (float)log10( Ecode );
    gcode0 = (float) pow(10, 0.05 * (Es_pred - Ei));

    /*-----------------------------------------------------------------*
     * decode normalized codebook gain
     *-----------------------------------------------------------------*/

    index = (short)get_next_indice( st, (nBits+1)>>1 );
    *gain_code = gain_dequant( index, G_CODE_MIN, G_CODE_MAX, (nBits+1)>>1 );
    *gain_code *= gcode0;
    *norm_gain_code = *gain_code / *gain_inov;

    return;
}


/*-------------------------------------------------*
 * gain_dec_gaus()
 *
 * Decoding of gains for Gaussian codebook
 *-------------------------------------------------*/

float gain_dec_gaus(            /* o  : quantized codebook gain         */
    const short index,          /* i  : quantization index              */
    const short bits,           /* i  : number of bits to quantize      */
    const float lowBound,       /* i  : lower bound of quantizer (dB)   */
    const float topBound,       /* i  : upper bound of quantizer (dB)   */
    const float gain_inov,      /* i  : unscaled innovation gain        */
    float *norm_gain_code /* o  : gain of normalized gaus. excit. */
)
{
    float gain, enr, stepSize;

    /*-----------------------------------------------------------------*
     * quantize linearly the log E
     *-----------------------------------------------------------------*/

    stepSize = (topBound - lowBound)/((float)(1<<bits));

    /*-----------------------------------------------------------------*
     * Gaussian codebook gain
     *-----------------------------------------------------------------*/

    enr = (float)index * stepSize + lowBound;         /* quantized codebook gain in dB */
    gain = (float)pow( 10.0f, enr/20.0f );

    *norm_gain_code = gain / gain_inov;

    return gain;
}
