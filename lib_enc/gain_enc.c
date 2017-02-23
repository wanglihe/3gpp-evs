/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * Local constants
 *-------------------------------------------------------------------*/

#define RANGE          64

/*---------------------------------------------------------------------*
 * Es_pred_enc()
 *
 * Calculation and quantization of average predicted innovation energy to be
 *---------------------------------------------------------------------*/

void Es_pred_enc(
    float *Es_pred,        /* o  : predicited scaled innovation energy       */
    int *Es_pred_indice, /* o  : indice corresponding to above parameter   */
    const short L_frame,         /* i  : length of the frame                       */
    const short L_subfr,         /* i  : length of the subframe                    */
    const float *res,            /* i  : residual signal                           */
    const float *voicing,        /* i  : normalized correlation in three 1/2frames */
    const short nb_bits,         /* i  : allocated number of bits                  */
    const short no_ltp           /* i  : no_ltp flag                               */
)
{
    short i, i_subfr, size;
    float tmp, dist, mean_ener_code, ener;
    float weight;
    const float *qua_table;

    /*----------------------------------------------------------*
     * calculate the average residual signal energy
     *----------------------------------------------------------*/
    if( L_frame == L_FRAME )
    {
        weight = 0.25f;
    }
    else  /* L_frame == L_FRAME16k */
    {
        weight = 0.2f;
    }

    mean_ener_code = 0.0f;
    for (i_subfr = 0; i_subfr < L_frame; i_subfr += L_subfr)
    {
        /* calculate the energy of residual signal */
        ener = sum2_f(&res[i_subfr], L_subfr) + 0.01f;
        ener = 10 * (float)log10( ener / ((float) L_subfr));
        if ((ener < 0) && !(no_ltp))
        {
            ener = 0;
        }

        /* update the average energy of residual signal */
        mean_ener_code += weight * ener;
    }

    if( !no_ltp )
    {
        /*----------------------------------------------------------*
         * subtract an estimate of adaptive codebook contribution
         *----------------------------------------------------------*/

        mean_ener_code -= 10.0f * (0.5f * voicing[0] + 0.5f * voicing[1]);

        /*----------------------------------------------------------*
         * quantize the average predicted innovation energy
         *----------------------------------------------------------*/
        switch ( nb_bits )
        {
        case 5:
        {
            qua_table = Es_pred_qua_5b;
            break;
        }
        case 4:
        {
            qua_table = Es_pred_qua_4b;
            break;
        }
        case 3:
        {
            qua_table = Es_pred_qua_3b;
            break;
        }
        default:
        {
            qua_table = Es_pred_qua_5b;
            break;
        }
        }
    }
    else
    {
        qua_table = Es_pred_qua_4b_no_ltp;
    }

    /* select codebook, size and number of bits */
    size = (short)pow2[nb_bits];

    /* find the nearest neighbour (codevector) */
    tmp = 1e30f;
    *Es_pred_indice=0;
    for (i=0; i<size; i++)
    {
        dist = (float)fabs(mean_ener_code - qua_table[i]);
        if (dist < tmp)
        {
            tmp = dist;
            *Es_pred = qua_table[i];
            *Es_pred_indice = i;
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * gain_enc_amr_wb()
 *
 * Quantization of pitch and codebook gains (used also in AMR-WB IO mode)
 * MA prediction is performed on the innovation energy (in dB with mean removed).
 * An initial predicted gain, gcode0, is first determined and the correction
 * factor     alpha = g_code / gcode0   is quantized.
 * The pitch gain and the correction factor are vector quantized and the
 * mean-squared weighted error criterion is used in the quantizer search.
 *-------------------------------------------------------------------*/

void gain_enc_amr_wb(
    Encoder_State *st,               /* i/o: encoder state structure      */
    const float *xn,               /* i  : target vector                                                   */
    const float *y1,               /* i  : zero-memory filtered adaptive excitation                        */
    const float *y2,               /* i  : zero-memory filtered algebraic codebook excitation              */
    const float *code,             /* i  : algebraic excitation                                            */
    const long core_brate,         /* i  : core bitrate                                                    */
    float *gain_pit,         /* i/o: pitch gain / Quantized pitch gain                               */
    float *gain_code,        /* o  : quantized codebook gain                                         */
    float *gain_inov,        /* o  : gain of the innovation (used for normalization)                 */
    float *norm_gain_code,   /* o  : norm. gain of the codebook excitation                           */
    float *g_corr,           /* i/o: correlations <y1,y1>, -2<xn,y1>,<y2,y2>, -2<xn,y2> and 2<y1,y2> */
    const short clip_gain,         /* i  : gain pitch clipping flag (1 = clipping)                         */
    float *past_qua_en       /* i/o: gain quantization memory (4 words)                              */
)
{
    short index, i, j, min_ind, size, nBits;
    float dist, dist_min, g_pitch, g_code, qua_en, gcode0;
    const float *p, *t_qua_gain;

    /*-----------------------------------------------------------------*
     * gain computation correlations
     * find raw innovation energy
     *-----------------------------------------------------------------*/

    E_corr_xy2( xn, y1, y2, g_corr,L_SUBFR);
    g_corr[2] += 0.01F;
    g_corr[3] -= 0.02F;
    g_corr[4] += 0.02F;
    *gain_inov = 1.0f / (float)sqrt( ( dotp( code, code, L_SUBFR ) + 0.01f ) / L_SUBFR );

    /*-----------------------------------------------------------------*
     * find the initial quantization pitch index
     * set gains search range
     *-----------------------------------------------------------------*/

    if ( core_brate < ACELP_12k65 )
    {
        t_qua_gain = t_qua_gain6b;
        nBits = 6;
        min_ind = 0;
        size = RANGE;
        if (clip_gain == 1)
        {
            size -= 16;                       /* limit pitch gain  to 1.0 */
        }
    }
    else
    {
        t_qua_gain = t_qua_gain7b;
        nBits = 7;
        p = t_qua_gain7b + RANGE;             /* pt at 1/4th of table */
        j = NB_QUA_GAIN7B - RANGE;
        if (clip_gain == 1)
        {
            j -= 27;                          /* limit pitch gain to 1.0 */
        }

        min_ind = 0;
        g_pitch = *gain_pit;


        for (i=0; i<j; i++, p+=2)
        {
            if (g_pitch > *p)
            {
                min_ind++;
            }
        }
        size = RANGE;
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
     * search the codebook
     *-----------------------------------------------------------------*/

    dist_min = 3.402823466e+38F;
    p = t_qua_gain + min_ind*2;

    index = 0;
    for (i = 0; i<size; i++)
    {
        g_pitch = *p++;                   /* pitch gain */
        g_code = gcode0 **p++;           /* code gain */
        dist = g_pitch*g_pitch * g_corr[0]
               + g_pitch         * g_corr[1]
               + g_code*g_code   * g_corr[2]
               + g_code          * g_corr[3]
               + g_pitch*g_code  * g_corr[4];
        if (dist < dist_min)
        {
            dist_min = dist;
            index = i;
        }
    }
    index = index + min_ind;
    *gain_pit  = t_qua_gain[index * 2];
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

    push_indice( st, IND_GAIN, index, nBits );

    *norm_gain_code = *gain_code / *gain_inov;

    return;
}

/*---------------------------------------------------------------------*
 * gain_enc_mless()
 *
 * Quantization of pitch and codebook gains without prediction (memory-less)
 * - an initial predicted gain, gcode0, is first determined based on
 *   the predicted average innovation energy
 * - a correction factor gamma = g_code / gcode0 is then vector quantized along with gain_pit
 * - the mean-squared weighted error criterion is used for codebook search
 *---------------------------------------------------------------------*/

void gain_enc_mless(
    Encoder_State *st,               /* i/o: encoder state structure      */
    const long  core_brate,        /* i  : core bitrate                                                    */
    const short L_frame,           /* i  : length of the frame                                             */
    const short coder_type,        /* i  : coding type                                                     */
    const short i_subfr,           /* i  : subframe index                                                  */
    const short tc_subfr,          /* i  : TC subframe index                                               */
    const float *xn,               /* i  : target vector                                                   */
    const float *y1,               /* i  : zero-memory filtered adaptive excitation                        */
    const float *y2,               /* i  : zero-memory filtered algebraic codebook excitation              */
    const float *code,             /* i  : algebraic excitation                                            */
    const float Es_pred,           /* i  : predicted scaled innovation energy                              */
    float *gain_pit,         /* o  : quantized pitch gain                                            */
    float *gain_code,        /* o  : quantized codebook gain                                         */
    float *gain_inov,        /* o  : gain of the innovation (used for normalization)                 */
    float *norm_gain_code,   /* o  : norm. gain of the codebook excitation                           */
    float *g_corr,           /* i/o: correlations <y1,y1>, -2<xn,y1>,<y2,y2>, -2<xn,y2> and 2<y1,y2> */
    const short clip_gain          /* i  : gain pitch clipping flag (1 = clipping)                         */
)
{
    short index, i, size, nBits;
    float dist, dist_min, g_pitch, g_code, gcode0, Ei, Ecode;
    short nBits2;
    float tmp1, tmp2;
    const float *p, *qua_table;

    /*-----------------------------------------------------------------*
     * calculate the rest of the correlation coefficients
     * c2 = <y2,y2>, c3 = -2<xn,y2>, c4 = 2<y1,y2>
     *-----------------------------------------------------------------*/

    E_corr_xy2( xn, y1, y2, g_corr,L_SUBFR);
    g_corr[2] += 0.01F;
    g_corr[3] -= 0.02F;
    g_corr[4] += 0.02F;

    /*-----------------------------------------------------------------*
     * calculate the unscaled innovation energy
     * calculate the predicted gain code
     *-----------------------------------------------------------------*/

    Ecode = (dotp( code, code, L_SUBFR) + 0.01f) / L_SUBFR;
    *gain_inov = 1.0f / (float)sqrt( Ecode );
    Ei = 10 * (float)log10( Ecode );
    gcode0 = (float) pow(10, 0.05 * (Es_pred - Ei));

    /*-----------------------------------------------------------------*
     * select the codebook, size and number of bits
     * set the gains searching range
     *-----------------------------------------------------------------*/

    if( L_frame == L_FRAME )
    {
        nBits = gain_bits_tbl[BIT_ALLOC_IDX(core_brate, coder_type, i_subfr, TC_SUBFR2IDX(tc_subfr))];
    }
    else  /* L_frame == L_FRAME16k */
    {
        nBits = gain_bits_16kHz_tbl[ BIT_ALLOC_IDX_16KHZ(core_brate, coder_type, i_subfr, TC_SUBFR2IDX_16KHZ(tc_subfr)) ];
    }

    if( (tc_subfr == 3*L_SUBFR && i_subfr == 3*L_SUBFR && L_frame == L_FRAME) ||
            (tc_subfr == 4*L_SUBFR && i_subfr == 4*L_SUBFR && L_frame == L_FRAME16k) )
    {
        /* in case of attack at the end of the frame, use scalar gain quantizers */
        tmp1 = (g_corr[0]*g_corr[2]) - (0.25f*g_corr[4]*g_corr[4]);
        tmp2 = -0.5f*g_corr[1]/tmp1;
        tmp1 = -0.5f*g_corr[3]/tmp1;

        *gain_pit =  (g_corr[2]*tmp2) - (0.5f*g_corr[4]*tmp1);
        *gain_code = (g_corr[0]*tmp1) - (0.5f*g_corr[4]*tmp2);

        *gain_pit = max(G_PITCH_MIN_TC192,min(*gain_pit,G_PITCH_MAX_TC192));

        /* set number of bits for two SQs */
        nBits2 = (nBits+1)>>1;
        nBits = nBits>>1;

        /* gain_pit Q */
        tmp1 = (G_PITCH_MAX_TC192 - G_PITCH_MIN_TC192) / ((1 << nBits) - 1);      /* set quantization step */
        index = usquant( *gain_pit, gain_pit, G_PITCH_MIN_TC192, tmp1, (1 << nBits) );
        push_indice( st, IND_GAIN_PIT, index, nBits );

        /* gain_code Q */
        *gain_code /= gcode0;
        index = gain_quant( gain_code, G_CODE_MIN_TC192, G_CODE_MAX_TC192, nBits2 );
        push_indice( st, IND_GAIN_CODE, index, nBits2 );
        *gain_code *= gcode0;
    }
    else
    {
        size = (short)pow2[nBits];

        switch ( nBits )
        {
        case 7:
        {
            qua_table = gain_qua_mless_7b;
            if ( clip_gain == 1 ) size -= 30;
            break;
        }
        case 6:
        {
            qua_table = gain_qua_mless_6b;
            if ( clip_gain == 1 ) size -= 14;
            break;
        }
        case 5:
        {
            qua_table = gain_qua_mless_5b;
            if ( clip_gain == 1 ) size -= 6;
            break;
        }
        default:
        {
            qua_table = gain_qua_mless_6b;
            if ( clip_gain == 1 ) size -= 14;
            break;
        }
        }

        /* in case of AVQ inactive, limit the gain_pit to 0.65 */
        if( clip_gain == 2 && nBits == 6 )
        {
            size -= 36;
            nBits--;
        }

        /*-----------------------------------------------------------------*
         * search for the best quantizer
         *-----------------------------------------------------------------*/

        p = qua_table;
        dist_min = 3.402823466e+38F;
        index = 0;
        for (i = 0; i<size; i++)
        {
            g_pitch = *p++;                   /* pitch gain */
            g_code = gcode0 **p++;           /* code gain */
            dist = g_pitch*g_pitch * g_corr[0]
                   + g_pitch         * g_corr[1]
                   + g_code*g_code   * g_corr[2]
                   + g_code          * g_corr[3]
                   + g_pitch*g_code  * g_corr[4];
            if (dist < dist_min)
            {
                dist_min = dist;
                index = i;
            }
        }
        *gain_pit  = qua_table[index*2];
        *gain_code = qua_table[index*2+1] * gcode0;

        push_indice( st, IND_GAIN, index, nBits );
    }

    *norm_gain_code = *gain_code / *gain_inov;

    return;
}


/*---------------------------------------------------------------------*
 * gain_enc_SQ()
 *
 * Scalar Quantization of pitch and codebook gains without prediction
 * - an initial predicted gain, gcode0, is first determined based on
 *   the predicted scaled innovation energy
 * - a correction factor gamma = g_code / gcode0 is then vector quantized
 *   along with gain_pit
 * - the mean-squared weighted error criterion is used for codebook search
 *---------------------------------------------------------------------*/

void gain_enc_SQ(
    Encoder_State *st,               /* i/o: encoder state structure      */
    const long  core_brate,        /* i  : core bitrate                                                    */
    const short coder_type,        /* i  : coding type                                                     */
    const short i_subfr,           /* i  : subframe index                                                  */
    const short tc_subfr,          /* i  : TC subframe index                                               */
    const float *xn,               /* i  : target vector                                                   */
    const float *yy1,              /* i  : zero-memory filtered adaptive excitation                        */
    const float *y2,               /* i  : zero-memory filtered algebraic codebook excitation              */
    const float *code,             /* i  : algebraic excitation                                            */
    const float Es_pred,           /* i  : predicted scaled innovation energy                              */
    float *gain_pit,         /* o  : quantized pitch gain                                            */
    float *gain_code,        /* o  : quantized codebook gain                                         */
    float *gain_inov,        /* o  : gain of the innovation (used for normalization)                 */
    float *norm_gain_code,   /* o  : norm. gain of the codebook excitation                           */
    float *g_corr,           /* i/o: correlations <y1,y1>, -2<xn,y1>,<y2,y2>, -2<xn,y2> and 2<y1,y2> */
    const short clip_gain          /* i  : gain pitch clipping flag (1 = clipping)                         */
)
{
    short index, nBits_pitch, nBits_code;
    float g_code, gcode0, Ei, Ecode, tmp1, tmp2;
    short tmp16;
    /*-----------------------------------------------------------------*
     * calculate the rest of the correlation coefficients
     * c2 = <y2,y2>, c3 = -2<xn,y2>, c4 = 2<y1,y2>
     *-----------------------------------------------------------------*/

    g_corr[1] *= -0.5;
    g_corr[2] = dotp( y2, y2, L_SUBFR )  + 0.01f;
    g_corr[3] = dotp( xn, y2, L_SUBFR )  - 0.02f;
    g_corr[4] = dotp( yy1, y2, L_SUBFR ) + 0.02f;

    /*-----------------------------------------------------------------*
     * calculate the unscaled innovation energy
     * calculate the predicted gain code
     * calculate optimal gains
     *-----------------------------------------------------------------*/

    Ecode = (dotp( code, code, L_SUBFR) + 0.01f) / L_SUBFR;
    *gain_inov = 1.0f / (float)sqrt( Ecode );
    Ei = 10 * (float)log10( Ecode );
    gcode0 = (float) pow(10, 0.05 * (Es_pred - Ei));

    tmp1 = (g_corr[0]*g_corr[2]) - (g_corr[4]*g_corr[4]);
    tmp2 = g_corr[1]/tmp1;
    tmp1 = g_corr[3]/tmp1;

    *gain_pit =  (g_corr[2]*tmp2) - (g_corr[4]*tmp1);
    *gain_code = (g_corr[0]*tmp1) - (g_corr[4]*tmp2);

    *gain_pit = max(G_PITCH_MIN,min(*gain_pit,G_PITCH_MAX));

    /*-----------------------------------------------------------------*
     * limit the pitch gain searching range (if indicated by clip_gain)
     *-----------------------------------------------------------------*/

    if( clip_gain == 1 && *gain_pit > 0.95f )
    {
        *gain_pit = 0.95f;
    }
    else if( clip_gain == 2 && *gain_pit > 0.65f )
    {
        *gain_pit = 0.65f;
    }

    /*-----------------------------------------------------------------*
     * search for the best quantized values
     *-----------------------------------------------------------------*/

    nBits_pitch = gain_bits_16kHz_tbl[ BIT_ALLOC_IDX_16KHZ(core_brate, coder_type, i_subfr, TC_SUBFR2IDX_16KHZ(tc_subfr)) ];
    nBits_code = (nBits_pitch+1)>>1;
    nBits_pitch = nBits_pitch>>1;

    tmp16 = div_s(1,((1 << (nBits_pitch)) - 1) );      /* Q15*/
    tmp1 = (float)mult_r((short)(G_PITCH_MAX*8192.0f+0.5f),tmp16)/8192.0f;
    index = usquant( *gain_pit, gain_pit, G_PITCH_MIN, tmp1, (1 << nBits_pitch) );
    push_indice( st, IND_GAIN_PIT, index, nBits_pitch );

    g_code = *gain_code / gcode0;
    index = gain_quant( &g_code, G_CODE_MIN, G_CODE_MAX, nBits_code );
    *gain_code = g_code * gcode0;
    push_indice( st, IND_GAIN_CODE, index, nBits_code);

    *norm_gain_code = *gain_code / *gain_inov;

    return;
}

/*-------------------------------------------------------------------*
 * gain_enc_gaus()
 *
 * Quantization of gain for Gaussian codebook
 *-------------------------------------------------------------------*/

short gain_enc_gaus(          /* o  : Return index of quantization  */
    float *gain,        /* i/o: Code gain to quantize         */
    const short bits,         /* i  : number of bits to quantize    */
    const float lowBound,     /* i  : lower bound of quantizer (dB) */
    const float topBound      /* i  : upper bound of quantizer (dB) */
)
{
    short index;
    float enr, stepSize;


    enr = (float)( 20.0 * log10( *gain + 0.001f ) );    /* codebook gain in dB  */

    /*-----------------------------------------------------------------*
     * quantize linearly the log E
     *-----------------------------------------------------------------*/

    stepSize = (topBound-lowBound)/((float)(1<<bits));
    index = (short)(((enr - lowBound)/stepSize)+0.5f);
    if( index >= (1<<bits) )
    {
        index = (1<<bits)-1;
    }

    if( index < 0 )
    {
        index = 0;
    }

    enr = (float)index*stepSize + lowBound;           /* quantized codebook gain in dB */
    *gain = (float)pow( 10.0f, enr/20.0f );           /* quantized codebook gain */
    return index;
}

/*-----------------------------------------------------------------*
 * gain_enc_tc()
 *
 * Search and quantization of gain_code for subframes (in the
 * beginning of frame) without pulses in TC - 3b coding.
 * In this case:
 * - gain_pit = 0
 * - gain_code - scalar quantization (no prediciton history used)
 *-----------------------------------------------------------------*/

void gain_enc_tc(
    Encoder_State *st,               /* i/o: encoder state structure      */
    const long  core_brate,        /* i  : core bitrate                                       */
    const short L_frame,           /* i  : length of the frame                                */
    const short i_subfr,           /* i  : subframe index                                     */
    const short tc_subfr,          /* i  : TC subframe index                                  */
    const float xn[],              /* i  : target vector                                      */
    const float y2[],              /* i  : zero-memory filtered algebraic codebook excitation */
    const float code[],            /* i  : algebraic excitation                               */
    const float Es_pred,           /* i  : predicted scaled innovation energy                 */
    float *gain_pit,         /* o  : Pitch gain / Quantized pitch gain                  */
    float *gain_code,        /* o  : quantized codebook gain                            */
    float *gain_inov,        /* o  : innovation gain                                    */
    float *norm_gain_code    /* o  : norm. gain of the codebook excitation                           */
)
{
    short i, index, nBits;
    float Ei, g_code, gcode0, Ecode;

    /*----------------------------------------------------------------*
     * get number of bits for gain quantization
     *----------------------------------------------------------------*/

    if( L_frame == L_FRAME )
    {
        nBits = gain_bits_tbl[BIT_ALLOC_IDX(core_brate, TRANSITION, i_subfr, TC_SUBFR2IDX(tc_subfr))];
    }
    else  /* L_frame == L_FRAME16k */
    {
        nBits = gain_bits_16kHz_tbl[BIT_ALLOC_IDX_16KHZ(core_brate, TRANSITION, i_subfr, TC_SUBFR2IDX_16KHZ(tc_subfr))];
    }

    /*----------------------------------------------------------------*
     * find the code pitch (for current subframe)
     *----------------------------------------------------------------*/

    *gain_code = dotp( xn, y2, L_SUBFR )/( dotp( y2, y2, L_SUBFR ) + 0.01f );

    /*----------------------------------------------------------------*
     * calculate the predicted gain code
     * decode codebook gain
     *----------------------------------------------------------------*/

    *gain_pit = 0;

    Ecode = (dotp( code, code, L_SUBFR) + 0.01f) / L_SUBFR;
    *gain_inov = 1.0f / (float)sqrt( Ecode );
    Ei = 10 * (float)log10( Ecode );
    gcode0 = (float) pow(10, 0.05 * (Es_pred - Ei));

    if( nBits > 3 )
    {
        g_code = *gain_code / gcode0;
        index = gain_quant( &g_code, G_CODE_MIN, G_CODE_MAX, nBits );
        *gain_code = g_code * gcode0;
        push_indice( st, IND_GAIN_CODE, index, nBits );
    }
    else
    {
        index = N_GAIN_CODE_TC-1;
        for( i=0; i < N_GAIN_CODE_TC-1; i++ )
        {
            if( *gain_code < ((tbl_gain_code_tc[i]+(tbl_gain_code_tc[i+1]-tbl_gain_code_tc[i])/2) * gcode0) )
            {
                index = i;
                break;
            }
        }

        /*----------------------------------------------------------------*
         * 3-bit -> 2-bit encoding
         *----------------------------------------------------------------*/

        if( nBits == 2 )
        {
            index /= 2;
            *gain_code = tbl_gain_code_tc[index*2] * gcode0;
            push_indice( st, IND_GAIN_CODE, index, nBits );
        }
        else
        {
            *gain_code = tbl_gain_code_tc[index] * gcode0;
            push_indice( st, IND_GAIN_CODE, index, nBits );
        }
    }

    *norm_gain_code = *gain_code / *gain_inov;

    return;
}

/*---------------------------------------------------------------------*
 * E_corr_xy2()
 *
 * Find the correlations between the target xn[], the filtered adaptive
 * codebook exc. y1[], and the filtered fixed codebook innovation y2[].
 *  ( <y2,y2> , -2<xn,y2> and 2<y1,y2> )   (stored in g_corr[2..4])
 *---------------------------------------------------------------------*/

void E_corr_xy2(
    const float xn[],   /* i  : target vector                          */
    const float y1[],   /* i  : filtered excitation components 1       */
    const float y2[],   /* i  : filtered excitation components 2       */
    float g_corr[],/* o  : correlations between x, y1, y2, y3, y4 */
    const short L_subfr  /* i : subframe size                            */
)
{


    g_corr[2] = dotp(y2, y2, L_subfr);
    g_corr[3] = -2.0f * dotp(xn, y2, L_subfr);
    g_corr[4] = 2.0f * dotp(y1, y2, L_subfr);

    return;
}


/*---------------------------------------------------------------------*
 * gain_enc_lbr()
 *
 * Quantization of pitch and codebook gains without prediction (memory-less)
 * in ACELP at 7.2 and 8.0 kbps
 * - the gain codebooks and gain estimation constants are different in each subframe
 * - the estimated gain, gcode0, is first determined based on
 *   classification and/or previous quantized gains (from previous subframes in the current frame)
 * - a correction factor gamma = g_code / gcode0 is then vector quantized
 *   along with gain_pit
 * - the mean-squared error criterion is used for codebook search
 *---------------------------------------------------------------------*/

void gain_enc_lbr(
    Encoder_State *st,               /* i/o: encoder state structure      */
    const long  core_brate,        /* i  : core bitrate                                                    */
    const short coder_type,        /* i  : coding type                                                     */
    const short i_subfr,           /* i  : subframe index                                                  */
    const float *xn,               /* i  : target vector                                                   */
    const float *y1,               /* i  : zero-memory filtered adaptive excitation                        */
    const float *y2,               /* i  : zero-memory filtered algebraic codebook excitation              */
    const float *code,             /* i  : algebraic excitation                                            */
    float *gain_pit,         /* o  : quantized pitch gain                                            */
    float *gain_code,        /* o  : quantized codebook gain                                         */
    float *gain_inov,        /* o  : gain of the innovation (used for normalization)                 */
    float *norm_gain_code,   /* o  : norm. gain of the codebook excitation                           */
    float *g_corr,           /* i/o: correlations <y1,y1>, -2<xn,y1>,<y2,y2>, -2<xn,y2> and 2<y1,y2> */
    float gains_mem[],       /* i/o: pitch gain and code gain from previous subframes                */
    const short clip_gain          /* i  : gain pitch clipping flag (1 = clipping)                         */
)
{
    short index = 0, i, size, nBits, n_pred, ctype;
    float dist, dist_min, g_pitch, g_code, gcode0, aux[10], Ecode;
    short rf_flag = 0;
    const float *p, *b, *cdbk = 0;

    /*-----------------------------------------------------------------*
     * calculate the rest of the correlation coefficients
     * c2 = <y2,y2>, c3 = -2<xn,y2>, c4 = 2<y1,y2>, c5* = <xn,xn>
     * c5* - not necessary to calculate
     *-----------------------------------------------------------------*/

    E_corr_xy2( xn, y1, y2, g_corr,L_SUBFR);
    g_corr[2] += 0.01F;
    g_corr[3] -= 0.02F;
    g_corr[4] += 0.02F;
    Ecode = ( dotp( code, code, L_SUBFR ) + 0.01f ) / L_SUBFR;
    *gain_inov = 1.0f / (float)sqrt(Ecode);

    /*-----------------------------------------------------------------*
     * select the codebook, size and number of bits
     * set the gains searching range
     *-----------------------------------------------------------------*/

    nBits = gain_bits_tbl[BIT_ALLOC_IDX(core_brate, coder_type, i_subfr, 0)];
    size = (short)pow2[nBits];

    /*-----------------------------------------------------------------*
     * calculate prediction of gcode
     * search for the best codeword
     *-----------------------------------------------------------------*/

    ctype = 2*(coder_type - 1);
    if (i_subfr == 0)
    {
        b = b_1sfr;
        n_pred = 2;

        switch ( nBits )
        {
        case 8:
        {
            cdbk = gp_gamma_1sfr_8b;
            if ( clip_gain == 1 ) size -= 60;
            break;
        }
        case 7:
        {
            cdbk = gp_gamma_1sfr_7b;
            if ( clip_gain == 1 ) size -= 27;
            break;
        }
        case 6:
        {
            cdbk = gp_gamma_1sfr_6b;
            if ( clip_gain == 1 ) size -= 10;
            break;
        }
        }

        /* calculate predicted gain */
        aux[0] = 1.0f;
        aux[1] = ctype;
        gcode0 = (float)pow(10, dotp(b, aux, n_pred) - 0.5f * (float)log10(Ecode));

        /* searching of codebook */
        p = cdbk;
        dist_min = 3.402823466e+38F;
        index = 0;
        for (i = 0; i<size; i++)
        {
            g_pitch = *p++;
            g_code = gcode0 **p++;
            dist = g_pitch*g_pitch * g_corr[0]
                   + g_pitch         * g_corr[1]
                   + g_code*g_code   * g_corr[2]
                   + g_code          * g_corr[3]
                   + g_pitch*g_code  * g_corr[4];

            if (dist < dist_min)
            {
                dist_min = dist;
                index = i;
            }
        }

        *gain_pit  = cdbk[index*2];
        *gain_code = cdbk[index*2+1] * gcode0;

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
            if ( clip_gain == 1 ) size -= 30;
            break;
        }
        case 6:
        {
            cdbk = gp_gamma_2sfr_6b;
            if ( clip_gain == 1 ) size -= 12;
            break;
        }
        }

        /* calculate predicted gain */
        aux[0] = 1.0f;
        aux[1] = ctype;
        aux[2] = (float)log10(gains_mem[0]);
        aux[3] = gains_mem[3];
        gcode0 = (float)pow(10, dotp(b, aux, n_pred));

        /* searching of codebook */
        p = cdbk;
        dist_min = 3.402823466e+38F;
        index = 0;
        for (i = 0; i<size; i++)
        {
            g_pitch = *p++;
            g_code = gcode0 **p++;
            dist = g_pitch*g_pitch * g_corr[0]
                   + g_pitch         * g_corr[1]
                   + g_code*g_code   * g_corr[2]
                   + g_code          * g_corr[3]
                   + g_pitch*g_code  * g_corr[4];

            if (dist < dist_min)
            {
                dist_min = dist;
                index = i;
            }
        }
        *gain_pit  = cdbk[index*2];
        *gain_code = cdbk[index*2+1] * gcode0;

        gains_mem[1] = *gain_code;
        gains_mem[4] = *gain_pit;
    }
    else if (i_subfr == 2*L_SUBFR)
    {
        if ( rf_flag == 1 )
        {
            gains_mem[1] = gains_mem[0];
            gains_mem[4] = gains_mem[3];
        }

        b = b_3sfr;
        n_pred = 6;

        cdbk = gp_gamma_3sfr_6b;
        if ( clip_gain == 1 )
        {
            size -= 11;
        }

        /* calculate predicted gain */
        aux[0] = 1.0f;
        aux[1] = ctype;
        aux[2] = (float)log10(gains_mem[0]);
        aux[3] = (float)log10(gains_mem[1]);
        aux[4] = gains_mem[3];
        aux[5] = gains_mem[4];
        gcode0 = (float)pow(10, dotp(b, aux, n_pred));

        /* searching of codebook */
        p = cdbk;
        dist_min = 3.402823466e+38F;
        index = 0;
        for (i = 0; i<size; i++)
        {
            g_pitch = *p++;
            g_code = gcode0 **p++;
            dist = g_pitch*g_pitch * g_corr[0]
                   + g_pitch         * g_corr[1]
                   + g_code*g_code   * g_corr[2]
                   + g_code          * g_corr[3]
                   + g_pitch*g_code  * g_corr[4];

            if (dist < dist_min)
            {
                dist_min = dist;
                index = i;
            }
        }
        *gain_pit  = cdbk[index*2];
        *gain_code = cdbk[index*2+1] * gcode0;

        gains_mem[2] = *gain_code;
        gains_mem[5] = *gain_pit;
    }
    else if (i_subfr == 3*L_SUBFR)
    {
        b = b_4sfr;
        n_pred = 8;

        cdbk = gp_gamma_4sfr_6b;
        if ( clip_gain == 1 )
        {
            size -= 11;
        }

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

        /* searching of codebook */
        p = cdbk;
        dist_min = 3.402823466e+38F;
        index = 0;
        for (i = 0; i<size; i++)
        {
            g_pitch = *p++;
            g_code = gcode0 **p++;
            dist = g_pitch*g_pitch * g_corr[0]
                   + g_pitch         * g_corr[1]
                   + g_code*g_code   * g_corr[2]
                   + g_code          * g_corr[3]
                   + g_pitch*g_code  * g_corr[4];

            if (dist < dist_min)
            {
                dist_min = dist;
                index = i;
            }
        }
        *gain_pit  = cdbk[index*2];
        *gain_code = cdbk[index*2+1] * gcode0;
    }

    *norm_gain_code = *gain_code / *gain_inov;

    push_indice( st, IND_GAIN, index, nBits );

    return;
}
