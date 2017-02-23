/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*-------------------------------------------------------------------*
 * Local functions
 *-------------------------------------------------------------------*/

static short VDQ_vec( float *Qvec_out, const float *mean_dic, const float *dic,
                      const short index, const short vec_en );

/*-------------------------------------------------------------------*
 * Comp_and_apply_gain()
 *
 * Compute and apply the quantized per band gain
 *-------------------------------------------------------------------*/

void Comp_and_apply_gain(
    float exc_diffQ[],       /* i/o: Quantized excitation              */
    float Ener_per_bd_iQ[],  /* o  : Target ener per band              */
    float Ener_per_bd_yQ[],  /* o  : Ener per band for norm vector     */
    short Mbands_gn,         /* i  : number of bands                   */
    const short ReUseGain          /* i  : Reuse the gain in Ener_per_bd_yQ  */
)
{
    short i, i_band;
    short StartBin, NB_Qbins;
    float y_gain;

    /* Recreate excitation for local synthesis and decoder */
    StartBin = 0;
    NB_Qbins  = 0;
    for( i_band = 0; i_band < Mbands_gn; i_band++ )
    {
        StartBin += NB_Qbins;
        NB_Qbins = mfreq_bindiv_loc[i_band];
        if( ReUseGain == 1 )
        {
            y_gain = Ener_per_bd_yQ[i_band];
        }
        else
        {
            y_gain = (float)pow(10, (Ener_per_bd_iQ[i_band]-Ener_per_bd_yQ[i_band]));
            Ener_per_bd_yQ[i_band] = y_gain;
        }

        for(i = StartBin ; i < NB_Qbins + StartBin ; i++)
        {
            exc_diffQ[i] *= y_gain;
        }
    }

    return;
}

/*------------------------------------------------------------------*
 * Ener_per_band_comp()
 *
 * Compute the energy per band in log domain for quantization purposes
 * Loops are decomposed to accomodate the PVQ quantization
 *------------------------------------------------------------------*/

void Ener_per_band_comp(
    const float exc_diff[],    /* i  : target signal               */
    float y_gain4[],     /* o  : Energy per band to quantize */
    const short Mband,         /* i  : Max band                    */
    const short Eflag          /* i  : flag of highest band        */
)
{
    float etmp;
    const float *pt;
    short i,j;

    pt  = exc_diff;
    for(j = 0; j < 2; j++)
    {
        y_gain4[j] = 0;
        etmp = 0.02f;

        pt = exc_diff + j*8;
        for(i = 0; i < 8; i++)
        {
            etmp += (*pt **pt);
            pt++;
        }

        /* normalized to 16 bins to easy the quantization */
        y_gain4[j] = (float)log10(sqrt(2*etmp));
    }

    for(j = 1; j < Mband-2; j++)
    {
        etmp = 0.01f;

        pt = exc_diff + j*16;
        for(i = 0; i < 16; i++)
        {
            etmp += (*pt **pt);
            pt++;
        }

        y_gain4[j+1] = (float)log10(sqrt(etmp));
    }

    if( Eflag == 1 )
    {
        etmp = 0.01f;

        pt = exc_diff + j*16;
        for(i = 0; i < 32; i++)
        {
            etmp += (*pt **pt);
            pt++;
        }

        y_gain4[j+1] = (float)log10(sqrt(etmp/2));
    }

    return;
}

/*-------------------------------------------------------------------*
 * gsc_gainQ()
 *
 * Quantization of the energy per band
 *-------------------------------------------------------------------*/

float gsc_gainQ(
    Encoder_State *st,         /* i/o: encoder state structure      */
    const float y_gain4[],     /* i  : Energy per band              */
    float y_gainQ[],     /* o  : quantized energy per band    */
    const long  core_brate,    /* i  : Core rate                    */
    const short coder_type,    /* i  : coding type                  */
    const short bwidth         /* i  : input signal bandwidth       */
)
{
    float y_gain_tmp[MBANDS_GN], y_gain_tmp2[MBANDS_GN];
    short i, idx_g = 0;
    float mean_4g[1], ftmp1;
    float Gain_off = 0.0f;
    short Mbands_gn = MBANDS_GN;
    float y_gain_tmp3[MBANDS_GN];

    mean_4g[0] = 0;

    if( (coder_type == AUDIO || coder_type == INACTIVE) && bwidth == NB )
    {
        ftmp1 =  mean(y_gain4, 10) - 0.6f;
        for(i = 0; i < Mbands_gn; i++)
        {
            if(y_gain4[i] < ftmp1)
            {
                y_gain_tmp2[i] = ftmp1;
            }
            else
            {
                y_gain_tmp2[i] = y_gain4[i];
            }
        }

        /* Quantized mean gain without clipping */
        mean_4g[0] = mean(y_gain_tmp2, 10);
        idx_g = (short)vquant(mean_4g, Gain_meanNB, mean_4g, Gain_mean_dicNB, 1, 64);
        push_indice( st, IND_MEAN_GAIN2, idx_g, 6 );

        for(i = 0; i < Mbands_gn; i++)
        {
            y_gain_tmp[i] = y_gain_tmp2[i] - mean_4g[0];
        }

        if( y_gain_tmp[9] < -0.3f )
        {
            y_gain_tmp[9] = -0.3f;
        }

        set_f(y_gain_tmp+10, 0.0f, MBANDS_GN-10);
        idx_g = (short)vquant(y_gain_tmp, Mean_dic_NB, y_gain_tmp, Gain_dic1_NB, 3, 64);
        push_indice( st, IND_Y_GAIN_TMP, idx_g, 6 );

        if(core_brate < ACELP_9k60)
        {
            idx_g = (short)vquant(y_gain_tmp+3, Mean_dic_NB+3, y_gain_tmp+3, Gain_dic2_NB, 3, 32);
            push_indice( st, IND_Y_GAIN_TMP, idx_g, 5 );

            idx_g = (short)vquant(y_gain_tmp+6, Mean_dic_NB+6, y_gain_tmp+6, Gain_dic3_NB, 4, 16);
            push_indice( st, IND_Y_GAIN_TMP, idx_g, 4 );
        }
        else
        {
            idx_g = (short)vquant(y_gain_tmp+3, Mean_dic_NB+3, y_gain_tmp+3, Gain_dic2_NBHR, 3, 64);
            push_indice( st, IND_Y_GAIN_TMP, idx_g, 6 );

            idx_g = (short)vquant(y_gain_tmp+6, Mean_dic_NB+6, y_gain_tmp+6, Gain_dic3_NBHR, 4, 128);
            push_indice( st, IND_Y_GAIN_TMP, idx_g, 7 );
        }

        if( core_brate <= ACELP_9k60 && coder_type == INACTIVE )
        {
            /* Some energy is needed in high band for stat_noise_uv_enc() to be functional in inactive speech */
            y_gain_tmp[10] = mean(y_gain_tmp+6, 3);
            y_gain_tmp[11] = mean(y_gain_tmp+7, 3);
            y_gain_tmp[12] = mean(y_gain_tmp+8, 3);
            y_gain_tmp[13] = mean(y_gain_tmp+9, 3);
            y_gain_tmp[14] = mean(y_gain_tmp+10, 3);
            y_gain_tmp[15] = mean(y_gain_tmp+11, 3);
        }
        else
        {
            set_f( y_gain_tmp + 10, 0, MBANDS_GN - 10 );
        }
    }
    else
    {
        ftmp1 =  mean(y_gain4, 16);
        for(i = 0; i < Mbands_gn; i++)
        {
            if(y_gain4[i] < ftmp1-0.6f)
            {
                y_gain_tmp2[i] = ftmp1-.6f;
            }
            else if(y_gain4[i] > ftmp1+0.6f)
            {
                y_gain_tmp2[i] = ftmp1+0.6f;
            }
            else
            {
                y_gain_tmp2[i] = y_gain4[i];
            }
        }

        mean_4g[0] = mean(y_gain_tmp2, 16);
        idx_g = (short)vquant(mean_4g, mean_m, mean_4g, mean_gain_dic, 1, 64);
        push_indice( st, IND_MEAN_GAIN2, idx_g, 6 );

        /* Subtraction of the average gain */
        for(i = 0; i < Mbands_gn; i++)
        {
            y_gain_tmp[i] = y_gain_tmp2[i] - mean_4g[0];
        }

        if( core_brate < ACELP_9k60 )
        {
            /* prediction and quantization of the average gain */

            /*--------------------------------------------------------------------------------------*
             * Quantization of the first 8 bands
             * Keep only 4 bands out of the last 8 bands
             *--------------------------------------------------------------------------------------*/

            mvr2r(y_gain_tmp, y_gain_tmp2, 8);

            y_gain_tmp2[8] = y_gain_tmp[8];
            y_gain_tmp2[9] = y_gain_tmp[10];
            y_gain_tmp2[10] = y_gain_tmp[12];
            y_gain_tmp2[11] = y_gain_tmp[14];

            idx_g = 0;
            idx_g = (short)vquant(y_gain_tmp2, YGain_mean_LR, y_gain_tmp2, YGain_dic1_LR, 3, 32);
            push_indice( st, IND_Y_GAIN_TMP, idx_g, 5 );

            idx_g = (short)vquant(y_gain_tmp2+3, YGain_mean_LR+3, y_gain_tmp2+3, YGain_dic2_LR, 4, 32);
            push_indice( st, IND_Y_GAIN_TMP, idx_g, 5 );

            /*----------------------------------------------------------------------*
             * Vector quantization of the first 8 bands + quantization of the 4 bands out of the last 8
             * Interpolation of the last 4 bands Q to create bands 8-16
             *----------------------------------------------------------------------*/

            idx_g = (short)vquant(y_gain_tmp2+7, YGain_mean_LR+7, y_gain_tmp2+7, YGain_dic3_LR, 5, 32);
            push_indice( st, IND_Y_GAIN_TMP, idx_g, 5 );

            set_f(y_gain_tmp2+12, 0, MBANDS_GN-12);

            /* Update to quantized vector */
            mvr2r(y_gain_tmp2, y_gain_tmp, 8);

            mvr2r(y_gain_tmp2+8, y_gain_tmp3, 4);
            set_f(y_gain_tmp+8, 0,8);
            fft_rel(y_gain_tmp2+8, 4, 2);

            mvr2r(y_gain_tmp2+8, y_gain_tmp+8, 3);
            y_gain_tmp[15] = y_gain_tmp2[11];
            ifft_rel(y_gain_tmp+8, 8, 3);

            for(i = 8; i < 16; i++)
            {
                y_gain_tmp[i] *=  1.41f;
            }

            y_gain_tmp[8] = y_gain_tmp3[0];
            y_gain_tmp[10]= y_gain_tmp3[1];
            y_gain_tmp[12]= y_gain_tmp3[2];
            y_gain_tmp[14]= y_gain_tmp3[3];
        }
        else
        {
            idx_g = (short)vquant(y_gain_tmp, YG_mean16, y_gain_tmp, YG_dicMR_1, 4, 64);
            push_indice( st, IND_Y_GAIN_TMP, idx_g, 6 );

            idx_g = (short)vquant(y_gain_tmp+4, YG_mean16+4, y_gain_tmp+4, YG_dicMR_2, 4, 32);
            push_indice( st, IND_Y_GAIN_TMP, idx_g, 5 );

            idx_g = (short)vquant(y_gain_tmp+8, YG_mean16+8, y_gain_tmp+8, YG_dicMR_3, 4, 32);
            push_indice( st, IND_Y_GAIN_TMP, idx_g, 5 );

            idx_g = (short)vquant(y_gain_tmp+12, YG_mean16+12, y_gain_tmp+12, YG_dicMR_4, 4, 16);
            push_indice( st, IND_Y_GAIN_TMP, idx_g, 4 );
        }
    }

    /* Gain adjustment to fit ACELP generic inactive coding gain at low rate */
    if( coder_type == INACTIVE )
    {
        Gain_off = 0.0f;
        {
            if(core_brate <= ACELP_7k20 )
            {
                Gain_off = 8.f;   /* 0 dB */
            }
            else if (core_brate <= ACELP_8k00 )
            {
                Gain_off = 6.6f;   /* ~-3.3 dB */
            }
            else if (core_brate <= ACELP_9k60 )
            {
                Gain_off = 4.8f;   /* ~-2.4 dB */
            }
            else if (core_brate <= ACELP_11k60 )
            {
                Gain_off = 3.5f;   /* ~-2.4 dB */
            }
            else if (core_brate <= ACELP_13k20 )
            {
                Gain_off = 3.0f;   /* ~-2.4 dB */
            }
        }

    }

    if( coder_type != INACTIVE )
    {
        for( i = 0; i < Mbands_gn; i++ )
        {
            y_gainQ[i] = y_gain_tmp[i] + mean_4g[0];
        }
    }
    else
    {
        /*mimic ACELP decay of energy for low rates*/
        for( i = 0; i < Mbands_gn; i++ )
        {
            y_gainQ[i] = y_gain_tmp[i]+mean_4g[0]-(i*(Gain_off/20.f)/((float) Mbands_gn));
        }
    }

    return mean_4g[0];
}

/*-------------------------------------------------------------------*
 * gsc_gaindec()
 *
 * Generic signal frequency band decoding and application
 *-------------------------------------------------------------------*/

float gsc_gaindec(                /* o  : average frequency gain    */
    Decoder_State *st,              /* i/o: decoder state structure   */
    float y_gainQ[],        /* o  : quantized gain per band   */
    const long  core_brate,       /* i  : core used                 */
    float old_y_gain[],     /* i/o: AR gain quantizer for low rate */
    const short coder_type,       /* i  : coding type               */
    const short bwidth            /* i  : input signal bandwidth    */
)
{
    short idx_g, i;
    float mean_4g;
    float Gain_off = 0.0;
    short Mbands_gn = MBANDS_GN;
    float y_gain_tmp3[MBANDS_GN];

    if( (coder_type == AUDIO || coder_type == INACTIVE) && bwidth == NB )
    {
        idx_g = (short) get_next_indice( st, 6 );
        VDQ_vec(&mean_4g, Gain_meanNB, Gain_mean_dicNB, idx_g, 1 );

        idx_g = (short) get_next_indice( st, 6 );
        VDQ_vec(y_gainQ, Mean_dic_NB, Gain_dic1_NB, idx_g, 3 );

        if(core_brate < ACELP_9k60)
        {
            idx_g = (short) get_next_indice( st, 5 );
            VDQ_vec(y_gainQ+3, Mean_dic_NB+3, Gain_dic2_NB, idx_g, 3 );

            idx_g = (short) get_next_indice( st, 4 );
            VDQ_vec(y_gainQ+6, Mean_dic_NB+6, Gain_dic3_NB, idx_g, 4 );
        }
        else
        {
            idx_g = (short) get_next_indice( st, 6 );
            VDQ_vec(y_gainQ+3, Mean_dic_NB+3, Gain_dic2_NBHR, idx_g, 3 );

            idx_g = (short) get_next_indice( st, 7 );
            VDQ_vec(y_gainQ+6, Mean_dic_NB+6, Gain_dic3_NBHR, idx_g, 4 );
        }

        if( core_brate <= ACELP_9k60 && coder_type == INACTIVE )
        {
            /* Some energy is needed in high band for stat_noise_uv_enc
              to be functional in inactive speech */
            y_gainQ[10] = mean(y_gainQ+6, 3);
            y_gainQ[11] = mean(y_gainQ+7, 3);
            y_gainQ[12] = mean(y_gainQ+8, 3);
            y_gainQ[13] = mean(y_gainQ+9, 3);
            y_gainQ[14] = mean(y_gainQ+10, 3);
            y_gainQ[15] = mean(y_gainQ+11, 3);
        }
        else
        {
            set_f( y_gainQ + 10, 0, MBANDS_GN-10 );
        }
    }
    else
    {
        idx_g = (short) get_next_indice( st, 6 );
        VDQ_vec(&mean_4g, mean_m, mean_gain_dic, idx_g, 1 );

        if( core_brate < ACELP_9k60 )
        {
            /*--------------------------------------------------------------------------------------*
             * UQ of the first 8 bands and half of the last 8 bands
             *--------------------------------------------------------------------------------------*/

            idx_g = (short) get_next_indice( st, 5 );
            VDQ_vec(y_gainQ, YGain_mean_LR, YGain_dic1_LR, idx_g, 3 );

            idx_g = (short) get_next_indice( st, 5 );
            VDQ_vec(y_gainQ+3, YGain_mean_LR+3, YGain_dic2_LR, idx_g, 4 );

            /*----------------------------------------------------------------------*
             * Interpolation of the last 4 Q bands to create bands 8-16
             * And scaling
             *----------------------------------------------------------------------*/

            idx_g = (short) get_next_indice( st, 5 );
            VDQ_vec(y_gainQ+7, YGain_mean_LR+7, YGain_dic3_LR, idx_g, 5 );

            mvr2r(y_gainQ + 8, y_gain_tmp3, 4);
            set_f(y_gainQ + 12, 0.0f, 4);

            fft_rel(y_gainQ + 8, 4, 2);

            y_gainQ[15] = y_gainQ[11];
            y_gainQ[11] = 0.0f;

            ifft_rel(y_gainQ + 8, 8, 3);

            for(i = 8; i < 16; i++)
            {
                y_gainQ[i] *=  1.41f;
            }

            /*----------------------------------------------------------------------*
             * Copy the true Q values in the specific bands
             *----------------------------------------------------------------------*/

            y_gainQ[8] = y_gain_tmp3[0];
            y_gainQ[10]= y_gain_tmp3[1];
            y_gainQ[12]= y_gain_tmp3[2];
            y_gainQ[14]= y_gain_tmp3[3];
        }
        else
        {
            idx_g = (short) get_next_indice( st, 6 );
            VDQ_vec(y_gainQ, YG_mean16, YG_dicMR_1, idx_g, 4 );

            idx_g = (short) get_next_indice( st, 5 );
            VDQ_vec(y_gainQ+4, YG_mean16+4, YG_dicMR_2, idx_g, 4 );

            idx_g = (short) get_next_indice( st, 5 );
            VDQ_vec(y_gainQ+8, YG_mean16+8, YG_dicMR_3, idx_g, 4 );

            idx_g = (short) get_next_indice( st, 4 );
            VDQ_vec(y_gainQ+12, YG_mean16+12, YG_dicMR_4, idx_g, 4);
        }
    }

    /* Gain adjustment to fit ACELP generic inactive coding gain at low rate */
    if( coder_type == INACTIVE )
    {
        Gain_off = 0.0f;
        {
            if(core_brate <= ACELP_7k20 )
            {
                Gain_off = 8.f;   /* 0 dB */
            }
            else if (core_brate <= ACELP_8k00 )
            {
                Gain_off = 6.6f;   /* ~-3.3 dB */
            }
            else if (core_brate <= ACELP_9k60 )
            {
                Gain_off = 4.8f;   /* ~-2.4 dB */
            }
            else if (core_brate <= ACELP_11k60 )
            {
                Gain_off = 3.5f;   /* ~-2.4 dB */
            }
            else if (core_brate <= ACELP_13k20 )
            {
                Gain_off = 3.0f;   /* ~-2.4 dB */
            }
        }

    }

    if( coder_type != INACTIVE )
    {
        for( i = 0; i < Mbands_gn; i++ )
        {
            old_y_gain[i] = y_gainQ[i];
            y_gainQ[i] += mean_4g;
        }
    }
    else
    {
        /*mimic ACELP decay of energy for low rates*/
        for( i = 0; i < Mbands_gn; i++ )
        {
            old_y_gain[i] = y_gainQ[i];
            y_gainQ[i] += mean_4g-i*(Gain_off/20.f)/((float) Mbands_gn);
        }
    }

    return mean_4g;
}

/*-------------------------------------------------------------------*
 * VDQ_vec()
 *
 * Return the dequantized vector of index
 *-------------------------------------------------------------------*/

static short VDQ_vec(
    float *Qvec_out,      /* o:  Quanitzed vector */
    const float *mean_dic,      /* i:  average codebook */
    const float *dic,           /* i:  codebook         */
    const short index,          /* i:  index of codebook*/
    const short vec_en          /* i:  vector length    */
)
{
    short i, j;

    j =  index*vec_en;
    for ( i = 0; i < vec_en; i++)
    {
        Qvec_out[i] = dic[j++];
    }

    for(i = 0; i < vec_en; i++)
    {
        Qvec_out[i] += mean_dic[i];
    }

    return index;
}
