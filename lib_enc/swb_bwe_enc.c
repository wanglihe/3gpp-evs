/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include "rom_enc.h"
#include "basop_util.h"
#include "basop_proto_func.h"

/*---------------------------------------------------------------------*
 * Local functions
 *---------------------------------------------------------------------*/
static short SWB_BWE_encoding( Encoder_State *st, const float *insig, const float *insig_lp, const float *insig_hp, const float *synth,
                               const float *yos, float *SWB_fenv, const float tilt_nb, const short st_offset, const short coder_type );
static void MSVQ_Interpol_Tran( float *SWB_env_energy, short *indice );
static void calculate_tonality( const float *org, const float *gen, float *SFM_org, float *SFM_gen, const short length );
static short WB_BWE_encoding( const short coder_type, const float *yos, float *WB_fenv, Encoder_State *st );
static void energy_control( Encoder_State *st, const short core, const short mode, const short coder_type, const float *org, const short offset, float *energy_factor );
static short decision_hq_generic_class (const float *coefs, const short hq_generic_offset );

/*-------------------------------------------------------------------*
 * wb_bwe_enc()
 *
 * WB BWE encoder
 *-------------------------------------------------------------------*/

void wb_bwe_enc(
    Encoder_State *st,                /* i/o: encoder state structure                  */
    const float *new_wb_speech,     /* i  : original input signal at 16kHz           */
    short coder_type          /* i  : coding type                              */
)
{
    float old_input[NS2SA(16000, DELAY_FD_BWE_ENC_NS + DELAY_FIR_RESAMPL_NS) + L_FRAME16k];
    float *new_input;                                /* pointer to original input signal         */
    float yorig[L_FRAME16k];                         /* MDCT spectrum of weighted original       */
    short mode = 0;
    float wtda_old_input[2*L_FRAME16k];
    float WB_fenv[SWB_FENV];
    short Sample_Delay_WB_BWE;

    if( st->total_brate == ACELP_13k20 )
    {
        /*---------------------------------------------------------------------*
         * Delay the original input signal to be synchronized with ACELP core synthesis
         *---------------------------------------------------------------------*/

        set_f( old_input, 0, NS2SA(16000, DELAY_FD_BWE_ENC_12k8_NS + DELAY_FIR_RESAMPL_NS) + L_FRAME16k );

        Sample_Delay_WB_BWE = NS2SA( 16000, DELAY_FD_BWE_ENC_12k8_NS );

        new_input = old_input + Sample_Delay_WB_BWE;
        mvr2r( st->old_input_wb, old_input, Sample_Delay_WB_BWE );
        mvr2r( new_wb_speech, new_input, L_FRAME16k );
        mvr2r( old_input + L_FRAME16k, st->old_input_wb, Sample_Delay_WB_BWE );

        /*---------------------------------------------------------------------*
         * WB BWE encoding
         *---------------------------------------------------------------------*/

        /* MDCT of the original signal */
        wtda( old_input, wtda_old_input, st->old_wtda_swb, ALDO_WINDOW, ALDO_WINDOW, L_FRAME16k );

        direct_transform( wtda_old_input, yorig, 0, L_FRAME16k );

        mode = WB_BWE_encoding( coder_type, yorig, WB_fenv, st );

        push_indice( st, IND_WB_CLASS, mode - 2, 1 );
    }


    st->prev_mode = mode;

    return;
}

/*-------------------------------------------------------------------*
 * get_normalize_spec()
 *
 *-------------------------------------------------------------------*/

static void get_normalize_spec(
    const short core,                   /* i  : core selected            */
    const short extl,                   /* i  : extension layer selected */
    const short mode,                   /* i  : SHB BWE class            */
    const short core_type,             /* i  : coding type              */
    const float *org,                   /* i  : input spectrum           */
    float *SWB_signal,            /* o  : output spectrum          */
    short *prev_L_swb_norm,       /* i  : previous norm. len       */
    const short offset                  /* i  : frequency offset         */
)
{
    short n_freq, L_swb_norm;
    float envelope[L_FRAME32k];
    short frq_end;

    set_f( SWB_signal, 0, HQ_GENERIC_HIGH0+offset );
    calc_normal_length(core, org, mode, extl, &L_swb_norm, prev_L_swb_norm);
    if( extl == SWB_BWE || extl == FB_BWE )
    {
        if( mode == HARMONIC )
        {
            mvr2r( org, &SWB_signal[240+offset], 240 );
            mvr2r( &org[128], &SWB_signal[480+offset], 80 );
        }
        else
        {
            mvr2r( &org[112], &SWB_signal[240+offset], 128 );
            mvr2r( &org[112], &SWB_signal[368+offset], 128 );
            mvr2r( &org[176], &SWB_signal[496+offset], 64 );
        }
        frq_end = 560+offset;
    }
    else if (extl == WB_BWE)
    {
        if( core_type == 0 )
        {
            mvr2r(&org[160], &SWB_signal[240], 80);
        }
        else
        {
            mvr2r(&org[80], &SWB_signal[240], 80);
        }
        frq_end = L_FRAME16k;
    }
    else
    {
        mvr2r( org+HQ_GENERIC_OFFSET, SWB_signal+HQ_GENERIC_HIGH0+offset, HQ_GENERIC_LEN0 );
        mvr2r( org+HQ_GENERIC_OFFSET, SWB_signal+HQ_GENERIC_HIGH1+offset, HQ_GENERIC_LEN0 );
        if ( offset == HQ_GENERIC_FOFFSET_24K4 )
        {
            mvr2r( org+HQ_GENERIC_LOW0, SWB_signal+HQ_GENERIC_HIGH2+offset, HQ_GENERIC_END_FREQ-HQ_GENERIC_HIGH2 );
        }
        frq_end = L_FRAME32k;
    }

    /* calculate envelope */
    calc_norm_envelop( SWB_signal, envelope, L_swb_norm, frq_end - offset, offset );

    /* Normalize with envelope */
    for( n_freq = swb_bwe_subband[0]+offset; n_freq<frq_end; n_freq++ )
    {
        SWB_signal[n_freq] /= envelope[n_freq];
    }

    return;
}

/*-------------------------------------------------------------------*
 * WB_BWE_fenv_q()
 *
 * Scalar quantizer routine
 *-------------------------------------------------------------------*/

static short WB_BWE_fenv_q(    /* o:   quantized gain index  */
    float *x,                  /* i/o: energy of WB envelop  */
    const float *cb,                 /* i:   quantizer codebook    */
    const short cb_length,           /* i:   length of codebook    */
    const short cb_dim               /* i:   dimension of codebook */
)
{
    short i, j, indx = 0;
    float dist, min_dist;
    const float *pit = cb;

    min_dist = FLT_MAX;
    for ( i = 0; i < cb_length; i++ )
    {
        dist = 0.0f;
        for(j=0; j<cb_dim; j++)
        {
            dist += (x[j]-(*pit)) * (x[j]-(*pit));
            pit++;
        }

        if( dist < min_dist)
        {
            min_dist = dist;
            indx = i;
        }
    }

    for(j=0; j<cb_dim; j++)
    {
        x[j] = cb[cb_dim*indx+j];
    }

    return (indx);
}

/*-------------------------------------------------------------------*
 * FD_BWE_class()
 *
 * classify signal of above 6.4kHz, can be used for WB/SWB switch
 *-------------------------------------------------------------------*/

static short FD_BWE_class(     /* o  : FD BWE class        */
    const float *fSpectrum,          /* i  : input spectrum      */
    const float fGain,               /* i  : global gain         */
    const float tilt_nb,             /* i  : BWE tilt            */
    Encoder_State *st                  /* i/o: Encoder structure   */
)
{
    short i, j, k, noise, sharpMod = 0;
    float peak, mean[20], mag;
    float sharpPeak;
    const float *input_hi = 0;
    float sharp;
    float gain_tmp = 0;
    short mode;
    float meanH, mean_d = 0;
    short sharplimit;
    short numsharp;
    short numharmonic;

    mode = NORMAL;
    k = 0;
    noise = 0;
    sharpPeak = 0;
    numsharp = 0;
    numharmonic = 4;
    sharplimit = 10;

    if ( st->extl == SWB_BWE || st->extl == FB_BWE )
    {
        input_hi = &fSpectrum[256];
        numsharp = NUM_SHARP;
        if ( ( st->last_extl == SWB_BWE && st->extl == SWB_BWE ) || ( st->last_extl == FB_BWE && st->extl == FB_BWE ) )
        {
            gain_tmp = fGain / (st->prev_global_gain + EPSILON);
            if (st->prev_mode == TRANSIENT)
            {
                numharmonic = numharmonic * 2;
            }
            else if (st->prev_mode == NORMAL || st->prev_mode == NOISE)
            {
                numharmonic = 3 * numharmonic / 2;
            }
        }
        else
        {
            gain_tmp = 1;
            if (st->prev_mode == HARMONIC)
            {
                numharmonic = numharmonic / 2;
                sharplimit = sharplimit / 2;
            }
            else
            {
                numharmonic = numharmonic * 2;
                sharplimit = sharplimit * 2;
            }
        }
    }
    else if (st->extl == WB_BWE)
    {
        input_hi = &fSpectrum[224];
        numsharp = NUM_SHARP / 3;
        if (st->prev_mode == HARMONIC)
        {
            numharmonic = numharmonic / 4;
        }
        else
        {
            numharmonic = numharmonic / 2;
        }
        if (st->last_extl != WB_BWE)
        {
            if (st->prev_mode == HARMONIC)
            {
                sharplimit = sharplimit / 2;
            }
            else
            {
                sharplimit = sharplimit * 2;
            }
        }
    }

    meanH = EPSILON;
    for(i = 0; i < numsharp; i ++)
    {
        peak = 0.0f;
        mean[i] = 0;

        for(j = 0; j < SHARP_WIDTH; j ++)
        {
            mag = (float) fabs(*input_hi);
            if (mag > peak)
            {
                peak = mag;
            }
            mean[i] += mag;
            input_hi ++;
        }
        meanH += mean[i];

        if(mean[i] != peak)
        {
            sharp = (float) (peak * (SHARP_WIDTH - 1) / (mean[i] - peak));
        }
        else
        {
            sharp = 0.0f;
        }

        if (sharp > 4.5 && peak > 8)
        {
            k += 1;
        }
        else if (sharp < 3.0)
        {
            noise += 1;
        }

        if (sharp > sharpPeak)
        {
            sharpPeak = sharp;
        }
    }

    if ( st->extl == SWB_BWE || st->extl == FB_BWE )
    {
        if(k >= numharmonic && gain_tmp > 0.5f && gain_tmp < 1.8f && sharpPeak > sharplimit)
        {
            sharpMod = 1;
        }
        else
        {
            sharpMod = 0;
        }
        meanH /= 288;
        mean_d = 0.0f;
        for(i=0; i<NUM_SHARP; i++)
        {
            mean_d += (float)fabs(mean[i]/32 - meanH);
        }
    }
    else if (st->extl == WB_BWE)
    {
        if (k >= numharmonic && sharpPeak > sharplimit)
        {
            sharpMod = 1;
        }
        else
        {
            sharpMod = 0;
        }
    }

    if (sharpMod && st->modeCount < 12)
    {
        st->modeCount++;
    }
    else if (sharpMod == 0 && st->modeCount > 0)
    {
        st->modeCount--;
    }

    if (st->modeCount >= 2)
    {
        sharpMod = 1;
    }

    if (sharpMod)
    {
        mode = HARMONIC;
    }
    else if ( st->extl == SWB_BWE || st->extl == FB_BWE )
    {
        if (noise > 4 && mean_d < 4.8f*meanH && tilt_nb < 5)
        {
            mode = NOISE;
        }
    }

    return mode;
}

/*-------------------------------------------------------------------*
 * WB_BWE_encoding()
 *
 * WB BWE main encoder
 *-------------------------------------------------------------------*/

static short WB_BWE_encoding(  /* o  : classification of wb signal            */
    const short coder_type,
    const float *yos,                /* i  : MDCT coefficients of weighted original */
    float *WB_fenv,            /* i/o: energy of WB envelope                  */
    Encoder_State *st                  /* i/o: Encoder structure                      */
)
{
    short mode;
    float global_gain;
    float energy;
    short i, n_coeff, n_band;
    short index;
    float energy_factor[4];

    /* Energy for the different bands and global energies */
    global_gain = EPSILON;

    for (i = 0, n_band = 0; i < 2; i++)
    {
        energy = EPSILON;
        for (n_coeff = swb_bwe_subband[n_band]; n_coeff < swb_bwe_subband[n_band+2]; n_coeff++)
        {
            energy += yos[n_coeff] * yos[n_coeff];
        }

        WB_fenv[i] = energy;
        n_band += 2;
        global_gain += energy;
    }

    mode = FD_BWE_class(yos, global_gain, 0, st);

    energy_control( st, ACELP_CORE, mode, coder_type, yos, 0, energy_factor );

    for (i = 0; i < 2; i++)
    {
        WB_fenv[i] = (float)(log10( WB_fenv[i]*energy_factor[i<<1]/40 )*FAC_LOG2);
    }

    index = WB_BWE_fenv_q( WB_fenv, F_2_5, 32, 2 );

    push_indice( st, IND_WB_FENV, index, 5 );


    return (mode);

}


/*-------------------------------------------------------------------*
 * swb_bwe_enc()
 *
 * SWB BWE encoder (only for 32kHz signals)
 *-------------------------------------------------------------------*/

void swb_bwe_enc(
    Encoder_State *st,                /* i/o: encoder state structure                 */
    const float *old_input_12k8,    /* i  : input signal @12.8kHz for SWB BWE       */
    const float *old_input_16k,     /* i  : input signal @16kHz for SWB BWE         */
    const float *old_syn_12k8_16k,  /* i  : ACELP core synthesis at 12.8kHz or 16kHz*/
    const float *new_swb_speech,    /* i  : original input signal at 32kHz          */
    const float *shb_speech,        /* i  : SHB target signal (6-14kHz) at 16kHz    */
    const short coder_type          /* i  : coding type                             */
)
{
    short i, inner_frame, idxGain;
    float *new_input;
    long  inner_Fs;
    float old_input[NS2SA(48000, DELAY_FD_BWE_ENC_NS + DELAY_FIR_RESAMPL_NS) + L_FRAME48k];
    float old_input_lp[L_FRAME16k], new_input_hp[L_FRAME16k];
    float yorig[L_FRAME48k];
    float wtda_old_input[2*L_FRAME48k];
    float SWB_fenv[SWB_FENV];
    float tilt_nb;
    short Sample_Delay_SWB_BWE, Sample_Delay_HP, Sample_Delay_LP;
    float ener_low, energy_fbe_fb, fb_ener_adjust, ener_adjust_quan;

    ener_adjust_quan = 0.0f;
    idxGain = 0;

    /*---------------------------------------------------------------------*
     * Delay the original input signal to be synchronized with ACELP core synthesis
     *---------------------------------------------------------------------*/

    if( st->extl == FB_BWE )
    {
        inner_frame = L_FRAME48k;
        inner_Fs = 48000;
    }
    else
    {
        inner_frame = L_FRAME32k;
        inner_Fs = 32000;
    }

    set_f( old_input, 0.0f, NS2SA(inner_Fs, DELAY_FD_BWE_ENC_12k8_NS + DELAY_FIR_RESAMPL_NS) + inner_frame );

    if( st->L_frame == L_FRAME )
    {
        Sample_Delay_SWB_BWE = NS2SA(inner_Fs, DELAY_FD_BWE_ENC_12k8_NS + DELAY_FIR_RESAMPL_NS);
        Sample_Delay_HP      = NS2SA(16000, ACELP_LOOK_NS + DELAY_FD_BWE_ENC_12k8_NS + DELAY_FIR_RESAMPL_NS - DELAY_CLDFB_NS);
        Sample_Delay_LP      = NS2SA(12800, ACELP_LOOK_NS + DELAY_FD_BWE_ENC_12k8_NS);

        mvr2r( st->old_input_lp, old_input_lp, Sample_Delay_LP );
        mvr2r( old_input_12k8 + L_INP_MEM, &old_input_lp[Sample_Delay_LP], L_FRAME-Sample_Delay_LP );
        mvr2r( old_input_12k8 + L_INP_MEM + L_FRAME - Sample_Delay_LP, st->old_input_lp, Sample_Delay_LP );
    }
    else
    {
        Sample_Delay_SWB_BWE = NS2SA(inner_Fs, DELAY_FD_BWE_ENC_16k_NS + DELAY_FIR_RESAMPL_NS);
        Sample_Delay_HP      = NS2SA(16000, ACELP_LOOK_NS + DELAY_FD_BWE_ENC_16k_NS + DELAY_FIR_RESAMPL_NS - DELAY_CLDFB_NS);
        Sample_Delay_LP      = NS2SA(16000, ACELP_LOOK_NS + DELAY_FD_BWE_ENC_16k_NS);

        mvr2r( st->old_input_lp, old_input_lp, Sample_Delay_LP );
        mvr2r( old_input_16k + L_INP_MEM, &old_input_lp[Sample_Delay_LP], L_FRAME16k-Sample_Delay_LP );
        mvr2r( old_input_16k + L_INP_MEM + L_FRAME16k - Sample_Delay_LP, st->old_input_lp, Sample_Delay_LP );
    }

    mvr2r( st->new_input_hp, new_input_hp, Sample_Delay_HP );
    mvr2r( shb_speech, new_input_hp + Sample_Delay_HP, L_FRAME16k - Sample_Delay_HP );
    mvr2r( shb_speech + L_FRAME16k - Sample_Delay_HP, st->new_input_hp, Sample_Delay_HP );

    new_input = old_input + Sample_Delay_SWB_BWE;
    mvr2r( st->old_input, old_input, Sample_Delay_SWB_BWE );
    mvr2r( new_swb_speech, new_input, inner_frame );
    mvr2r( old_input + inner_frame, st->old_input, Sample_Delay_SWB_BWE );

    /*----------------------------------------------------------------------*
     * Calculate tilt of the input signal and the ACELP core synthesis
     *----------------------------------------------------------------------*/

    calc_tilt_bwe( old_input_lp, &tilt_nb, L_FRAME );

    /*---------------------------------------------------------------------*
     * SWB BWE encoding
     * FB BWE encoding
     *---------------------------------------------------------------------*/

    /* windowing of the original input signal */
    wtda( old_input, wtda_old_input, st->old_wtda_swb, ALDO_WINDOW,ALDO_WINDOW,inner_frame );

    /* DCT of the original input signal */
    direct_transform( wtda_old_input, yorig, 0, inner_frame );

    /* FB BWE encoding */
    if ( st->extl == FB_BWE )
    {
        energy_fbe_fb = sum2_f( yorig + FB_BAND_BEGIN, FB_BAND_WIDTH ) + EPSILON;
        ener_low = EPSILON;
        for( i=FB_BAND_BEGIN - FB_BAND_WIDTH; i<FB_BAND_BEGIN; i++)
        {
            ener_low += yorig[i] * yorig[i];
        }

        fb_ener_adjust = (float)sqrt(energy_fbe_fb/ener_low);
        fb_ener_adjust = min(fb_ener_adjust, FB_MAX_GAIN_VAR);
        idxGain = (short)usquant(fb_ener_adjust, &ener_adjust_quan, FB_GAIN_QLOW, FB_GAIN_QDELTA, 1<<NUM_BITS_FB_FRAMEGAIN);
    }

    /* SWB BWE encoding */
    if (st->L_frame == L_FRAME16k)
    {
        SWB_BWE_encoding( st, old_input, old_input_lp, new_input_hp, old_syn_12k8_16k, yorig, SWB_fenv, tilt_nb, 80, coder_type );
    }
    else
    {
        SWB_BWE_encoding( st, old_input, old_input_lp, new_input_hp, old_syn_12k8_16k, yorig, SWB_fenv, tilt_nb, 6, coder_type );
    }

    /* write FB BWE frame gain to the bitstream */
    if( st->extl == FB_BWE )
    {
        push_indice( st, IND_FB_SLOPE, idxGain, NUM_BITS_FB_FRAMEGAIN );
    }


    return;
}

/*-------------------------------------------------------------------*
 * Freq_weights()
 *
 *-------------------------------------------------------------------*/
static
void freq_weights(
    const float Band_Ener[],            /* i  : Band energy                 */
    const float f_weighting[],          /* i  : weigting coefs.             */
    float w[],                    /* o  : Freq. weighting             */
    const short Nb                      /* i  : Number of bands             */
)
{
    short i;
    float tmp, w1[SWB_FENV], w2[SWB_FENV];
    float min_b, max_b;

    /* Find Max band energy */
    min_b = Band_Ener[0];
    max_b = Band_Ener[0];
    for( i=1; i<Nb; i++ )
    {
        if( Band_Ener[i] < min_b )
        {
            min_b = Band_Ener[i];
        }
        if( Band_Ener[i] > max_b )
        {
            max_b = Band_Ener[i];
        }
    }

    /* Find weighting function */
    tmp = 1.f/(max_b-min_b);
    for( i=0; i<Nb; i++ )
    {
        w1[i] = (Band_Ener[i]-min_b)*tmp + 1.f;         /*1<= var <=2 */
        w2[i] = f_weighting[i];                         /*1~0.75*/
        w[i] = w1[i]* w2[i];
    }

    return;
}


/*-------------------------------------------------------------------*
 * VQwithCand_w()
 *
 *-------------------------------------------------------------------*/
static
void vqWithCand_w(
    const float *x,                        /* i  : input vector                                */
    const float *E_ROM_dico,            /* i  : codebook                                    */
    const short dim,                    /* i  : codebook dimension                          */
    const short E_ROM_dico_size,        /* i  : codebook size                               */
    short *index,                    /* o  : survivors indices                           */
    const short surv,                    /* i  : survivor number                             */
    float dist_min[],             /* o  : minimum distortion                          */
    const float *w,                     /* i  : weighting                                   */
    const short flag                    /* i  : flag indicationg weighted distortion metric */
)
{
    short i, j, k, l;
    const float *p_E_ROM_dico;
    float dist, temp1;

    if( flag )
    {
        set_f( dist_min, 3.402823466e+38F, surv );  /* FLT_MAX */

        for( i = 0; i < surv; i++ )
        {
            index[i] = i;
        }

        p_E_ROM_dico = E_ROM_dico;

        for( i = 0; i < E_ROM_dico_size; i++ )
        {
            dist = x[0] - *p_E_ROM_dico++;
            dist *= (dist * w[0]);

            for( j = 1; j < dim; j++ )
            {
                temp1 = x[j] - *p_E_ROM_dico++;
                dist += temp1 * temp1 * w[j];
            }

            for( k = 0; k < surv; k++ )
            {
                if( dist < dist_min[k] )
                {
                    for( l = surv - 1; l > k; l-- )
                    {
                        dist_min[l] = dist_min[l - 1];
                        index[l] = index[l - 1];
                    }
                    dist_min[k] = dist;
                    index[k] = i;
                    break;
                }
            }
        }
    }
    else
    {
        set_f( dist_min, 3.402823466e+38F, surv );  /* FLT_MAX */

        for (i = 0; i < surv; i++)
        {
            index[i] = i;
        }

        p_E_ROM_dico = E_ROM_dico;

        for( i = 0; i < E_ROM_dico_size; i++ )
        {
            dist = x[0] - *p_E_ROM_dico++;
            dist *= dist;

            for( j = 1; j < dim; j++ )
            {
                temp1 = x[j] - *p_E_ROM_dico++;
                dist += temp1 * temp1;
            }

            for( k = 0; k < surv; k++ )
            {
                if( dist < dist_min[k] )
                {
                    for( l = surv - 1; l > k; l-- )
                    {
                        dist_min[l] = dist_min[l-1];
                        index[l] = index[l-1];
                    }
                    dist_min[k] = dist;
                    index[k] = i;
                    break;
                }
            }
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * vqSimple_w()
 *
 *-------------------------------------------------------------------*/
static
short vqSimple_w(
    const float *x,                 /* i  : input for quantizer                         */
    float *y,                 /* i  : quantized value                             */
    const float *cb,                /* i  : codebooks                                   */
    const float *w,                 /* i  : weight                                      */
    const short dim,                /* i  : dimension                                   */
    const short l,                  /* i  : number of candidates                        */
    const short flag                /* i  : flag indicationg weighted distortion metric */
)
{
    short i, j, index;
    const float *cbP;
    float dist_min, dist, temp;

    index = 0;
    dist_min = FLT_MAX;
    cbP = cb;

    if( flag )
    {
        for( i = 0; i < l; i++ )
        {
            dist = x[0] - *cbP++;
            dist *= (dist * w[0]);
            for( j = 1; j < dim; j++ )
            {
                temp = x[j] - *cbP++;
                dist += temp * temp * w[j];
            }

            if( dist < dist_min )
            {
                dist_min = dist;
                index = i;
            }
        }
    }
    else
    {
        for( i = 0; i < l; i++ )
        {
            dist = x[0] - *cbP++;
            dist *= dist;
            for( j = 1; j < dim; j++ )
            {
                temp = x[j] - *cbP++;
                dist += temp * temp;
            }

            if( dist < dist_min )
            {
                dist_min = dist;
                index = i;
            }
        }
    }

    /* Reading the selected vector */
    mvr2r( &cb[index * dim], y, dim );

    return(index);
}

/*-------------------------------------------------------------------*
 * MSVQ_Interpol_Tran()
 *
 *-------------------------------------------------------------------*/

static void MSVQ_Interpol_Tran(
    float *SWB_env_energy,  /* i/o  : (original/quantized) energy   */
    short *indice           /* o    : quantized index               */
)
{
    short k, n_band, candInd[N_CAND_TR], ind_tmp[2];
    float env_temp11[SWB_FENV_TRANS/2], env_temp12[SWB_FENV_TRANS/2];
    float dist, minDist, tmp_q;
    float quant_tmp[SWB_FENV_TRANS], quant_tmp2[SWB_FENV_TRANS];
    float distCand[N_CAND_TR], quant_select[SWB_FENV_TRANS];

    /* Extract target vector */
    for( n_band = 0; n_band < DIM_TR1; n_band++ )
    {
        env_temp11[n_band] = SWB_env_energy[2*n_band];
        env_temp12[n_band] = SWB_env_energy[2*n_band+1];
    }

    vqWithCand_w( env_temp11, Env_TR_Cdbk1, DIM_TR1, N_CB_TR1, candInd, N_CAND_TR, distCand, NULL, 0 );

    minDist = FLT_MAX;
    for( k=0; k<N_CAND_TR; k++ )
    {
        for( n_band = 0; n_band < DIM_TR1; n_band++ )
        {
            quant_tmp[n_band] = Env_TR_Cdbk1[candInd[k] * DIM_TR1 + n_band];
        }

        for( n_band = 0; n_band < DIM_TR2-1; n_band++ )
        {
            quant_tmp2[n_band] = env_temp12[n_band] - ((quant_tmp[n_band]+quant_tmp[n_band+1])/2.f);
        }
        quant_tmp2[n_band] = env_temp12[n_band] - quant_tmp[n_band];
        ind_tmp[0] = vqSimple_w( quant_tmp2, quant_tmp2, Env_TR_Cdbk2, NULL, DIM_TR2, N_CB_TR2, 0 );

        for( n_band = 0; n_band < DIM_TR1; n_band++ )
        {
            quant_select[n_band*2] = quant_tmp[n_band];
        }

        for( n_band = 0; n_band < DIM_TR2-1; n_band++ )
        {
            quant_select[n_band*2+1] = ((quant_tmp[n_band]+quant_tmp[n_band+1])/2.f) + quant_tmp2[n_band];
        }
        quant_select[n_band*2+1] = quant_tmp[n_band]+quant_tmp2[n_band];

        dist = 0.f;
        for( n_band = 0; n_band < SWB_FENV_TRANS; n_band++ )
        {
            tmp_q = SWB_env_energy[n_band] - quant_select[n_band];
            dist += tmp_q*tmp_q;
        }

        /* Check optimal candidate */
        if( dist < minDist )
        {
            minDist = dist;
            indice[0] = candInd[k];
            indice[1] = ind_tmp[0];
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * MSVQ_Interpol()
 *
 *-------------------------------------------------------------------*/
static
void msvq_interpol(
    float *SWB_env_energy,    /* i/o: (original/quantized) energy */
    float *w_env,             /* i/o: weighting coffecients       */
    short *indice             /* o  : quantized index             */
)
{
    short k, n_band, candInd[N_CAND], ind_tmp[4];
    float dist, minDist, tmp_q;
    float env_temp11[SWB_FENV/2], env_temp12[SWB_FENV/2];
    float quant_tmp[SWB_FENV], quant_tmp1[SWB_FENV], quant_tmp2[SWB_FENV],distCand[N_CAND];
    float quant_select[SWB_FENV], w_env11[SWB_FENV/2], w_env12[SWB_FENV/2];
    float synth_energy[SWB_FENV];

    /* Extract target vector */
    for(n_band = 0; n_band < DIM11; n_band++)
    {
        env_temp11[n_band] = SWB_env_energy[2*n_band];
        env_temp12[n_band] = SWB_env_energy[2*n_band+1];

        w_env11[n_band] = w_env[2*n_band];
        w_env12[n_band] = w_env[2*n_band+1];
    }

    vqWithCand_w( env_temp11, EnvCdbk11, DIM11, N_CB11, candInd, N_CAND, distCand, w_env11, 1 );

    minDist = FLT_MAX;
    for( k=0; k<N_CAND; k++ )
    {
        for( n_band = 0; n_band < DIM11; n_band++ )
        {
            quant_tmp1[n_band] = EnvCdbk11[candInd[k] * DIM11 + n_band];
            quant_tmp2[n_band] = env_temp11[n_band] - quant_tmp1[n_band];
        }

        ind_tmp[0] = vqSimple_w( quant_tmp2, quant_tmp2, EnvCdbk1st, w_env11, DIM1ST, N_CB1ST, 1 );
        ind_tmp[1] = vqSimple_w( quant_tmp2+DIM1ST, quant_tmp2+DIM1ST, EnvCdbk2nd, w_env11+DIM1ST, DIM2ND, N_CB2ND, 1 );

        /* Extract vector for odd position */
        for( n_band = 0; n_band < DIM11; n_band++ )
        {
            quant_tmp[n_band] = quant_tmp1[n_band] + quant_tmp2[n_band];
        }

        for( n_band = 0; n_band < DIM12-1; n_band++ )
        {
            quant_tmp2[n_band] = env_temp12[n_band] - ((quant_tmp[n_band]+quant_tmp[n_band+1])/2.f);
        }

        quant_tmp2[n_band] = env_temp12[n_band]-quant_tmp[n_band];

        ind_tmp[2] = vqSimple_w( quant_tmp2, quant_tmp2, EnvCdbk3rd, w_env12, DIM3RD, N_CB3RD, 1 );
        ind_tmp[3] = vqSimple_w( quant_tmp2+DIM3RD, quant_tmp2+DIM3RD, EnvCdbk4th, w_env12+DIM3RD, DIM4TH, N_CB4TH, 1 );

        for( n_band = 0; n_band < DIM11; n_band++ )
        {
            quant_select[n_band*2] = quant_tmp[n_band];
        }

        for( n_band = 0; n_band < DIM12-1; n_band++ )
        {
            quant_select[n_band*2+1] = ((quant_tmp[n_band]+quant_tmp[n_band+1])/2.f) + quant_tmp2[n_band];
        }
        quant_select[n_band*2+1] = quant_tmp[n_band] + quant_tmp2[n_band];

        dist = 0.f;
        for( n_band = 0; n_band < SWB_FENV; n_band++ )
        {
            tmp_q = SWB_env_energy[n_band] - quant_select[n_band];
            tmp_q = tmp_q*tmp_q;
            dist += tmp_q*w_env[n_band];
        }

        /* Check optimal candidate */
        if( dist < minDist )
        {
            minDist = dist;
            mvr2r( quant_select, synth_energy, SWB_FENV );

            indice[0] = candInd[k];
            indice[1] = ind_tmp[0];
            indice[2] = ind_tmp[1];
            indice[3] = ind_tmp[2];
            indice[4] = ind_tmp[3];
        }
    }

    mvr2r( synth_energy, SWB_env_energy, SWB_FENV );

    return;
}

static
void msvq_interpol_2(
    float *hq_generic_fenv,   /* i/o: (original/quantized) energy */
    const float *w_env,             /* i  : weighting coffecients       */
    short *indice,            /* o  : quantized index             */
    const short nenv                /* i  : the number of envelopes     */
)
{
    short k, n_band, candInd[N_CAND], ind_tmp[4];
    float dist, minDist, tmp_q;
    float env_temp11[SWB_FENV/2], env_temp12[SWB_FENV/2];
    float quant_tmp[SWB_FENV], quant_tmp1[SWB_FENV], quant_tmp2[SWB_FENV],distCand[N_CAND];
    float quant_select[SWB_FENV], w_env11[SWB_FENV/2], w_env12[SWB_FENV/2];
    float synth_energy[SWB_FENV];

    /* Extract target vector */
    for(n_band = 0; n_band < DIM11-1; n_band++)
    {
        env_temp11[n_band] = hq_generic_fenv[2*n_band];
        w_env11[n_band] = w_env[2*n_band];
    }
    env_temp11[DIM11-1] = hq_generic_fenv[2*(DIM11-2)+1];
    w_env11[DIM11-1] = w_env[2*(DIM11-2)+1];

    env_temp12[0] = hq_generic_fenv[0];
    w_env12[0] = w_env[0];
    for(n_band = 1; n_band < DIM11-1; n_band++)
    {
        env_temp12[n_band] = hq_generic_fenv[2*n_band-1];
        w_env12[n_band] = w_env[2*n_band-1];
    }

    vqWithCand_w( env_temp11, EnvCdbk11, DIM11, N_CB11, candInd, N_CAND, distCand, w_env11, 1 );

    minDist = FLT_MAX;
    for( k=0; k<N_CAND; k++ )
    {
        for( n_band = 0; n_band < DIM11; n_band++ )
        {
            quant_tmp1[n_band] = EnvCdbk11[candInd[k] * DIM11 + n_band];
            quant_tmp2[n_band] = env_temp11[n_band] - quant_tmp1[n_band];
        }

        ind_tmp[0] = vqSimple_w( quant_tmp2, quant_tmp2, EnvCdbk1st, w_env11, DIM1ST, N_CB1ST, 1 );
        ind_tmp[1] = vqSimple_w( quant_tmp2+DIM1ST, quant_tmp2+DIM1ST, EnvCdbk2nd, w_env11+DIM1ST, DIM2ND, N_CB2ND, 1 );

        /* Extract vector for odd position */
        for( n_band = 0; n_band < DIM11; n_band++ )
        {
            quant_tmp[n_band] = quant_tmp1[n_band] + quant_tmp2[n_band];
        }

        quant_tmp2[0] = env_temp12[0] - quant_tmp[0];
        for( n_band = 1; n_band < DIM12-1; n_band++ )
        {
            quant_tmp2[n_band] = env_temp12[n_band] - ((quant_tmp[n_band-1]+quant_tmp[n_band])/2.f);
        }

        ind_tmp[2] = vqSimple_w( quant_tmp2, quant_tmp2, EnvCdbk3rd, w_env12, DIM3RD, N_CB3RD, 1 );
        ind_tmp[3] = vqSimple_w( quant_tmp2+DIM3RD, quant_tmp2+DIM3RD, EnvCdbk3rd, w_env12+DIM3RD, DIM3RD, N_CB3RD, 1 );

        for( n_band = 0; n_band < DIM12-1; n_band++ )
        {
            quant_select[n_band*2] = quant_tmp[n_band];
        }
        quant_select[11] = quant_tmp[DIM12-1];

        quant_select[0] += quant_tmp2[0];
        for( n_band = 1; n_band < DIM12-1; n_band++ )
        {
            quant_select[n_band*2-1] = ((quant_tmp[n_band-1]+quant_tmp[n_band])/2.f) + quant_tmp2[n_band];
        }

        dist = 0.f;
        for( n_band = 0; n_band < SWB_FENV-2; n_band++ )
        {
            tmp_q = hq_generic_fenv[n_band] - quant_select[n_band];
            tmp_q = tmp_q*tmp_q;
            dist += tmp_q*w_env[n_band];
        }

        /* Check optimal candidate */
        if( dist < minDist )
        {
            minDist = dist;
            mvr2r( quant_select, synth_energy, SWB_FENV-2 );
            synth_energy[SWB_FENV-2] = 0;
            synth_energy[SWB_FENV-1] = 0;

            indice[0] = candInd[k];
            indice[1] = ind_tmp[0];
            indice[2] = ind_tmp[1];
            indice[3] = ind_tmp[2];
            indice[4] = ind_tmp[3];
        }
    }

    mvr2r( synth_energy, hq_generic_fenv, nenv );

    return;
}

/*-------------------------------------------------------------------*
 * calculate_tonality()
 *
 * Calculate tonality
 *-------------------------------------------------------------------*/

static void calculate_tonality(
    const float *org,               /* i  : MDCT coefficients of original           */
    const float *gen,               /* i  : MDCT coefficients of generated signal   */
    float *SFM_org,           /* o  : Spectral Flatness results               */
    float *SFM_gen,           /* o  : Spectral Flatness results               */
    const short length              /* i  : length for calculating tonality         */
)
{
    short n_coeff;
    float am_org, am_gen, gm_org, gm_gen;
    float inv_len, max, mult;
    float org_spec[80], gen_spec[80];

    /* to reduce dynamic range of original spectrum */
    max = EPSILON;
    for( n_coeff=0; n_coeff<length; n_coeff++ )
    {
        org_spec[n_coeff] = (float)fabs(org[n_coeff]);

        if( max<org_spec[n_coeff] )
        {
            max = org_spec[n_coeff];
        }
    }
    mult = 25.f/max;

    for( n_coeff=0; n_coeff<length; n_coeff++ )
    {
        org_spec[n_coeff] *= mult;
    }

    max = EPSILON;
    for( n_coeff=0; n_coeff<length; n_coeff++ )
    {
        gen_spec[n_coeff] = (float)fabs(gen[n_coeff]);

        if( max < gen_spec[n_coeff] )
        {
            max = gen_spec[n_coeff];
        }
    }
    mult = 25.f/max;

    for( n_coeff=0; n_coeff<length; n_coeff++ )
    {
        gen_spec[n_coeff] *= mult;
    }

    inv_len = 1.f/(float)length;

    am_org = EPSILON;
    am_gen = EPSILON;
    gm_org = 1.f;
    gm_gen = 1.f;

    for( n_coeff = 0; n_coeff<length; n_coeff++ )
    {
        am_org += org_spec[n_coeff];
        am_gen += gen_spec[n_coeff];
        gm_org *= org_spec[n_coeff];
        gm_gen *= gen_spec[n_coeff];
    }

    *SFM_org = 10.f * ((float)log10(am_org*inv_len)-inv_len*(float)log10(gm_org));
    *SFM_org = max( 0.0001f, min(*SFM_org, 5.993f) );
    *SFM_gen = 10.f * ((float)log10(am_gen*inv_len)-inv_len*(float)log10(gm_gen));
    *SFM_gen = max( 0.0001f,min(*SFM_gen, 5.993f) );

    return;
}

/*-------------------------------------------------------------------*
 * energy_control()
 *
 *-------------------------------------------------------------------*/

static void energy_control(
    Encoder_State *st,                /* i/o: encoder structure   */
    const short core,               /* i  : core                */
    const short mode,               /* i  : SHB BWE class       */
    const short coder_type,         /* i  : SHB BWE class       */
    const float *org,               /* i  : input spectrum      */
    const short offset,             /* i  : frequency offset    */
    float *energy_factor      /* o  : energy factor       */
)
{
    short n_band;
    float gamma;
    short core_type;
    float SWB_signal[L_FRAME32k], SFM_org[SWB_FENV], SFM_gen[SWB_FENV];
    short max_band=SWB_FENV,band_step=1;

    if( core == ACELP_CORE )
    {
        gamma = 0.35f;
        if( coder_type != AUDIO && st->total_brate <= ACELP_8k00 )
        {
            core_type = 0;
        }
        else
        {
            core_type = 1;
        }

        get_normalize_spec( core, st->extl, mode, core_type, org, SWB_signal, &(st->prev_L_swb_norm1), offset );

        if ( st->extl == WB_BWE)
        {
            max_band = 4;
            band_step = 2;
        }
    }
    else  /* HQ core */
    {
        gamma = 0.55f;
        get_normalize_spec( core, -1, mode, -1, org, SWB_signal, &(st->prev_L_swb_norm1), offset );

        if ( offset == HQ_GENERIC_FOFFSET_32K )
        {
            max_band = 12;
        }
    }

    for( n_band=0; n_band<max_band; n_band+=band_step )
    {
        calculate_tonality( org+swb_bwe_subband[n_band]+offset, SWB_signal+swb_bwe_subband[n_band]+offset,
                            &SFM_org[n_band], &SFM_gen[n_band], swb_bwe_subband[n_band+band_step]-swb_bwe_subband[n_band] );

        if( SFM_gen[n_band] < 0.75*SFM_org[n_band] )
        {
            energy_factor[n_band] = (SFM_gen[n_band]/SFM_org[n_band]);

            if( energy_factor[n_band] < gamma )
            {
                energy_factor[n_band] = gamma;
            }
        }
        else
        {
            energy_factor[n_band] = 1.0f;
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * SWB_BWE_encoding()
 *
 * SWB BWE encoder
 *-------------------------------------------------------------------*/

static short SWB_BWE_encoding(
    Encoder_State *st,                 /* i/o: Encoder state structure */
    const float *insig,              /* i  : delayed original input signal at 32kHz */
    const float *insig_lp,           /* i  : delayed original lowband input signal at 16kHz */
    const float *insig_hp,           /* i  : delayed original highband input signal at 16kHz */
    const float *synth,              /* i  : delayed ACELP core synthesis at 12.8kHz  */
    const float *yos,                /* i  : MDCT coefficients of the windowed original input signal at 32kHz */
    float *SWB_fenv,           /* o  : frequency-domain quantized BWE envelope */
    const float tilt_nb,             /* i  : SWB tilt */
    const short st_offset,           /* i  : start frequency offset for BWE envelope */
    const short coder_type           /* i  : coding type                              */
)
{
    short IsTransient, mode;
    short index;
    float SWB_tenv_tmp[SWB_TENV];
    float SWB_tenv[SWB_TENV];
    float global_gain;
    float energy;
    float max;
    short i, n_coeff, n_band, pos, indice[6];
    float tilt, WB_tenv_orig, WB_tenv_syn, Rat_tenv;
    float energy_factor[SWB_FENV], w_env[SWB_FENV];
    short L;
    short IsTransient_LF;

    if(st->L_frame == L_FRAME )
    {
        L = L_SUBFR;
    }
    else
    {
        L = L_SUBFR16k;
    }

    /* HF transient detect */
    IsTransient = detect_transient( insig_hp, st, L_FRAME16k, coder_type );

    /* LF transient detect */
    IsTransient_LF = 0;
    for ( n_band = 0; n_band < 4; n_band++ )
    {
        energy = EPSILON;
        for ( i = 0; i < L; i++ )
        {
            energy += insig_lp[i + n_band*L] * insig_lp[i + n_band*L];
        }

        if( energy > 5.5f * st->EnergyLF )
        {
            IsTransient_LF = 1;
        }

        st->EnergyLF = energy;
    }
    calc_tilt_bwe(insig, &tilt, L_FRAME32k);
    if( IsTransient == 1 && (tilt > 8.0 || st->clas > 1) )
    {
        IsTransient = 0;
        st->TransientHangOver = 0;
    }

    if( IsTransient == 1 )
    {
        mode = IsTransient;
        push_indice( st, IND_SWB_CLASS, mode, 2 );

        /* Energy for the different bands and global energies */
        global_gain = 0;
        for (n_band = 0; n_band < SWB_FENV_TRANS; n_band++)
        {
            energy = EPSILON;
            for (n_coeff = swb_bwe_trans_subband[n_band]+st_offset; n_coeff < swb_bwe_trans_subband[n_band+1]+st_offset; n_coeff++)
            {
                energy += yos[n_coeff] * yos[n_coeff];
            }
            global_gain += energy;
            SWB_fenv[n_band] = energy;
        }
        global_gain *= 0.5f;

        for (n_band = 0; n_band < SWB_FENV_TRANS; n_band++)
        {
            SWB_fenv[n_band] = 10.0f * (float)log10( SWB_fenv[n_band]/swb_bwe_trans_subband_width[n_band] ) - Mean_env_tr[n_band];
        }

        WB_tenv_orig = EPSILON;
        WB_tenv_syn = EPSILON;
        for(n_band = 0; n_band < SWB_TENV; n_band++)
        {
            SWB_tenv[n_band] = EPSILON;

            for(i = 0; i < L_SUBFR16k; i++)
            {
                SWB_tenv[n_band] += insig_hp[i + n_band*L_SUBFR16k] * insig_hp[i + n_band*L_SUBFR16k];
            }

            for(i=0; i<L; i++)
            {
                WB_tenv_syn += synth[i + n_band*L] * synth[i + n_band*L];
                WB_tenv_orig += insig_lp[i + n_band*L] * insig_lp[i + n_band*L];
            }

            SWB_tenv[n_band] = (float)(sqrt(SWB_tenv[n_band]*INV_L_SUBFR16k));
        }

        Rat_tenv = (float)sqrt(WB_tenv_syn / WB_tenv_orig);

        if(Rat_tenv < 0.5)
        {
            Rat_tenv *= 1.2f;
        }
        else if (Rat_tenv > 1)
        {
            Rat_tenv = 1.0f;
        }

        for(n_band = 0; n_band < SWB_TENV; n_band++)
        {
            SWB_tenv[n_band] *= Rat_tenv;
        }

        max = SWB_tenv[0];
        pos = 0;
        for(n_band = 1; n_band < SWB_TENV; n_band++)
        {
            if(SWB_tenv[n_band] > max)
            {
                max = SWB_tenv[n_band];
                pos = n_band;
            }
        }

        max = SWB_tenv[0];
        for(n_band = 1; n_band < SWB_TENV; n_band++)
        {
            if(SWB_tenv[n_band] > 5.0f*SWB_tenv[n_band-1])
            {
                break;
            }
        }

        if(n_band < SWB_TENV)
        {
            energy = 0.0f;
            for(n_band = (pos+1); n_band < SWB_TENV; n_band++)
            {
                energy += SWB_tenv[n_band];
            }
            if(pos == SWB_TENV-1)
            {
                energy = 0.0f;
            }
            else
            {
                energy /= (SWB_TENV-pos-1);
            }

            for(n_band = 0; n_band < pos; n_band++)
            {
                SWB_tenv[n_band] *= 0.5f;
            }

            SWB_tenv[pos] *= 1.005f;
            if(energy < SWB_tenv[pos])
            {
                for(n_band = pos+1; n_band < SWB_TENV; n_band++)
                {
                    SWB_tenv[n_band] *= 0.9f;
                }
            }
        }
        else
        {
            for(n_band = 1; n_band < SWB_TENV; n_band++)
            {
                if(SWB_tenv[n_band-1] > SWB_tenv[n_band])
                {
                    SWB_tenv[n_band-1] = 0.5f*(SWB_tenv[n_band-1]+SWB_tenv[n_band]);
                }
                else
                {
                    SWB_tenv[n_band] = 0.5f*(SWB_tenv[n_band-1]+SWB_tenv[n_band]);
                }
            }

            for(n_band = 0; n_band < SWB_TENV; n_band++)
            {
                SWB_tenv[n_band] *= 0.9f;
            }
        }

        if(IsTransient_LF == 0 && coder_type == INACTIVE && st->TransientHangOver == 1)
        {
            for(n_band = 0; n_band < SWB_TENV; n_band++)
            {
                SWB_tenv[n_band] *= 0.5f;
            }
            for(n_band = 0; n_band < SWB_FENV_TRANS; n_band++)
            {
                SWB_fenv[n_band] *= 0.05f;
            }
        }
        else
        {
            SWB_fenv[2] *= 0.1f;
            SWB_fenv[3] *= 0.05f;
        }

        for(n_band = 0; n_band < SWB_TENV; n_band++)
        {
            SWB_tenv_tmp[n_band] = (float) log10( SWB_tenv[n_band] + EPSILON ) * FAC_LOG2;
            if (SWB_tenv_tmp[n_band] > 15)
            {
                index = 15;
            }
            else if (SWB_tenv_tmp[n_band] < 0)
            {
                index = 0;
            }
            else
            {
                index = (short)(SWB_tenv_tmp[n_band]+0.5f);
            }

            push_indice( st, IND_SWB_TENV, index, 4 );

        }

        MSVQ_Interpol_Tran(SWB_fenv, indice);

        push_indice( st, IND_SWB_FENV, indice[0], 7 );
        push_indice( st, IND_SWB_FENV, indice[1], 6 );

    }
    else
    {
        /* Energy for the different bands and global energies */
        global_gain = 0;
        for (n_band = 0; n_band < SWB_FENV; n_band++)
        {
            energy = EPSILON;
            for (n_coeff = swb_bwe_subband[n_band]+st_offset; n_coeff < swb_bwe_subband[n_band+1]+st_offset; n_coeff++)
            {
                energy += yos[n_coeff] * yos[n_coeff];
            }

            if (n_band<SWB_FENV-2)
            {
                global_gain += energy;
            }
            SWB_fenv[n_band] = energy;
        }

        global_gain *= 0.5f;

        mode = FD_BWE_class(yos, global_gain, tilt_nb, st);

        push_indice( st, IND_SWB_CLASS, mode, 2 );

        energy_control( st, ACELP_CORE, mode, -1, yos, st_offset, energy_factor );

        for (n_band = 0; n_band < SWB_FENV; n_band++)
        {
            SWB_fenv[n_band] *= energy_factor[n_band];
            SWB_fenv[n_band] = 10.0f * (float)log10( SWB_fenv[n_band]*swb_inv_bwe_subband_width[n_band] );
        }

        freq_weights(SWB_fenv, w_NOR, w_env, SWB_FENV);

        for (n_band = 0; n_band < SWB_FENV; n_band++)
        {
            SWB_fenv[n_band] -= Mean_env[n_band];
        }

        /* Energy VQ */
        msvq_interpol( SWB_fenv, w_env, indice );

        push_indice( st, IND_SWB_FENV, indice[0], 5 );
        push_indice( st, IND_SWB_FENV, indice[1], 7 );
        push_indice( st, IND_SWB_FENV, indice[2], 6 );
        push_indice( st, IND_SWB_FENV, indice[3], 5 );
        push_indice( st, IND_SWB_FENV, indice[4], 6 );

    }

    st->prev_mode = mode;
    st->prev_global_gain = global_gain;

    return mode;
}

static short decision_hq_generic_class (
    const float *coefs,           /* i: original MDCT spectrum                      */
    const short hq_generic_offset /* i: frequency offset of high frequency spectrum */
)
{
    short i, k;
    float p, a, e;
    float p2a;
    float avgp2a;
    short nband;

    if ( hq_generic_offset == HQ_GENERIC_FOFFSET_24K4 )
    {
        nband = 10;
    }
    else
    {
        nband = 8;
    }

    avgp2a = 0.f;
    for (k = 0; k < nband; k++)
    {
        a = 0.0f;
        p = 0.0f;
        for (i = swb_bwe_subband[k]+hq_generic_offset; i <swb_bwe_subband[k+1]+hq_generic_offset; i++)
        {
            e = coefs[i] * coefs[i];

            if (e > p)
            {
                p = e;
            }

            a += e;
        }

        if (a > 0.0f)
        {
            a *= swb_inv_bwe_subband_width[k];
            p2a = 10.0f * (float) log10 (p / a);
            avgp2a += p2a;
        }
    }

    avgp2a /= (float)(nband);

    if ( avgp2a > 8.6f )
    {
        return HQ_GENERIC_EXC1;
    }
    else
    {
        return HQ_GENERIC_EXC0;
    }
}

/*-------------------------------------------------------------------*
 * hq_generic_hf_encoding()
 *
 *-------------------------------------------------------------------*/

void hq_generic_hf_encoding(
    const float *coefs,                     /* i  : MDCT coefficients of weighted original      */
    float *hq_generic_fenv,           /* i/o: energy of SWB envelope                      */
    const short hq_generic_offset,          /* i  : frequency offset for extracting energy      */
    Encoder_State *st,                        /* i/o: encoder state structure                     */
    short *hq_generic_exc_clas        /* o  : HF excitation class                        */
)
{
    short n_coeff, n_band;
    float energy;
    float energy_factor[SWB_FENV], w_env[SWB_FENV];
    short indice[HQ_GENERIC_NVQIDX];
    short nenv;

    if ( hq_generic_offset <= HQ_GENERIC_FOFFSET_24K4 )
    {
        nenv = SWB_FENV;
    }
    else
    {
        nenv = SWB_FENV-2;
    }

    for( n_band = 0; n_band < nenv; n_band++ )
    {
        energy = EPSILON;
        for( n_coeff = swb_bwe_subband[n_band]+hq_generic_offset; n_coeff < swb_bwe_subband[n_band+1] + hq_generic_offset; n_coeff++ )
        {
            energy += coefs[n_coeff] * coefs[n_coeff];
        }

        hq_generic_fenv[n_band] = energy;
    }

    if ( st->bwidth == FB )
    {
        for( n_band = 0; n_band < DIM_FB; n_band++ )
        {
            energy = EPSILON;
            for( n_coeff = fb_bwe_subband[n_band]; n_coeff < fb_bwe_subband[n_band+1]; n_coeff++ )
            {
                energy += coefs[n_coeff] * coefs[n_coeff];
            }

            hq_generic_fenv[n_band+nenv] = energy;
        }
    }

    energy_control( st, HQ_CORE, -1, -1, coefs, hq_generic_offset, energy_factor );

    if ( st->hq_generic_speech_class == 1 )
    {
        push_indice( st, IND_HQ_SWB_EXC_SP_CLAS, 1, 1 );
        *hq_generic_exc_clas = HQ_GENERIC_SP_EXC;
    }
    else
    {
        *hq_generic_exc_clas = decision_hq_generic_class(coefs, hq_generic_offset);
        push_indice( st, IND_HQ_SWB_EXC_SP_CLAS, 0, 1 );
        push_indice( st, IND_HQ_SWB_EXC_CLAS, *hq_generic_exc_clas, 1 );
    }

    for( n_band = 0; n_band < nenv; n_band++ )
    {
        hq_generic_fenv[n_band] *= energy_factor[n_band];
        hq_generic_fenv[n_band] = 10.0f * (float)log10( hq_generic_fenv[n_band]*swb_inv_bwe_subband_width[n_band] );
    }

    if ( st->bwidth == FB )
    {
        for( n_band = 0; n_band < DIM_FB; n_band++ )
        {
            hq_generic_fenv[n_band+nenv] = 10.0f * (float)log10( hq_generic_fenv[n_band+nenv]*fb_inv_bwe_subband_width[n_band] );
        }
    }

    freq_weights( hq_generic_fenv, w_NOR, w_env, nenv );

    for( n_band = 0; n_band < nenv; n_band++ )
    {
        hq_generic_fenv[n_band] -= Mean_env[n_band];
    }

    if ( st->bwidth == FB )
    {
        for( n_band = 0; n_band < DIM_FB; n_band++ )
        {
            hq_generic_fenv[n_band+nenv] -= Mean_env_fb[n_band];
        }
    }

    /* Energy VQ */
    if ( hq_generic_offset <= HQ_GENERIC_FOFFSET_24K4 )
    {
        msvq_interpol( hq_generic_fenv, w_env, indice );
    }
    else
    {
        msvq_interpol_2( hq_generic_fenv, w_env, indice, nenv );
    }

    if ( st->bwidth == FB )
    {
        indice[5] = vqSimple_w( hq_generic_fenv+nenv, hq_generic_fenv+nenv, EnvCdbkFB, NULL, DIM_FB, N_CB_FB, 0 );
    }

    push_indice( st, IND_SWB_FENV_HQ, indice[0], 5 );
    push_indice( st, IND_SWB_FENV_HQ, indice[1], 7 );
    push_indice( st, IND_SWB_FENV_HQ, indice[2], 6 );
    push_indice( st, IND_SWB_FENV_HQ, indice[3], 5 );

    if ( hq_generic_offset <= HQ_GENERIC_FOFFSET_24K4 )
    {
        push_indice( st, IND_SWB_FENV_HQ, indice[4], 6 );
    }
    else
    {
        push_indice( st, IND_SWB_FENV_HQ, indice[4], 5 );
    }

    if ( st->bwidth == FB )
    {
        push_indice( st, IND_FB_FENV_HQ, indice[5], 5 );
    }

    for( n_band = 0; n_band < nenv; n_band++ )
    {
        Word16 tmp,frac,exp;
        Word32 L_tmp;
        tmp = add((short)(hq_generic_fenv[n_band]*256),(short)(Mean_env[n_band]*256)); /*Q8 */

        L_tmp = L_mult(tmp, 21771); /* 0.166096 in Q17 -> Q26 */
        L_tmp = L_shr(L_tmp, 10); /* From Q26 to Q16 */
        frac = L_Extract_lc(L_tmp, &exp); /* Extract exponent of L_tmp */

        tmp = extract_l(Pow2(13, frac));/* Put 13 as exponent so that */
        /* output of Pow2() will be: */
        /* 16384 < Pow2() <= 32767 */
        exp = sub(exp, 13);
        tmp = shl(tmp, add(exp,1)); /*Q1 */
        hq_generic_fenv[n_band] = (float)tmp*0.5f;; /*Q1 */
    }

    if ( st->bwidth == FB )
    {
        for( n_band = 0; n_band < DIM_FB; n_band++ )
        {
            Word16 tmp,frac,exp;
            Word32 L_tmp;

            tmp = add((short)(hq_generic_fenv[n_band + nenv]*128),(short)(Mean_env_fb[n_band]*128)); /*Q7 */
            L_tmp = L_mult(tmp, 21771); /* 0.166096 in Q17 -> Q25 */
            L_tmp = L_shr(L_tmp, 9); /* From Q25 to Q16 */
            frac = L_Extract_lc(L_tmp, &exp); /* Extract exponent of L_tmp */

            tmp = extract_l(Pow2(13, frac));/* Put 13 as exponent so that */
            /* output of Pow2() will be: */
            /* 16384 < Pow2() <= 32767 */
            exp = sub(exp, 13);
            tmp = shl(tmp, add(exp, 1));/*Q1 */
            hq_generic_fenv[add(n_band,nenv)] = (float)tmp*0.5f;
        }
    }

    return;
}
