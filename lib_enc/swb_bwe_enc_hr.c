/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"

/*-------------------------------------------------------------------*
 * en_band_quant()
 *
 * Quantize the band envelop
 *-------------------------------------------------------------------*/

static short en_band_quant( /* o  : quantization index              */
    float *en_band,   /* i/o: (un)quantized envelope value    */
    const float *env_code,  /* i  : envelope codebook               */
    const short N           /* i  : codebook dimension              */
)
{
    float maxerr, err;
    short i, j, ind;

    maxerr = FLT_MAX;
    ind = 0;

    for( i = 0; i < N; i++ )
    {
        err = FLT_MIN;
        for (j = 0; j < 2; j++)
        {
            err += (en_band[j] - env_code[i*2+j]) * (en_band[j] - env_code[i*2+j]);
        }
        if (err < maxerr)
        {
            maxerr = err;
            ind = i;
        }
    }

    en_band[0] = env_code[2*ind];
    en_band[1] = env_code[2*ind+1];

    return( ind );
}

/*-------------------------------------------------------------------*
 * swb_bwe_enc_hr()
 *
 * HR SWB BWE encoder
 *-------------------------------------------------------------------*/

void swb_bwe_enc_hr(
    Encoder_State *st,            /* i/o: encoder state structure     */
    const float *new_input,     /* i  : input signal                */
    const short input_frame,    /* i  : frame length                */
    const short coder_type,     /* i  : coding type                 */
    const short unbits          /* i  : number of core unused bits  */
)
{
    short i, j, k, nBits, nBits_total, nBits_block, Nsv, Nsv2, width_noncoded;
    short is_transient, pos;
    int   x_norm[NSV_MAX*(WIDTH_BAND+1)], x_norm1[NSV_MAX*(WIDTH_BAND+1)];
    float t_audio[L_FRAME48k], t_audio_tmp[L_FRAME48k];
    float gain, gain2, en_band[N_BANDS_BWE_HR];
    short ind1, ind2;
    short nq[NSV_MAX], nq2[NSV_MAX];
    float tmpF, min_env;
    float en_noncoded;

    /*---------------------------------------------------------------------*
     * initializations
     *---------------------------------------------------------------------*/

    ind2 = 0;       /* only to suppress warnings */
    Nsv2 = 0;       /* only to suppress warnings */
    gain2 = 0 ;     /* only to suppress warnings */
    en_noncoded = 0;/* only to suppress warnings */

    /* reset memories in case that last frame was a different technology */
    if( st->last_core == HQ_CORE || st->last_extl != st->extl )
    {
        set_f( st->old_wtda_swb, 0, L_FRAME48k );
    }

    /* calculate SWB BWE bit-budget (extension layer bit-rate + AVQ unused bits from the core layer) */
    nBits = (short)(st->extl_brate) / 50 + unbits;
    nBits_total = nBits;


    /*---------------------------------------------------------------------*
     * detect transient frames
     *---------------------------------------------------------------------*/

    is_transient = detect_transient( new_input, st, input_frame, coder_type );
    push_indice( st, IND_HR_IS_TRANSIENT, is_transient, 1 );

    /*---------------------------------------------------------------------*
     * OLA and MDCT
     *---------------------------------------------------------------------*/

    wtda( new_input, t_audio_tmp, st->old_wtda_swb, ALDO_WINDOW, ALDO_WINDOW, input_frame );

    direct_transform( t_audio_tmp, t_audio, is_transient, input_frame );

    if( is_transient )
    {
        nBits = -1;     /* is_transient flag */
        nBits_block = nBits_total / NUM_TIME_SWITCHING_BLOCKS;
        nBits += nBits_total % NUM_TIME_SWITCHING_BLOCKS;

        /* set width of noncoded (blind estimated) spectrum */
        if( st->extl == SWB_BWE_HIGHRATE )
        {
            width_noncoded = L_FRAME32k/NUM_TIME_SWITCHING_BLOCKS - NUM_TRANS_END_FREQ_COEF;
        }
        else  /* st->extl == FB_BWE_HIGHRATE */
        {
            width_noncoded = (2*END_FREQ_BWE_FULL_FB/50)/NUM_TIME_SWITCHING_BLOCKS - NUM_TRANS_END_FREQ_COEF;
        }

        /*---------------------------------------------------------------------*
         * transient frames: processing in blocks (subframes)
         *---------------------------------------------------------------------*/

        for( k = 0; k < NUM_TIME_SWITCHING_BLOCKS; k++ )
        {
            nBits += nBits_block;

            /* compute energy of noncoded (14.4-20kHz) spectrum */
            if( st->extl == FB_BWE_HIGHRATE )
            {
                en_noncoded = sum2_f( t_audio + k*input_frame/NUM_TIME_SWITCHING_BLOCKS + NUM_TRANS_END_FREQ_COEF, width_noncoded ) + 0.001f;
                en_noncoded = (float)sqrt( en_noncoded / width_noncoded );
            }

            /* keep only frequencies in interest */
            set_f( t_audio + k*input_frame/NUM_TIME_SWITCHING_BLOCKS, 0, NUM_TRANS_START_FREQ_COEF );
            set_f( t_audio + k*input_frame/NUM_TIME_SWITCHING_BLOCKS + L_FRAME32k/NUM_TIME_SWITCHING_BLOCKS, 0, (input_frame-L_FRAME32k)/NUM_TIME_SWITCHING_BLOCKS );

            /*---------------------------------------------------------------------*
             * global gain coding
             *---------------------------------------------------------------------*/

            /* compute and quantize global energy */
            gain = sum2_f( t_audio + k*input_frame/NUM_TIME_SWITCHING_BLOCKS + NUM_TRANS_START_FREQ_COEF, WIDTH_TRANS_FREQ_COEF*N_BANDS_TRANS_BWE_HR ) + 0.001f;
            gain = (float)sqrt( gain ) / (WIDTH_TRANS_FREQ_COEF*N_BANDS_TRANS_BWE_HR);
            ind1 = gain_quant( &gain, MIN_GLOB_GAIN_BWE_HR, MAX_GLOB_GAIN_BWE_HR, NBITS_GLOB_GAIN_BWE_HR );

            push_indice( st, IND_HR_GAIN, ind1, NBITS_GLOB_GAIN_BWE_HR );
            nBits -= NBITS_GLOB_GAIN_BWE_HR;

            /* normalization with global gain  */
            tmpF = 1/gain;
            for( i=0; i<WIDTH_TRANS_FREQ_COEF*N_BANDS_TRANS_BWE_HR; i++ )
            {
                t_audio[NUM_TRANS_START_FREQ_COEF + k*input_frame/NUM_TIME_SWITCHING_BLOCKS + i] *= tmpF;
            }

            /*---------------------------------------------------------------------*
             * envelope coding
             *---------------------------------------------------------------------*/

            /* compute energy per band */
            for( i=0; i<N_BANDS_TRANS_BWE_HR; i++ )
            {
                en_band[i] = sum2_f( t_audio + k*input_frame/NUM_TIME_SWITCHING_BLOCKS + NUM_TRANS_START_FREQ_COEF + i*WIDTH_TRANS_FREQ_COEF, WIDTH_TRANS_FREQ_COEF ) + 0.001f;
                en_band[i] = (float)sqrt(en_band[i]/(WIDTH_TRANS_FREQ_COEF));
            }

            /* Q energy per band */
            if( k == 0 )
            {
                ind1 = en_band_quant( en_band, swb_hr_env_code3, NUM_ENVLOPE_CODE_HR_TR );
                push_indice( st, IND_HR_ENVELOPE, ind1, NBITS_ENVELOPE_BWE_HR_TR );
                nBits -= NBITS_ENVELOPE_BWE_HR_TR;
                ind2 = ind1;
            }
            else
            {
                if( ind2 < NUM_ENVLOPE_CODE_HR_TR2 )
                {
                    ind1 = en_band_quant( en_band, swb_hr_env_code3, NUM_ENVLOPE_CODE_HR_TR2 );
                }
                else
                {
                    ind1 = en_band_quant( en_band, swb_hr_env_code3 + (NUM_ENVLOPE_CODE_HR_TR2*2), NUM_ENVLOPE_CODE_HR_TR2 );
                }

                push_indice( st, IND_HR_ENVELOPE, ind1, NBITS_ENVELOPE_BWE_HR_TR - 1 );
                nBits -= (NBITS_ENVELOPE_BWE_HR_TR - 1);
            }

            /* normalize spectrum per bands */
            for( i = 0; i < N_BANDS_TRANS_BWE_HR; i++ )
            {
                tmpF = 1/en_band[i];
                for( j = 0; j < WIDTH_TRANS_FREQ_COEF; j++ )
                {
                    t_audio[k*input_frame/NUM_TIME_SWITCHING_BLOCKS + NUM_TRANS_START_FREQ_COEF + i*WIDTH_TRANS_FREQ_COEF + j] *= tmpF;
                }
            }

            /*---------------------------------------------------------------------*
             * estimate energy of noncoded spectrum (14.4-20kHz)
             *---------------------------------------------------------------------*/

            if( st->extl == SWB_BWE_HIGHRATE )
            {
                en_noncoded = en_band[N_BANDS_TRANS_BWE_HR-1];
            }
            else  /* st->extl == FB_BWE_HIGHRATE */
            {
                en_noncoded /= (gain * en_band[N_BANDS_TRANS_BWE_HR-1]);

                ind1 = 0;
                if( en_noncoded < BWE_HR_TRANS_EN_LIMIT1 )
                {
                    ind1 = 1;
                    en_noncoded = en_band[N_BANDS_TRANS_BWE_HR-1] * BWE_HR_TRANS_EN_LIMIT1;
                }
                else if( en_noncoded < BWE_HR_TRANS_EN_LIMIT2 )
                {
                    ind1 = 2;
                    en_noncoded = en_band[N_BANDS_TRANS_BWE_HR-1] * BWE_HR_TRANS_EN_LIMIT2;
                }
                else if( en_noncoded < BWE_HR_TRANS_EN_LIMIT3 )
                {
                    ind1 = 3;
                    en_noncoded = en_band[N_BANDS_TRANS_BWE_HR-1] * BWE_HR_TRANS_EN_LIMIT3;
                }
                else
                {
                    en_noncoded = en_band[N_BANDS_TRANS_BWE_HR-1];
                }
                push_indice( st, IND_HR_HF_GAIN, ind1, NBITS_HF_GAIN_BWE_HR );
                nBits -= NBITS_HF_GAIN_BWE_HR;
            }

            /*---------------------------------------------------------------------*
             * AVQ coding (quantize normalized spectrum)
             *---------------------------------------------------------------------*/

            Nsv = (NUM_TRANS_END_FREQ_COEF - NUM_TRANS_START_FREQ_COEF) / WIDTH_BAND;
            AVQ_cod( t_audio + k*input_frame/NUM_TIME_SWITCHING_BLOCKS + NUM_TRANS_START_FREQ_COEF, x_norm, nBits, Nsv );
            AVQ_encmux( st, st->extl, x_norm, &nBits, Nsv, nq );

        }
    }
    else /* !is_transient */
    {
        nBits--;        /* is_transient flag */

        /*---------------------------------------------------------------------*
         * processing of normal (non-transient) frames
         *---------------------------------------------------------------------*/

        /* set width of noncoded (blind estimated) spectrum */
        if( st->extl == SWB_BWE_HIGHRATE )
        {
            width_noncoded = L_FRAME32k - NUM_NONTRANS_END_FREQ_COEF;
        }
        else  /* st->extl == FB_BWE_HIGHRATE */
        {
            width_noncoded = 2*END_FREQ_BWE_FULL_FB/50 - NUM_NONTRANS_END_FREQ_COEF;
        }

        /* compute energy of noncoded (14.4-20kHz) spectrum */
        if( st->extl == FB_BWE_HIGHRATE )
        {
            en_noncoded = sum2_f( t_audio + NUM_NONTRANS_END_FREQ_COEF, width_noncoded ) + 0.001f;
            en_noncoded = (float)sqrt( en_noncoded / width_noncoded );
        }

        /* keep only frequencies in interest */
        set_f( t_audio, 0, NUM_NONTRANS_START_FREQ_COEF );
        set_f( t_audio + NUM_NONTRANS_END_FREQ_COEF, 0, input_frame-NUM_NONTRANS_END_FREQ_COEF );

        /*---------------------------------------------------------------------*
         * global gain coding
         *---------------------------------------------------------------------*/

        /* compute and quantize global gain */
        gain = sum2_f( t_audio + NUM_NONTRANS_START_FREQ_COEF, WIDTH_NONTRANS_FREQ_COEF*N_BANDS_BWE_HR ) + 0.001f;
        gain = (float)sqrt(gain) / (WIDTH_NONTRANS_FREQ_COEF*N_BANDS_BWE_HR);
        ind1 = gain_quant( &gain, MIN_GLOB_GAIN_BWE_HR, MAX_GLOB_GAIN_BWE_HR, NBITS_GLOB_GAIN_BWE_HR );

        push_indice( st, IND_HR_GAIN, ind1, NBITS_GLOB_GAIN_BWE_HR );
        nBits -= NBITS_GLOB_GAIN_BWE_HR;

        /* normalization with global gain */
        tmpF = 1/gain;
        for( i=0; i<WIDTH_NONTRANS_FREQ_COEF*N_BANDS_BWE_HR; i++ )
        {
            t_audio[NUM_NONTRANS_START_FREQ_COEF + i] *= tmpF;
        }

        /*---------------------------------------------------------------------*
         * envelope coding
         *---------------------------------------------------------------------*/

        /* compute energy per band */
        for( i=0; i<N_BANDS_BWE_HR; i++ )
        {
            en_band[i] = sum2_f( t_audio + NUM_NONTRANS_START_FREQ_COEF + i*WIDTH_NONTRANS_FREQ_COEF, WIDTH_NONTRANS_FREQ_COEF ) + 0.001f;
            en_band[i] = (float)sqrt( en_band[i]/WIDTH_NONTRANS_FREQ_COEF );
        }

        /* Q energy per band */
        ind1 = en_band_quant( en_band, swb_hr_env_code1, NUM_ENVLOPE_CODE_HR1 );
        ind2 = en_band_quant( en_band + 2, swb_hr_env_code2, NUM_ENVLOPE_CODE_HR2 );

        push_indice( st, IND_HR_ENVELOPE, ind1, NBITS_ENVELOPE_BWE_HR1 );
        push_indice( st, IND_HR_ENVELOPE, ind2, NBITS_ENVELOPE_BWE_HR2 );

        nBits -= NBITS_ENVELOPE_BWE_HR1 + NBITS_ENVELOPE_BWE_HR2;

        /* normalize spectrum per bands */
        for( i=0; i<N_BANDS_BWE_HR; i++ )
        {
            tmpF = 1/en_band[i];
            for( j=0; j<WIDTH_NONTRANS_FREQ_COEF; j++ )
            {
                t_audio[NUM_NONTRANS_START_FREQ_COEF + i*WIDTH_NONTRANS_FREQ_COEF + j] *= tmpF;
            }
        }

        /*---------------------------------------------------------------------*
         * choose sub-bands to be quantized
         *---------------------------------------------------------------------*/

        /* find the subband with the min envelope */
        pos = minimum( en_band, N_BANDS_BWE_HR, &min_env );

        /* decide the spectrum to be quantized */
        if( nBits_total > NBITS_THRESH_BWE_HR )
        {
            i = NUM_NONTRANS_END_FREQ_COEF - NUM_NONTRANS_START_FREQ_COEF;
            mvr2r( t_audio + NUM_NONTRANS_START_FREQ_COEF, t_audio_tmp, NUM_NONTRANS_END_FREQ_COEF - NUM_NONTRANS_START_FREQ_COEF );
        }
        else
        {
            /* reorder the spectrum */
            ind1 = (pos*64 + pos/2 * WIDTH_BAND);
            mvr2r( t_audio + NUM_NONTRANS_START_FREQ_COEF, t_audio_tmp, ind1 );

            ind2 = ((pos+1)*64 + (pos+1)/2 * WIDTH_BAND);
            mvr2r( t_audio + NUM_NONTRANS_START_FREQ_COEF + ind2, t_audio_tmp + ind1, NUM_NONTRANS_END_FREQ_COEF - NUM_NONTRANS_START_FREQ_COEF - ind2 );

            i = ind1 + NUM_NONTRANS_END_FREQ_COEF - NUM_NONTRANS_START_FREQ_COEF - ind2;
        }

        /*---------------------------------------------------------------------*
         * estimate energy of noncoded spectrum (14.4-20kHz)
         *---------------------------------------------------------------------*/

        if( st->extl == SWB_BWE_HIGHRATE )
        {
            en_noncoded = 0.5f * min_env;
        }
        else  /* st->extl == FB_BWE_HIGHRATE */
        {
            en_noncoded /= (gain * min_env);

            ind1 = 0;
            if( en_noncoded < BWE_HR_NONTRANS_EN_LIMIT1 )
            {
                ind1 = 1;
                en_noncoded = 0.5f * min_env * BWE_HR_NONTRANS_EN_LIMIT1;
            }
            else if( en_noncoded > BWE_HR_NONTRANS_EN_LIMIT2 )
            {
                ind1 = 2;
                en_noncoded = min_env * BWE_HR_NONTRANS_EN_LIMIT2;
            }
            else if( en_noncoded > BWE_HR_NONTRANS_EN_LIMIT3 )
            {
                ind1 = 3;
                en_noncoded = min_env * BWE_HR_NONTRANS_EN_LIMIT3;
            }
            else
            {
                en_noncoded = 0.5f * min_env;
            }

            push_indice( st, IND_HR_HF_GAIN, ind1, NBITS_HF_GAIN_BWE_HR );
            nBits -= NBITS_HF_GAIN_BWE_HR;
        }

        /*---------------------------------------------------------------------*
         * AVQ coding (quantize normalized spectrum)
         *---------------------------------------------------------------------*/

        Nsv = i / WIDTH_BAND;
        AVQ_cod( t_audio_tmp, x_norm, nBits, Nsv );
        AVQ_encmux( st, st->extl, x_norm, &nBits, Nsv, nq );

        /*---------------------------------------------------------------------*
         * second stage coding
         *---------------------------------------------------------------------*/

        if( nBits >= 9 + NBITS_GLOB_GAIN_BWE_HR && sum_s( nq, Nsv) > 0 )
        {
            /* select spectrum of the second stage coding */
            k = 0;
            for( i=0; i<Nsv; i++ )
            {
                if( nq[i] == 0 )
                {
                    for( j=0; j<WIDTH_BAND; j++ )
                    {
                        t_audio[k++] = t_audio_tmp[i*WIDTH_BAND + j];
                    }
                }
            }

            for( i=0; i<Nsv; i++ )
            {
                if( nq[i] != 0 )
                {
                    for( j=0; j<WIDTH_BAND; j++ )
                    {
                        t_audio[k++] = t_audio_tmp[i*WIDTH_BAND + j] - x_norm[i*WIDTH_BAND + j];
                    }
                }
            }

            /* calculate the number of subbands according to the rest bits */
            if( nBits > 396 )
            {
                Nsv2 = 33;
            }
            else
            {
                Nsv2 = nBits/12;
            }

            /* second stage global gain estimation and coding */
            gain2 = sum2_f( t_audio, Nsv2*WIDTH_BAND ) + 0.001f;
            gain2 = (float)(16*sqrt( gain2 / (Nsv2*WIDTH_BAND) ));
            ind1 = gain_quant( &gain2, MIN_GLOB_GAIN_BWE_HR, MAX_GLOB_GAIN_BWE_HR, NBITS_GLOB_GAIN_BWE_HR );

            push_indice( st, IND_HR_GAIN, ind1, NBITS_GLOB_GAIN_BWE_HR );
            nBits -= NBITS_GLOB_GAIN_BWE_HR;

            /* normalize with global gain */
            gain2 *= 0.0625f;                 /* 1/16 */
            tmpF = 1/gain2;
            for( i=0; i<Nsv2*WIDTH_BAND; i++ )
            {
                t_audio[i] *= tmpF;
            }

            set_s( nq2, 0, Nsv );

            AVQ_cod( t_audio, x_norm1, nBits, Nsv2 );
            AVQ_encmux( st, st->extl, x_norm1, &nBits, Nsv2, nq2 );
        }

    }

    /* write unused bits */
    while( nBits > 0 )
    {
        i = min( nBits, 16 );
        push_indice( st, IND_UNUSED, 0, i );
        nBits -= i;
    }


    return;
}
