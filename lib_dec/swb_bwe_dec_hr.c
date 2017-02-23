/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"

/*-------------------------------------------------------------------*
 * swb_bwe_dec_hr()
 *
 * HR SWB BWE decoder
 *-------------------------------------------------------------------*/

void swb_bwe_dec_hr(
    Decoder_State *st,              /* i/o: decoder state structure     */
    const float *syn_12k8_16k,    /* i  : ACELP core synthesis @16kHz */
    float *hb_synth,        /* o  : SHB synthesis               */
    const short output_frame,     /* i  : frame length                */
    const short unbits,           /* i  : number of core unused bits  */
    const float pitch_buf[]       /* i  : pitch buffer                */
)
{
    short i, j, k, nBits, nBits_total, nBits_block, Nsv, Nsv2, width_noncoded;
    short is_transient, tmpS, incr, IsTransient, pos;
    int   x_norm[NSV_MAX*(WIDTH_BAND+1)], x_norm1[NSV_MAX*(WIDTH_BAND+1)];
    float t_audio[L_FRAME48k], t_audio_tmp[L_FRAME48k];
    float gain = 0.0f, gain2, en_band[N_BANDS_BWE_HR];
    short ind1, ind2;
    float EnergyLT, Energy;
    short nq[NSV_MAX], nq2[NSV_MAX], nq_tmp[NSV_MAX];
    float tilt_wb, min_env, max_env;
    float tmpF, tmp2, pitch, alpha;
    float en_noncoded;
    float env = 0.0f;

    /*---------------------------------------------------------------------*
     * initializations
     *---------------------------------------------------------------------*/

    set_f( t_audio, 0.0f, output_frame );

    st->bwe_highrate_seed = (short)((short)(pitch_buf[0]*64.0f)*(short)(pitch_buf[3]*64.0f));

    ind2 = 0;       /* only to suppress warnings */
    Nsv2 = 0;       /* only to suppress warnings */
    gain2 = 0 ;     /* only to suppress warnings */

    /* reset memories in case that last frame was a different technology */
    if( st->last_core == HQ_CORE || st->last_extl != st->extl )
    {
        set_f( st->old_wtda_swb, 0, L_FRAME48k );
    }

    /* calculate SWB BWE bit-budget */
    nBits = (short)(st->extl_brate) / 50 + unbits;
    nBits_total = nBits;

    /*---------------------------------------------------------------------*
     * calculate tilt of the core synthesis
     *---------------------------------------------------------------------*/
    calc_tilt_bwe( syn_12k8_16k, &tilt_wb, L_FRAME16k );
    pitch = sum_f( pitch_buf, NB_SUBFR16k ) + EPSILON;


    /*---------------------------------------------------------------------*
     * FEC, or good frame decoding
     *---------------------------------------------------------------------*/

    if( st->bfi )
    {
        is_transient = st->old_is_transient_hr_bwe;

        /* Replication of the last spectrum, with an attenuation */
        if( (st->clas_dec == VOICED_CLAS || st->clas_dec == INACTIVE_CLAS) && st->nbLostCmpt <= 3 )
        {
            alpha = 0.8f;
        }
        else if( is_transient )
        {
            alpha = 0.15f;
        }
        else
        {
            alpha = 0.3f;
        }

        if( is_transient )
        {
            /* set BWE spectrum length */
            if( output_frame == L_FRAME32k )
            {
                tmpS = L_FRAME32k/NUM_TIME_SWITCHING_BLOCKS - NUM_TRANS_START_FREQ_COEF;
            }
            else  /* output_frame == L_FRAME48k */
            {
                tmpS = (2*END_FREQ_BWE_FULL_FB/50)/NUM_TIME_SWITCHING_BLOCKS - NUM_TRANS_START_FREQ_COEF;
            }

            /* reconstruct */
            for( k=0; k<NUM_TIME_SWITCHING_BLOCKS; k++ )
            {
                for( i=0; i<tmpS; i++ )
                {
                    t_audio[NUM_TRANS_START_FREQ_COEF + k*output_frame/NUM_TIME_SWITCHING_BLOCKS + i] = alpha * st->t_audio_prev[i + k*tmpS];
                }
                /* save transform coefficients for the next frame (needed in case of frame erasures) */
                mvr2r( t_audio + NUM_TRANS_START_FREQ_COEF + k*output_frame/NUM_TIME_SWITCHING_BLOCKS, st->t_audio_prev + k*tmpS, tmpS );
            }
        }
        else
        {
            /* set BWE spectrum length */
            if( output_frame == L_FRAME32k )
            {
                tmpS = L_FRAME32k - NUM_NONTRANS_START_FREQ_COEF;
            }
            else  /* output_frame == L_FRAME48k */
            {
                tmpS = 2*END_FREQ_BWE_FULL_FB/50 - NUM_NONTRANS_START_FREQ_COEF;
            }

            /* reconstruct */
            for( i=0; i<tmpS; i++ )
            {
                t_audio[NUM_NONTRANS_START_FREQ_COEF + i] = alpha * st->t_audio_prev[i];
            }
            /* Save transform coefficients for the next frame (needed in case of frame erasures) */
            mvr2r( t_audio + NUM_NONTRANS_START_FREQ_COEF, st->t_audio_prev, tmpS );
        }
        st->mem_EnergyLT *= alpha;
        gain = (float)(2.0f*sqrt(st->mem_EnergyLT/output_frame));
        env = 1.0f;
    }
    else
    {
        /*---------------------------------------------------------------------*
         * get transient frame flag
         *---------------------------------------------------------------------*/

        is_transient = (short) get_next_indice( st, 1 );

        if( is_transient )
        {
            nBits = -1;     /* is_transient flag */
            nBits_block = nBits_total / NUM_TIME_SWITCHING_BLOCKS;
            nBits += nBits_total % NUM_TIME_SWITCHING_BLOCKS;

            /* set width of noncoded (blind estimated) spectrum */
            if( st->extl == SWB_BWE_HIGHRATE || output_frame == L_FRAME32k )
            {
                width_noncoded = L_FRAME32k/NUM_TIME_SWITCHING_BLOCKS - NUM_TRANS_END_FREQ_COEF;
                tmpS = L_FRAME32k/NUM_TIME_SWITCHING_BLOCKS - NUM_TRANS_END_FREQ_COEF_EFF;
            }
            else  /* st->extl == FB_BWE_HIGHRATE */
            {
                width_noncoded = (2*END_FREQ_BWE_FULL_FB/50)/NUM_TIME_SWITCHING_BLOCKS - NUM_TRANS_END_FREQ_COEF;
                tmpS = (2*END_FREQ_BWE_FULL_FB/50)/NUM_TIME_SWITCHING_BLOCKS - NUM_TRANS_END_FREQ_COEF_EFF;
            }

            /*---------------------------------------------------------------------*
             * transient frames: processing in blocks (subframes)
             *---------------------------------------------------------------------*/

            for( k = 0; k < NUM_TIME_SWITCHING_BLOCKS; k++ )
            {
                nBits += nBits_block;

                /*---------------------------------------------------------------------*
                 * global gain and envelope decoding
                 *---------------------------------------------------------------------*/

                /* get global gain */
                ind1 = (short) get_next_indice( st, NBITS_GLOB_GAIN_BWE_HR );
                gain = gain_dequant( ind1, MIN_GLOB_GAIN_BWE_HR, MAX_GLOB_GAIN_BWE_HR, NBITS_GLOB_GAIN_BWE_HR );
                nBits -= NBITS_GLOB_GAIN_BWE_HR;

                /* get energy per band */
                if( k == 0 )
                {
                    ind1 = (short) get_next_indice( st, NBITS_ENVELOPE_BWE_HR_TR );
                    ind2 = ind1;
                    nBits -= NBITS_ENVELOPE_BWE_HR_TR;
                }
                else
                {
                    if( ind2 < 8 )
                    {
                        ind1 = (short) get_next_indice( st, NBITS_ENVELOPE_BWE_HR_TR - 1 );
                    }
                    else
                    {
                        ind1 = (short) get_next_indice( st, NBITS_ENVELOPE_BWE_HR_TR - 1 ) + NUM_ENVLOPE_CODE_HR_TR2;
                    }
                    nBits -= (NBITS_ENVELOPE_BWE_HR_TR - 1);
                }

                en_band[0] = swb_hr_env_code3[2 * ind1];
                en_band[1] = swb_hr_env_code3[2 * ind1 + 1];
                env = 0.5f*(en_band[0] + en_band[1]);

                /*---------------------------------------------------------------------*
                 * estimate energy of noncoded spectrum (14.4-20kHz)
                 *---------------------------------------------------------------------*/

                en_noncoded = en_band[N_BANDS_TRANS_BWE_HR-1];

                if( st->extl == FB_BWE_HIGHRATE )
                {
                    ind1 = (short) get_next_indice( st, NBITS_HF_GAIN_BWE_HR );
                    nBits -= NBITS_HF_GAIN_BWE_HR;

                    if( ind1 == 1 )
                    {
                        en_noncoded *= BWE_HR_TRANS_EN_LIMIT1;
                    }
                    else if( ind1 == 2 )
                    {
                        en_noncoded *= BWE_HR_TRANS_EN_LIMIT2;
                    }
                    else if( ind1 == 3 )
                    {
                        en_noncoded *= BWE_HR_TRANS_EN_LIMIT3;
                    }
                }

                /*---------------------------------------------------------------------*
                 * AVQ decoding (dequantize normalized spectrum)
                 *---------------------------------------------------------------------*/

                Nsv = (NUM_TRANS_END_FREQ_COEF - NUM_TRANS_START_FREQ_COEF) / WIDTH_BAND;
                AVQ_demuxdec( st, x_norm, &nBits, Nsv, nq );

                for( i=0; i < Nsv*WIDTH_BAND; i++ )
                {
                    t_audio[k*output_frame/NUM_TIME_SWITCHING_BLOCKS + NUM_TRANS_START_FREQ_COEF + i] = (float) (x_norm[i]);
                }

                /* apply noise-fill */
                swb_hr_noise_fill( is_transient, NUM_TRANS_START_FREQ_COEF, NUM_TRANS_END_FREQ_COEF, tilt_wb, pitch, nq, Nsv, &st->bwe_highrate_seed,
                                   t_audio + NUM_TRANS_START_FREQ_COEF + k*output_frame/NUM_TIME_SWITCHING_BLOCKS );

                /*---------------------------------------------------------------------*
                 * reconstruction
                 *---------------------------------------------------------------------*/

                /* reconstruct 14-16(20) kHz spectrum */
                for( j = 0; j < tmpS; j++ )
                {
                    t_audio[k*output_frame/NUM_TIME_SWITCHING_BLOCKS + NUM_TRANS_END_FREQ_COEF_EFF + j] = 0.5f*t_audio[k*output_frame/NUM_TIME_SWITCHING_BLOCKS + NUM_TRANS_END_FREQ_COEF_EFF - tmpS + j];
                }

                /* envelope denormalization */
                for( i=0; i<N_BANDS_TRANS_BWE_HR; i++ )
                {
                    for( j=0; j<WIDTH_TRANS_FREQ_COEF; j++ )
                    {
                        t_audio[k*output_frame/NUM_TIME_SWITCHING_BLOCKS + NUM_TRANS_START_FREQ_COEF + i*WIDTH_TRANS_FREQ_COEF + j] *= en_band[i];
                    }
                }

                /* envelope denormalization of 14.4-16(20) kHz spectrum */
                for( j = tmpS-width_noncoded; j < tmpS; j++ )
                {
                    t_audio[k*output_frame/NUM_TIME_SWITCHING_BLOCKS + NUM_TRANS_END_FREQ_COEF_EFF + j] *= en_noncoded;
                }

                /* overlap region */
                if( output_frame == L_FRAME48k )
                {
                    for( i=0; i<NSV_OVERLAP*WIDTH_BAND/NUM_TIME_SWITCHING_BLOCKS; i++ )
                    {
                        t_audio[k*output_frame/NUM_TIME_SWITCHING_BLOCKS + NUM_TRANS_START_FREQ_COEF + i] *= overlap_coefs_48kHz[i*4];
                    }
                }
                else
                {
                    for( i=0; i<NSV_OVERLAP*WIDTH_BAND/NUM_TIME_SWITCHING_BLOCKS; i++ )
                    {
                        t_audio[k*output_frame/NUM_TIME_SWITCHING_BLOCKS + NUM_TRANS_START_FREQ_COEF + i] *= overlap_coefs[i*4];
                    }
                }

                /* apply global gain */
                for( i=0; i<WIDTH_TRANS_FREQ_COEF*N_BANDS_TRANS_BWE_HR + width_noncoded; i++ )
                {
                    t_audio[NUM_TRANS_START_FREQ_COEF + k*output_frame/NUM_TIME_SWITCHING_BLOCKS + i] *= gain;
                }

                /* save transform coefficients for the next frame (needed in case of frame erasures) */
                if( output_frame == L_FRAME32k )
                {
                    mvr2r( t_audio + NUM_TRANS_START_FREQ_COEF + k*output_frame/NUM_TIME_SWITCHING_BLOCKS, st->t_audio_prev + k*(L_FRAME32k/NUM_TIME_SWITCHING_BLOCKS - NUM_TRANS_START_FREQ_COEF), L_FRAME32k/NUM_TIME_SWITCHING_BLOCKS - NUM_TRANS_START_FREQ_COEF );
                }
                else  /* output_frame == L_FRAME48k */
                {
                    mvr2r( t_audio + NUM_TRANS_START_FREQ_COEF + k*output_frame/NUM_TIME_SWITCHING_BLOCKS, st->t_audio_prev + k*((2*END_FREQ_BWE_FULL_FB/50)/NUM_TIME_SWITCHING_BLOCKS - NUM_TRANS_START_FREQ_COEF), (2*END_FREQ_BWE_FULL_FB/50)/NUM_TIME_SWITCHING_BLOCKS - NUM_TRANS_START_FREQ_COEF );
                }

                /* attenuate HFs in case of band-width switching */
                if( st->bws_cnt1 > 0 )
                {
                    if( output_frame == L_FRAME32k )
                    {
                        j = L_FRAME32k/NUM_TIME_SWITCHING_BLOCKS - NUM_TRANS_START_FREQ_COEF;
                    }
                    else  /* output_frame == L_FRAME48k */
                    {
                        j = (2*END_FREQ_BWE_FULL_FB/50)/NUM_TIME_SWITCHING_BLOCKS - NUM_TRANS_START_FREQ_COEF;
                    }

                    for( i=0; i<j; i++ )
                    {
                        t_audio[NUM_TRANS_START_FREQ_COEF + k*output_frame/NUM_TIME_SWITCHING_BLOCKS + i] *= (float)st->bws_cnt1 / (float)N_NS2W_FRAMES;
                    }
                }
            }
        }
        else    /* !is_transient */
        {
            /* subtract one bit for is_transient flag */
            nBits--;

            /*---------------------------------------------------------------------*
             * global gain and envelope decoding
             *---------------------------------------------------------------------*/

            /* get global gain */
            ind1 = (short) get_next_indice( st, NBITS_GLOB_GAIN_BWE_HR );
            gain = gain_dequant( ind1, MIN_GLOB_GAIN_BWE_HR, MAX_GLOB_GAIN_BWE_HR, NBITS_GLOB_GAIN_BWE_HR );

            /* get energy per band */
            ind1 = (short) get_next_indice( st, NBITS_ENVELOPE_BWE_HR1 );
            ind2 = (short) get_next_indice( st, NBITS_ENVELOPE_BWE_HR2 );

            en_band[0] = swb_hr_env_code1[2 * ind1];
            en_band[1] = swb_hr_env_code1[2 * ind1 + 1];
            en_band[2] = swb_hr_env_code2[2 * ind2];
            en_band[3] = swb_hr_env_code2[2 * ind2 + 1];
            env = 0.25f*(en_band[0] + en_band[1] + en_band[2] + en_band[3]);

            /*---------------------------------------------------------------------*
             * choose sub-bands to be dequantized
             *---------------------------------------------------------------------*/

            /* find the subband with the min envelope */
            pos = 0;
            min_env = en_band[0];
            max_env = en_band[0];
            for (j = 1; j < N_BANDS_BWE_HR; j++)
            {
                if(en_band[j] < min_env)
                {
                    pos = j;
                    min_env = en_band[j];
                }
                if(en_band[j] > max_env)
                {
                    max_env = en_band[j];
                }
            }

            /* decide the spectrum to be dequantized */
            if( nBits_total > NBITS_THRESH_BWE_HR )
            {
                i = NUM_NONTRANS_END_FREQ_COEF - NUM_NONTRANS_START_FREQ_COEF;
            }
            else
            {
                i = NUM_NONTRANS_END_FREQ_COEF - NUM_NONTRANS_START_FREQ_COEF - 64 - 8*(pos%2);
            }

            nBits -= NBITS_GLOB_GAIN_BWE_HR + NBITS_ENVELOPE_BWE_HR1 + NBITS_ENVELOPE_BWE_HR2;

            /*---------------------------------------------------------------------*
             * estimate energy of noncoded spectrum (14.4-20kHz)
             *---------------------------------------------------------------------*/

            en_noncoded = 0.5f * min_env;

            if( st->extl == FB_BWE_HIGHRATE )
            {
                ind1 = (short) get_next_indice( st, NBITS_HF_GAIN_BWE_HR );
                nBits -= NBITS_HF_GAIN_BWE_HR;

                if( ind1 == 1 )
                {
                    en_noncoded *= BWE_HR_NONTRANS_EN_LIMIT1;
                }
                else if( ind1 == 2 )
                {
                    en_noncoded *= 2.0f * BWE_HR_NONTRANS_EN_LIMIT2;
                }
                else if( ind1 == 3 )
                {
                    en_noncoded *= 2.0f * BWE_HR_NONTRANS_EN_LIMIT3;
                }
            }

            /*---------------------------------------------------------------------*
             * AVQ decoding (dequantize normalized spectrum)
             *---------------------------------------------------------------------*/

            Nsv = i / WIDTH_BAND;
            AVQ_demuxdec( st, x_norm, &nBits, Nsv, nq );

            /*---------------------------------------------------------------------*
             * second stage decoding
             *---------------------------------------------------------------------*/

            if( nBits >= 9 + NBITS_GLOB_GAIN_BWE_HR && sum_s( nq, Nsv) > 0 )
            {
                ind1 = (short) get_next_indice( st, NBITS_GLOB_GAIN_BWE_HR );
                gain2 = gain_dequant( ind1, MIN_GLOB_GAIN_BWE_HR, MAX_GLOB_GAIN_BWE_HR, NBITS_GLOB_GAIN_BWE_HR );
                gain2 *= 0.0625f;

                /* calculate the number of subbands according to the rest bits */
                if( nBits > 396 )
                {
                    Nsv2 = 33;
                }
                else
                {
                    Nsv2 = nBits/12;
                }

                nBits -= NBITS_GLOB_GAIN_BWE_HR;
                AVQ_demuxdec( st, x_norm1, &nBits, Nsv2, nq2 );
            }

            /*---------------------------------------------------------------------*
             * dequantization
             *---------------------------------------------------------------------*/

            for( i=0; i<Nsv*WIDTH_BAND; i++ )
            {
                t_audio_tmp[i] = (float) (x_norm[i]);
            }

            mvs2s( nq, nq_tmp, Nsv );
            if( Nsv2 > Nsv )
            {
                /* Safety check, happens rarely */
                set_s( nq_tmp + Nsv, 0, Nsv2 - Nsv );
            }

            k = 0;
            incr = 0;
            for( i=0; i<Nsv; i++ )
            {
                if( nq[i] == 0 && incr < Nsv2 )
                {
                    for( j=0; j<WIDTH_BAND; j++ )
                    {
                        t_audio_tmp[i*WIDTH_BAND + j] = gain2 * x_norm1[k++];
                    }
                    nq[i] += nq2[incr++];
                }
            }

            for( i=0; incr<Nsv2; i++ )
            {
                /* safety check, happens rarely */
                if( i >= Nsv2 )
                {
                    break;
                }

                if( nq_tmp[i] != 0 )
                {
                    for( j=0; j<WIDTH_BAND; j++ )
                    {
                        t_audio_tmp[i*WIDTH_BAND + j] += gain2 * x_norm1[k++];
                    }
                    nq[i] += nq2[incr++];
                }
            }

            /*---------------------------------------------------------------------*
             * reorder the decoded spectrum
             *---------------------------------------------------------------------*/

            if( nBits_total > NBITS_THRESH_BWE_HR )
            {
                mvr2r( t_audio_tmp, t_audio + NUM_NONTRANS_START_FREQ_COEF, NUM_NONTRANS_END_FREQ_COEF - NUM_NONTRANS_START_FREQ_COEF );
            }
            else
            {
                ind1 = (pos*64 + pos/2 * WIDTH_BAND);
                mvr2r( t_audio_tmp, t_audio + NUM_NONTRANS_START_FREQ_COEF, ind1 );

                ind2 = ((pos+1)*64 + (pos+1)/2 * WIDTH_BAND);
                mvr2r( t_audio_tmp + ind1, t_audio + NUM_NONTRANS_START_FREQ_COEF + ind2, NUM_NONTRANS_END_FREQ_COEF - NUM_NONTRANS_START_FREQ_COEF - ind2 );

                /* reconstruct non-encoded subband */
                if( pos == 3 )
                {
                    mvr2r( t_audio + NUM_NONTRANS_START_FREQ_COEF + 128, t_audio + NUM_NONTRANS_START_FREQ_COEF + 200, 72 );

                    mvs2s( nq + 16, nq + 25, 9 );
                }
                else
                {
                    pos %= 2;
                    mvr2r( t_audio + NUM_NONTRANS_START_FREQ_COEF + ind2, t_audio + NUM_NONTRANS_START_FREQ_COEF + ind1, 64 + pos*WIDTH_BAND );
                    ind1 /= WIDTH_BAND;
                    ind2 /= WIDTH_BAND;

                    j = 0;
                    for( i=Nsv-1; i>=ind1; i-- )
                    {
                        nq[33 - j++] = nq[i];
                    }

                    mvs2s( nq + ind2, nq + ind1, WIDTH_BAND+pos );
                }
            }

            /* apply noise-fill */
            if( nBits < 200 )
            {
                swb_hr_noise_fill( is_transient, NUM_NONTRANS_START_FREQ_COEF, NUM_NONTRANS_END_FREQ_COEF, tilt_wb,
                                   pitch, nq, Nsv, &st->bwe_highrate_seed, t_audio + NUM_NONTRANS_START_FREQ_COEF );
            }

            /*---------------------------------------------------------------------*
             * reconstruction
             *---------------------------------------------------------------------*/

            /* smoothing 12.6-12.8kHz */
            if( pos == 3 && nBits_total <= 400 )
            {
                tmpF = sum2_f( t_audio + NUM_NONTRANS_START_FREQ_COEF + 200 - WIDTH_BAND, WIDTH_BAND ) + EPSILON;
                tmp2 = sum2_f( t_audio + NUM_NONTRANS_START_FREQ_COEF + 200, WIDTH_BAND ) + EPSILON;
                tmpF = (float)sqrt( tmpF/tmp2 );
                for( i=0; i<WIDTH_BAND; i++ )
                {
                    t_audio[NUM_NONTRANS_START_FREQ_COEF + 200 + i] *= ((1.0f - i/(float)WIDTH_BAND)*tmpF + i/(float)WIDTH_BAND);
                }
            }

            /* reconstruct 14.4-16(20) kHz spectrum */
            if( st->extl == SWB_BWE_HIGHRATE || output_frame == L_FRAME32k  )
            {
                width_noncoded = L_FRAME32k - NUM_NONTRANS_END_FREQ_COEF;
            }
            else  /* st->extl == FB_BWE_HIGHRATE */
            {
                width_noncoded = 2*END_FREQ_BWE_FULL_FB/50 - NUM_NONTRANS_END_FREQ_COEF;
            }
            mvr2r( t_audio + NUM_NONTRANS_END_FREQ_COEF - width_noncoded, t_audio + NUM_NONTRANS_END_FREQ_COEF, width_noncoded );

            /* smoothing 14.4-14.8kHz */
            tmpF = sum2_f( t_audio + NUM_NONTRANS_END_FREQ_COEF - WIDTH_BAND, WIDTH_BAND ) + EPSILON;
            tmp2 = sum2_f( t_audio + NUM_NONTRANS_END_FREQ_COEF, WIDTH_BAND ) + EPSILON;
            tmpF = (float)sqrt( tmpF/tmp2 );
            for( i=0; i<WIDTH_BAND; i++ )
            {
                t_audio[NUM_NONTRANS_END_FREQ_COEF + i] *= tmpF;
            }

            /* envelope denormalization */
            for( i=0; i<N_BANDS_BWE_HR; i++ )
            {
                for( j=0; j<WIDTH_NONTRANS_FREQ_COEF; j++ )
                {
                    t_audio[NUM_NONTRANS_START_FREQ_COEF + i*WIDTH_NONTRANS_FREQ_COEF + j] *= en_band[i];
                }
            }

            /* equalize 14.4-16(20) kHz spectrum */
            tmpF = max_env / min_env;
            if( st->extl == SWB_BWE_HIGHRATE || tmpF < 2.2f )
            {
                for( j=0; j<WIDTH_BAND; j++ )
                {
                    t_audio[NUM_NONTRANS_END_FREQ_COEF + j] *= ((1.0f - j/(float)WIDTH_BAND)*en_band[3] + (j/(float)WIDTH_BAND)*en_noncoded);
                }

                for( j=WIDTH_BAND; j<width_noncoded; j++ )
                {
                    t_audio[NUM_NONTRANS_END_FREQ_COEF + j] *= en_noncoded;
                }
            }
            else
            {
                if( output_frame == L_FRAME48k )
                {
                    tmpS = width_noncoded - 2*WIDTH_NONTRANS_FREQ_COEF;
                }
                else
                {
                    tmpS = L_FRAME32k - NUM_NONTRANS_END_FREQ_COEF;
                }

                k = 0;
                for( j=0; j<tmpS; j++ )
                {
                    t_audio[NUM_NONTRANS_END_FREQ_COEF + j] *= 2.2f * en_noncoded * (1-(float)k/(float)160);
                    k++;
                }

                k = 0;
                for( ; j<width_noncoded; j++ )
                {
                    t_audio[NUM_NONTRANS_END_FREQ_COEF + j] *= 0.65f * en_noncoded * (1-(float)k/(float)320);
                    k++;
                }
            }

            /* overlap region */
            if( output_frame == L_FRAME48k )
            {
                for( i=0; i<NSV_OVERLAP*WIDTH_BAND; i++ )
                {
                    t_audio[NUM_NONTRANS_START_FREQ_COEF + i] *= overlap_coefs_48kHz[i];
                }
            }
            else
            {
                for( i=0; i<NSV_OVERLAP*WIDTH_BAND; i++ )
                {
                    t_audio[NUM_NONTRANS_START_FREQ_COEF + i] *= overlap_coefs[i];
                }
            }

            /* apply global gain */
            if( nBits_total <= NBITS_THRESH_BWE_HR )
            {
                gain *= 0.85f;
            }

            for( i=0; i<WIDTH_NONTRANS_FREQ_COEF*N_BANDS_BWE_HR + width_noncoded; i++ )
            {
                t_audio[NUM_NONTRANS_START_FREQ_COEF + i] *= gain;
            }

            /* save transform coefficients for the next frame (needed in case of frame erasures) */
            if( output_frame == L_FRAME32k )
            {
                mvr2r( t_audio + NUM_NONTRANS_START_FREQ_COEF, st->t_audio_prev, L_FRAME32k - NUM_NONTRANS_START_FREQ_COEF );
            }
            else  /* output_frame == L_FRAME48k */
            {
                mvr2r( t_audio + NUM_NONTRANS_START_FREQ_COEF, st->t_audio_prev, 2*END_FREQ_BWE_FULL_FB/50 - NUM_NONTRANS_START_FREQ_COEF );
            }

            /* attenuate HFs in case of band-width switching */
            if( st->bws_cnt1 > 0 )
            {
                if( output_frame == L_FRAME32k )
                {
                    j = L_FRAME32k - NUM_NONTRANS_START_FREQ_COEF;
                }
                else  /* output_frame == L_FRAME48k */
                {
                    j = 2*END_FREQ_BWE_FULL_FB/50 - NUM_NONTRANS_START_FREQ_COEF;
                }

                for( i=0; i<j; i++ )
                {
                    t_audio[NUM_NONTRANS_START_FREQ_COEF+i] *= (float)st->bws_cnt1 / (float)N_NS2W_FRAMES;
                }
            }
        }
    }

    st->prev_ener_shb = gain*env;
    for(i=0; i<SWB_FENV; i++)
    {
        st->prev_SWB_fenv[i] = gain*env;
    }

    /*---------------------------------------------------------------------*
     * iOLA and iMDCT
     *---------------------------------------------------------------------*/

    inverse_transform( t_audio, t_audio_tmp, is_transient, output_frame, output_frame );
    window_ola( t_audio_tmp, hb_synth, st->old_wtda_swb, output_frame, ALDO_WINDOW, ALDO_WINDOW, 0,0,0 );

    /*---------------------------------------------------------------------*
     * final adjustments
     *---------------------------------------------------------------------*/

    if( !st->bfi )
    {
        IsTransient = 0;
        EnergyLT = st->mem_EnergyLT;
        pos = 0;
        for( j=0; j<4; j++ )
        {
            Energy = sum2_f( hb_synth + j*(output_frame/4), output_frame/4 ) + EPSILON;
            if( Energy > 12.5f * EnergyLT )
            {
                IsTransient = 1;
                pos = j;
            }

            EnergyLT = 0.75f*EnergyLT + 0.25f*Energy;
        }

        if( IsTransient == 1 && pos > 0 && tilt_wb < 3.0f && pitch > 500 )
        {
            Nsv = pos*(output_frame/4);
            Energy = sum2_f( hb_synth, Nsv ) + EPSILON;
            if( st->last_extl != st->extl )
            {
                st->mem_EnergyLT = Energy;
            }
            gain = (float)sqrt( pos*st->mem_EnergyLT/Energy );

            gain *= 0.2f;
            for( i=0; i<Nsv; i++ )
            {
                hb_synth[i] *= gain;
            }

            alpha = (float)WIDTH_BAND/output_frame;
            for ( i = 0; i < output_frame/WIDTH_BAND; i++ )
            {
                hb_synth[i+Nsv] *= ((1.0f - i*alpha)*gain + i*alpha);
            }
        }

        st->mem_EnergyLT = EnergyLT;
        st->old_is_transient_hr_bwe = is_transient;
    }

    /* post-processing in case of TD/FD switching */
    if( st->last_core == HQ_CORE || st->last_extl != st->extl )
    {
        if( tilt_wb < 3.0f )
        {
            gain = td_postprocess( hb_synth, output_frame, st->last_extl );

            for( i=0; i<output_frame; i++ )
            {
                st->old_wtda_swb[i] *= gain;
            }

            tmpS = L_FRAME32k - NUM_NONTRANS_START_FREQ_COEF;
            if( output_frame == L_FRAME48k )
            {
                tmpS = 2*END_FREQ_BWE_FULL_FB/50 - NUM_NONTRANS_START_FREQ_COEF;
            }

            for( i=0; i<tmpS; i++ )
            {
                st->t_audio_prev[i] *= gain;
            }
        }
    }

    return;
}
