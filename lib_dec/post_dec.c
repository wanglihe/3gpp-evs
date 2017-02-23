/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"


/*---------------------------------------------------------------------*
 * Function prototypes
 *---------------------------------------------------------------------*/

static void bass_pf_1sf_delay( float *syn, const int *T_sf, const float *gainT_sf, const short l_frame,
                               const short l_subfr, float *bpf_noise_buf, int *gain_factor_param,
                               const short disable_bpf, float *mem_deemph_err, float *lp_ener );


/*---------------------------------------------------------------------*
 * post_decoder()
 *
 * Perform post-processing
 *---------------------------------------------------------------------*/

void post_decoder(
    Decoder_State *st,            /* i/o: decoder memory state pointer    */
    const short coder_type,     /* i  : coder type                      */
    float synth_buf[],
    const float pit_gain[],
    const int   pitch[],
    float signal_out[],
    float *bpf_noise_buf
)
{
    short L_frame, nb_subfr, i;
    float *synth2;
    float *synth;
    short pfstat_on_previous;
    int pitch_gain_adjust[NB_SUBFR16k];
    float synth_buf2[NBPSF_PIT_MAX+L_FRAME_MAX+M];
    long bitrate;
    float A[M+1];
    float pitch_buf[NB_SUBFR16k];
    float tmp;
    short L_subfr;

    L_frame = st->L_frame;
    nb_subfr = st->nb_subfr;
    bitrate = (st->core_brate > SID_2k40 )?st->total_brate:st->last_active_brate;
    pfstat_on_previous = st->pfstat.on;
    st->pfstat.on = 0;
    set_i( pitch_gain_adjust, st->bpf_gain_param, nb_subfr );
    synth = synth_buf + st->old_synth_len;
    synth2 = synth_buf2 + NBPSF_PIT_MAX;
    mvr2r( st->pst_old_syn, synth_buf2, NBPSF_PIT_MAX );

    if( st->tcxonly )
    {
        /* High bitrates (48kbps and above), high sampling rates (25.6kHz and above) */

        mvr2r( synth, synth2, L_frame );

        if( pfstat_on_previous )
        {
            /* Past frame was low-bitrate with formant post-filter */
            lsp2a_stab( st->lsp_old, A, M );
            mvr2r( st->pfstat.mem_pf_in+L_SYN_MEM-M, synth-M, M );
            L_subfr = st->L_frame/st->nb_subfr;
            residu( A, M, synth, synth_buf, L_subfr );
            syn_filt ( A, M, synth_buf, synth2, L_subfr, st->pfstat.mem_stp+L_SYN_MEM-M, 0 );
            scale_st ( synth, synth2, &st->pfstat.gain_prec, L_subfr, -1 );
            blend_subfr2(synth2+L_subfr/2, synth+L_subfr/2, synth2+L_subfr/2);
        }
    }
    else
    {
        /* Low bitrates (32kbps and below), low sampling rates (12.8kHz and 16kHz) */
        if( st->last_bwidth == NB )
        {
            /* NB Post-filter (pitch+formant post-filter) */
            mvr2r( synth, synth_buf, L_frame );
            tmp = synth[-1];
            preemph( synth_buf, st->preemph_fac, L_frame, &tmp );

            tmp = 0.0f;
            for( i=0; i< nb_subfr; i ++ )
            {
                pitch_buf[i] = pitch[i];
            }

            if( pfstat_on_previous == 0 )
            {
                st->pfstat.reset = 1;
            }

            if( st->bwidth == NB )
            {
                st->pfstat.on = 1;
                nb_post_filt( L_frame, L_SUBFR, &(st->pfstat), &tmp, 0, synth_buf, st->mem_Aq, pitch_buf, GENERIC, st->BER_detect,
                              st->lp_noise>LP_NOISE_THRESH?1:((st->core != ACELP_CORE)||(coder_type==UNVOICED)) );
            }
            else
            {
                st->pfstat.on = 0;
                nb_post_filt( L_frame, L_SUBFR, &(st->pfstat), &tmp, 0, synth_buf, st->mem_Aq, pitch_buf, AUDIO, st->BER_detect,
                              st->lp_noise>LP_NOISE_THRESH?1:((st->core != ACELP_CORE)||(coder_type==UNVOICED)) );
            }


            mvr2r( synth_buf, synth2, L_frame );
            tmp = synth2[-1];
            deemph( synth2, st->preemph_fac, L_frame, &tmp );
        }
        else
        {
            /* Formant Post-filter */
            if( pfstat_on_previous == 0 )
            {
                st->pfstat.reset = 1;
            }

            if( st->bwidth >= WB )
            {
                st->pfstat.on = 1;
                formant_post_filt( &(st->pfstat), synth, st->mem_Aq, synth2, L_frame, L_SUBFR, st->lp_noise, bitrate, 0 );
            }
            else
            {
                st->pfstat.on = 0;
                formant_post_filt( &(st->pfstat), synth, st->mem_Aq, synth2, L_frame, L_SUBFR, st->lp_noise, bitrate, 1 );
            }
        }

        /*Bass Post-filter */
        bass_pf_1sf_delay( synth2, pitch, pit_gain, L_frame, L_SUBFR, bpf_noise_buf, pitch_gain_adjust,
                           (st->lp_noise>LP_NOISE_THRESH && st->narrowBand)?1:0, &(st->pst_mem_deemp_err), &(st->pst_lp_ener) );
    }

    /* Output */
    mvr2r( synth2, signal_out, L_frame );

    /* Update synth2 memory */
    mvr2r( synth_buf2 + L_frame, st->pst_old_syn, NBPSF_PIT_MAX );

    return;
}


/*---------------------------------------------------------------------*
 * bass_pf_1sf_delay()
 *
 * Perform low-frequency postfiltering
 *---------------------------------------------------------------------*/

static void bass_pf_1sf_delay(
    float *syn,                 /* i  : synthesis to postfilter                     */
    const int   *T_sf,                /* i  : Pitch period for all subframes (T_sf[4])    */
    const float *gainT_sf,            /* i  : Pitch gain for all subframes (gainT_sf[4])  */
    const short L_frame,              /* i  : frame length (multiple of l_subfr)          */
    const short L_subfr_in,           /* i  : sub-frame length (80/64)                    */
    float *bpf_noise_buf,       /* i  : harmoninc filtered signal                   */
    int   *gain_factor_param,   /* i  : gain factor param 0-> minimum BPF, 3-> full BPF  */
    const short disable_bpf,          /* i  : flag to disable BPF                         */
    float *mem_deemph_err,      /* i/o: Error deemphasis memory                     */
    float *lp_ener              /* i/o: long_term error signal energy               */
)
{
    short i, sf, i_subfr, T, lg, L_subfr;
    float tmp, corr, ener, gain;
    float noise_buf[(2*L_SUBFR)], *noise_in;
    float error[L_SUBFR];
    float ener2;

    noise_in = noise_buf;

    sf = 0;
    L_subfr = L_subfr_in;

    for( i_subfr = 0; i_subfr < L_frame; i_subfr += L_subfr, sf++ )
    {
        if( i_subfr == 0 )
        {
            L_subfr = L_subfr_in;
        }
        else if( i_subfr == L_frame )
        {
            L_subfr  = 0;
        }
        else
        {
            L_subfr = L_subfr_in;
        }

        T = T_sf[sf];
        gain = gainT_sf[sf];

        if (gain > 1.0f) gain = 1.0f;
        if (gain < 0.0f) gain = 0.0f;

        lg = L_frame - T - i_subfr;
        if (lg < 0) lg = 0;
        if (lg > L_subfr) lg = L_subfr;

        if( !disable_bpf && gain > 0 )
        {
            corr = 0.01f;
            ener = 0.01f;

            for( i=0; i<lg; i++ )
            {
                corr += syn[i+i_subfr] * (0.5f*syn[i+i_subfr-T] + 0.5f*syn[i+i_subfr+T]);
                ener += (0.5f*syn[i+i_subfr-T] + 0.5f*syn[i+i_subfr+T])*(0.5f*syn[i+i_subfr-T] + 0.5f*syn[i+i_subfr+T]);
            }

            for( i=lg; i<L_subfr; i++ )
            {
                corr += syn[i+i_subfr]*syn[i+i_subfr-T];
                ener += syn[i+i_subfr-T]*syn[i+i_subfr-T];
            }
            gain = corr/ener;

            if( gain > 1.f )
            {
                gain = 1.0f;
            }
            else if( gain<0.f )
            {
                gain = 0.f;
            }

            ener2 = 0.01f;
            for( i=0; i<lg; i++ )
            {
                error[i] = gain * (syn[i+i_subfr] - 0.5f*syn[i+i_subfr-T] - 0.5f*syn[i+i_subfr+T]);
                error[i] = error[i] + 0.9f* *mem_deemph_err;
                *mem_deemph_err = error[i];
                ener2 += error[i]*error[i];
            }

            for( i=lg; i<L_subfr; i++ )
            {
                error[i] = 0.5f*gain * (syn[i+i_subfr] - syn[i+i_subfr-T]);
                error[i] = error[i] + 0.9f* *mem_deemph_err;
                *mem_deemph_err = error[i];
                ener2 += error[i]*error[i];
            }

            ener2 = (float) (10.f*log10(ener2));
            *lp_ener = (float)(0.99f* *lp_ener + 0.01f*ener2);
            ener2 = (float)pow(10.f,0.1f* *lp_ener);
            tmp = 0.5f*corr/(ener+ener2);

            if( tmp > 0.5f )
            {
                tmp = 0.5f;
            }
            else if( tmp < 0.f )
            {
                tmp = 0.0f;
            }

            /*Adjust gain*/
            /* full gain = gainLTP*0.5*/
            /* adaptive gain = gainLTP*0.5*gain_factor*0.5*/
            tmp *= max(0.5f*gain_factor_param[sf],0.125f);

            /* calculate noise based on voiced pitch */
            for( i=0; i<lg; i++ )
            {
                noise_in[i] = tmp * (syn[i+i_subfr] - 0.5f*syn[i+i_subfr-T] - 0.5f*syn[i+i_subfr+T]);
            }

            for( i=lg; i<L_subfr; i++ )
            {
                noise_in[i] = tmp * (syn[i+i_subfr] - syn[i+i_subfr-T]);
                /*It simulates an extrapolation of the buffer syn: syn[i+i_subfr+T]=syn[i+i_subfr]
                 * -> reduce nrg of noise_in and avoid too much post-filtering*/
                /*noise_in[i] = tmp * (syn[i+i_subfr] - 0.5f*syn[i+i_subfr-T] - 0.5f*syn[i+i_subfr]);*/
                /*->noise_in[i] = tmp * 0.5f * (syn[i+i_subfr] - syn[i+i_subfr-T]);*/
                noise_in[i] *= 0.5f;
            }
        }
        else
        {
            set_zero( noise_in, L_subfr );
        }

        /* copy bpf noise signal to buffer */
        mvr2r( noise_in, bpf_noise_buf + i_subfr, L_subfr );

    }

    return;
}
/*---------------------------------------------------------------------*
 * cldfb_synth_set_bandsToZero()
 *
 *
 *---------------------------------------------------------------------*/

void cldfb_synth_set_bandsToZero(
    Decoder_State *st,
    float **rAnalysis,
    float **iAnalysis,
    const short nTimeSlots
)
{
    float nrg_bwddec, nrg_band[CLDFB_NO_CHANNELS_MAX], thr_bwddwc, max_nrg, realQ1, imagQ1;
    short flag, offset, WBcnt, i, k, update_perc;
    float perc_detect, perc_miss;

    realQ1 = 0.0f;
    imagQ1 = 0.0f;

    set_f( nrg_band, 0.0f, CLDFB_NO_CHANNELS_MAX );
    max_nrg = 0.0f;

    offset = 250;
    WBcnt = 20;
    perc_miss = 0.83f;
    perc_detect = 0.93f;

    if(st->VAD==1)
    {
        st->active_frame_cnt_bwddec++;
        st->total_frame_cnt_bwddec++;
        if(st->active_frame_cnt_bwddec > 99)
        {
            st->active_frame_cnt_bwddec = 100;
        }
        if(st->total_frame_cnt_bwddec > 500)
        {
            st->total_frame_cnt_bwddec = 500;
        }

        for (i = 0; i < (st->cldfbSyn->no_channels - st->cldfbSyn->bandsToZero); i++)
        {
            nrg_bwddec = 0.0f;
            for (k = 0; k < nTimeSlots; k++)
            {
                realQ1 = rAnalysis[k][i];
                imagQ1 = iAnalysis[k][i];
                nrg_bwddec += (realQ1*realQ1);
                nrg_bwddec += (imagQ1*imagQ1);
            }
            nrg_band[i] = (nrg_bwddec);
            if( (nrg_band[i] > max_nrg) && (i > 11) )
            {
                max_nrg = nrg_band[i];
            }
        }
        for(; i < st->cldfbSyn->no_channels; i++)
        {
            nrg_band[i] = 0;
        }

        nrg_bwddec = 0;
        for(i = 2; i < 9; i++)
        {
            nrg_bwddec += (nrg_band[i]/7.0f);
        }

        thr_bwddwc = (nrg_bwddec/512.0f);

        st->avg_nrg_LT = 0.98999f*st->avg_nrg_LT + 0.009979f*thr_bwddwc;
        update_perc = 1;
        if(st->ini_frame >= 25 && thr_bwddwc < st->avg_nrg_LT*0.005f)
        {
            update_perc = 0;
        }

        flag = 1;
        if(max_nrg >= thr_bwddwc)
        {
            flag = 0;
        }

        for(i = 0; i < WBcnt-1; i++)
        {
            st->flag_buffer[i] = st->flag_buffer[i+1];
        }
        st->flag_buffer[WBcnt-1] = flag;

        /*long term percentage*/
        if(update_perc == 1)
        {
            st->perc_bwddec += (flag - st->perc_bwddec)/st->active_frame_cnt_bwddec;
        }
        if((st->total_frame_cnt_bwddec > offset) && (st->active_frame_cnt_bwddec > 50) )
        {
            if( (st->perc_bwddec >= perc_detect || (st->perc_bwddec >= perc_miss && st->last_flag_filter_NB)) && (sum_s(st->flag_buffer, WBcnt) != 0)) /*decision hysterysis*/
            {
                st->cldfbSyn->bandsToZero = ( st->cldfbSyn->no_channels - 10 );
                st->last_flag_filter_NB = 1; /*VAD processing must be dependent on hysterysis, as if hysterysis fails, but threshold passes, we dont want next vad frames to have NB only*/
            }
            else
            {
                st->last_flag_filter_NB = 0;
            }
        }
        else
        {
            st->last_flag_filter_NB = 0;
        }
        if(sum_s(st->flag_buffer, WBcnt) == 0)
        {
            st->perc_bwddec = 0.0f;
            st->active_frame_cnt_bwddec = 0;
            st->total_frame_cnt_bwddec = 0;
            st->last_flag_filter_NB = 0;
        }
    }
    else
    {
        if(st->last_flag_filter_NB == 1)
        {
            st->cldfbSyn->bandsToZero = st->last_active_bandsToZero_bwdec;
        }
        st->total_frame_cnt_bwddec++;
        if(st->total_frame_cnt_bwddec > 500)
        {
            st->total_frame_cnt_bwddec = 500;
        }
    }

    st->last_active_bandsToZero_bwdec = st->cldfbSyn->bandsToZero;

    return;
}
