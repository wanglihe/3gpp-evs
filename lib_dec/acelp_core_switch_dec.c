/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*---------------------------------------------------------------------*
 * Local functions
 *---------------------------------------------------------------------*/

static void decod_gen_voic_core_switch( Decoder_State *st, const short L_frame, const short sharpFlag, const float *Aq,
                                        const short coder_type, float *exc, const long  core_brate );

/*-------------------------------------------------------------------*
 * acelp_core_switch_dec()
 *
 * ACELP core decoder in the first ACELP->HQ switching frame
 *-------------------------------------------------------------------*/

void acelp_core_switch_dec(
    Decoder_State *st,                /* i/o: decoder state structure         */
    float *synth_subfr_out,   /* o  : synthesized ACELP subframe      */
    float *tmp_synth_bwe,     /* o  : synthesized ACELP subframe BWE  */
    const short output_frame,       /* i  : input frame legth               */
    const short core_switching_flag,/* i  : core switching flag             */
    float *mem_synth          /* o  : synthesis to overlap            */
)
{
    short i, delta, L_frame_for_cs, decode_bwe;
    short d1m, ind1, fdelay, gapsize;
    long cbrate;
    float synth_intFreq[2*L_SUBFR];
    float old_exc[L_EXC_DEC], *exc;
    float tmp_mem2[2*L_FILT48k], gain;
    float hb_synth_tmp[NS2SA(48000, 10000000L)];
    const float *hp_filter;
    float Aq[2*(M+1)];
    float bpf_error_signal[2*L_SUBFR];
    float *realBuffer[CLDFB_NO_COL_MAX_SWITCH], *imagBuffer[CLDFB_NO_COL_MAX_SWITCH];
    float realBufferTmp[CLDFB_NO_COL_MAX_SWITCH][CLDFB_NO_CHANNELS_MAX], imagBufferTmp[CLDFB_NO_COL_MAX_SWITCH][CLDFB_NO_CHANNELS_MAX];

    /*----------------------------------------------------------------*
     * Initializations
     *----------------------------------------------------------------*/

    /* open CLDFB buffer up to CLDFB_NO_CHANNELS_MAX bands for 48kHz */
    for( i=0; i<CLDFB_NO_COL_MAX_SWITCH; i++ )
    {
        set_f( realBufferTmp[i], 0, CLDFB_NO_CHANNELS_MAX );
        set_f( imagBufferTmp[i], 0, CLDFB_NO_CHANNELS_MAX );
        realBuffer[i] = realBufferTmp[i];
        imagBuffer[i] = imagBufferTmp[i];
    }

    d1m = 0;
    gain = 0;

    mvr2r( st->old_Aq_12_8, Aq, M+1 );
    mvr2r( st->old_Aq_12_8, Aq + (M+1), M+1 );

    set_f(mem_synth, 0, NS2SA(16000, DELAY_CLDFB_NS)+2 );

    /* set multiplication factor according to the sampling rate */
    delta = 1;
    if( output_frame == L_FRAME32k )
    {
        delta = 2;
    }
    else if( output_frame == L_FRAME48k )
    {
        delta = 3;
    }

    /*----------------------------------------------------------------*
     * set switching frame bit-rate
     *----------------------------------------------------------------*/

    if( core_switching_flag && st->last_L_frame == st->last_L_frame_ori && (st->last_core == ACELP_CORE || st->last_core == AMR_WB_CORE) )
    {
        exc = old_exc + L_EXC_MEM_DEC;
        mvr2r(st->old_exc, old_exc, L_EXC_MEM_DEC );

        if( st->last_L_frame == L_FRAME )
        {
            if( st->core_brate > ACELP_24k40 )
            {
                cbrate = ACELP_24k40;
            }
            else
            {
                cbrate = st->core_brate;
            }

            L_frame_for_cs = L_FRAME;
        }
        else
        {
            if( st->core_brate <= ACELP_8k00 )
            {
                cbrate = ACELP_8k00;
            }
            else if( st->core_brate <= ACELP_14k80 )
            {
                cbrate = ACELP_14k80;
            }
            else
            {
                cbrate = min(st->core_brate, ACELP_22k60 );
            }

            L_frame_for_cs = L_FRAME16k;
        }

        /*----------------------------------------------------------------*
         * Excitation decoding
         *----------------------------------------------------------------*/

        decod_gen_voic_core_switch( st, L_frame_for_cs, 0, Aq, GENERIC, exc, cbrate );

        /*----------------------------------------------------------------*
         * synthesis, deemphasis, postprocessing and resampling
         *----------------------------------------------------------------*/

        /* synthesis and deemphasis */
        syn_12k8( 2*L_SUBFR, Aq, exc, synth_intFreq, st->mem_syn2, 1 );

        if( st->pfstat.on && st->last_bwidth == NB )
        {
            float tmp_noise, pitch_buf_tmp[2];
            tmp_noise = 0;

            for( i=0; i<2; i++ )
            {
                pitch_buf_tmp[i] = L_SUBFR;
            }

            nb_post_filt( 2*L_SUBFR, L_SUBFR, &(st->pfstat), &tmp_noise, 0, synth_intFreq, Aq, pitch_buf_tmp, AUDIO, st->BER_detect, 0 );
        }

        if( L_frame_for_cs == L_FRAME )
        {
            deemph( synth_intFreq, PREEMPH_FAC, 2*L_SUBFR, &st->mem_deemph );
        }
        else
        {
            deemph( synth_intFreq, PREEMPH_FAC_16k, 2*L_SUBFR, &st->mem_deemph );
        }

        AGC_dec( synth_intFreq, st->agc_mem2, 2*L_SUBFR );

        if( st->pfstat.on && st->last_bwidth != NB )
        {
            mvr2r( st->pfstat.mem_pf_in+L_SYN_MEM-M, bpf_error_signal, M ); /*bpf_error_signal used as temporary buffer*/
            mvr2r( synth_intFreq, bpf_error_signal+M, L_SUBFR );
            residu( Aq, M, bpf_error_signal+M,old_exc, L_SUBFR );
            syn_filt( Aq, M, old_exc, bpf_error_signal, L_SUBFR, st->pfstat.mem_stp+L_SYN_MEM-M, 0 );
            scale_st( synth_intFreq, bpf_error_signal, &st->pfstat.gain_prec, L_SUBFR, -1 );
            mvr2r( bpf_error_signal, synth_intFreq, L_SUBFR/2 );
            blend_subfr2( bpf_error_signal + L_SUBFR/2, synth_intFreq + L_SUBFR/2, synth_intFreq + L_SUBFR/2 );
        }
        st->pfstat.on = 0;

        if( st->flag_cna )
        {
            generate_masking_noise( synth_intFreq, st->hFdCngDec->hFdCngCom, 2*L_SUBFR, 0 );
        }

        /*----------------------------------------------------------------*
         * Resample to the output sampling rate (8/16/32/48 kHz)
         * Bass post-filter
         *----------------------------------------------------------------*/

        /* bass post-filter */
        bass_psfilter( st->Opt_AMR_WB, synth_intFreq, 2*L_SUBFR, NULL, st->pst_old_syn,
                       &st->pst_mem_deemp_err, &st->pst_lp_ener, st->bpf_off, st->stab_fac, &st->stab_fac_smooth,
                       st->mem_mean_pit, st->Track_on_hist, st->vibrato_hist, &st->psf_att, GENERIC, bpf_error_signal );

        /* CLDFB analysis of the synthesis at internal sampling rate */
        cldfb_save_memory( st->cldfbAna );
        cldfbAnalysis( synth_intFreq, realBuffer, imagBuffer, NS2SA(50*L_frame_for_cs, SWITCH_GAP_LENGTH_NS+DELAY_CLDFB_NS), st->cldfbAna );
        cldfb_restore_memory( st->cldfbAna );

        /* CLDFB analysis and add the BPF error signal */
        cldfb_save_memory( st->cldfbBPF );
        addBassPostFilter( bpf_error_signal, st->bpf_off?0:NS2SA(50*L_frame_for_cs, SWITCH_GAP_LENGTH_NS+DELAY_CLDFB_NS), realBuffer, imagBuffer, st->cldfbBPF );
        cldfb_restore_memory( st->cldfbBPF );

        /* CLDFB synthesis of the combined signal */
        cldfb_save_memory( st->cldfbSyn );
        cldfbSynthesis( realBuffer, imagBuffer, synth_subfr_out, NS2SA(st->output_Fs, SWITCH_GAP_LENGTH_NS+DELAY_CLDFB_NS), st->cldfbSyn );
        cldfb_restore_memory( st->cldfbSyn );

        mvr2r( synth_intFreq + NS2SA(L_frame_for_cs * 50, SWITCH_GAP_LENGTH_NS-DELAY_CLDFB_NS)-2, mem_synth, NS2SA(L_frame_for_cs * 50, DELAY_CLDFB_NS)+2 );     /* need for switching (-2 is due to 0 delay filtering) */

        /*----------------------------------------------------------------*
         * BWE decoding
         *----------------------------------------------------------------*/

        decode_bwe = 0;
        if( !((inner_frame_tbl[st->bwidth] == L_FRAME16k && st->last_L_frame == L_FRAME16k) || inner_frame_tbl[st->bwidth] == L_FRAME8k) )
        {
            /* Decoding of BWE */
            d1m = (short)get_next_indice(st, AUDIODELAYBITS);
            ind1 = (short)get_next_indice(st, NOOFGAINBITS1);
            gain = usdequant( ind1, MINVALUEOFFIRSTGAIN, DELTAOFFIRSTGAIN );
            decode_bwe = 1;
        }



        if( decode_bwe && !((output_frame == L_FRAME16k && st->last_L_frame == L_FRAME16k) || output_frame == L_FRAME8k) )
        {
            set_f( tmp_mem2, 0, 2*L_FILT48k );

            hp_filter = hp16000_48000;
            fdelay = 48;
            if( st->output_Fs == 16000 )
            {
                if( st->last_L_frame == L_FRAME )
                {
                    hp_filter = hp12800_16000;
                    fdelay = 20;
                }
            }
            else if( st->output_Fs == 32000 )
            {

                if( st->last_L_frame == L_FRAME )
                {
                    hp_filter = hp12800_32000;
                    fdelay = 40;
                }
                else
                {
                    hp_filter = hp16000_32000;
                    fdelay = 32;
                }
            }
            else if( st->last_L_frame == L_FRAME )
            {
                hp_filter = hp12800_48000;
                fdelay = 60;
            }

            /* safety check in case of bit errors */
            i = MAX_D1M_16k;
            if( st->last_L_frame == L_FRAME )
            {
                i = MAX_D1M_12k8;
            }

            if( d1m >= i )
            {
                d1m = i - 1;
                gain = 0; /* force muting */
                st->BER_detect = 1;
            }

            i = NS2SA(st->output_Fs, FRAME_SIZE_NS-ACELP_LOOK_NS-DELAY_BWE_TOTAL_NS);
            mvr2r( st->old_synth_sw, hb_synth_tmp, i );
            set_f( hb_synth_tmp + i, 0, NS2SA(st->output_Fs, 10000000L) - i );
            fir( hb_synth_tmp, hp_filter, hb_synth_tmp, tmp_mem2, output_frame>>1, fdelay, 0 ); /* put the 40 past samples into the memory */
            set_f( tmp_synth_bwe, 0, SWITCH_MAX_GAP );

            gapsize = delta * (NS2SA(16000,SWITCH_GAP_LENGTH_NS));
            for( i=0; i<gapsize; i++ )
            {
                tmp_synth_bwe[i] = hb_synth_tmp[i+d1m*delta+fdelay] * gain;
            }
        }
        else
        {
            set_f( tmp_synth_bwe, 0, SWITCH_MAX_GAP );
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * acelp_core_switch_dec_bfi()
 *
 * ACELP core decoder in the first ACELP->HQ switching frame in case of BAD frame
 *-------------------------------------------------------------------*/

void acelp_core_switch_dec_bfi(
    Decoder_State *st,                  /* i/o: decoder state structure  */
    float synth_out[],          /* o  : synthesis                */
    const short coder_type            /* i  : coder type               */
)
{
    short i;
    float old_exc[L_EXC_DEC], *exc;                     /* excitation signal buffer              */
    float syn[L_FRAME16k];                              /* synthesis signal buffer               */
    float lsf_new[M];                                   /* LSFs at the end of the frame          */
    float lsp_new[M];                                   /* LSPs at the end of the frame          */
    float Aq[NB_SUBFR16k*(M+1)];                        /* A(q)   quantized for the 4 subframes  */
    float old_exc2[L_FRAME16k + L_EXC_MEM], *exc2;      /* total excitation buffer               */
    float tmp_noise;                                    /* Long term temporary noise energy      */
    float FEC_pitch;                                    /* FEC pitch                             */
    float old_bwe_exc[((PIT16k_MAX + (L_FRAME16k+1) + L_SUBFR16k) * 2)]; /* excitation buffer    */
    float *bwe_exc;                                     /* Excitation for SWB TBE                */
    float tmp_float[NBPSF_PIT_MAX];
    float tmp_float2[M];
    float tmp_float3;
    float tmp_float4[L_TRACK_HIST];
    short tmp_float5[L_TRACK_HIST];
    short tmp_float6[L_TRACK_HIST];
    float tmp_float7;
    float voice_factors[NB_SUBFR16k];
    float pitch_buf[NB_SUBFR16k];
    float *realBuffer[CLDFB_NO_COL_MAX_SWITCH_BFI], *imagBuffer[CLDFB_NO_COL_MAX_SWITCH_BFI];
    float realBufferTmp[CLDFB_NO_COL_MAX_SWITCH_BFI][CLDFB_NO_CHANNELS_MAX], imagBufferTmp[CLDFB_NO_COL_MAX_SWITCH_BFI][CLDFB_NO_CHANNELS_MAX];

    /*----------------------------------------------------------------*
     * Initialization
     *----------------------------------------------------------------*/

    /* initialize CLDFB buffer up to CLDFB_NO_CHANNELS_MAX bands for 48kHz */
    for( i=0; i<CLDFB_NO_COL_MAX_SWITCH_BFI; i++ )
    {
        set_f( realBufferTmp[i], 0, CLDFB_NO_CHANNELS_MAX );
        set_f( imagBufferTmp[i], 0, CLDFB_NO_CHANNELS_MAX );
        realBuffer[i] = realBufferTmp[i];
        imagBuffer[i] = imagBufferTmp[i];
    }

    st->bpf_off = 1;
    st->clas_dec = st->last_good;
    tmp_noise = 0.0f;

    mvr2r( st->old_exc, old_exc, L_EXC_MEM_DEC );
    exc = old_exc + L_EXC_MEM_DEC;
    mvr2r( st->old_exc2, old_exc2, L_EXC_MEM );
    exc2 = old_exc2 + L_EXC_MEM;
    mvr2r( st->old_bwe_exc, old_bwe_exc, PIT16k_MAX * 2);
    bwe_exc = old_bwe_exc + PIT16k_MAX * 2;
    st->GSC_noisy_speech = 0;
    st->relax_prev_lsf_interp = 0;

    /* SC-VBR */
    if( st->last_nelp_mode_dec == 1 )
    {
        st->nelp_mode_dec = 1;
    }

    mvr2r( st->mem_AR, tmp_float, M );
    mvr2r( st->mem_MA, tmp_float2, M );

    /* LSF estimation and A(z) calculation */
    lsf_dec_bfi( MODE1, lsf_new, st->lsf_old, st->lsf_adaptive_mean, NULL, st->mem_MA, st->mem_AR,
                 st->stab_fac, st->last_coder_type, st->L_frame,  st->last_good,
                 st->nbLostCmpt, 0, NULL, NULL, NULL, st->Last_GSC_pit_band_idx, st->Opt_AMR_WB, st->bwidth );

    FEC_lsf2lsp_interp( st, st->L_frame, Aq, lsf_new, lsp_new );

    mvr2r( tmp_float, st->mem_AR, M );
    mvr2r( tmp_float2, st->mem_MA, M );

    /*----------------------------------------------------------------*
     * Excitation decoding
     *----------------------------------------------------------------*/

    if( st->nelp_mode_dec == 1 )
    {
        /* SC-VBR */
        decod_nelp( st, coder_type, &tmp_noise, pitch_buf, exc, exc2, voice_factors, bwe_exc, st->bfi, tmp_float2 );
        FEC_pitch = pitch_buf[3];
    }
    else
    {
        tmp_float[0] = st->bfi_pitch;
        tmp_float[1] = (float)st->bfi_pitch_frame;
        tmp_float[2] = st->lp_gainp;
        tmp_float[3] = st->lp_gainc;
        tmp_float[4] = st->upd_cnt;
        tmp_float[5] = st->seed;

        /* calculation of excitation signal */
        FEC_exc_estim( st, st->L_frame, exc, exc2, syn /* dummy buffer */, pitch_buf, voice_factors, &FEC_pitch, bwe_exc, lsf_new, &tmp_noise );

        st->seed = (short)tmp_float[5];
        st->bfi_pitch = tmp_float[0];
        st->bfi_pitch_frame = (short)tmp_float[1];
        st->lp_gainp = tmp_float[2];
        st->lp_gainc = tmp_float[3];
        st->upd_cnt = (short)tmp_float[4];
    }

    /*------------------------------------------------------------------*
     * Synthesis
     *-----------------------------------------------------------------*/

    mvr2r( st->mem_syn2, tmp_float, M );
    syn_12k8( st->L_frame, Aq, exc2, syn, tmp_float, 1 );

    tmp_float[0] = st->enr_old;
    fer_energy( st->L_frame, st->last_good, syn, FEC_pitch, tmp_float, st->L_frame );

    /*------------------------------------------------------------------*
     * Perform fixed deemphasis through 1/(1 - g*z^-1)
     *-----------------------------------------------------------------*/

    mvr2r( &(st->mem_deemph), tmp_float, 1 );
    if( st->L_frame == L_FRAME )
    {
        deemph( syn, PREEMPH_FAC, L_FRAME, tmp_float );
    }
    else
    {
        deemph( syn, PREEMPH_FAC_16k, L_FRAME16k, tmp_float );
    }

    /*----------------------------------------------------------------*
     * Bass post-filter
     *----------------------------------------------------------------*/

    st->bpf_off = 1;
    mvr2r( st->pst_old_syn, tmp_float, NBPSF_PIT_MAX );
    tmp_float3 = st->stab_fac_smooth;
    mvr2r( st->mem_mean_pit, tmp_float4, L_TRACK_HIST );
    mvs2s( st->Track_on_hist, tmp_float5, L_TRACK_HIST );
    mvs2s( st->vibrato_hist, tmp_float6, L_TRACK_HIST );
    tmp_float7 = st->psf_att;

    /* apply bass post-filter */
    bass_psfilter( st->Opt_AMR_WB, syn, st->L_frame, pitch_buf, tmp_float, &st->pst_mem_deemp_err, &st->pst_lp_ener, st->bpf_off,
                   st->stab_fac, &tmp_float3, tmp_float4, tmp_float5, tmp_float6, &tmp_float7, coder_type, old_exc /* tmp buffer */ );

    /*----------------------------------------------------------------*
     * Resampling to the output sampling frequency
     *----------------------------------------------------------------*/

    /* CLDFB analysis of the synthesis at internal sampling rate */
    cldfb_save_memory( st->cldfbAna );
    cldfbAnalysis ( syn, realBuffer, imagBuffer, st->L_frame/2, st->cldfbAna );
    cldfb_restore_memory( st->cldfbAna );

    /* CLDFB synthesis of the combined signal */
    cldfb_save_memory( st->cldfbSyn );
    cldfbSynthesis ( realBuffer, imagBuffer, synth_out, st->output_Fs*0.01, st->cldfbSyn );
    cldfb_restore_memory( st->cldfbSyn );

    return;
}


/*-------------------------------------------------------------------*
 * decod_gen_voic_core_switch()
 *
 * Decode excitation signal in teh first ACELP->HQ switching frame
 *-------------------------------------------------------------------*/

static void decod_gen_voic_core_switch(
    Decoder_State *st,            /* i/o: decoder static memory       */
    const short L_frame,          /* i  : length of the frame         */
    const short sharpFlag,        /* i  : flag for formant sharpening */
    const float *Aq,              /* i  : LP filter coefficient       */
    const short coder_type,       /* i  : coding type                 */
    float *exc,             /* i/o: adapt. excitation exc       */
    const long  core_brate        /* i  : switching frame bit-rate    */
)
{
    short T0, T0_frac, T0_min, T0_max;/* integer pitch variables                          */
    float gain_pit;               /* pitch gain                                           */
    float gain_code;              /* gain/normalized gain of the algebraic excitation     */
    float norm_gain_code;         /* normalized gain of the algebraic excitation          */
    float gain_inov;              /* Innovation gain                                      */
    float voice_fac;              /* voicing factor                                       */
    float code[L_SUBFR];          /* algebraic codevector                                 */
    float pitch;                  /* pointer to floating pitch                            */
    short i;                      /* tmp variables                                        */
    short pitch_limit_flag;
    float tmpF;

    /* initializations */
    if( L_frame == L_FRAME )
    {
        T0_max = PIT_MAX;
        T0_min = PIT_MIN;
    }
    else /* L_frame == L_FRAME16k */
    {
        T0_max = PIT16k_MAX;
        T0_min = PIT16k_MIN;
    }

    /*----------------------------------------------------------------------*
    * Decode pitch lag
    *----------------------------------------------------------------------*/

    pitch = pit_decode( st, core_brate, 0, L_frame, 0, coder_type, &pitch_limit_flag, &T0, &T0_frac, &T0_min, &T0_max, L_SUBFR );

    /*--------------------------------------------------------------*
    * Find the adaptive codebook vector.
    *--------------------------------------------------------------*/

    pred_lt4( &exc[0], &exc[0], T0, T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP );

    /*--------------------------------------------------------------*
    * LP filtering of the adaptive excitation
    *--------------------------------------------------------------*/

    lp_filt_exc_dec( st, MODE1, core_brate,0, coder_type, 0, L_SUBFR, L_frame, 0, exc );

    /*--------------------------------------------------------------*
    * Innovation decoding
    *--------------------------------------------------------------*/

    inov_decode( st, core_brate, 0, L_frame, coder_type, sharpFlag, 0, -1, Aq, st->tilt_code, pitch, code );

    /*--------------------------------------------------------------*
    * Gain decoding
    * Estimate spectrum tilt and voicing
    *--------------------------------------------------------------*/

    if( L_frame == L_FRAME )
    {
        gain_dec_mless( st, core_brate, L_frame, TRANSITION, 0, -1, code, st->old_Es_pred, &gain_pit, &gain_code, &gain_inov, &norm_gain_code );
    }
    else
    {
        gain_dec_mless( st, core_brate, L_frame, coder_type, 0, -1, code, st->old_Es_pred, &gain_pit, &gain_code, &gain_inov, &norm_gain_code );
    }

    st->tilt_code = est_tilt( exc, gain_pit, code, gain_code, &voice_fac, L_SUBFR, 0 );

    /*----------------------------------------------------------------------*
    * Find the total excitation
    *----------------------------------------------------------------------*/

    if( st->prev_bfi )
    {
        gain_code = min( gain_code, 0.5f*gain_code+0.5f*st->lp_gainc );
    }

    for( i = 0; i < L_SUBFR; i++ )
    {
        tmpF = gain_pit * exc[i];
        exc[i] = tmpF + gain_code * code[i];
    }

    /*-----------------------------------------------------------------*
    * long term prediction on the 2nd sub frame
    *-----------------------------------------------------------------*/

    pred_lt4( &exc[L_SUBFR], &exc[L_SUBFR], T0, T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP );

    for( i = 0; i < L_SUBFR; i++ )
    {
        exc[i+L_SUBFR] *= gain_pit;
    }

    return;
}
