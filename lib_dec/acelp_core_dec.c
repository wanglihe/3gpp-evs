/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * acelp_core_dec()
 *
 * ACELP core decoder
 *-------------------------------------------------------------------*/

void acelp_core_dec(
    Decoder_State *st,                  /* i/o: decoder state structure         */
    float synth[],                      /* o  : synthesis                       */
    float bwe_exc_extended[],           /* i/o: bandwidth extended excitation   */
    float *voice_factors,               /* o  : voicing factors                 */
    float old_syn_12k8_16k[],           /* o  : intermediate ACELP synthesis at 12.8kHz or 16kHz to be used by SWB BWE */
    short coder_type,                   /* i  : coder type                      */
    short sharpFlag,                    /* i  : formant sharpening flag         */
    float pitch_buf[NB_SUBFR16k],       /* o  : floating pitch for each subframe*/
    short *unbits,                      /* o  : number of unused bits           */
    short *sid_bw                       /* o  : 0-NB/WB, 1-SWB SID              */
)
{
    float old_exc[L_EXC_DEC], *exc;                     /* excitation signal buffer              */
    float syn_tmp[L_FRAME16k+L_SUBFR], *syn;            /* synthesis signal buffer               */
    short output_frame;                                 /* frame length at output sampling freq. */
    float lsf_new[M];                                   /* LSFs at the end of the frame          */
    float lsp_new[M];                                   /* LSPs at the end of the frame          */
    float lsp_mid[M];                                   /* LSPs in the middle of the frame       */
    float Aq[NB_SUBFR16k*(M+1)];                        /* A(q)   quantized for the 4 subframes  */
    float old_exc2[L_FRAME16k + L_EXC_MEM], *exc2;      /* total excitation buffer               */
    float mem_tmp[M];                                   /* temporary synthesis filter memory     */
    float enr_q;                                        /* E information for FER protection      */
    float tmp_noise;                                    /* Long term temporary noise energy      */
    float Es_pred;                                      /* predicted scaled innov. energy        */
    float FEC_pitch;                                    /* FEC pitch                             */
    float old_bwe_exc[((PIT16k_MAX + (L_FRAME16k+1) + L_SUBFR16k) * 2)]; /* excitation buffer    */
    float *bwe_exc;                                     /* Excitation for SWB TBE                */
    short i, int_fs;
    short tc_subfr;
    short allow_cn_step;
    float temp_buf[L_FRAME16k + L_SYN_MEM];
    short last_pulse_pos;
    short T0_tmp;
    short do_WI;
    float dct_buffer[DCT_L_POST];
    float exc_buffer[DCT_L_POST];
    float dct_exc_tmp[L_FRAME16k];
    float bpf_error_signal[L_FRAME16k];
    short nb_bits;                                     /* number of bits                        */
    int   indice;                                      /* parameter indices to write            */
    float gain_buf[NB_SUBFR16k];
    float q_env[20];
    float exc3[L_FRAME16k];
    float syn1_tmp[L_FRAME16k+2], *syn1;
    float *realBuffer[CLDFB_NO_COL_MAX], *imagBuffer[CLDFB_NO_COL_MAX];
    float realBufferTmp[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX];
    float imagBufferTmp[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX];
    short LSF_Q_prediction;  /* o  : LSF prediction mode                 */
    float tmpF;

    /* open CLDFB buffer up to CLDFB_NO_CHANNELS_MAX bands for 48kHz */
    for( i=0; i<CLDFB_NO_COL_MAX; i++ )
    {
        set_f( realBufferTmp[i], 0, CLDFB_NO_CHANNELS_MAX );
        set_f( imagBufferTmp[i], 0, CLDFB_NO_CHANNELS_MAX );
        realBuffer[i] = realBufferTmp[i];
        imagBuffer[i] = imagBufferTmp[i];
    }


    /*----------------------------------------------------------------*
     * Initialization
     *----------------------------------------------------------------*/

    LSF_Q_prediction = -1;
    set_f( syn_tmp, 0, L_SUBFR );
    syn = syn_tmp + L_SUBFR;
    syn1_tmp[0] = 0;
    syn1_tmp[1] = 0;
    syn1 = syn1_tmp+2;

    output_frame = (short)(st->output_Fs / 50);

    st->bpf_off = 0;
    if( st->last_core == HQ_CORE )
    {
        /* in case of HQ->ACELP switching, do not apply BPF */
        st->bpf_off = 1;

        /* in case of core switching, rest post-filter memories */
        st->pfstat.on = 0;

        /* reset the GSC pre echo energy threshold in case of switching */
        st->Last_frame_ener = (float)MAX_32;
    }

    if( st->prev_bfi > 0 )
    {
        /* reset the GSC pre echo energy threshold in case of FEC */
        st->Last_frame_ener = (float)MAX_32;
    }

    st->clas_dec = st->last_good;
    enr_q = 0.0f;
    Es_pred = 0.0f;
    tmp_noise = 0.0f;

    mvr2r(st->old_exc, old_exc, L_EXC_MEM_DEC );
    exc = old_exc + L_EXC_MEM_DEC;
    mvr2r( st->old_exc2, old_exc2, L_EXC_MEM );
    exc2 = old_exc2 + L_EXC_MEM;
    mvr2r( st->old_bwe_exc, old_bwe_exc, PIT16k_MAX * 2);
    bwe_exc = old_bwe_exc + PIT16k_MAX * 2;

    last_pulse_pos = 0;
    do_WI = 0;
    st->GSC_noisy_speech = 0;
    st->relax_prev_lsf_interp = 0;
    set_zero( gain_buf, NB_SUBFR16k );

    if( st->L_frame == L_FRAME )
    {
        st->gamma = GAMMA1;
        st->preemph_fac = PREEMPH_FAC;
        int_fs = INT_FS_12k8;
    }
    else
    {
        st->gamma = GAMMA16k;
        st->preemph_fac = PREEMPH_FAC_16k;
        int_fs = INT_FS_16k;
    }

    /* reset post-filter in case of switching */
    if( st->pfstat.on == 0 )
    {
        st->pfstat.reset = 1;
    }

    /*----------------------------------------------------------------*
     * Updates in case of internal sampling rate switching
     *----------------------------------------------------------------*/

    if( st->last_L_frame != st->L_frame && st->last_core != HQ_CORE )
    {
        if( st->pfstat.on != 0 )
        {
            short mem_syn_r_size_old, mem_syn_r_size_new;

            mem_syn_r_size_old = (short)(1.25*st->last_L_frame/20.f);
            mem_syn_r_size_new = (short)(1.25*st->L_frame/20.f);
            lerp( st->pfstat.mem_stp+L_SYN_MEM-mem_syn_r_size_old, st->pfstat.mem_stp+L_SYN_MEM-mem_syn_r_size_new, mem_syn_r_size_new, mem_syn_r_size_old );
            lerp( st->pfstat.mem_pf_in+L_SYN_MEM-mem_syn_r_size_old, st->pfstat.mem_pf_in+L_SYN_MEM-mem_syn_r_size_new, mem_syn_r_size_new, mem_syn_r_size_old );
        }

        /* convert quantized LSP vector */
        st->rate_switching_reset = lsp_convert_poly( st->lsp_old, st->L_frame, 0 );

        /* convert old quantized LSF vector */
        lsp2lsf( st->lsp_old, st->lsf_old, M, int_fs );

        /* FEC - update adaptive LSF mean vector */
        mvr2r( st->lsf_old, st->lsfoldbfi1, M );
        mvr2r( st->lsf_old, st->lsfoldbfi0, M );
        mvr2r( st->lsf_old, st->lsf_adaptive_mean, M );

        /* Reset LPC mem */
        if( st->sr_core == 16000 )
        {
            mvr2r( GEWB2_Ave, st->mem_AR, M );
        }
        else
        {
            mvr2r( GEWB_Ave, st->mem_AR, M );
        }
        set_zero( st->mem_MA, M );

        /* update synthesis filter memories */
        synth_mem_updt2( st->L_frame, st->last_L_frame, st->old_exc, st->mem_syn_r, st->mem_syn2, NULL, DEC );
        mvr2r( st->old_exc, old_exc, L_EXC_MEM_DEC );
        mvr2r( st->mem_syn2, st->mem_syn1, M );
        mvr2r( st->mem_syn2, st->mem_syn3, M );

    }

    /* update buffer of old subframe pitch values */
    if( st->last_L_frame != st->L_frame )
    {
        if( st->L_frame == L_FRAME )
        {
            if( st->last_L_frame == L_FRAME32k )
            {
                tmpF = (float)12800/(float)32000;
            }
            else if( st->last_L_frame == 512 )
            {
                tmpF = (float)12800/(float)25600;
            }
            else /* st->last_L_frame == L_FRAME16k */
            {
                tmpF = (float)12800/(float)16000;
            }

            for( i=NB_SUBFR16k-NB_SUBFR; i<NB_SUBFR16k; i++ )
            {
                st->old_pitch_buf[i-1] = tmpF * st->old_pitch_buf[i];
            }

            for( i=2*NB_SUBFR16k-NB_SUBFR; i<2*NB_SUBFR16k; i++ )
            {
                st->old_pitch_buf[i-2] = tmpF * st->old_pitch_buf[i];
            }
        }
        else
        {

            if( st->last_L_frame == L_FRAME32k )
            {
                tmpF = (float)16000/(float)32000;
            }
            else if( st->last_L_frame == 512 )
            {
                tmpF = (float)16000/(float)25600;
            }
            else /* st->last_L_frame == L_FRAME12k8 */
            {
                tmpF = (float)16000/(float)12800;
            }

            for( i=2*NB_SUBFR-1; i>=NB_SUBFR; i-- )
            {
                st->old_pitch_buf[i+2] = tmpF * st->old_pitch_buf[i];
            }
            st->old_pitch_buf[NB_SUBFR+1] = st->old_pitch_buf[NB_SUBFR+2];

            for( i=NB_SUBFR-1; i>=0; i-- )
            {
                st->old_pitch_buf[i+1] = tmpF * st->old_pitch_buf[i];
            }
            st->old_pitch_buf[0] = st->old_pitch_buf[1];
        }
    }

    if( st->bfi_pitch_frame != st->L_frame )
    {
        if( st->L_frame == L_FRAME )
        {
            if( st->bfi_pitch_frame == L_FRAME32k )
            {
                tmpF = (float)12800/(float)32000;
            }
            else if( st->bfi_pitch_frame == 512 )
            {
                tmpF = (float)12800/(float)25600;
            }
            else /* st->bfi_pitch_frame == L_FRAME16k */
            {
                tmpF = (float)12800/(float)16000;
            }
            st->bfi_pitch *= tmpF;
            st->bfi_pitch_frame = L_FRAME;
        }
        else
        {
            if( st->bfi_pitch_frame == L_FRAME32k )
            {
                tmpF = (float)16000/(float)32000;
            }
            else if( st->bfi_pitch_frame == 512 )
            {
                tmpF = (float)16000/(float)25600;
            }
            else /* st->bfi_pitch_frame == L_FRAME12k8 */
            {
                tmpF = (float)16000/(float)12800;
            }
            st->bfi_pitch *= tmpF;
            st->bfi_pitch_frame = L_FRAME16k;
        }
    }

    if( st->last_bwidth == NB && st->bwidth != NB && st->ini_frame != 0 )
    {
        st->rate_switching_reset = 1;
    }

    /*----------------------------------------------------------------------*
     * GOOD frame
     *----------------------------------------------------------------------*/

    if( !st->bfi )
    {

        /*----------------------------------------------------------------*
         * Decoding of TC subframe clasification
         *----------------------------------------------------------------*/

        tc_subfr = -1;
        if( coder_type == TRANSITION )
        {
            tc_subfr = tc_classif( st, st->L_frame );
        }

        /*----------------------------------------------------------------*
         * Decoding of inactive CNG frames
         *----------------------------------------------------------------*/

        if ( st->core_brate == FRAME_NO_DATA || st->core_brate == SID_2k40 )
        {
            /* decode CNG parameters */
            if( st->cng_type == LP_CNG )
            {
                CNG_dec( st, st->L_frame, Aq, st->core_brate, lsp_new, lsf_new, &allow_cn_step, sid_bw, q_env );

                /* comfort noise generation */
                CNG_exc( st->core_brate, st->L_frame, &st->Enew, &st->cng_seed, exc, exc2, &st->lp_ener, st->last_core_brate,
                         &st->first_CNG, &(st->cng_ener_seed), bwe_exc, allow_cn_step, &st->last_allow_cn_step, st->num_ho,
                         q_env, st->lp_env, st->old_env, st->exc_mem, st->exc_mem1, sid_bw, &st->cng_ener_seed1, exc3 ,st->Opt_AMR_WB );
            }
            else
            {
                if( st->core_brate == SID_2k40 )
                {
                    FdCng_decodeSID( st );
                    *sid_bw = 0;
                }

                generate_comfort_noise_dec( NULL, NULL, st );

                FdCng_exc( st->hFdCngDec->hFdCngCom, &st->CNG_mode, st->L_frame, st->lsp_old, st->first_CNG, st->lspCNG, Aq, lsp_new,lsf_new, exc, exc2, bwe_exc );

                mvr2r( exc2, exc3, st->L_frame );
            }

            /* update past excitation signals for LD music post-filter */
            mvr2r( st->dct_post_old_exc + L_FRAME, st->dct_post_old_exc, DCT_L_POST-L_FRAME-OFFSET2 );
            mvr2r( exc2, st->dct_post_old_exc + (DCT_L_POST-L_FRAME-OFFSET2), L_FRAME );

            /* synthesis at 12.8kHz sampling rate */
            syn_12k8( st->L_frame, Aq, exc2, syn, st->mem_syn2, 1 );
            syn_12k8( st->L_frame, Aq, exc3, syn1, st->mem_syn3, 1 );

            /* reset the decoder */
            CNG_reset_dec( st, pitch_buf, voice_factors );

            /* update st->mem_syn1 for ACELP core switching */
            mvr2r( st->mem_syn3, st->mem_syn1, M );

            /* update old synthesis for classification */
            mvr2r( syn1 + st->L_frame - L_SYN_MEM_CLAS_ESTIM, st->mem_syn_clas_estim, L_SYN_MEM_CLAS_ESTIM );

            /* Update music post processing values */
            /* Filter energies update */
            for( i = 0; i < DCT_L_POST; i++ )
            {
                st->filt_lfE[i] = 0.3f + 0.7f*st->filt_lfE[i];
            }

            /* save and delay synthesis to be used by SWB BWE */
            save_old_syn( st->L_frame, syn1, old_syn_12k8_16k, st->old_syn_12k8_16k, st->preemph_fac, &st->mem_deemph_old_syn );
        }
        /*----------------------------------------------------------------*
         * Decoding of all other frames
         *----------------------------------------------------------------*/

        else
        {
            /*-----------------------------------------------------------------*
             * After CNG period, use the most up-to-date LSPs
             *-----------------------------------------------------------------*/

            if( st->last_core_brate == FRAME_NO_DATA || st->last_core_brate == SID_2k40 )
            {
                mvr2r( st->lspCNG, st->lsp_old, M );

                lsp2lsf( st->lspCNG, st->lsf_old, M, int_fs );
            }

            /*-----------------------------------------------------------------*
             * Reset higher ACELP pre-quantizer in case of switching
             *-----------------------------------------------------------------*/

            if( !st->use_acelp_preq )
            {
                st->mem_preemp_preQ = 0.0f;
                st->last_nq_preQ = 0;
            }
            st->use_acelp_preq = 0;

            /*-----------------------------------------------------------------*
             * LSF de-quantization and interpolation
             *-----------------------------------------------------------------*/

            lsf_dec( st, tc_subfr, st->L_frame, coder_type, st->bwidth, Aq, &LSF_Q_prediction, lsf_new, lsp_new, lsp_mid );

            /*-----------------------------------------------------------------*
             * FEC - first good frame after lost frame(s) (possibility to correct the ACB)
             *-----------------------------------------------------------------*/

            if( st->core_brate >= ACELP_11k60 )
            {
                last_pulse_pos = 0;

                /* decode the last glottal pulse position */
                T0_tmp = FEC_pos_dec( st, coder_type, st->last_good, &last_pulse_pos, &st->clas_dec, &enr_q, st->core_brate );

                if( st->last_core != HQ_CORE || (st->last_core == HQ_CORE && st->last_con_tcx) )
                {
                    if( st->clas_dec == SIN_ONSET && last_pulse_pos != 0 && st->prev_bfi == 1 )
                    {
                        FEC_SinOnset( old_exc+L_EXC_MEM_DEC-L_EXC_MEM, last_pulse_pos, T0_tmp, enr_q, Aq, st->L_frame);
                    }
                    else if( (coder_type == GENERIC || coder_type == VOICED) && last_pulse_pos != 0 && st->old_bfi_cnt == 1  && output_frame == L_FRAME16k )
                    {
                        do_WI = FEC_enhACB(st->L_frame, st->last_L_frame, old_exc+L_EXC_MEM_DEC-L_EXC_MEM, T0_tmp, last_pulse_pos, st->bfi_pitch );
                    }
                }
            }

            /*------------------------------------------------------------*
             * In case of first frame after an erasure and transition from voiced to unvoiced or inactive
             * redo the LPC interpolation
             *------------------------------------------------------------*/

            if( st->stab_fac == 0 && st->old_bfi_cnt > 0 && st->clas_dec != VOICED_CLAS && st->clas_dec != ONSET && st->relax_prev_lsf_interp == 0 )
            {
                int_lsp4( st->L_frame, st->lsp_old, lsp_mid, lsp_new, Aq, M, 2 );
            }

            /*---------------------------------------------------------------*
             * Decoding of the scaled predicted innovation energy
             *---------------------------------------------------------------*/

            if( ( coder_type != UNVOICED && coder_type != AUDIO && coder_type != INACTIVE && !(st->core_brate <= ACELP_8k00 && coder_type != TRANSITION) )
                    || (coder_type == INACTIVE && st->total_brate >= ACELP_32k) )
            {
                nb_bits = Es_pred_bits_tbl[BIT_ALLOC_IDX(st->core_brate, coder_type, -1, -1)];
                indice = (short)get_next_indice( st, nb_bits );
                Es_pred_dec( &Es_pred, indice, nb_bits, 0 );
            }

            /*------------------------------------------------------------*
             * Decode excitation according to coding type
             *------------------------------------------------------------*/

            if( st->nelp_mode_dec )
            {
                /* SC-VBR - NELP frames */
                decod_nelp( st, coder_type, &tmp_noise, pitch_buf, exc, exc2, voice_factors, bwe_exc, st->bfi, gain_buf );
            }
            else if( coder_type == UNVOICED )
            {
                /* UNVOICED frames */
                decod_unvoiced( st, Aq, coder_type, &tmp_noise, pitch_buf, voice_factors, exc, exc2, bwe_exc, gain_buf );
            }
            else if( st->ppp_mode_dec )
            {
                /* SC-VBR - PPP frames */
                decod_ppp( st, Aq, pitch_buf, exc, exc2, voice_factors, bwe_exc, gain_buf, st->bfi );
            }
            else if( coder_type == TRANSITION )
            {
                decod_tran( st, st->L_frame, tc_subfr, Aq, coder_type, Es_pred, pitch_buf, voice_factors, exc, exc2, bwe_exc, unbits, sharpFlag, gain_buf );
            }
            else if( coder_type == AUDIO || ( coder_type == INACTIVE && st->core_brate <= ACELP_24k40 ) )
            {
                /* AUDIO and INACTIVE frames (coded by GSC technology) */
                decod_audio( st, dct_exc_tmp, Aq, coder_type, &tmp_noise, pitch_buf, voice_factors, exc, exc2, bwe_exc, lsf_new, gain_buf );
            }
            else
            {
                /* GENERIC, VOICED and INACTIVE frames (coded by AVQ technology) */
                decod_gen_voic( st, st->L_frame, sharpFlag, Aq, coder_type, Es_pred, do_WI, pitch_buf, voice_factors, exc, exc2, bwe_exc, unbits, gain_buf );
            }

            /* synthesis for ACELP core switching and SWB BWE */
            syn_12k8( st->L_frame, Aq, exc, temp_buf, st->mem_syn1, 1 );

            /* save and delay synthesis to be used by SWB BWE */
            save_old_syn( st->L_frame, temp_buf, old_syn_12k8_16k, st->old_syn_12k8_16k, st->preemph_fac, &st->mem_deemph_old_syn );

            /*-----------------------------------------------------------------*
             * Apply energy matching when switching to inactive frames
             *-----------------------------------------------------------------*/

            inact_switch_ematch( exc2, dct_exc_tmp, st->lt_ener_per_band, coder_type, st->L_frame, st->core_brate, st->bfi, st->last_core, st->last_codec_mode );

            /*------------------------------------------------------------*
             * Decode information and modify the excitation signal of stationary unvoiced frames
             *------------------------------------------------------------*/

            if( st->nelp_mode_dec != 1 )
            {
                stat_noise_uv_dec( st, coder_type, lsp_new, lsp_mid, Aq, exc2 );
            }

            /*------------------------------------------------------------*
             * Save filter memory in case the synthesis is redone after scaling
             * Synthesis at 12k8 Hz sampling rate
             *------------------------------------------------------------*/

            /* update past excitation signals for LD music post-filter */
            mvr2r( st->dct_post_old_exc + L_FRAME, st->dct_post_old_exc, DCT_L_POST-L_FRAME-OFFSET2 );
            mvr2r( exc2, st->dct_post_old_exc + (DCT_L_POST-L_FRAME-OFFSET2), L_FRAME );
            mvr2r( st->dct_post_old_exc, exc_buffer, DCT_L_POST-OFFSET2 );

            if( coder_type == AUDIO && !st->GSC_noisy_speech )
            {
                /* Extrapolation of the last future part, windowing and high resolution DCT transform */
                Prep_music_postP( exc_buffer, dct_buffer, st->filt_lfE, st->last_core, pitch_buf, st->LDm_enh_lp_gbin );

                /* LD music post-filter */
                LD_music_post_filter( dct_buffer, dct_buffer, st->core_brate, &st->LDm_last_music_flag,
                                      st->LDm_thres, &st->LDm_nb_thr_1, &st->LDm_nb_thr_3, st->LDm_lt_diff_etot,
                                      &st->LDm_mem_etot, st->LDm_enh_min_ns_gain, st->LDm_bckr_noise, st->LDm_enh_lf_EO,
                                      st->LDm_enh_lp_gbin, st->filt_lfE, &st->last_nonfull_music, AUDIO, st->last_coder_type );

                /* Inverse DCT transform, retrieval of the aligned excitation, re-synthesis */
                mvr2r( st->mem_syn2, mem_tmp, M );
                Post_music_postP( dct_buffer, exc_buffer, exc2, st->mem_syn2, st->mem_syn2, Aq, syn );
            }
            else
            {
                /* Core synthesis at 12.8kHz or 16kHz */
                mvr2r( st->mem_syn2, mem_tmp, M );
                syn_12k8( st->L_frame, Aq, exc2, syn, st->mem_syn2, 1 );

                for( i = 0; i < DCT_L_POST; i++ )
                {
                    st->filt_lfE[i] = 0.3f + 0.7f * st->filt_lfE[i];
                }
            }

            /*------------------------------------------------------------*
             * FEC - Estimate the classification information
             *------------------------------------------------------------*/

            FEC_clas_estim( syn, pitch_buf, st->L_frame, coder_type, st->codec_mode, st->mem_syn_clas_estim, &st->clas_dec,
                            &st->lp_ener_bfi, st->core_brate, st->Opt_AMR_WB, &st->decision_hyst, NULL, NULL, NULL,
                            NULL, NULL, NULL, temp_buf, 0, 0, 0,
                            0, 0, 0, st->last_core_brate );

            /*------------------------------------------------------------*
             * FEC - Estimate pitch
             *------------------------------------------------------------*/

            FEC_pitch_estim( st->Opt_AMR_WB, st->last_core, st->L_frame, st->clas_dec, st->last_good, pitch_buf, st->old_pitch_buf,
                             &st->bfi_pitch, &st->bfi_pitch_frame, &st->upd_cnt, coder_type );

            /*------------------------------------------------------------*
             * FEC - Smooth the speech energy evolution when recovering after a BAD frame
             * (smoothing is performed in the excitation domain and signal is resynthesized after)
             *------------------------------------------------------------*/

            FEC_scale_syn( st->L_frame, st->clas_dec, st->last_good, syn, pitch_buf, st->enr_old, enr_q, coder_type,
                           LSF_Q_prediction, &st->scaling_flag, &st->lp_ener_FEC_av, &st->lp_ener_FEC_max, st->bfi,
                           st->total_brate, st->prev_bfi, st->last_core_brate, exc, exc2, Aq, &st->old_enr_LP, mem_tmp, st->mem_syn2,
                           st->last_con_tcx && (st->L_frameTCX_past != st->L_frame) && (st->last_core != 0), 0 );

            /* estimate the pitch-synchronous speech energy per sample to be used when normal operation recovers */
            if( (st->total_brate == ACELP_7k20) || (st->total_brate == ACELP_8k00) )
            {
                fer_energy( st->L_frame, st->clas_dec, syn, pitch_buf[((st->L_frame)>>6)-1], &st->enr_old, st->L_frame );
            }
        }

    } /* End of GOOD FRAME */
    /*----------------------------------------------------------------*
     * BAD frame
     *----------------------------------------------------------------*/

    else
    {
        /* SC-VBR */
        if( st->last_nelp_mode_dec == 1 )
        {
            st->nelp_mode_dec = 1;
        }

        /* long burst frame erasures */
        if( st->nbLostCmpt > 5 && st->clas_dec >= VOICED_CLAS && st->clas_dec < INACTIVE_CLAS )
        {
            st->last_good = VOICED_TRANSITION;
        }

        /* LSF estimation and A(z) calculation */
        lsf_dec_bfi( MODE1, lsf_new, st->lsf_old, st->lsf_adaptive_mean, NULL, st->mem_MA, st->mem_AR,
                     st->stab_fac, st->last_coder_type, st->L_frame,  st->last_good,
                     st->nbLostCmpt, 0, NULL, NULL, NULL, st->Last_GSC_pit_band_idx, st->Opt_AMR_WB, st->bwidth );

        FEC_lsf2lsp_interp( st, st->L_frame, Aq, lsf_new, lsp_new );

        if( st->nelp_mode_dec == 1 )
        {
            /* SC-VBR */
            decod_nelp( st, coder_type, &tmp_noise, pitch_buf, exc, exc2, voice_factors, bwe_exc, st->bfi, gain_buf );
            FEC_pitch = pitch_buf[3];
        }
        else
        {
            /* calculation of excitation signal */
            FEC_exc_estim( st, st->L_frame, exc, exc2, dct_exc_tmp, pitch_buf, voice_factors, &FEC_pitch, bwe_exc, lsf_new, &tmp_noise );

            tmp_noise = st->lp_gainc;

            /* SC-VBR */
            st->prev_gain_pit_dec = st->lp_gainp;
        }

        /* synthesis for ACELP core switching and SWB BWE */
        syn_12k8( st->L_frame, Aq, exc, temp_buf, st->mem_syn1, 1 );

        /* save and delay synthesis to be used by SWB BWE */
        save_old_syn( st->L_frame, temp_buf, old_syn_12k8_16k, st->old_syn_12k8_16k, st->preemph_fac, &st->mem_deemph_old_syn );

        /* Apply energy matching when switching to inactive frames */
        inact_switch_ematch( exc2, dct_exc_tmp, st->lt_ener_per_band, coder_type, st->L_frame, st->core_brate, st->bfi, st->last_core, st->last_codec_mode );

        /* update past excitation signals for LD music post-filter */
        mvr2r( st->dct_post_old_exc + L_FRAME, st->dct_post_old_exc, DCT_L_POST-L_FRAME-OFFSET2 );
        mvr2r( exc2, st->dct_post_old_exc + (DCT_L_POST-L_FRAME-OFFSET2), L_FRAME );

        /* synthesis at 12k8 Hz sampling rate */
        if( (st->total_brate == ACELP_7k20) || (st->total_brate == ACELP_8k00) )
        {
            mvr2r( st->mem_syn2, mem_tmp, M );
        }
        syn_12k8( st->L_frame, Aq, exc2, syn, st->mem_syn2, 1 );

        /* update buffer for classifier */
        mvr2r( exc2 + st->L_frame - L_EXC_MEM, st->old_exc2, L_EXC_MEM );
        mvr2r( syn + st->L_frame - L_EXC_MEM, st->old_syn2, L_EXC_MEM );
        mvr2r( syn + st->L_frame - L_SYN_MEM_CLAS_ESTIM, st->mem_syn_clas_estim, L_SYN_MEM_CLAS_ESTIM );

        /* Update music post processing values */
        /* Filter energies update */
        for( i = 0; i < DCT_L_POST; i++ )
        {
            st->filt_lfE[i] = 0.3f + 0.7f*st->filt_lfE[i];
        }
        /* Update circular buffer, keep last energy difference unchanged */
        for (i = 1; i<MAX_LT; i++)
        {
            st->LDm_lt_diff_etot[i-1] = st->LDm_lt_diff_etot[i];
        }

        /*------------------------------------------------------------*
         * FEC - Smooth the speech energy evolution when recovering after a BAD frame
         * (smoothing is performed in the excitation domain and signal is resynthesized after)
         *------------------------------------------------------------*/

        if( (st->total_brate == ACELP_7k20) || (st->total_brate == ACELP_8k00) )
        {
            FEC_scale_syn( st->L_frame, st->clas_dec, st->last_good, syn, pitch_buf, st->enr_old, enr_q, coder_type,
                           LSF_Q_prediction, &st->scaling_flag, &st->lp_ener_FEC_av, &st->lp_ener_FEC_max, st->bfi,
                           st->total_brate, st->prev_bfi, st->last_core_brate, exc, exc2, Aq, &st->old_enr_LP, mem_tmp, st->mem_syn2,
                           st->last_con_tcx && (st->L_frameTCX_past != st->L_frame) && (st->last_core != 0), 0 );
        }

        /* estimate the pitch-synchronous speech energy per sample to be used when normal operation recovers */
        fer_energy( st->L_frame, st->last_good, syn, FEC_pitch, &st->enr_old, st->L_frame );

        if( st->nelp_mode_dec !=1 )
        {
            /* modify the excitation signal of stationary unvoiced frames */
            stat_noise_uv_mod( coder_type, 0, st->lsp_old, lsp_new, lsp_new, Aq, exc2, 1, &st->ge_sm, &st->uv_count, &st->act_count,
                               st->lspold_s, &st->noimix_seed, &st->min_alpha, &st->exc_pe, st->core_brate, st->bwidth );
        }
        /* SC-VBR */
        st->FadeScale = st->FadeScale*0.75;
    }


    if( st->L_frame == L_FRAME )
    {
        mvr2r( Aq+2*(M+1), st->cur_sub_Aq, (M+1) );
    }
    else
    {
        mvr2r( Aq+3*(M+1), st->cur_sub_Aq, (M+1) );
    }

    /*--------------------------------------------------------*
     * Apply NB postfilter in case of 8kHz output
     *--------------------------------------------------------*/

    if( st->last_bwidth == NB )
    {
        if( st->bwidth == NB )
        {
            st->pfstat.on = 1;
            nb_post_filt( st->L_frame, L_SUBFR, &(st->pfstat), &st->psf_lp_noise, tmp_noise, syn, Aq, pitch_buf, coder_type, st->BER_detect, 0 );
        }
        else
        {
            st->pfstat.on = 0;
            nb_post_filt( st->L_frame, L_SUBFR, &(st->pfstat), &st->psf_lp_noise, tmp_noise, syn, Aq, pitch_buf, AUDIO, st->BER_detect, 0 );
        }
    }
    else
    {
        st->psf_lp_noise = st->lp_noise;
    }

    /*------------------------------------------------------------------*
     * Perform fixed deemphasis through 1/(1 - g*z^-1)
     *-----------------------------------------------------------------*/

    /* update old synthesis buffer - needed for ACELP internal sampling rate switching */
    mvr2r( syn + st->L_frame - L_SYN_MEM, st->mem_syn_r, L_SYN_MEM );
    deemph( syn, st->preemph_fac, st->L_frame, &(st->mem_deemph) );
    AGC_dec(syn, st->agc_mem2, st->L_frame);

    mvr2r( syn+st->L_frame/2, st->old_syn_Overl, st->L_frame/2 );
    mvr2r( syn+st->L_frame-M-1, st->syn, M+1 );


    /*------------------------------------------------------------------*
     * Formant post-filter
     *-----------------------------------------------------------------*/

    if( st->last_bwidth>=WB && st->core_brate > ACELP_24k40 && st->core_brate <= ACELP_32k )
    {
        mvr2r( syn, temp_buf + L_SYN_MEM, L_FRAME16k );

        st->pfstat.on = 1;
        formant_post_filt( &(st->pfstat), temp_buf + L_SYN_MEM, Aq, syn, L_FRAME16k, L_SUBFR, st->lp_noise, st->total_brate, 0 );
    }
    else if( st->last_bwidth >= WB )
    {
        if( st->pfstat.on )
        {
            mvr2r( st->pfstat.mem_pf_in+L_SYN_MEM-M, temp_buf, M );
            mvr2r( syn, temp_buf+M, L_SUBFR );
            residu ( Aq, M, temp_buf+M,temp_buf+M+L_SUBFR, L_SUBFR );
            syn_filt ( Aq, M, temp_buf+M+L_SUBFR, temp_buf, L_SUBFR, st->pfstat.mem_stp+L_SYN_MEM-M, 0 );
            scale_st ( syn, temp_buf, &st->pfstat.gain_prec, L_SUBFR, -1 );
            mvr2r( temp_buf, syn, L_SUBFR/2 );
            blend_subfr2( temp_buf + L_SUBFR/2, syn + L_SUBFR/2, syn + L_SUBFR/2 );

        }
        st->pfstat.on = 0;
    }

    /*----------------------------------------------------------------*
     * Comfort noise addition
     *----------------------------------------------------------------*/

    if( st->flag_cna || (st->cng_type == FD_CNG && st->total_brate <= ACELP_32k) || (st->cng_type == LP_CNG && st->core_brate <= SID_2k40) )
    {
        /*VAD only for non inactive frame*/
        st->VAD = st->VAD && (coder_type != INACTIVE);

        /*Noisy speech detector*/
        noisy_speech_detection( st->VAD, syn, st->hFdCngDec->hFdCngCom->frameSize, st->hFdCngDec->msNoiseEst, st->hFdCngDec->psize_shaping,
                                st->hFdCngDec->nFFTpart_shaping, &(st->hFdCngDec->lp_noise), &(st->hFdCngDec->lp_speech), &(st->hFdCngDec->hFdCngCom->flag_noisy_speech) );

        st->hFdCngDec->hFdCngCom->likelihood_noisy_speech = 0.99f*st->hFdCngDec->hFdCngCom->likelihood_noisy_speech + 0.01f*(float)st->hFdCngDec->hFdCngCom->flag_noisy_speech;
        st->lp_noise = st->hFdCngDec->lp_noise;

        /*Noise estimate*/
        ApplyFdCng( syn, realBuffer, imagBuffer, st->hFdCngDec, st->m_frame_type, st, 0,
                    ( coder_type == AUDIO && !st->GSC_noisy_speech ) );

        /* CNA: Generate additional comfort noise to mask potential coding artefacts */
        if( st->flag_cna && coder_type != AUDIO )
        {
            generate_masking_noise( syn, st->hFdCngDec->hFdCngCom, st->hFdCngDec->hFdCngCom->frameSize, 0 );
        }
        else if( st->flag_cna && coder_type == AUDIO && st->last_core == ACELP_CORE && st->last_coder_type != AUDIO )
        {
            v_multc( st->hFdCngDec->hFdCngCom->olapBufferSynth2+5*st->hFdCngDec->hFdCngCom->frameSize/4, (float)(st->hFdCngDec->hFdCngCom->fftlen/2), temp_buf, st->hFdCngDec->hFdCngCom->frameSize/2);
            v_add( temp_buf, syn, syn, st->hFdCngDec->hFdCngCom->frameSize/2);
        }
    }

    if( st->flag_cna == 0 && st->L_frame == L_FRAME16k && st->last_flag_cna == 1 && ( (st->last_core == ACELP_CORE && st->last_coder_type != AUDIO) || st->last_core == AMR_WB_CORE) )
    {
        v_multc( st->hFdCngDec->hFdCngCom->olapBufferSynth2+5*st->L_frame/4, 256.f, temp_buf, st->L_frame/2 );
        v_add( temp_buf, syn, syn, st->L_frame/2 );
    }

    if( st->flag_cna == 0 || coder_type == AUDIO )
    {
        set_f( st->hFdCngDec->hFdCngCom->olapBufferSynth2, 0.f, st->hFdCngDec->hFdCngCom->fftlen );
    }

    /*----------------------------------------------------------------*
     * Resample to the output sampling rate (8/16/32/48 kHz)
     * Bass post-filter
     *----------------------------------------------------------------*/

    /* check if the CLDFB works on the right sample rate */
    if( (st->cldfbAna->no_channels * st->cldfbAna->no_col) != st->L_frame )
    {
        resampleCldfb( st->cldfbAna, st->L_frame*50 );
        resampleCldfb( st->cldfbBPF, st->L_frame*50 );

        if( st->ini_frame > 0 )
        {
            st->cldfbSyn->bandsToZero = st->cldfbSyn->no_channels - st->cldfbAna->no_channels;
        }
    }

    if( st->L_frame != st->last_L_frame && st->last_codec_mode != MODE2 )
    {
        if( st->L_frame == L_FRAME )
        {
            retro_interp5_4( st->pst_old_syn );
        }
        else if( st->L_frame == L_FRAME16k )
        {
            retro_interp4_5( syn, st->pst_old_syn );
        }
    }

    /* bass post-filter */
    bass_psfilter( st->Opt_AMR_WB, syn, st->L_frame, pitch_buf, st->pst_old_syn, &st->pst_mem_deemp_err, &st->pst_lp_ener, st->bpf_off, st->stab_fac,
                   &st->stab_fac_smooth, st->mem_mean_pit, st->Track_on_hist, st->vibrato_hist, &st->psf_att, coder_type, bpf_error_signal );

    /* analysis of the synthesis at internal sampling rate */
    cldfbAnalysis( syn, realBuffer, imagBuffer, -1, st->cldfbAna );

    /* analysis and add the BPF error signal */
    addBassPostFilter( bpf_error_signal, st->bpf_off?0:-1, realBuffer, imagBuffer, st->cldfbBPF );

    /* set output mask for upsampling */
    if( st->bwidth == NB )
    {
        /* set NB mask for upsampling */
        st->cldfbSyn->bandsToZero = st->cldfbSyn->no_channels - 10;
    }
    else if( st->cldfbSyn->bandsToZero != st->cldfbSyn->no_channels - st->cldfbAna->no_channels )
    {
        /* in case of BW switching, re-init to default */
        st->cldfbSyn->bandsToZero = st->cldfbSyn->no_channels - st->cldfbAna->no_channels;
    }

    /*WB/SWB-FD_CNG*/
    if( ( st->core_brate == FRAME_NO_DATA || st->core_brate == SID_2k40 ) && ( st->cng_type == FD_CNG ) && ( st->hFdCngDec->hFdCngCom->numCoreBands < st->cldfbSyn->no_channels ) )
    {
        generate_comfort_noise_dec_hf( realBuffer,imagBuffer, st );

        if( st->hFdCngDec->hFdCngCom->regularStopBand < st->cldfbSyn->no_channels )
        {
            st->cldfbSyn->bandsToZero = st->cldfbSyn->no_channels - st->hFdCngDec->hFdCngCom->regularStopBand;
        }
        else
        {
            st->cldfbSyn->bandsToZero = 0;
        }
    }

    /* synthesis of the combined signal */
    cldfbSynthesis( realBuffer, imagBuffer, synth, -1, st->cldfbSyn );

    /* save synthesis - needed in case of core switching */
    mvr2r( synth, st->previoussynth, output_frame );

    /*-----------------------------------------------------------------*
     * Bandwidth extension 6kHz-7kHz (only for 16kHz input signals)
     *-----------------------------------------------------------------*/

    if( (st->L_frame == L_FRAME && st->bwidth != NB && output_frame >= L_FRAME16k &&
            ( st->extl == -1 || st->extl == SWB_CNG || (st->extl == WB_BWE && st->extl_brate == 0 && coder_type != AUDIO) ) ) )
    {
        hf_synth( st->core_brate, output_frame, Aq, exc2, syn, synth, &st->seed2, st->mem_hp400,
                  st->mem_syn_hf, st->mem_hf, st->delay_syn_hf, st->mem_hp_interp );
    }
    else
    {
        hf_synth_reset( &st->seed2, st->mem_hf, st->mem_syn_hf, st->mem_hp400, st->mem_hp_interp, st->delay_syn_hf );
    }

    /*-----------------------------------------------------------------*
     * Populate parameters for SWB TBE
     *-----------------------------------------------------------------*/

    if( ( !st->bfi && st->prev_bfi) || (st->last_vbr_hw_BWE_disable_dec == 1 && st->vbr_hw_BWE_disable_dec == 0) || ((st->extl == SWB_TBE || st->extl == WB_TBE || st->extl == FB_TBE)&&st->last_extl != SWB_TBE && st->last_extl != WB_TBE && st->last_extl != FB_TBE) )
    {
        st->bwe_non_lin_prev_scale = 0.0f;
        set_f( st->old_bwe_exc_extended, 0.0f, NL_BUFF_OFFSET );
    }

    if( !st->ppp_mode_dec )
    {
        non_linearity( bwe_exc, bwe_exc_extended, st->old_bwe_exc_extended, L_FRAME32k, &st->bwe_non_lin_prev_scale, coder_type, voice_factors, st->L_frame );
    }

    if( st->core_brate == FRAME_NO_DATA || st->core_brate == SID_2k40 )
    {
        st->bwe_non_lin_prev_scale = 0.0f;
    }

    /*----------------------------------------------------------------------*
     * Updates
     *----------------------------------------------------------------------*/

    updt_dec( st, st->L_frame, coder_type, old_exc, pitch_buf, Es_pred, Aq, lsf_new, lsp_new, voice_factors, old_bwe_exc, gain_buf );

    if( st->core_brate > SID_2k40 )
    {
        /* update CNG parameters in active frames */
        cng_params_upd( lsp_new, exc, st->L_frame, &st->ho_circ_ptr, st->ho_ener_circ,
                        &st->ho_circ_size, st->ho_lsp_circ, DEC, st->ho_env_circ, NULL,
                        NULL, NULL, st->last_active_brate );

        /* Set 16k LSP flag for CNG buffer */
        st->ho_16k_lsp[st->ho_circ_ptr] = (st->L_frame == L_FRAME ? 0 : 1 );
    }

    return;
}
