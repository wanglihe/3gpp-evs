/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*------------------------------------------------------------------*
 * amr_wb_dec()
 *
 * AMR-WB decoder
 *------------------------------------------------------------------*/

void amr_wb_dec(
    Decoder_State *st,              /* i/o: decoder state structure */
    float *output           /* o  : synthesis output        */
)
{
    short i;
    short vad_flag;
    short coder_type;
    float ftmp;
    float old_exc[L_EXC_DEC], *exc;     /* excitation signal buffer              */
    float syn_tmp[L_FRAME+2], *syn;     /* synthesis signal buffer               */
    float synth_out[L_FRAME48k];        /* synthesis output                      */
    short output_frame;                 /* frame length at output sampling freq. */
    float lsf_new[M];                   /* LSFs at the end of the frame          */
    float lsp_new[M];                   /* LSPs at the end of the frame          */
    float Aq[NB_SUBFR*(M+1)];           /* A(q) quantized for the 4 subframes    */
    float exc2[L_FRAME];                 /* total excitation buffer               */
    float mem_tmp[M];                   /* temporary synthesis filter memory     */
    float pitch_buf[NB_SUBFR];          /* floating pitch for each subframe      */
    float enr_q;                        /* E information for FER protection      */
    float tmp_noise;                    /* Long term temporary noise energy      */
    float FEC_pitch;                    /* FEC pitch                             */
    float dummy_buf[L_FRAME32k];        /* dummy buffer - no usage               */
    short allow_cn_step;
    short locattack, amr_io_class;
    short tmps;
    float xsp_tmp[M];
    float tmp_buffer[L_FRAME48k];
    float dct_buffer[DCT_L_POST];
    float frame_e;
    float exc_buffer[DCT_L_POST];
    float class_para;
    float voice_factors[NB_SUBFR];
    short hf_gain[NB_SUBFR];
    short delay_comp;
    short nZeros;
    float tmp;
    short last_core_ori;
    short sid_bw = 0;
    float bpf_error_signal[L_FRAME];
    float ng_ener;
    float gain_buf[NB_SUBFR16k];
    float q_env[20];
    float exc3[L_FRAME];
    float *realBuffer[CLDFB_NO_COL_MAX], *imagBuffer[CLDFB_NO_COL_MAX];
    float realBufferTmp[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX], imagBufferTmp[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX];
    short flag_cna;
    short waveadj_rec = 0;


    /*------------------------------------------------------------------*
     * Initialization
     *------------------------------------------------------------------*/

    syn_tmp[0] = 0;
    syn_tmp[1] = 0;
    syn = syn_tmp+2;

    for( i=0; i<CLDFB_NO_COL_MAX; i++ )
    {
        set_f( realBufferTmp[i], 0, CLDFB_NO_CHANNELS_MAX );
        set_f( imagBufferTmp[i], 0, CLDFB_NO_CHANNELS_MAX );
        realBuffer[i] = realBufferTmp[i];
        imagBuffer[i] = imagBufferTmp[i];
    }

    set_f( gain_buf, 0, NB_SUBFR16k );

    st->use_partial_copy = 0;
    st->rf_flag = 0;
    st->rf_flag_last = 0;

    st->L_frame = L_FRAME;
    st->nb_subfr = NB_SUBFR;
    st->core = AMR_WB_CORE;
    st->core_brate = st->total_brate;
    st->extl = -1;
    st->bwidth = WB;
    coder_type = GENERIC;
    output_frame = (short)(st->output_Fs / 50);                   /* frame length of the input signal */

    st->bpf_off = 0;
    if( st->last_core == HQ_CORE )
    {
        st->bpf_off = 1;
        st->pfstat.on = 0;
    }
    st->igf = 0;

    st->sr_core     = st->L_frame*50;
    st->fscale_old  = st->fscale;
    st->fscale      = sr2fscale(st->sr_core);

    /* Initialization in case that the first frame is the good received AMR-WB (IO) frame */
    if( st->ini_frame == 0 )
    {
        st->last_core = AMR_WB_CORE;
        mvr2r( mean_isf_amr_wb, st->lsf_old, M );
        isf2isp( st->lsf_old, st->lsp_old, M, INT_FS_12k8 );
    }

    /* Updates in case of EVS primary mode -> AMR-WB IO mode switching */
    if( st->last_core != AMR_WB_CORE )
    {
        updt_IO_switch_dec( output_frame, st );
    }

    /* Updates in case of HQ -> AMR-WB IO switching */
    core_switching_pre_dec( st, output_frame );

    last_core_ori = st->last_core;
    set_s( hf_gain, 0, NB_SUBFR );

    enr_q = 0.0f;
    tmp_noise = 0.0f;
    amr_io_class = UNVOICED_CLAS;

    mvr2r( st->old_exc, old_exc, L_EXC_MEM_DEC );
    exc = old_exc + L_EXC_MEM_DEC;

    /* reset post-filter in case of switching */
    if( st->pfstat.on == 0 )
    {
        st->pfstat.reset = 1;
    }

    if( st->bfi )
    {
        st->nbLostCmpt++;
    }
    else
    {
        st->nbLostCmpt = 0;
    }

    /* PLC: [TCX: Fade-out-recovery]
     * PLC: overlapping part needs to be attenuated for first good frame */
    if (!st->bfi && st->prev_bfi && (st->last_codec_mode == MODE2) && (st->last_core_bfi == TCX_20_CORE || st->last_core_bfi == TCX_10_CORE))
    {
        v_multc( st->old_out, st->plcInfo.recovery_gain, st->old_out, st->L_frameTCX );
    }

    /*-----------------------------------------------------------------*
     * switching from ACELP@16k core to AMR-WB IO mode
     *-----------------------------------------------------------------*/

    st->rate_switching_reset = 0;

    if( st->last_core != AMR_WB_CORE && st->last_L_frame == L_FRAME16k  && st->last_core != HQ_CORE )
    {
        /* in case of switching, do not apply BPF */
        st->bpf_off = 1;
        if( st->pfstat.on != 0 )
        {
            short mem_syn_r_size_old, mem_syn_r_size_new;

            mem_syn_r_size_old = (short)(1.25*st->last_L_frame/20.f);
            mem_syn_r_size_new = (short)(1.25*st->L_frame/20.f);
            lerp( st->pfstat.mem_stp+L_SYN_MEM-mem_syn_r_size_old, st->pfstat.mem_stp+L_SYN_MEM-mem_syn_r_size_new, mem_syn_r_size_new, mem_syn_r_size_old );
            lerp( st->pfstat.mem_pf_in+L_SYN_MEM-mem_syn_r_size_old, st->pfstat.mem_pf_in+L_SYN_MEM-mem_syn_r_size_new, mem_syn_r_size_new, mem_syn_r_size_old );
        }

        /* convert old quantized LSP vector */
        st->rate_switching_reset = lsp_convert_poly( st->lsp_old, L_FRAME, 1 );

        /* convert old quantized LSF vector */
        lsp2lsf( st->lsp_old, st->lsf_old, M, INT_FS_12k8 );

        /* FEC - update adaptive LSF mean vector */
        mvr2r( st->lsf_old, st->lsfoldbfi1, M );
        mvr2r( st->lsf_old, st->lsfoldbfi0, M );
        mvr2r( st->lsf_old, st->lsf_adaptive_mean, M );

        /* Reset LPC mem */
        mvr2r( GEWB_Ave, st->mem_AR, M );
        set_zero( st->mem_MA, M );

        /* update synthesis filter memories */
        synth_mem_updt2( L_FRAME, st->last_L_frame, st->old_exc, st->mem_syn_r, st->mem_syn2, NULL, DEC );
        mvr2r( st->old_exc, old_exc, L_EXC_MEM_DEC );
        mvr2r( st->mem_syn2, st->mem_syn1, M );
        mvr2r( st->mem_syn2, st->mem_syn3, M );

        /* LSP -> ISP */
        mvr2r( stable_ISP, xsp_tmp, M );
        lsp2isp( st->lsp_old, st->lsp_old, xsp_tmp, M );

    }

    /* update buffer of old subframe pitch values */
    if( st->last_L_frame != L_FRAME )
    {
        if( st->last_L_frame == L_FRAME32k )
        {
            tmp = (float)12800/(float)32000;
        }
        else if( st->last_L_frame == 512 )
        {
            tmp = (float)12800/(float)25600;
        }
        else /* st->last_L_frame == L_FRAME16k */
        {
            tmp = (float)12800/(float)16000;
        }

        for( i=NB_SUBFR16k-NB_SUBFR; i<NB_SUBFR16k; i++ )
        {
            st->old_pitch_buf[i-1] = tmp * st->old_pitch_buf[i];
        }

        for( i=2*NB_SUBFR16k-NB_SUBFR; i<2*NB_SUBFR16k; i++ )
        {
            st->old_pitch_buf[i-2] = tmp * st->old_pitch_buf[i];
        }
    }

    if( st->bfi_pitch_frame != L_FRAME )
    {
        if( st->bfi_pitch_frame == L_FRAME32k )
        {
            tmp = (float)12800/(float)32000;
        }
        else if( st->bfi_pitch_frame == 512 )
        {
            tmp = (float)12800/(float)25600;
        }
        else /* st->bfi_pitch_frame == L_FRAME16k */
        {
            tmp = (float)12800/(float)16000;
        }

        st->bfi_pitch *= tmp;
        st->bfi_pitch_frame = L_FRAME;
    }


    if( st->last_core != AMR_WB_CORE )
    {
        /* reset the unvoiced/audio signal improvement memories */
        isp2a( st->lsp_old, st->old_Aq, M );
        mvr2r( st->old_Aq, st->old_Aq + (M+1), M+1 );
        mvr2r( st->old_Aq, st->old_Aq + 2*(M+1), M+1 );
        mvr2r( st->old_Aq, st->old_Aq + 3*(M+1), M+1 );
    }

    if( st->last_bwidth == NB && st->ini_frame != 0 )
    {
        st->rate_switching_reset = 1;
    }

    /*----------------------------------------------------------------------*
     * GOOD frame
     *----------------------------------------------------------------------*/

    if( !st->bfi )
    {
        /*----------------------------------------------------------------*
         * Processing of FRAME_NO_DATA frames
         * Decoding of SID frames
         *----------------------------------------------------------------*/

        if ( st->core_brate == FRAME_NO_DATA || st->core_brate == SID_1k75 )
        {
            /* decode CNG parameters */
            CNG_dec( st, L_FRAME, Aq, st->core_brate, lsp_new, lsf_new, &allow_cn_step, &sid_bw, q_env );

            /* comfort noise generation */
            CNG_exc( st->core_brate, L_FRAME, &st->Enew, &st->cng_seed, exc, exc2, &st->lp_ener, st->last_core_brate,
                     &st->first_CNG, &st->cng_ener_seed, dummy_buf, allow_cn_step, &st->last_allow_cn_step, st->num_ho,
                     q_env, st->lp_env, st->old_env, st->exc_mem, st->exc_mem1, &sid_bw, &st->cng_ener_seed1, exc3, st->Opt_AMR_WB );

            set_f( voice_factors, 1.0f, NB_SUBFR );
            class_para = 0.0f;

            if ( st->first_CNG == 0 )
            {
                st->first_CNG = 1;
            }

            /* update past excitation signals for LD music post-filter */
            mvr2r( st->dct_post_old_exc + L_FRAME, st->dct_post_old_exc, DCT_L_POST-L_FRAME-OFFSET2 );
            mvr2r( exc2, st->dct_post_old_exc + (DCT_L_POST-L_FRAME-OFFSET2), L_FRAME );

            /* synthesis at 12k8 Hz sampling rate */
            syn_12k8( L_FRAME, Aq, exc2, syn, st->mem_syn2, 1 );
            syn_12k8( st->L_frame, Aq, exc2, dummy_buf, st->mem_syn3, 1 );

            /* reset the decoder */
            CNG_reset_dec( st, pitch_buf, dummy_buf+L_FRAME );

            /* update st->mem_syn1 for ACELP core switching */
            mvr2r( st->mem_syn3, st->mem_syn1, M );
            if( output_frame != L_FRAME8k )
            {
                frame_e = 10.0f * (float) log10( dotp( syn, syn, L_FRAME ) / (float) L_FRAME );
                st->psf_lp_noise = 0.99f * st->psf_lp_noise + 0.01f * frame_e;
            }

            /* update old synthesis for classification */
            mvr2r( syn + L_FRAME - L_SYN_MEM_CLAS_ESTIM, st->mem_syn_clas_estim, L_SYN_MEM_CLAS_ESTIM );

            /* Update music post processing values */
            /* Filter energies update */
            for( i = 0; i < DCT_L_POST; i++ )
            {
                st->filt_lfE[i] = 0.3f + 0.7f*st->filt_lfE[i];
            }

            vad_flag = 0;
        }

        /*----------------------------------------------------------------*
         * Decoding of all other frames
         *----------------------------------------------------------------*/

        else
        {
            /*-----------------------------------------------------------------*
             * After CNG period, use the most up-to-date ISPs
             *-----------------------------------------------------------------*/

            if( st->last_core_brate == FRAME_NO_DATA || st->last_core_brate == SID_1k75 )
            {
                mvr2r( st->lspCNG, st->lsp_old, M );
                isp2isf( st->lspCNG, st->lsf_old, M, INT_FS_12k8 );
                set_f( old_exc, 0, L_EXC_MEM_DEC );
            }

            /*------------------------------------------------------------*
             * Extracts VAD information from the bitstream in AMR-WB IO mode
             *------------------------------------------------------------*/

            vad_flag = (short)get_next_indice( st, 1 );

            if( vad_flag == 0 )
            {
                coder_type = INACTIVE;
            }
            else
            {
                coder_type = GENERIC;
            }

            /*-----------------------------------------------------------------*
             * ISF de-quantization and interpolation
             *-----------------------------------------------------------------*/

            isf_dec_amr_wb( st, Aq, lsf_new, lsp_new );

            /*------------------------------------------------------------*
             * Decode excitation
             *------------------------------------------------------------*/

            decod_amr_wb( st, Aq, pitch_buf, exc, exc2, hf_gain, voice_factors, gain_buf );

            /* synthesis for ACELP core switching and SWB BWE */
            syn_12k8( L_FRAME, Aq, exc, tmp_buffer, st->mem_syn1, 1 );

            /*------------------------------------------------------------*
             * Update long-term energies for FEC
             * Update ISP vector for CNG
             *------------------------------------------------------------*/

            if( coder_type == INACTIVE )
            {
                if( st->unv_cnt > 20 )
                {
                    ftmp = st->lp_gainc * st->lp_gainc;
                    st->lp_ener = 0.7f * st->lp_ener + 0.3f * ftmp;
                    for( i=0; i<M; i++ )
                    {
                        st->lspCNG[i] = (float)(0.9f * st->lspCNG[i] + 0.1f * lsp_new[i]);
                    }
                }
                else
                {
                    st->unv_cnt++;
                }
            }
            else
            {
                st->unv_cnt = 0;
            }

            /*------------------------------------------------------------*
             * Save filter memory in case the synthesis is redone after scaling
             * Core synthesis at 12k8 Hz
             *------------------------------------------------------------*/

            mvr2r( st->mem_syn2, mem_tmp, M );
            syn_12k8( L_FRAME, Aq, exc2, syn, st->mem_syn2, 1 );

            /*------------------------------------------------------------*
             * FEC - Estimate the classification information
             *------------------------------------------------------------*/

            FEC_clas_estim( syn, pitch_buf, st->L_frame, coder_type, st->codec_mode, st->mem_syn_clas_estim, &st->clas_dec,
                            &st->lp_ener_bfi, st->core_brate, st->Opt_AMR_WB, &st->decision_hyst, &locattack, &st->UV_cnt,
                            &st->LT_UV_cnt, &st->Last_ener, &amr_io_class, st->lt_diff_etot, &class_para, 0, 0, 0,
                            0, 0, 0, st->last_core_brate );

            /* update past excitation signals for LD music post-filter */
            mvr2r( st->dct_post_old_exc + L_FRAME, st->dct_post_old_exc, DCT_L_POST-L_FRAME-OFFSET2 );
            mvr2r( exc2, st->dct_post_old_exc + (DCT_L_POST-L_FRAME-OFFSET2), L_FRAME );
            mvr2r( st->dct_post_old_exc, exc_buffer, DCT_L_POST-OFFSET2 );

            if( output_frame != L_FRAME8k )
            {
                if ( coder_type == INACTIVE )
                {
                    frame_energy( L_FRAME, pitch_buf, syn, 0.0f, &frame_e );
                    st->psf_lp_noise = 0.99f * st->psf_lp_noise + 0.01f * frame_e;
                }
            }

            if( amr_io_class != UNVOICED_CLAS && coder_type != INACTIVE && st->psf_lp_noise < 15.0f )
            {
                short tmp_coder_type = AUDIO;
                if( st->last_coder_type == INACTIVE || st->last_coder_type == UNVOICED )
                {
                    tmp_coder_type = INACTIVE;
                }
                /* Extrapolation of the last future part, windowing and high resolution DCT transform */
                Prep_music_postP( exc_buffer, dct_buffer, st->filt_lfE, st->last_core, pitch_buf, st->LDm_enh_lp_gbin );

                /* LD music post-filter */
                LD_music_post_filter( dct_buffer, dct_buffer, st->core_brate, &st->LDm_last_music_flag,
                                      st->LDm_thres, &st->LDm_nb_thr_1, &st->LDm_nb_thr_3, st->LDm_lt_diff_etot,
                                      &st->LDm_mem_etot, st->LDm_enh_min_ns_gain, st->LDm_bckr_noise, st->LDm_enh_lf_EO,
                                      st->LDm_enh_lp_gbin, st->filt_lfE, &st->last_nonfull_music, -1, tmp_coder_type );

                /* Inverse DCT transform, retrieval of the aligned excitation, re-synthesis */
                Post_music_postP( dct_buffer, exc_buffer, exc2, mem_tmp, st->mem_syn2, Aq, syn );
            }
            else
            {
                /*------------------------------------------------------------*
                 * Improvement for unvoiced and audio signals
                 *------------------------------------------------------------*/

                improv_amr_wb_gs( amr_io_class, coder_type, st->core_brate, &st->seed_tcx, st->old_Aq, st->mem_syn2, st->lt_voice_fac,
                                  locattack, Aq, exc2, mem_tmp, syn, pitch_buf, st->Last_ener, st->rate_switching_reset, st->last_coder_type );

                for( i = 0; i < DCT_L_POST; i++ )
                {
                    st->filt_lfE[i] = 0.3f + 0.7f * st->filt_lfE[i] ;
                }
            }

            /*------------------------------------------------------------*
             * FEC - Estimate pitch
             *------------------------------------------------------------*/

            FEC_pitch_estim( 1, st->last_core, L_FRAME, st->clas_dec, st->last_good, pitch_buf, st->old_pitch_buf,
                             &st->bfi_pitch, &st->bfi_pitch_frame, &st->upd_cnt, GENERIC );

            /*------------------------------------------------------------*
             * FEC - Smooth the speech energy evolution when recovering after a BAD frame
             * (smoothing is performed in the excitation domain and signal is resynthesized after)
             *------------------------------------------------------------*/

            FEC_scale_syn( L_FRAME, st->clas_dec, st->last_good, syn, pitch_buf, st->enr_old, enr_q, -1,
                           MOVING_AVERAGE, &st->scaling_flag, &st->lp_ener_FEC_av, &st->lp_ener_FEC_max, st->bfi,
                           st->total_brate, st->prev_bfi, st->last_core_brate, exc, exc2, Aq, &st->old_enr_LP, mem_tmp, st->mem_syn2,
                           st->last_con_tcx && (st->L_frameTCX_past != st->L_frame) && (st->last_core != 0), 0 );

            /* estimate the pitch-synchronous speech energy per sample to be used when normal operation recovers */
            fer_energy( L_FRAME, st->clas_dec, syn, pitch_buf[3], &st->enr_old, L_FRAME );
        }

    } /* End of GOOD FRAME */

    /*----------------------------------------------------------------*
     * BAD frame
     *----------------------------------------------------------------*/

    else
    {
        /* long burst frame erasures */
        if( st->nbLostCmpt > 5 && st->clas_dec >= VOICED_CLAS )
        {
            st->last_good = VOICED_TRANSITION;
        }

        vad_flag = st->last_vad;
        amr_io_class = st->last_good;
        class_para = 0.0f;

        /* LSF estimation and A(z) calculation */
        lsf_dec_bfi( MODE1, lsf_new, st->lsf_old, st->lsf_adaptive_mean, NULL, st->mem_MA, st->mem_AR,
                     st->stab_fac, st->last_coder_type, L_FRAME,  st->last_good,
                     st->nbLostCmpt, 0, NULL, NULL, NULL, st->Last_GSC_pit_band_idx, st->Opt_AMR_WB, st->bwidth );

        FEC_lsf2lsp_interp( st, L_FRAME, Aq, lsf_new, lsp_new );

        /* calculation of excitation signal */
        FEC_exc_estim( st, L_FRAME, exc, exc2, tmp_buffer, pitch_buf, voice_factors, &FEC_pitch, dummy_buf, lsf_new, &tmp_noise );

        /* synthesis for ACELP core switching and SWB BWE */
        syn_12k8( L_FRAME, Aq, exc, tmp_buffer, st->mem_syn1, 1 );

        /* update past excitation signals */
        mvr2r( st->dct_post_old_exc + L_FRAME, st->dct_post_old_exc, DCT_L_POST-L_FRAME-OFFSET2 );
        mvr2r( exc2, st->dct_post_old_exc + (DCT_L_POST-L_FRAME-OFFSET2), L_FRAME );

        /* Update music post processing values */
        /* Update circular buffer, keep last energy difference unchanged */
        for (i = 1; i<MAX_LT; i++)
        {
            st->LDm_lt_diff_etot[i-1] = st->LDm_lt_diff_etot[i];
        }
        /* Filter energies update */
        for( i = 0; i < DCT_L_POST; i++ )
        {
            st->filt_lfE[i] = 0.3f + 0.7f * st->filt_lfE[i];
        }

        /* synthesis at 12k8 Hz sampling rate */
        mvr2r( st->mem_syn2, mem_tmp, M );
        syn_12k8( L_FRAME, Aq, exc2, syn, st->mem_syn2, 1 );

        /* update old synthesis for classification */
        mvr2r( syn + L_FRAME - L_SYN_MEM_CLAS_ESTIM, st->mem_syn_clas_estim, L_SYN_MEM_CLAS_ESTIM );


        /*------------------------------------------------------------*
         * FEC - Smooth the speech energy evolution when recovering after a BAD frame
         * (smoothing is performed in the excitation domain and signal is resynthesized after)
         *------------------------------------------------------------*/

        FEC_scale_syn( L_FRAME, st->clas_dec, st->last_good, syn, pitch_buf, st->enr_old, enr_q, -1,
                       MOVING_AVERAGE, &st->scaling_flag, &st->lp_ener_FEC_av, &st->lp_ener_FEC_max, st->bfi,
                       st->total_brate, st->prev_bfi, st->last_core_brate, exc, exc2, Aq, &st->old_enr_LP, mem_tmp, st->mem_syn2,
                       0,0 );

        /* estimate the pitch-synchronous speech energy per sample to be used when normal operation recovers */
        fer_energy( L_FRAME, st->last_good, syn, FEC_pitch, &st->enr_old, L_FRAME );
    }

    /*--------------------------------------------------------*
     * NB post-filter
     *--------------------------------------------------------*/

    if( output_frame == L_FRAME8k || st->last_bwidth == NB )
    {
        if( output_frame == L_FRAME8k )
        {
            st->pfstat.on = 1;
            nb_post_filt( L_FRAME, L_SUBFR, &(st->pfstat), &st->psf_lp_noise, tmp_noise, syn, Aq, pitch_buf, coder_type, st->BER_detect, 0 );
        }
        else
        {
            st->pfstat.on = 0;
            nb_post_filt( L_FRAME, L_SUBFR, &(st->pfstat), &st->psf_lp_noise, tmp_noise, syn, Aq, pitch_buf, AUDIO, st->BER_detect, 0);
        }
    }

    /*------------------------------------------------------------------*
     * Perform fixed deemphasis through 1/(1 - g*z^-1)
     *-----------------------------------------------------------------*/

    /* update old synthesis buffer - needed for ACELP internal sampling rate switching */
    mvr2r( syn + L_FRAME - L_SYN_MEM, st->mem_syn_r, L_SYN_MEM );


    deemph( syn, PREEMPH_FAC, L_FRAME, &(st->mem_deemph) );
    AGC_dec( syn, st->agc_mem2, L_FRAME );
    mvr2r( syn+L_FRAME/2, st->old_syn_Overl, L_FRAME/2 );
    mvr2r( syn+L_FRAME-M-1, st->syn, M+1 );

    /*------------------------------------------------------------------*
     * Formant post-filter
     *-----------------------------------------------------------------*/

    mvr2r( syn, tmp_buffer + L_SYN_MEM, L_FRAME );

    if( output_frame != L_FRAME8k && st->last_bwidth != NB )
    {
        st->pfstat.on = 1;
        formant_post_filt( &(st->pfstat), tmp_buffer + L_SYN_MEM, Aq, syn, L_FRAME, L_SUBFR, st->psf_lp_noise, st->total_brate, amr_io_class == AUDIO_CLAS );
    }

    /*----------------------------------------------------------------*
     * Comfort Noise Addition
     *----------------------------------------------------------------*/

    flag_cna = 0;
    if( (st->psf_lp_noise >= 15.f) || (coder_type == INACTIVE) )
    {
        /*VAD only for non inactive frame*/
        st->VAD = (st->VAD && (coder_type != INACTIVE));

        ApplyFdCng( syn, NULL, NULL, st->hFdCngDec, st->m_frame_type, st, 0, 0 );

        st->hFdCngDec->hFdCngCom->frame_type_previous = st->m_frame_type;

        /*Noisy speech detector*/
        noisy_speech_detection( st->VAD, syn, st->hFdCngDec->hFdCngCom->frameSize, st->hFdCngDec->msNoiseEst, st->hFdCngDec->psize_shaping,
                                st->hFdCngDec->nFFTpart_shaping, &(st->hFdCngDec->lp_noise), &(st->hFdCngDec->lp_speech), &(st->hFdCngDec->hFdCngCom->flag_noisy_speech));

        st->hFdCngDec->hFdCngCom->likelihood_noisy_speech = 0.99f*st->hFdCngDec->hFdCngCom->likelihood_noisy_speech + 0.01f*(float)st->hFdCngDec->hFdCngCom->flag_noisy_speech;
        st->lp_noise = st->hFdCngDec->lp_noise;

        if( st->flag_cna && (st->psf_lp_noise >= 15.f) )
        {
            flag_cna = 1;
            generate_masking_noise( syn, st->hFdCngDec->hFdCngCom, st->hFdCngDec->hFdCngCom->frameSize, AMR_WB_CORE );
        }
        else if ( st->flag_cna )
        {
            generate_masking_noise_update_seed(st->hFdCngDec->hFdCngCom);
        }
    }
    else if ( st->flag_cna )
    {
        generate_masking_noise_update_seed(st->hFdCngDec->hFdCngCom);
    }

    if( flag_cna == 0 )
    {
        if( st->last_flag_cna == 1 && ( (st->last_core == ACELP_CORE && st->last_coder_type != AUDIO) || st->last_core == AMR_WB_CORE) )
        {
            v_multc( st->hFdCngDec->hFdCngCom->olapBufferSynth2+5*L_FRAME/4, 256.f, tmp_buffer, L_FRAME/2 );
            v_add( tmp_buffer, syn, syn, L_FRAME/2 );
        }
        set_f( st->hFdCngDec->hFdCngCom->olapBufferSynth2, 0.f, L_FRAME*2 );
    }

    /*----------------------------------------------------------------*
     * Change the sampling frequency to 8/16/32 kHz
     * Bass post-filter
     *----------------------------------------------------------------*/

    /* check if the cldfb works on the right sample rate */
    if( (st->cldfbAna->no_channels * st->cldfbAna->no_col) != L_FRAME )
    {
        /* resample to ACELP internal sampling rate */
        resampleCldfb (st->cldfbAna, INT_FS_12k8);
        resampleCldfb (st->cldfbBPF, INT_FS_12k8);

        if( st->ini_frame > 0 )
        {
            st->cldfbSyn->bandsToZero = st->cldfbSyn->no_channels - st->cldfbAna->no_channels;
        }
    }

    /* bass post-filter */
    bass_psfilter( st->Opt_AMR_WB, syn, L_FRAME, pitch_buf, st->pst_old_syn,
                   &st->pst_mem_deemp_err, &st->pst_lp_ener, st->bpf_off, st->stab_fac, &st->stab_fac_smooth,
                   st->mem_mean_pit, st->Track_on_hist, st->vibrato_hist, &st->psf_att, GENERIC, bpf_error_signal );

    /* analysis of the synthesis at internal sampling rate */
    cldfbAnalysis( syn, realBuffer, imagBuffer, -1, st->cldfbAna );

    /* analysis and add the BPF error signal */
    addBassPostFilter( bpf_error_signal, st->bpf_off?0:-1, realBuffer, imagBuffer, st->cldfbBPF );

    if( st->cldfbSyn->bandsToZero != st->cldfbSyn->no_channels - st->cldfbAna->no_channels )
    {
        /* in case of BW switching, re-init to default */
        st->cldfbSyn->bandsToZero = st->cldfbSyn->no_channels-st->cldfbAna->no_channels;
    }

    cldfb_synth_set_bandsToZero( st, realBuffer, imagBuffer, CLDFB_NO_COL_MAX );

    /* synthesis of the combined signal */
    cldfbSynthesis( realBuffer, imagBuffer, synth_out, -1, st->cldfbSyn );

    /* save synthesis - needed in case of core switching */
    mvr2r( synth_out, st->previoussynth, output_frame );

    /*--------------------------------------------------------*
     * calculate the average frame energy
     *--------------------------------------------------------*/

    fer_energy( L_FRAME, st->clas_dec, syn, pitch_buf[3], &ng_ener, L_FRAME );

    /*--------------------------------------------------------*
     * optimized for NO_S@-26dBov with street noise @ SNR=25dB
     *--------------------------------------------------------*/

    ng_ener = 10.0f * (float)log10(ng_ener + 0.01f) - 90.3087f + 15;
    st->ng_ener_ST = 0.7f * st->ng_ener_ST + 0.3f * ng_ener;

    /*-----------------------------------------------------------------*
     * Bandwidth extension 6kHz-8kHz
     *-----------------------------------------------------------------*/
    if( output_frame >= L_FRAME16k && ((st->cldfbSyn->bandsToZero -  st->cldfbSyn->no_channels + 10 ) != 0 || st->last_flag_filter_NB != 1) )
    {
        hf_synth_amr_wb( st->core_brate, output_frame, Aq, exc2, syn, st->mem_syn_hf, st->delay_syn_hf, &st->prev_r, &st->fmerit_w_sm, &amr_io_class, st->mem_hp_interp,
                         synth_out, class_para, hf_gain, voice_factors, pitch_buf, st->ng_ener_ST, lsf_new, &st->frame_count, &st->ne_min, &st->fmerit_m_sm, &st->voice_fac_amr_wb_hf,
                         &st->unvoicing, &st->unvoicing_sm, &st->unvoicing_flag, &st->voicing_flag, &st->start_band_old, &st->OptCrit_old );
    }
    else
    {
        hf_synth_amr_wb_reset( &st->seed2, st->mem_syn_hf, st->mem_hp_interp, &st->prev_r, &st->fmerit_w_sm, st->delay_syn_hf,
                               &st->frame_count, &st->ne_min, &st->fmerit_m_sm, &st->voice_fac_amr_wb_hf, &st->unvoicing,
                               &st->unvoicing_sm, &st->unvoicing_flag, &st->voicing_flag, &st->start_band_old, &st->OptCrit_old );
    }
    /*----------------------------------------------------------------------*
     * Updates
     *----------------------------------------------------------------------*/

    updt_dec( st, L_FRAME, coder_type, old_exc, pitch_buf, 0, Aq, lsf_new, lsp_new, voice_factors, dummy_buf, gain_buf );

    /* update old_Aq[] - needed in improv_amr_wb_gs_fx() */
    mvr2r( Aq, st->old_Aq, NB_SUBFR * (M+1) );

    if( !st->bfi
            && st->prev_bfi
            && st->last_total_brate >= HQ_48k
            && st->last_codec_mode == MODE2
            && (st->last_core_bfi == TCX_20_CORE || st->last_core_bfi == TCX_10_CORE)
            && st->plcInfo.concealment_method == TCX_NONTONAL
            && st->plcInfo.nbLostCmpt < 4 )
    {
        waveadj_rec = 1;
    }

    /* update main codec parameters */
    st->last_core = st->core;
    st->last_extl = -1;
    st->last_L_frame = L_FRAME;
    st->last_core_brate = st->core_brate;
    st->last_codec_mode = st->codec_mode;
    st->last_bwidth = WB;
    st->last_total_brate_ber = st->total_brate;
    if( !st->bfi )
    {
        st->last_total_brate = st->total_brate;
        st->last_good = st->clas_dec;
    }
    st->last_vad = vad_flag;
    st->last_flag_cna = flag_cna;

    /*----------------------------------------------------------------*
     * Overlap of ACELP synthesis with old MDCT memory
     *----------------------------------------------------------------*/

    if( st->bfi )
    {
        /* calculate another loss frame to fill gap in case of switching frame loss */
        acelp_core_switch_dec_bfi( st, st->fer_samples, coder_type );
    }

    delay_comp = NS2SA(st->output_Fs, DELAY_CLDFB_NS);
    if( last_core_ori == HQ_CORE )
    {
        nZeros = (short)(NS2SA(st->output_Fs,N_ZERO_MDCT_NS));

        if( st->prev_bfi && st->HqVoicing )
        {
            mvr2r( st->fer_samples, st->old_out+nZeros, NS2SA(st->output_Fs,6000000) );
        }

        /* copy the HQ/ACELP delay synchroniation buffer to the beginning of ACELP frame */
        mvr2r( st->delay_buf_out, synth_out, delay_comp );

        tmp = 1.0f/(float)NS2SA(st->output_Fs,6000000);
        for( i=0; i<NS2SA(st->output_Fs,6000000); i++ )
        {
            synth_out[i+delay_comp] = (1-tmp*(float)i)*st->old_out[i+nZeros] + tmp*(float)i*synth_out[i+delay_comp];
        }
    }

    st->prev_bfi = st->bfi;
    st->last_con_tcx = st->con_tcx;

    if( st->core_brate > SID_1k75 )
    {
        st->last_active_brate = st->total_brate;
    }

    if ( st->core_brate > SID_1k75 && st->first_CNG )
    {
        if( st->act_cnt >= BUF_DEC_RATE )
        {
            st->act_cnt = 0;
        }

        st->act_cnt++;

        if( st->act_cnt == BUF_DEC_RATE && st->ho_hist_size > 0 )
        {
            st->ho_hist_size--;
        }

        if( ++(st->act_cnt2) >= MIN_ACT_CNG_UPD )
        {
            st->act_cnt2 = MIN_ACT_CNG_UPD;
        }
    }

    st->prev_bws_cnt = 0;
    st->bws_cnt = 0;
    st->bws_cnt1 = 0;

    /*----------------------------------------------------------------*
     * HP filtering
     * Final synthesis output
     *----------------------------------------------------------------*/

    /* Delay ACELP synthesis by DELAY_BWE_TOTAL_NS - DELAY_CLDFB_NS delay */
    if ( output_frame >= L_FRAME16k )
    {
        tmps = NS2SA(st->output_Fs, DELAY_BWE_TOTAL_NS - DELAY_CLDFB_NS);
        mvr2r( synth_out, tmp_buffer, output_frame );
        mvr2r( st->prev_synth_buffer, synth_out, tmps );
        mvr2r( tmp_buffer, synth_out + tmps, output_frame - tmps );
        mvr2r( tmp_buffer + output_frame - tmps, st->prev_synth_buffer, tmps );
    }

    if (waveadj_rec)
    {
        tmps = 0;
        if( output_frame >= L_FRAME16k )
        {
            tmps = NS2SA(st->output_Fs, DELAY_BWE_TOTAL_NS);
        }
        waveform_adj2(st->tonalMDCTconceal.secondLastPcmOut,
                      synth_out+tmps,
                      st->plcInfo.data_noise,
                      &st->plcInfo.outx_new_n1,
                      &st->plcInfo.nsapp_gain,
                      &st->plcInfo.nsapp_gain_n,
                      &st->plcInfo.recovery_gain, st->plcInfo.step_concealgain,
                      st->plcInfo.Pitch, st->plcInfo.FrameSize,
                      tmps, st->plcInfo.nbLostCmpt + 1, st->bfi );
    }

    /* HP filter */
    hp20( synth_out, output_frame, st->mem_hp20_out, st->output_Fs );

    /* save synthesis for core switching */
    mvr2r( synth_out + NS2SA(st->output_Fs,ACELP_LOOK_NS+DELAY_BWE_TOTAL_NS), st->old_synth_sw, NS2SA(st->output_Fs,FRAME_SIZE_NS-ACELP_LOOK_NS-DELAY_BWE_TOTAL_NS) );

    /* TCX-LTP Postfilter: used in AMR-WB IO to update memories and to avoid discontinuities when the past frame was TCX */
    tcx_ltp_post( st->tcxltp, ACELP_CORE, output_frame, st->L_frame_past, 0, synth_out, NULL, NS2SA( st->output_Fs, TCXLTP_DELAY_NS ),
                  0, 0, 0, 0.f, &st->tcxltp_pitch_int_post_prev, &st->tcxltp_pitch_fr_post_prev, &st->tcxltp_gain_post_prev,
                  &st->tcxltp_filt_idx_prev, st->pit_res_max, &st->pit_res_max_past, 0.f, 0, st->tcxltp_mem_in, st->tcxltp_mem_out, st->total_brate );

    /* final output of synthesis signal */
    mvr2r( synth_out, output, output_frame );




    return;
}
