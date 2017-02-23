/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*--------------------------------------------------------------------------*
 * evs_dec()
 *
 * Principal decoder routine
 *--------------------------------------------------------------------------*/

void evs_dec(
    Decoder_State *st,                    /* i/o: Decoder state structure   */
    float *output,                /* o  : output synthesis signal   */
    frameMode frameMode               /* i  : Decoder frame mode        */
)
{
    short i, output_frame, coder_type;
    short sharpFlag;
    float synth[L_FRAME48k + HQ_DELTA_MAX*HQ_DELAY_COMP];
    float hb_synth[L_FRAME48k];
    float tmp_buffer[L_FRAME48k];
    short tmps, incr;
    float bwe_exc_extended[L_FRAME32k+NL_BUFF_OFFSET];
    float voice_factors[NB_SUBFR16k];
    float fb_exc[L_FRAME16k];
    short core_switching_flag;
    float old_syn_12k8_16k[L_FRAME16k];
    float tmp;
    float pitch_buf[NB_SUBFR16k];
    short unbits;
    short hq_core_type;
    short post_hq_delay;
    short sid_bw;
    short delay_comp, delta, delay_tdbwe;
    float tmpF;
    short concealWholeFrame;    /* status after decoding */
    short concealWholeFrameTmp;
    float pcmbufFB[L_FRAME_MAX];


    /*------------------------------------------------------------------*
     * Initialization
     *-----------------------------------------------------------------*/

    delay_tdbwe=0;
    sid_bw = -1;
    concealWholeFrameTmp = -1;
    if( !st->bfi )
    {
        st->extl = -1;
    }

    output_frame = (short)(st->output_Fs / 50);

    core_switching_flag = 0;
    sharpFlag = 0;
    unbits = 0;

    st->use_partial_copy = 0;
    st->rf_flag = 0;

    if( st->bfi == 1 )
    {
        hq_core_type = st->last_hq_core_type;
        coder_type = st->last_coder_type;
    }
    else
    {
        hq_core_type = -1;
        coder_type = INACTIVE;
    }

    /* PLC: [TCX: Fade-out-recovery]
     * PLC: overlapping part needs to be attenuated for first good frame */
    if (!st->bfi && st->prev_bfi && (st->last_codec_mode == MODE2) && (st->last_core_bfi == TCX_20_CORE || st->last_core_bfi == TCX_10_CORE))
    {
        v_multc( st->old_out, st->plcInfo.recovery_gain, st->old_out, st->L_frameTCX );
        v_multc( st->old_outLB, st->plcInfo.recovery_gain, st->old_outLB, st->L_frame );

        if( !st->tcx_cfg.last_aldo )
        {
            v_multc( st->syn_OverlFB, st->plcInfo.recovery_gain, st->syn_OverlFB, st->tcx_cfg.tcx_mdct_window_lengthFB );
            v_multc( st->syn_Overl, st->plcInfo.recovery_gain, st->syn_Overl, st->tcx_cfg.tcx_mdct_window_length );
        }
    }

    set_f( voice_factors, 0.f, NB_SUBFR16k );
    set_f( hb_synth, 0.0f, output_frame );

    st->rate_switching_reset = 0;

    /*----------------------------------------------------------------*
     * Updates in case of AMR-WB IO mode -> EVS primary mode switching
     *----------------------------------------------------------------*/

    if( st->last_core == AMR_WB_CORE )
    {
        updt_IO_switch_dec( output_frame, st );
    }

    if( frameMode != FRAMEMODE_MISSING )  /* frame mode normal or future frame */
    {
        getPartialCopyInfo(st, &coder_type, &sharpFlag);

        frameMode = st->bfi;
    }

    if( st->rf_frame_type == RF_NO_DATA && st->use_partial_copy )
    {
        /* the partial copy is a RF FRAME_NO_DATA frame and should follow the concealment path*/
        st->bfi = 1;
        st->codec_mode = st->last_codec_mode;
        frameMode = FRAMEMODE_MISSING;
        st->use_partial_copy = 0;
    }

    /* if previous frame was concealed via ACELP, drop TCX partial copy info and continue ACELP concealment */
    if( st->use_partial_copy && st->core == TCX_20_CORE && st->prev_bfi && st->last_core == ACELP_CORE )
    {
        st->bfi = 1;
        st->codec_mode = st->last_codec_mode;
        frameMode = FRAMEMODE_MISSING;
        st->use_partial_copy = 0;
        st->core = ACELP_CORE;
    }

    /*------------------------------------------------------------------*
     * Decoding
     *-----------------------------------------------------------------*/

    if( st->codec_mode == MODE1 )
    {
        /*------------------------------------------------------------------*
         * Decision matrix (selection of technologies)
         *-----------------------------------------------------------------*/

        /* decision matrix (selection of technologies) */
        if( st->bfi != 1 )
        {
            decision_matrix_dec( st, &coder_type, &sharpFlag, &hq_core_type, &core_switching_flag );

            if( st->bfi != 1 )
            {
                st->sr_core = 50*st->L_frame;
                st->fscale_old  = st->fscale;
                st->fscale      = sr2fscale(st->sr_core);
            }
            else
            {
                frameMode = FRAMEMODE_MISSING;
            }
        }
    }


    if( st->codec_mode == MODE1 )
    {
        /*------------------------------------------------------------------*
         * Initialization
         *-----------------------------------------------------------------*/

        if( st->bfi == 1 )
        {
            st->nbLostCmpt++;
        }
        else
        {
            st->nbLostCmpt = 0;
        }
        st->enablePlcWaveadjust = 0;

        /*---------------------------------------------------------------------*
         * Detect bandwidth switching
         *---------------------------------------------------------------------*/

        bandwidth_switching_detect( st );

        /*---------------------------------------------------------------------*
         * Preprocessing (preparing) for ACELP/HQ core switching
         *---------------------------------------------------------------------*/

        core_switching_pre_dec( st, output_frame );

        /*---------------------------------------------------------------------*
         * ACELP core decoding
         * HQ core decoding
         *---------------------------------------------------------------------*/

        if ( st->core == ACELP_CORE )
        {
            /* ACELP core decoder */
            acelp_core_dec( st, synth, bwe_exc_extended, voice_factors, old_syn_12k8_16k, coder_type, sharpFlag, pitch_buf, &unbits, &sid_bw );
        }
        else
        {
            /* HQ core decoder */
            hq_core_dec( st, synth, output_frame, hq_core_type, core_switching_flag );
        }

        /*---------------------------------------------------------------------*
         * Postprocessing for ACELP/HQ core switching
         *---------------------------------------------------------------------*/

        core_switching_post_dec( st, synth, output_frame, core_switching_flag, coder_type );

        /*---------------------------------------------------------------------*
         * Pre-processing for bandwidth switching
         *---------------------------------------------------------------------*/

        bw_switching_pre_proc( st, old_syn_12k8_16k );

        /*---------------------------------------------------------------------*
         * WB TBE decoding
         * WB BWE decoding
         *---------------------------------------------------------------------*/


        if( st->extl == WB_TBE )
        {
            /* WB TBE decoder */
            wb_tbe_dec( st, coder_type, bwe_exc_extended, voice_factors, hb_synth );
        }

        if( st->extl == WB_BWE && st->bws_cnt == 0 )
        {
            /* WB BWE decoder */
            wb_bwe_dec( synth, hb_synth, output_frame, st, coder_type, voice_factors, pitch_buf );
        }

        /*---------------------------------------------------------------------*
         * SWB(FB) TBE decoding
         * SWB(FB) BWE decoding
         *---------------------------------------------------------------------*/

        if( st->extl == SWB_TBE || st->extl == FB_TBE
                ||( coder_type != AUDIO && coder_type != INACTIVE && st->core_brate >= SID_2k40 && st->core == ACELP_CORE
                    && st->output_Fs >= 32000 && st->bwidth > NB && st->bws_cnt > 0 && !st->ppp_mode_dec && !( st->nelp_mode_dec == 1 && st->bfi == 1 ) ) )

        {
            /* SWB TBE decoder */
            swb_tbe_dec( st, coder_type, bwe_exc_extended, voice_factors,old_syn_12k8_16k, fb_exc, hb_synth, pitch_buf );

            /* FB TBE decoder */
            if( output_frame == L_FRAME48k && st->extl == FB_TBE )
            {
                fb_tbe_dec( st, fb_exc, hb_synth );
            }
        }
        else if( st->extl == SWB_BWE || st->extl == FB_BWE || (st->output_Fs >= 32000 && st->core == ACELP_CORE && st->bwidth > NB
                 && st->bws_cnt > 0 && !st->ppp_mode_dec && !( st->nelp_mode_dec == 1 && st->bfi == 1 ) ) )

        {
            /* SWB BWE decoder */
            swb_bwe_dec( st, synth, hb_synth, output_frame, coder_type );
        }
        else if( st->extl == SWB_BWE_HIGHRATE || st->extl == FB_BWE_HIGHRATE )
        {
            swb_bwe_dec_hr( st, old_syn_12k8_16k, hb_synth, output_frame, unbits, pitch_buf );
        }

        /*---------------------------------------------------------------------*
         * FEC - recovery after lost HQ core (smoothing of the BWE component)
         *---------------------------------------------------------------------*/

        if( st->prev_bfi && st->last_core == HQ_CORE && st->extl != -1 )
        {
            tmp = FRAC_BWE_SMOOTH/output_frame;

            for (i = 0; i < output_frame/FRAC_BWE_SMOOTH; i++)
            {
                hb_synth[i] *= (i*tmp);
            }
        }

        /*---------------------------------------------------------------------*
         * SWB CNG
         *---------------------------------------------------------------------*/

        if( output_frame >= L_FRAME32k )
        {
            /* SHB CNG decoder */
            swb_CNG_dec( st, synth, hb_synth, sid_bw );
        }


        /*----------------------------------------------------------------*
         * Delay ACELP core synthesis to be synchronized with the components of bandwidth extension layers
         *----------------------------------------------------------------*/

        if ( output_frame >= L_FRAME16k )
        {
            tmps = NS2SA(st->output_Fs, DELAY_BWE_TOTAL_NS - DELAY_CLDFB_NS);
            mvr2r( synth, tmp_buffer, output_frame );
            mvr2r( st->prev_synth_buffer, synth, tmps );
            mvr2r( tmp_buffer, synth + tmps, output_frame - tmps );
            mvr2r( tmp_buffer + output_frame - tmps, st->prev_synth_buffer, tmps );
        }

        if( st->core == ACELP_CORE
                && !st->bfi && st->prev_bfi
                && st->last_total_brate >= HQ_48k
                && st->last_codec_mode == MODE2
                && (st->last_core_bfi == TCX_20_CORE || st->last_core_bfi == TCX_10_CORE)
                && st->plcInfo.concealment_method == TCX_NONTONAL
                && st->plcInfo.nbLostCmpt < 4 )
        {
            tmps = 0;
            if( output_frame >= L_FRAME16k )
            {
                tmps = NS2SA(st->output_Fs, DELAY_BWE_TOTAL_NS);
            }

            waveform_adj2( st->tonalMDCTconceal.secondLastPcmOut,
                           synth+tmps,
                           st->plcInfo.data_noise,
                           &st->plcInfo.outx_new_n1,
                           &st->plcInfo.nsapp_gain,
                           &st->plcInfo.nsapp_gain_n,
                           &st->plcInfo.recovery_gain,
                           st->plcInfo.step_concealgain,
                           st->plcInfo.Pitch,
                           st->plcInfo.FrameSize,
                           tmps,
                           st->plcInfo.nbLostCmpt + 1,
                           st->bfi );
            st->plcInfo.Pitch = 0;
        }

        /*----------------------------------------------------------------*
         * Addition of BWE components to the ACELP core synthesis
         *----------------------------------------------------------------*/

        if( st->extl != -1 || (st->bws_cnt > 0 && st->core == ACELP_CORE) )
        {
            /* Calculate an additional delay of extension layer components to be synchronized with ACELP synthesis */
            if (st->L_frame == L_FRAME )
            {
                /* TBE on top of ACELP@12.8kHz */
                tmps = NS2SA( st->output_Fs, MAX_DELAY_TBE_NS - DELAY_SWB_TBE_12k8_NS );
            }
            else
            {
                if( st->extl == SWB_BWE_HIGHRATE || st->extl == FB_BWE_HIGHRATE )
                {
                    /* HR SWB BWE on top of ACELP@16kHz */
                    tmps = NS2SA(st->output_Fs, DELAY_BWE_TOTAL_NS);
                }
                else
                {
                    /* TBE on top of ACELP@16kHz */
                    tmps = NS2SA( st->output_Fs, MAX_DELAY_TBE_NS - DELAY_SWB_TBE_16k_NS );
                }
            }

            /* Smooth transitions when switching between different technologies */
            if( (st->extl != st->last_extl || (st->extl == st->last_extl && (st->core ^ st->last_core) == HQ_CORE)) && !(st->extl == SWB_CNG && st->last_extl == SWB_TBE) )
            {
                /* switching between BWE and TBE technologies */
                incr = (short) ( L_FRAME / (tmps + 0.5f) );
                for (i=0; i<tmps; i++)
                {
                    hb_synth[i] *= sin_table256[i * incr];
                }

                set_f( st->hb_prev_synth_buffer, 0.0f, tmps );
            }
            else if( tmps < st->old_bwe_delay )
            {
                /* the previous frame was TBE on top of ACELP@16kHz and the current frame is TBE on top of ACELP@12.8kHz */
                incr = (short) ( L_FRAME / (tmps + 0.5f) );
                for (i=0; i<tmps; i++)
                {
                    tmp_buffer[i] = st->hb_prev_synth_buffer[i] * sin_table256[255 - i * incr] +
                                    st->hb_prev_synth_buffer[st->old_bwe_delay - 1 - i] * sin_table256[i * incr];
                }

                mvr2r( tmp_buffer, st->hb_prev_synth_buffer, tmps );
            }
            else if( tmps > st->old_bwe_delay )
            {
                /* the previous frame was TBE on top of ACELP@12.8kHz and the current frame is TBE on top of ACELP@16kHz */
                incr = (short)( L_FRAME  / (st->old_bwe_delay + 0.5f) );
                for( i=0; i<st->old_bwe_delay; i++ )
                {
                    tmp_buffer[i] = st->hb_prev_synth_buffer[i] * sin_table256[255 - i * incr];
                }

                for( ; i<tmps; i++ )
                {
                    tmp_buffer[i] = 0.0f;
                }

                for( i=0; i<st->old_bwe_delay; i++ )
                {
                    tmp_buffer[tmps - 1 - i] += st->hb_prev_synth_buffer[st->old_bwe_delay - 1 - i] * sin_table256[i * incr];
                }

                mvr2r( tmp_buffer, st->hb_prev_synth_buffer, tmps );
            }

            /* Delay hb_synth */
            mvr2r( hb_synth, tmp_buffer, output_frame );
            mvr2r( st->hb_prev_synth_buffer, hb_synth, tmps );
            mvr2r( tmp_buffer, hb_synth + tmps, output_frame - tmps );
            mvr2r( tmp_buffer + output_frame - tmps, st->hb_prev_synth_buffer, tmps );

            st->old_bwe_delay = tmps;
            if( ( st->ppp_mode_dec || ( st->nelp_mode_dec == 1 && st->bfi == 1 ) ) && st->L_frame == st->last_L_frame && (st->bws_cnt > 1 || st->last_extl != -1) )

            {
                mvr2r( st->old_hb_synth, hb_synth, output_frame );
            }
            else
            {
                mvr2r(  hb_synth, st->old_hb_synth, output_frame );
            }

            /* Add the delayed hb_synth component to the delayed ACELP synthesis */
            v_add( synth, hb_synth, synth, output_frame );

            /* SWB CNG/DTX - calculate SHB energy */
            if ( output_frame >= L_FRAME32k && st->extl > SWB_CNG )
            {
                st->last_shb_ener = 0.001f;
                for ( i=0; i<output_frame; i++ )
                {
                    st->last_shb_ener += hb_synth[i] * hb_synth[i];
                }
                st->last_shb_ener /= (float)output_frame;
                st->last_shb_ener = 10 * (float)log10(st->last_shb_ener);
            }
        }

        /* TCX-LTP Postfilter: used in Mode 1 to update memories and to avoid discontinuities when the past frame was TCX */
        tcx_ltp_post( st->tcxltp, ACELP_CORE, output_frame, st->L_frame_past, 0, synth, NULL, NS2SA( st->output_Fs, TCXLTP_DELAY_NS ),
                      0, 0, 0, 0.f, &st->tcxltp_pitch_int_post_prev, &st->tcxltp_pitch_fr_post_prev, &st->tcxltp_gain_post_prev,
                      &st->tcxltp_filt_idx_prev, st->pit_res_max, &st->pit_res_max_past, 0.f, 0, st->tcxltp_mem_in, st->tcxltp_mem_out, st->total_brate );

        /* final output of synthesis signal */
        mvr2r( synth, output, output_frame );
    }
    else    /* Mode 2 */
    {

        /* -------------------------------------------------------------- *
         * Mode 2 concealment
         * -------------------------------------------------------------- */

        concealWholeFrame = 0;

        if( frameMode == FRAMEMODE_NORMAL )
        {
            st->m_decodeMode = DEC_NO_FRAM_LOSS;
        }

        if( frameMode == FRAMEMODE_MISSING )
        {
            if( st->use_partial_copy && st->rf_frame_type >= RF_TCXFD && st->rf_frame_type <= RF_TCXTD2 )
            {
                st->m_decodeMode = DEC_NO_FRAM_LOSS;
            }
            else
            {
                st->m_decodeMode = DEC_CONCEALMENT_EXT;
            }
        }

        if( st->m_decodeMode == DEC_CONCEALMENT_EXT )
        {
            concealWholeFrame = 1;
        }

        /* -------------------------------------------------------------- *
         * Decode core
         * -------------------------------------------------------------- */

        dec_acelp_tcx_frame( st, &coder_type, &concealWholeFrame, output,
                             st->p_bpf_noise_buf, pcmbufFB, bwe_exc_extended, voice_factors, pitch_buf );

        concealWholeFrameTmp = concealWholeFrame;
        if(st->bfi)
        {
            frameMode = FRAMEMODE_MISSING;
        }

        if( st->igf )
        {
            /* TBE for Mode 2 interface */
            if( (st->bfi == 0 || st->last_core == ACELP_CORE) && st->core == ACELP_CORE )
            {
                switch( st->bwidth )
                {
                case WB:
                    st->extl = WB_TBE;
                    st->extl_brate = WB_TBE_0k35;
                    break;
                case SWB:
                    st->extl = SWB_TBE;
                    st->extl_brate = SWB_TBE_1k6;
                    break;

                case FB:
                    st->extl = FB_TBE;
                    st->extl_brate = FB_TBE_1k8;
                    break;
                }
            }
            else
            {
                st->extl = IGF_BWE;
                st->extl_brate = 0;
            }

            if( st->output_Fs == 8000 || ( st->output_Fs == 16000 && st->L_frame == L_FRAME16k ) )
            {
                st->extl = -1;
            }

            st->core_brate = st->total_brate - st->extl_brate;

            st->bws_cnt = 0;
            st->bws_cnt1 = 0;
            st->tilt_wb = 0;

            if( st->m_frame_type == ACTIVE_FRAME )
            {
                if( (st->bfi == 0 || st->last_core == ACELP_CORE) && st->core == ACELP_CORE )
                {
                    if( st->extl == WB_TBE )
                    {
                        wb_tbe_dec( st, coder_type, bwe_exc_extended, voice_factors, hb_synth );
                    }
                    else if( st->extl == SWB_TBE || st->extl == FB_TBE )
                    {
                        /* SWB TBE decoder */
                        swb_tbe_dec( st, coder_type, bwe_exc_extended, voice_factors, st->old_core_synth, fb_exc, hb_synth, pitch_buf );

                        if( st->extl == FB_TBE && output_frame == L_FRAME48k )
                        {
                            fb_tbe_dec( st, fb_exc, hb_synth );
                        }
                    }
                    mvr2r( hb_synth, st->old_hb_synth, output_frame );
                }
                else
                {
                    if( st->last_core == ACELP_CORE )
                    {
                        if( ( st->bwidth == SWB || st->bwidth == FB ) &&
                                (( st->last_extl == SWB_TBE || st->last_extl == FB_TBE) && st->last_codec_mode == MODE2 ) )
                        {
                            GenTransition( st->syn_overlap, st->old_tbe_synth, 2*NS2SA(st->output_Fs, DELAY_BWE_TOTAL_NS), hb_synth,
                                           st->genSHBsynth_Hilbert_Mem, st->genSHBsynth_state_lsyn_filt_shb_local, &(st->syn_dm_phase),
                                           st->output_Fs, st->int_3_over_2_tbemem_dec, st->rf_flag, st->total_brate );
                        }
                        else if ( st->bwidth == WB && st->last_extl == WB_TBE )
                        {
                            GenTransition_WB( st->syn_overlap, st->old_tbe_synth, 2*NS2SA(st->output_Fs, DELAY_BWE_TOTAL_NS), hb_synth,
                                              st->state_lsyn_filt_shb, st->state_lsyn_filt_dwn_shb, st->output_Fs, st->mem_resamp_HB );
                        }

                        TBEreset_dec( st, st->bwidth );
                    }
                    else if( st->last_codec_mode == MODE1 )
                    {
                        swb_tbe_reset( st->mem_csfilt, st->mem_genSHBexc_filt_down_shb, st->state_lpc_syn,
                                       st->syn_overlap, st->state_syn_shbexc, &(st->tbe_demph), &(st->tbe_premph), st->mem_stp_swb, &(st->gain_prec_swb) );
                        if( sub(st->extl, FB_TBE) == 0 )
                        {
                            set_f( st->fb_state_lpc_syn, 0, LPC_SHB_ORDER );
                            st->fb_tbe_demph = 0;
                            fb_tbe_reset_synth( st->fbbwe_hpf_mem, &st->prev_fbbwe_ratio );
                        }
                        swb_tbe_reset_synth( st->genSHBsynth_Hilbert_Mem, st->genSHBsynth_state_lsyn_filt_shb_local );
                    }
                }
            }
        }


        if( st->m_frame_type != ACTIVE_FRAME )
        {
            st->extl = -1;
            st->extl_brate = 0;
        }

        /* -------------------------------------------------------------- *
         * Postprocessing
         * -------------------------------------------------------------- */

        {
            float *realBuffer[CLDFB_NO_COL_MAX], *imagBuffer[CLDFB_NO_COL_MAX];
            float realBufferTmp[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX], imagBufferTmp[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX];

            for( i=0; i<CLDFB_NO_COL_MAX; i++ )
            {
                set_f( realBufferTmp[i], 0, CLDFB_NO_CHANNELS_MAX );
                set_f( imagBufferTmp[i], 0, CLDFB_NO_CHANNELS_MAX );
                realBuffer[i] = realBufferTmp[i];
                imagBuffer[i] = imagBufferTmp[i];
            }

            st->cldfbSyn->nab = min( st->cldfbAna->no_channels, st->cldfbSyn->no_channels );
            st->cldfbAna->nab = 0;

            if( st->hFdCngDec != NULL && (st->sr_core == 8000 || st->sr_core == 12800 || st->sr_core == 16000) && st->total_brate <= ACELP_32k )
            {
                /* -------------------------------------------------------------- *
                 * In CLDFB domain:
                 *   - perform noise estimation during active frames
                 *   - do CNG during inactive frames
                 * -------------------------------------------------------------- */

                noisy_speech_detection( st->VAD && st->m_frame_type==ACTIVE_FRAME, output, st->hFdCngDec->hFdCngCom->frameSize,
                                        st->hFdCngDec->msNoiseEst, st->hFdCngDec->psize_shaping, st->hFdCngDec->nFFTpart_shaping,
                                        &(st->hFdCngDec->lp_noise), &(st->hFdCngDec->lp_speech), &(st->hFdCngDec->hFdCngCom->flag_noisy_speech) );

                st->hFdCngDec->hFdCngCom->likelihood_noisy_speech = 0.99f*st->hFdCngDec->hFdCngCom->likelihood_noisy_speech + 0.01f*(float)st->hFdCngDec->hFdCngCom->flag_noisy_speech;

                st->lp_noise = st->hFdCngDec->lp_noise;

                ApplyFdCng( output, realBuffer, imagBuffer, st->hFdCngDec, st->m_frame_type, st, concealWholeFrame, 0 );

                /* Generate additional comfort noise to mask potential coding artefacts */
                if( st->m_frame_type == ACTIVE_FRAME && st->flag_cna )
                {
                    generate_masking_noise( output, st->hFdCngDec->hFdCngCom, st->hFdCngDec->hFdCngCom->frameSize, 0 );
                }
            }

            if( st->flag_cna == 0 && st->L_frame == L_FRAME16k && st->last_flag_cna == 1 && ( (st->last_core == ACELP_CORE && st->last_coder_type != AUDIO) || st->last_core == TCX_20_CORE || st->last_core == AMR_WB_CORE ) )
            {
                v_multc( st->hFdCngDec->hFdCngCom->olapBufferSynth2+5*st->L_frame/4, 256.f, tmp_buffer, st->L_frame/2 );
                v_add( tmp_buffer, output, output, st->L_frame/2 );
            }

            if( st->m_frame_type == ACTIVE_FRAME )
            {
                cldfbAnalysis( output, realBuffer, imagBuffer, -1, st->cldfbAna );
            }
            else
            {
                float timeDomainBuffer[L_FRAME16k];
                float A[M+1];

                mvr2r( st->hFdCngDec->hFdCngCom->timeDomainBuffer, timeDomainBuffer, st->L_frame );
                mvr2r( st->hFdCngDec->hFdCngCom->A_cng, A, M+1 );

                update_decoder_LPD_cng( st, coder_type, timeDomainBuffer, A, st->p_bpf_noise_buf );

                /* Generate additional comfort noise to mask potential coding artefacts */
                if( st->flag_cna )
                {
                    generate_masking_noise( timeDomainBuffer, st->hFdCngDec->hFdCngCom, st->hFdCngDec->hFdCngCom->frameSize, 0 );
                }
                else if( st->L_frame == L_FRAME16k && st->last_flag_cna == 1 && ( (st->last_core == ACELP_CORE && st->last_coder_type != AUDIO) || st->last_core == TCX_20_CORE || st->last_core == AMR_WB_CORE ) )
                {
                    v_multc( st->hFdCngDec->hFdCngCom->olapBufferSynth2+5*st->L_frame/4, 256.f, tmp_buffer, st->L_frame/2 );
                    v_add( tmp_buffer, timeDomainBuffer, timeDomainBuffer, st->L_frame/2 );
                }

                /* check if the CLDFB works on the right sample rate */
                if( (st->cldfbAna->no_channels * st->cldfbAna->no_col) != st->L_frame )
                {
                    resampleCldfb (st->cldfbAna, (st->L_frame * 50));
                    resampleCldfb (st->cldfbBPF, (st->L_frame * 50));

                }

                st->cldfbSyn->bandsToZero = 0;
                if( st->bwidth == NB && st->cldfbSyn->no_channels > 10 )
                {
                    st->cldfbSyn->bandsToZero = st->cldfbSyn->no_channels - 10;
                }
                else if( st->hFdCngDec->hFdCngCom->regularStopBand < st->cldfbSyn->no_channels )
                {
                    st->cldfbSyn->bandsToZero = st->cldfbSyn->no_channels - st->hFdCngDec->hFdCngCom->regularStopBand;
                }
                cldfbAnalysis( timeDomainBuffer, realBuffer, imagBuffer, -1, st->cldfbAna );
            }

            if( st->flag_cna == 0 )
            {
                set_f( st->hFdCngDec->hFdCngCom->olapBufferSynth2, 0.f, st->hFdCngDec->hFdCngCom->fftlen );
            }

            if( st->p_bpf_noise_buf )
            {
                addBassPostFilter( st->p_bpf_noise_buf, -1, realBuffer, imagBuffer, st->cldfbBPF );
            }

            if (st->output_Fs > 8000)
            {
                calcGainTemp_TBE( realBuffer, imagBuffer, st->tecDec.loBuffer, 0,
                                  st->cldfbAna->no_col, st->cldfbAna->no_channels, st->tecDec.pGainTemp, st->tec_flag );
            }

            /* set high band buffers to zero. Covering the current frame and the overlap area. */
            if( st->m_frame_type == ACTIVE_FRAME )
            {
                for( i = 0; i < CLDFB_NO_COL_MAX; i++ )
                {
                    set_f( &realBuffer[i][st->cldfbSyn->nab], 0.f, st->cldfbSyn->no_channels - st->cldfbSyn->nab );
                    set_f( &imagBuffer[i][st->cldfbSyn->nab], 0.f, st->cldfbSyn->no_channels - st->cldfbSyn->nab );
                }
            }

            cldfbSynthesis( realBuffer, imagBuffer, output, -1, st->cldfbSyn );

            /* set multiplication factor according to the sampling rate */
            delay_comp = NS2SA(st->output_Fs, DELAY_CLDFB_NS);
            delay_tdbwe= NS2SA(st->output_Fs, DELAY_BWE_TOTAL_NS- DELAY_CLDFB_NS);

            /* MODE1 MDCT to ACELP 2 transition */
            if( st->last_codec_mode == MODE1 && st->last_core_bfi > ACELP_CORE )
            {
                mvr2r( st->delay_buf_out, output, delay_comp );  /* copy the HQ/ACELP delay synchronization buffer at the beginning of ACELP frame */

                if( st->core == ACELP_CORE )
                {
                    tmpF = 1.0f/(float)NS2SA(st->output_Fs, 3000000);
                    if(st->prev_bfi && st->HqVoicing )
                    {
                        mvr2r( st->fer_samples, st->old_out+NS2SA(st->output_Fs, N_ZERO_MDCT_NS), NS2SA(st->output_Fs,3000000) );
                    }
                    for( i=0; i<NS2SA(st->output_Fs, 3000000); i++ )
                    {
                        output[i+delay_comp] = (1-tmpF*(float)i)*st->old_out[i+NS2SA(st->output_Fs, N_ZERO_MDCT_NS)] + tmpF*(float)i*output[i+delay_comp];
                    }
                }
                else
                {
                    if( st->output_Fs == 8000 )
                    {
                        mvr2r( st->delay_buf_out, st->FBTCXdelayBuf, delay_comp );
                    }
                    else
                    {
                        mvr2r( st->prev_synth_buffer, st->FBTCXdelayBuf, delay_tdbwe );
                        mvr2r( st->delay_buf_out, st->FBTCXdelayBuf + delay_tdbwe, delay_comp );
                    }
                }
            }

            /* set delay compensation between HQ synthesis and ACELP synthesis */
            if( st->core == ACELP_CORE && !st->con_tcx )
            {
                set_f( st->delay_buf_out, 0, delay_comp );
                mvr2r( output, st->previoussynth, output_frame );
            }
            else
            {
                mvr2r( st->old_synthFB+st->old_synth_lenFB-delay_comp, st->delay_buf_out, delay_comp );

                if( st->output_Fs == 8000 )
                {
                    mvr2r( st->FBTCXdelayBuf, st->previoussynth, delay_comp );
                }
                else
                {
                    mvr2r( st->FBTCXdelayBuf + delay_tdbwe, st->previoussynth, delay_comp );
                }

                mvr2r( pcmbufFB, st->previoussynth + delay_comp, output_frame - delay_comp );
            }

            /* Delay compensation for TBE */
            if( output_frame >= L_FRAME16k )
            {
                mvr2r( output, tmp_buffer, output_frame );
                mvr2r( st->prev_synth_buffer, output, delay_tdbwe );
                mvr2r( tmp_buffer, output + delay_tdbwe, output_frame - delay_tdbwe );
                mvr2r( tmp_buffer + output_frame - delay_tdbwe, st->prev_synth_buffer, delay_tdbwe );
            }

            if( st->igf && st->m_frame_type == ACTIVE_FRAME )
            {
                if( !st->bfi && st->core == ACELP_CORE && (st->tec_flag || st->tfa_flag) && (st->output_Fs > 8000) )
                {
                    procTecTfa_TBE( hb_synth, st->tecDec.pGainTemp, st->tfa_flag, st->last_core, (int)(output_frame / N_TEC_TFA_SUBFR), st->tec_flag == 2 ? 1 : 0 );
                }

                if( (((!st->bfi || st->last_core == ACELP_CORE) && st->core == ACELP_CORE) ||
                        (st->last_core == ACELP_CORE && st->bwidth != NB && st->last_codec_mode == MODE2))
                        && (st->output_Fs > 8000) )
                {
                    /* Add the delayed hb_synth component to the delayed core synthesis */
                    v_add( output, hb_synth, output, output_frame );
                }
            }
        }

        /* set delay */
        if( st->output_Fs == 8000 )
        {
            tmps = NS2SA( st->output_Fs, DELAY_CLDFB_NS );
        }
        else
        {
            tmps = NS2SA( st->output_Fs, DELAY_BWE_TOTAL_NS );
        }
        delta = NS2SA( st->output_Fs, TCXLTP_DELAY_NS );

        /* TCX/ACELP/HQ-CORE->TCX */
        if( (st->bfi && st->last_core > ACELP_CORE) || st->core > ACELP_CORE )
        {
            /* TCX / HQ-CORE / TD-TCX-PLC -> TCX / TD-TCX-PLC */
            if( st->last_core_bfi > ACELP_CORE || (st->bfi && st->last_core > ACELP_CORE) || (st->prev_bfi && st->last_con_tcx))
            {
                mvr2r( st->FBTCXdelayBuf, output, tmps );
                mvr2r( pcmbufFB, output + tmps, st->L_frameTCX - tmps );
            }
            /* ACELP -> TCX */
            else
            {
                /*cross-fading between LB-TCX and FB-TCX over 2.3125ms*/
                for( i = 0; i < tmps; i++ )
                {
                    output[i+tmps] = (output[i+tmps] * (tmps-i) + pcmbufFB[i] * i) / tmps;
                }
                mvr2r( pcmbufFB+tmps, output + 2*tmps, st->L_frameTCX - 2*tmps );
            }

            mvr2r( pcmbufFB + st->L_frameTCX - tmps, st->FBTCXdelayBuf, tmps );

            if( st->bfi && st->last_core > ACELP_CORE )
            {
                if( st->output_Fs == 8000 )
                {
                    mvr2r(st->FBTCXdelayBuf, st->delay_buf_out, NS2SA(st->output_Fs, DELAY_CLDFB_NS));
                }
                else
                {
                    mvr2r( st->FBTCXdelayBuf, st->prev_synth_buffer, NS2SA(st->output_Fs, DELAY_BWE_TOTAL_NS - DELAY_CLDFB_NS) );
                    mvr2r( st->FBTCXdelayBuf + NS2SA(st->output_Fs, DELAY_BWE_TOTAL_NS - DELAY_CLDFB_NS), st->delay_buf_out, NS2SA(st->output_Fs, DELAY_CLDFB_NS) );
                }
            }
        }
        /* TCX/TD TCX PLC->ACELP */
        else if( st->last_codec_mode == MODE2 && st->last_core > ACELP_CORE )
        {
            mvr2r( st->FBTCXdelayBuf, output, delta );
            for( i = delta; i < tmps; i++ )
            {
                output[i] = (output[i] * (i-delta) + st->FBTCXdelayBuf[i] * (tmps-i)) / (tmps-delta);
            }
        }

        tcx_ltp_post( st->tcxltp, st->core, output_frame, st->L_frame, NS2SA( st->output_Fs, ACELP_LOOK_NS ) + tmps, output,
                      st->FBTCXdelayBuf, delta, st->bfi, st->tcxltp_pitch_int, st->tcxltp_pitch_fr, st->tcxltp_gain,
                      &st->tcxltp_pitch_int_post_prev, &st->tcxltp_pitch_fr_post_prev, &st->tcxltp_gain_post_prev, &st->tcxltp_filt_idx_prev,
                      st->pit_res_max, &st->pit_res_max_past, st->damping, st->total_brate >= HQ_96k, st->tcxltp_mem_in, st->tcxltp_mem_out, st->total_brate );
    }   /* end of Mode 2 */

    /*----------------------------------------------------------------*
     * Save synthesis for HQ FEC
     *----------------------------------------------------------------*/

    post_hq_delay = NS2SA( st->output_Fs, POST_HQ_DELAY_NS );
    if( st->codec_mode == MODE1 )
    {
        mvr2r( st->synth_history+output_frame, st->synth_history, output_frame-post_hq_delay+NS2SA( st->output_Fs, PH_ECU_MEM_NS ));
        mvr2r( output, st->old_synthFB+output_frame-post_hq_delay, output_frame );
        /* reset the remaining buffer, which is read in TCX concealment the necessary samples to fill
           this buffer are not available for all cases, the impact on the output is limited */
        set_f( st->old_synthFB+2*output_frame-post_hq_delay, 0.f, post_hq_delay );

        if ( output_frame >= L_FRAME16k )
        {
            mvr2r( st->prev_synth_buffer, st->old_synthFB+2*output_frame-NS2SA(st->output_Fs, DELAY_BWE_TOTAL_NS), NS2SA(st->output_Fs, DELAY_BWE_TOTAL_NS - DELAY_CLDFB_NS));
        }

        if (st->core != ACELP_CORE)
        {
            if ( output_frame >= L_FRAME16k )
            {
                mvr2r( synth+output_frame, st->old_synthFB+2*output_frame-NS2SA(st->output_Fs, DELAY_CLDFB_NS), NS2SA(st->output_Fs, DELAY_CLDFB_NS));
                mvr2r( st->old_out+NS2SA(st->output_Fs, N_ZERO_MDCT_NS), st->old_synthFB+2*output_frame, NS2SA(st->output_Fs, PH_ECU_LOOKAHEAD_NS));
            }
            else
            {
                mvr2r( synth+output_frame, st->old_synthFB+2*output_frame-NS2SA(st->output_Fs, DELAY_BWE_TOTAL_NS), NS2SA(st->output_Fs, DELAY_CLDFB_NS));
                mvr2r( st->old_out+NS2SA(st->output_Fs, N_ZERO_MDCT_NS), st->old_synthFB+2*output_frame-NS2SA(st->output_Fs, DELAY_BWE_TOTAL_NS - DELAY_CLDFB_NS), NS2SA(st->output_Fs, PH_ECU_LOOKAHEAD_NS));
            }
        }
    }

    /*----------------------------------------------------------------*
     * HP filtering
     *----------------------------------------------------------------*/

    hp20( output, output_frame, st->mem_hp20_out, st->output_Fs );


    /*--------------------------------------------------------*
     * Updates
     *--------------------------------------------------------*/

    if( st->last_is_cng == 0 && st->codec_mode == MODE2 )
    {
        st->bfi = 0;
        if( st->use_partial_copy && st->rf_frame_type >= RF_TCXFD && st->rf_frame_type <= RF_TCXTD2 )
        {
            if( frameMode == FRAMEMODE_MISSING )
            {
                st->bfi = 1;
            }
        }
        else if( st->m_decodeMode == DEC_CONCEALMENT_EXT )
        {
            st->bfi = 1;
        }

        updt_dec_common( st, -1, output );
    }
    else
    {
        if( st->codec_mode == MODE2 )
        {
            st->bfi = 0;
        }
        updt_dec_common( st, hq_core_type, output );
    }

    if( st->codec_mode == MODE2 )
    {
        if( st->use_partial_copy && st->rf_frame_type == RF_NELP )
        {
            st->last_nelp_mode_dec = 1;
        }
        else
        {
            st->last_nelp_mode_dec = 0;
        }
    }

    st->prev_use_partial_copy = st->use_partial_copy;

    st->prev_tilt_code_dec = 0.0f;
    for( i=0; i<NB_SUBFR; i++ )
    {
        st->prev_tilt_code_dec += st->tilt_code_dec[i]*0.25f;
    }

    if( st->core == HQ_CORE )
    {
        st->prev_coder_type = GENERIC;
    }
    else
    {
        st->prev_coder_type = coder_type;
    }

    if( st->core_brate > SID_2k40 && st->first_CNG == 1 )
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

    if( st->core_brate <= SID_2k40 && st->first_CNG == 0 && st->cng_type == LP_CNG )
    {
        st->first_CNG = 1;
    }

    /* update bandwidth switching parameters */
    if( st->codec_mode == MODE1 )
    {
        updt_bw_switching( st, output, inner_frame_tbl );
    }
    else
    {
        st->last_bwidth = st->bwidth;
    }

    /* synchronisation of CNG seeds */
    if( st->bfi || (st->core_brate != FRAME_NO_DATA && st->core_brate != SID_2k40) )
    {
        own_random( &(st->cng_seed) );
        own_random( &(st->cng_ener_seed) );
    }

    if( st->enablePlcWaveadjust && !concealWholeFrameTmp )
    {
        /* update the parameters used in waveform adjustment */
        concealment_update2( output, &st->plcInfo, st->L_frameTCX );
    }

    st->last_total_brate_ber = st->total_brate;
    if( !st->bfi )
    {
        st->last_total_brate = st->total_brate;
    }

    st->last_flag_cna = st->flag_cna;
    st->hFdCngDec->hFdCngCom->frame_type_previous = st->m_frame_type;

    st->prev_last_core = st->last_core;
    st->prev_bws_cnt = st->bws_cnt;



    return;
}
