/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"
#include "cnst.h"

/*-------------------------------------------------------------------*
 * updt_dec()
 *
 * Common updates (all frame types)
 *-------------------------------------------------------------------*/

void updt_dec(
    Decoder_State *st,             /* i/o: state structure                          */
    const short L_frame,           /* i  : length of the frame                      */
    const short coder_type,        /* i  : coding type                              */
    const float *old_exc,          /* i  : buffer of excitation                     */
    const float *pitch_buf,        /* i  : floating pitch values for each subframe  */
    const float Es_pred,           /* i  : predicited scaled innovation energy      */
    const float *Aq,               /* i  : A(z) quantized for all subframes         */
    const float *lsf_new,          /* i  : current frame LSF vector                 */
    const float *lsp_new,          /* i  : current frame LSP vector                 */
    const float voice_factors[],   /* i  : voicing factors                          */
    const float *old_bwe_exc,      /* i  : buffer of excitation                     */
    const float *gain_buf
)
{
    short i;
    short tmp_seed;

    /* update old excitation buffer */
    mvr2r( &old_exc[L_frame], st->old_exc, L_EXC_MEM_DEC);
    if( !st->Opt_AMR_WB )
    {
        mvr2r( &old_bwe_exc[L_FRAME32k], st->old_bwe_exc, PIT16k_MAX * 2 );
    }

    /* update old LSP and LSF vector */
    mvr2r( lsf_new, st->lsf_old, M );
    mvr2r( lsp_new, st->lsp_old, M );

    /* update last coding type */
    st->last_coder_type = coder_type;
    if ( coder_type == INACTIVE || (st->bpf_off == 1 && coder_type != AUDIO && coder_type != TRANSITION) )
    {
        /* overwrite previous coding type to help FEC */
        st->last_coder_type = UNVOICED;
    }

    if( (coder_type != AUDIO || st->Last_GSC_noisy_speech_flag != 0) && st->Last_GSC_pit_band_idx > 0 )
    {
        st->Last_GSC_pit_band_idx = 0;    /*The temporal contribution of the GSC is meaningless after 1 frame lost for inactive & unvoiced content */
    }
    /* this ensures that st->last_coder_type is never set to INACTIVE in case of AVQ inactive because the FEC does not distinguish between GSC inactive and AVQ inactive */

    if ( coder_type == INACTIVE && st->total_brate > ACELP_24k40 )
    {
        st->last_coder_type = GENERIC;
    }

    if( st->Opt_AMR_WB && coder_type == INACTIVE && st->core_brate != SID_1k75 && st->core_brate != FRAME_NO_DATA )
    {
        /* overwrite previous coding type to help FEC */
        st->last_coder_type = UNVOICED;
        st->last_voice_factor = voice_factors[NB_SUBFR-1];
    }

    if( !st->Opt_AMR_WB )
    {
        /* update voicing factor of TBE to help FEC */
        if(st->L_frame == L_FRAME )
        {
            st->last_voice_factor = voice_factors[NB_SUBFR-1];
        }
        else  /* L_frame == L_FRAME16k */
        {
            st->last_voice_factor = voice_factors[NB_SUBFR16k-1];
        }
    }

    if ( coder_type != AUDIO && coder_type != INACTIVE )
    {
        st->noise_lev = NOISE_LEVEL_SP3;
        set_f( st->old_y_gain, 0.0f, MBANDS_GN );

        for( i = 0; i < L_FRAME; i++ )
        {
            tmp_seed = st->seed_tcx;
            st->Last_GSC_spectrum[i] = own_random( &tmp_seed ) / 32768.0f;
        }
    }

    /* update last GSC SWB speech flag for FEC */
    st->Last_GSC_noisy_speech_flag = st->GSC_noisy_speech;

    /* update counter for FEC pitch estimate */
    st->upd_cnt++;
    if( st->upd_cnt > MAX_UPD_CNT )
    {
        st->upd_cnt = MAX_UPD_CNT;
    }

    mvr2r( &st->old_pitch_buf[L_frame/L_SUBFR], st->old_pitch_buf, L_frame/L_SUBFR );
    mvr2r( pitch_buf, &st->old_pitch_buf[L_frame/L_SUBFR], L_frame/L_SUBFR );
    mvr2r( &st->mem_pitch_gain[2], &st->mem_pitch_gain[L_frame/L_SUBFR+2], L_frame/L_SUBFR );

    if( L_frame == L_FRAME )
    {
        st->mem_pitch_gain[2] = gain_buf[3];
        st->mem_pitch_gain[3] = gain_buf[2];
        st->mem_pitch_gain[4] = gain_buf[1];
        st->mem_pitch_gain[5] = gain_buf[0];
    }
    else
    {
        st->mem_pitch_gain[2] = gain_buf[4];
        st->mem_pitch_gain[3] = gain_buf[3];
        st->mem_pitch_gain[4] = gain_buf[2];
        st->mem_pitch_gain[5] = gain_buf[1];
        st->mem_pitch_gain[6] = gain_buf[0];
    }


    /* FEC - update adaptive LSF mean vector */
    mvr2r( st->lsfoldbfi0, st->lsfoldbfi1, M );
    mvr2r( lsf_new, st->lsfoldbfi0, M );

    /* update of pitch and voicing information for HQ FEC */
    if ( st->last_core != HQ_CORE )
    {
        if( !st->Opt_AMR_WB && coder_type == UNVOICED )
        {
            st->HqVoicing = 0;
        }
        else
        {
            st->HqVoicing = 1;
        }
    }

    /* SC-VBR */
    st->old_ppp_mode = st->last_ppp_mode_dec;
    st->last_ppp_mode_dec = st->ppp_mode_dec;
    st->last_nelp_mode_dec = st->nelp_mode_dec;
    st->last_vbr_hw_BWE_disable_dec = st->vbr_hw_BWE_disable_dec;

    /* core switching updates */
    mvr2r( &Aq[(st->L_frame/L_SUBFR-1)*(M+1)], st->old_Aq_12_8, M+1 );
    st->old_Es_pred = Es_pred;

    return;
}

/*-------------------------------------------------------------------*
 * updt_IO_switch()
 *
 * Common updates for AMR-WB IO mode and EVS primary mode switching
 *-------------------------------------------------------------------*/

void updt_IO_switch_dec(
    const short output_frame,   /* i  : output frame length         */
    Decoder_State *st             /* i/o: state structure             */
)
{
    float xsp_tmp[M];

    if( st->last_core == AMR_WB_CORE )      /* switching to EVS primary mode */
    {
        /* AMR-WB IO mode uses ISF(ISP), but EVS primary mode LSF(LSP) */
        mvr2r( stable_LSP, xsp_tmp, M );
        isf2lsf( st->lsf_old, st->lsf_old, xsp_tmp, M, INT_FS_12k8 );
        mvr2r( stable_LSP, xsp_tmp, M );
        isp2lsp( st->lsp_old, st->lsp_old, xsp_tmp, M );

        /* AMR-WB IO mode uses ISF(ISP), but EVS primary mode uses LSF(LSP) */
        mvr2r( stable_LSP, xsp_tmp, M );
        isp2lsp( st->lspCNG, st->lspCNG, xsp_tmp, M );
        st->old_enr_index = min( (short)((float)st->old_enr_index / STEP_AMR_WB_SID * STEP_SID), 127 );

        /* reset TD BWE buffers */
        set_f( st->old_bwe_exc, 0.0f, PIT16k_MAX * 2 );
        set_f( st->old_bwe_exc_extended, 0.0f, NL_BUFF_OFFSET );
        st->bwe_non_lin_prev_scale = 0.0;
        st->last_voice_factor = 0.0f;

        wb_tbe_extras_reset( st->mem_genSHBexc_filt_down_wb2, st->mem_genSHBexc_filt_down_wb3 );
        wb_tbe_extras_reset_synth( st->state_lsyn_filt_shb, st->state_lsyn_filt_dwn_shb, st->mem_resamp_HB );

        if( output_frame >= L_FRAME32k )
        {
            swb_tbe_reset( st->mem_csfilt, st->mem_genSHBexc_filt_down_shb, st->state_lpc_syn,
                           st->syn_overlap, st->state_syn_shbexc, &st->tbe_demph, &st->tbe_premph, st->mem_stp_swb,&(st->gain_prec_swb) );

            /* reset GainShape delay for SWB TBE FEC */
            set_f( st->GainShape_Delay, 0, NUM_SHB_SUBFR/2 );

            swb_tbe_reset_synth( st->genSHBsynth_Hilbert_Mem, st->genSHBsynth_state_lsyn_filt_shb_local );
        }

        if( output_frame == L_FRAME48k )
        {
            st->prev_fb_ener_adjust = 0.0f;
            set_f(st->fb_state_lpc_syn, 0, LPC_SHB_ORDER);
            st->fb_tbe_demph = 0;
            fb_tbe_reset_synth( st->fbbwe_hpf_mem, &st->prev_fbbwe_ratio );
        }

        /* reset FD BWE buffers */
        st->prev_mode = NORMAL;
        st->prev_Energy = 0.0f;
        st->prev_Energy_wb = 0.0f;
        st->prev_L_swb_norm = 8;
        st->prev_frica_flag = 0;
        set_f( st->mem_imdct, 0, L_FRAME48k );
        st->prev_td_energy = 0.0f;
        st->prev_weight = 0.2f;
        set_f( st->old_wtda_swb, 0, L_FRAME48k );

        /* HQ core buffers */
        set_f( st->delay_buf_out, 0, HQ_DELTA_MAX*HQ_DELAY_COMP );

        /* reset the unvoiced/audio signal improvement  memories */
        st->seed_tcx = 15687;
        st->UV_cnt = 30;
        st->LT_UV_cnt = 60.0f;

        st->use_acelp_preq = 0;
        if(st->last_flag_filter_NB == 1)
            st->cldfbSyn->bandsToZero = 0;
        st->last_active_bandsToZero_bwdec = 0;
        st->flag_NB_bwddec = 0;
        st->perc_bwddec = 0.0f;
        st->last_flag_filter_NB = 0;
        st->active_frame_cnt_bwddec = 0;
        set_s(st->flag_buffer, 0, 20);
    }
    else                                    /* switching to AMR-WB IO mode */
    {
        /* ISF Q memories */
        set_f(st->mem_MA, 0, M );

        /* AMR-WB IO mode uses ISF(ISP), but EVS primary mode LSF(LSP) */
        mvr2r( stable_ISP, xsp_tmp, M );
        lsf2isf( st->lsf_old, st->lsf_old, xsp_tmp, M, INT_FS_12k8 );
        mvr2r( stable_ISP, xsp_tmp, M );
        lsp2isp( st->lsp_old, st->lsp_old, xsp_tmp, M );

        /* AMR-WB IO mode uses ISF(ISP), but EVS primary mode LSF(LSP) */
        mvr2r( stable_ISP, xsp_tmp, M );
        lsp2isp( st->lspCNG, st->lspCNG, xsp_tmp, M );
        st->old_enr_index = min( (short)((float)st->old_enr_index / STEP_SID * STEP_AMR_WB_SID), 63 );

        /* gain quantization memory */
        set_f(st->past_qua_en, -14.0f, GAIN_PRED_ORDER );

        /* HF synthesis memories */
        st->ng_ener_ST = -51.0f;

        hf_synth_amr_wb_reset( &st->seed2, st->mem_syn_hf, st->mem_hp_interp, &st->prev_r, &st->fmerit_w_sm, st->delay_syn_hf,
                               &st->frame_count, &st->ne_min, &st->fmerit_m_sm, &st->voice_fac_amr_wb_hf, &st->unvoicing,
                               &st->unvoicing_sm, &st->unvoicing_flag, &st->voicing_flag, &st->start_band_old, &st->OptCrit_old );

        /* reset the unvoiced/audio signal improvement memories */
        st->seed_tcx = 15687;
        st->UV_cnt = 30;
        st->LT_UV_cnt = 60.0f;
        st->Last_ener = 0.0f;
        st->lt_voice_fac = 0.0f;

        st->psf_lp_noise = st->lp_noise;

        /* reset VBR signalling */
        st->last_ppp_mode_dec = 0;
        st->last_nelp_mode_dec = 0;
        st->ppp_mode_dec = 0;
        st->nelp_mode_dec = 0;
    }

    /* CNG - reset */
    st->ho_hist_size = 0;

    /* ISF Q memories */
    mvr2r( UVWB_Ave, st->mem_AR, M );

    /* FEC - update adaptive LSF mean vector */
    mvr2r( st->lsf_old, st->lsfoldbfi0, M );
    mvr2r( st->lsf_old, st->lsfoldbfi1, M );
    mvr2r( st->lsf_old, st->lsf_adaptive_mean, M );

    return;
}

/*-------------------------------------------------------------------*
 * updt_bw_switching()
 *
 * Updates for BW switching
 *-------------------------------------------------------------------*/

void updt_bw_switching(
    Decoder_State *st,                   /* i/o: decoder state structure                  */
    const float *synth,                /* i  : float synthesis signal                   */
    const short *inner_frame_tbl       /* i  : HQ inner_frame signalisation table       */
)
{
    if( st->output_Fs == 32000 && st->bwidth == SWB )
    {
        calc_tilt_bwe( synth, &(st->tilt_swb), L_FRAME32k );
    }

    st->prev_enerLH = st->enerLH;
    st->prev_enerLL = st->enerLL;
    st->last_bwidth = st->bwidth;

    if( st->core == ACELP_CORE )
    {
        if( st->bwidth == WB && st->bws_cnt == 0 )
        {
            st->last_inner_frame = L_FRAME16k;
        }
        else
        {
            st->last_inner_frame = L_FRAME32k;
        }

        if(st->prev_mode == HARMONIC)
        {
            st->prev_weight1 = 0.2f;
        }
        else
        {
            st->prev_weight1 = 0.5f;
        }
    }
    else
    {
        if( st->last_inner_frame >= L_FRAME16k && inner_frame_tbl[st->bwidth] <= L_FRAME16k && st->bws_cnt > 0 && st->bws_cnt < N_WS2N_FRAMES )
        {
            st->last_inner_frame = st->last_inner_frame;
        }
        else
        {
            st->last_inner_frame = inner_frame_tbl[st->bwidth];
        }

        if(inner_frame_tbl[st->bwidth] >= L_FRAME32k || st->core_brate <= HQ_16k40)
        {
            if(st->prev_hqswb_clas == HQ_HARMONIC || st->prev_hqswb_clas == HQ_HVQ)
            {
                st->prev_weight1 = 0.2f;
            }
            else
            {
                st->prev_weight1 = 0.5f;
            }
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * updt_dec_common()
 *
 * Common updates for MODE1 and MODE2
 *-------------------------------------------------------------------*/

void updt_dec_common(
    Decoder_State *st,            /* i/o: decoder state structure     */
    const short hq_core_type,   /* i  : HQ core type                */
    const float *synth          /* i  : decoded synthesis           */
)
{

    st->last_codec_mode = st->codec_mode;
    st->last_extl = st->extl;
    st->last_L_frame = st->L_frame;

    st->prev_old_bfi = st->prev_bfi;
    st->prev_bfi = st->bfi;
    st->old_bfi_cnt = st->nbLostCmpt;
    st->last_con_tcx = st->con_tcx;
    st->con_tcx = 0;

    if( st->use_partial_copy )
    {
        st->prev_rf_frame_type = st->rf_frame_type;
    }
    else
    {
        st->prev_rf_frame_type = INACTIVE;
    }

    st->tcxConceal_recalc_exc = 0;
    if( (st->rf_frame_type >= RF_TCXFD && st->rf_frame_type <= RF_TCXTD2 && st->use_partial_copy && st->bfi) || !st->bfi )
    {
        if( st->bfi && st->last_good <= UNVOICED_TRANSITION && st->clas_dec > UNVOICED_TRANSITION && st->last_con_tcx )
        {
            st->tcxConceal_recalc_exc = 1;
        }
        st->last_good = st->clas_dec;
    }

    if ( st->m_frame_type == ACTIVE_FRAME && ( !st->bfi || st->use_partial_copy ) )
    {
        st->rf_flag_last = st->rf_flag;
    }

    if( st->codec_mode == MODE1 )
    {
        if( !st->bfi && st->core_brate > SID_2k40 )
        {
            st->last_active_brate = st->total_brate;
        }

        st->last_core = st->core;
        st->last_hq_core_type = hq_core_type;
    }
    else if( st->codec_mode == MODE2 )
    {
        if( !st->bfi && st->last_is_cng == 0 )
        {
            st->last_active_brate = st->total_brate;
        }

        if( st->m_frame_type != ACTIVE_FRAME )
        {
            st->last_is_cng = 1;
        }

        if( !st->bfi )
        {
            st->last_core = st->core;
        }
        st->last_core_bfi = st->core; /* also required for clean channel decoding */
    }

    st->last_core_brate = st->core_brate;

    /* save synthesis for core switching */
    mvr2r( synth + NS2SA(st->output_Fs, ACELP_LOOK_NS+DELAY_BWE_TOTAL_NS), st->old_synth_sw, NS2SA(st->output_Fs, FRAME_SIZE_NS-ACELP_LOOK_NS-DELAY_BWE_TOTAL_NS) );

    if( (st->core_brate <= SID_2k40 && st->cng_type == FD_CNG) || (st->tcxonly && st->codec_mode == MODE2) )
    {
        /* reset LP memories */
        set_zero( st->mem_MA, M );
        if(st->sr_core == 16000)
        {
            mvr2r( GEWB2_Ave, st->mem_AR, M );
        }
        else
        {
            mvr2r( GEWB_Ave, st->mem_AR, M );
        }
    }

    return;
}

/*-------------------------------------------------------------------*
* update_decoder_LPD_cng()
*
*
*--------------------------------------------------------------------*/

void update_decoder_LPD_cng(
    Decoder_State *st,
    const short coder_type,
    float *timeDomainBuffer,
    float *A,
    float *bpf_noise_buf
)
{
    short  i;
    float lsp[M], lsf[M], pitch[NB_SUBFR16k];
    float *synth, synth_buf[M+1+L_FRAME_MAX+L_FRAME_MAX/2];
    float tmp;
    float buf_synth[OLD_SYNTH_SIZE_DEC+L_FRAME_MAX+M];
    int   pf_pitch[NB_SUBFR16k];
    float pf_gain[NB_SUBFR16k];

    /* LPC -> LSP/lsp */
    a2lsp_stab( A, lsp, st->lsp_old );

    /* LSP/lsp -> LSF/lsf */
    if( st->L_frame == L_FRAME16k )
    {
        lsp2lsf( lsp, lsf, M, INT_FS_16k );
    }
    else
    {
        lsp2lsf( lsp, lsf, M, INT_FS_12k8 );
    }

    mvr2r( st->old_synth, buf_synth, st->old_synth_len );
    mvr2r( timeDomainBuffer, buf_synth+st->old_synth_len, st->L_frame );

    /* Update synth memory */
    synth = synth_buf + (1+M);
    mvr2r( st->syn, synth_buf, 1+M );
    mvr2r( timeDomainBuffer, synth, st->L_frame );
    mvr2r( synth+st->L_frame-(1+M), st->syn, 1+M );
    mvr2r( st->old_synth+st->L_frame, st->old_synth, st->old_synth_len-st->L_frame );
    mvr2r( synth, st->old_synth+st->old_synth_len-st->L_frame, st->L_frame );

    mvr2r( synth+st->L_frame-(st->L_frame/2), st->old_syn_Overl, st->L_frame/2 );

    st->tcxltp_last_gain_unmodified   = 0.0f;

    /* Update pre-synth memory */
    tmp = synth[-(1+M)];
    preemph( synth-M, st->preemph_fac, M+st->L_frame, &tmp );
    mvr2r( synth+st->L_frame-M, st->mem_syn2, M );
    mvr2r( synth+st->L_frame-L_SYN_MEM, st->mem_syn_r, L_SYN_MEM );

    /* Update excitation memory */
    assert(st->L_frame < L_EXC_MEM_DEC);
    mvr2r( st->old_exc+st->L_frame, st->old_exc, L_EXC_MEM_DEC-st->L_frame );
    residu( A, M,synth, st->old_exc+L_EXC_MEM_DEC-st->L_frame, st->L_frame );

    /* Update LPC-related memories */
    mvr2r( lsp, st->lsp_old, M );
    mvr2r( lsf, st->lsf_old, M );
    mvr2r( lsp, st->lspold_uw, M );
    mvr2r( lsf, st->lsfold_uw, M );

    st->envWeighted = 0;
    mvr2r( A, st->old_Aq_12_8, M+1 );
    st->old_Es_pred = 0;

    /* Reset acelp memories */
    set_zero( st->dispMem, 8 );
    st->tilt_code = TILT_CODE;
    st->gc_threshold = 0.0f;

    /* Update ace/tcx mode */
    st->core = ACELP_CORE;
    st->last_is_cng = 1;

    /* Reset TCX overlap */
    st->tcx_cfg.tcx_curr_overlap_mode = st->tcx_cfg.tcx_last_overlap_mode = ALDO_WINDOW;

    /* For BBWE and Postfilter */
    mvr2r( A, &(st->mem_Aq[0]), M+1 );
    mvr2r( A, &(st->mem_Aq[(M+1)]), M+1 );
    mvr2r( A, &(st->mem_Aq[2*(M+1)]), M+1 );
    mvr2r( A, &(st->mem_Aq[3*(M+1)]), M+1 );
    if( st->L_frame == L_FRAME16k )
    {
        mvr2r( A, &(st->mem_Aq[4*(M+1)]), M+1 );
    }

    /* Update for concealment */
    st->nbLostCmpt = 0;
    st->prev_old_bfi = 0;

    for (i=0; i<M; i++)
    {
        st->lsf_adaptive_mean[i] = ( st->lsfoldbfi1[i]+ st->lsfoldbfi0[i] + lsf[i] )/3;
        st->lsfoldbfi1[i] = st->lsfoldbfi0[i];
        st->lsfoldbfi0[i] = lsf[i];
    }

    set_f( pitch, (float)L_SUBFR, NB_SUBFR16k );

    FEC_clas_estim( synth, pitch, st->L_frame, UNVOICED, st->codec_mode, st->mem_syn_clas_estim, &(st->clas_dec),
                    &st->lp_ener_bfi, st->core_brate, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -1.0f, st->narrowBand, 0,
                    0, st->preemph_fac, st->tcxonly, st->last_core_brate );

    /* Postfiltering */
    pf_pitch[0] = pf_pitch[1] = pf_pitch[2] = pf_pitch[3] = pf_pitch[4] = L_SUBFR;
    pf_gain[0]  = pf_gain[1]  = pf_gain[2]  = pf_gain[3]  = pf_gain[4] = 0.f;
    st->bpf_gain_param = 0;

    post_decoder( st, coder_type, buf_synth, pf_gain, pf_pitch, timeDomainBuffer, bpf_noise_buf );

    return;
}
