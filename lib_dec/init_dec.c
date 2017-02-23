/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"
#include <string.h>


/*----------------------------------------------------------------------*
 * init_decoder()
 *
 * Initialization of static variables for the decoder
 *----------------------------------------------------------------------*/

void init_decoder(
    Decoder_State *st  /* o:   Decoder static variables structure */
)
{
    short i, j;

    /*-----------------------------------------------------------------*
     * ACELP core parameters
     *-----------------------------------------------------------------*/

    st->codec_mode = MODE1;
    st->last_codec_mode = MODE1;
    st->core = ACELP_CORE;
    st->L_frame = L_FRAME;
    st->extl = -1;
    st->total_brate = ACELP_8k00;
    st->last_total_brate = -1;
    st->last_total_brate_ber = -1;
    st->core_brate = ACELP_8k00;
    st->ini_frame = 0;
    st->bwidth = NB;
    st->extl_brate = 0;

    st->last_coder_type = GENERIC;
    st->last_L_frame = st->L_frame;
    st->last_core_brate = st->core_brate;
    st->last_core = -1;
    st->prev_last_core = -1;
    st->last_hq_core_type = -1;
    st->last_extl = st->extl;

    /* LSF initilaizations */
    mvr2r( GEWB_Ave, st->mem_AR, M );
    init_lvq( st->offset_scale1, st->offset_scale2, st->offset_scale1_p, st->offset_scale2_p, st->no_scales, st->no_scales_p );

    set_f( st->mem_MA, 0, M );

    set_f( st->dispMem, 0, 8 );

    /* AMR-WB IO HF synth init */
    hf_synth_amr_wb_init( &st->prev_r, &st->fmerit_w_sm, st->mem_syn_hf, &st->frame_count, &st->ne_min, &st->fmerit_m_sm, &st->voice_fac_amr_wb_hf,
                          &st->unvoicing, &st->unvoicing_sm, &st->unvoicing_flag, &st->voicing_flag, &st->start_band_old, &st->OptCrit_old );

    hf_synth_init( st->mem_hp400, st->mem_hf );
    set_f( st->mem_hp_interp, 0, INTERP_3_1_MEM_LEN );
    set_f( st->delay_syn_hf, 0, NS2SA(16000,DELAY_CLDFB_NS) );

    st->tilt_code = 0.0f;
    st->gc_threshold = 0.0f;
    st->last_good = UNVOICED_CLAS;
    st->clas_dec = UNVOICED_CLAS;

    st->lp_gainp = 0.0f;
    st->lp_gainc = 0.0f;

    set_f( st->old_exc, 0, L_EXC_MEM_DEC );

    /* AVQ pre-quantizer memory */
    st->mem_preemp_preQ = 0.0f;
    st->last_nq_preQ = 0;
    st->use_acelp_preq = 0;

    st->mem_deemph = 0.0f;

    set_f( st->mem_syn1, 0, M );
    st->mem_deemph_old_syn = 0.0f;
    set_f( st->mem_syn2, 0, M );
    st->stab_fac = 0.0f;
    st->stab_fac_smooth = 0.0f;
    set_f( st->agc_mem2, 0, 2 );
    set_f( st->mem_hp20_out, 0.0f, 4 );
    set_f( st->mem_syn3, 0, M );

    for (i=0; i<GAIN_PRED_ORDER; i++)
    {
        st->past_qua_en[i] = -14.0f;   /* gain quantization memory (used in AMR-WB IO mode) */
    }

    mvr2r( GEWB_Ave, st->lsf_old, M );
    lsf2lsp( st->lsf_old, st->lsp_old, M, INT_FS_12k8 );

    st->mid_lsf_int = 0;
    st->safety_net = 0;

    /* parameters for AC mode (GSC) */
    st->seed_tcx = 15687;
    st->GSC_noisy_speech = 0;
    st->Last_GSC_noisy_speech_flag = 0;
    st->cor_strong_limit = 1;
    set_f( st->old_y_gain, 0, MBANDS_GN );
    st->noise_lev = NOISE_LEVEL_SP0;
    set_f( st->Last_GSC_spectrum, 0.0f, L_FRAME );
    st->Last_GSC_pit_band_idx = 0;

    set_f( st->lt_ener_per_band, 1.0f, MBANDS_GN );
    set_f( st->last_exc_dct_in, 0, L_FRAME );
    st->last_ener = 0.0f;
    set_s( st->last_bitallocation_band, 0, 6 );

    /* NB post-filter */
    Init_post_filter( &(st->pfstat) );
    st->psf_lp_noise = 0.0f;

    /* FEC */
    st->scaling_flag = 0;
    st->lp_ener_FEC_av = 5.0e5f;
    st->lp_ener_FEC_max = 5.0e5f;
    st->prev_bfi = 0;
    st->lp_ener_bfi = 60.0f;
    st->old_enr_LP = 0.0f;
    st->lp_ener = 0.0f;
    st->enr_old = 0.0f;
    st->bfi_pitch = (float)L_SUBFR;
    st->bfi_pitch_frame = L_FRAME;
    set_f( st->mem_syn_clas_estim, 0.0f, L_SYN_MEM_CLAS_ESTIM );
    st->last_con_tcx = 0;

    for (i=0; i<2*NB_SUBFR16k; i++)
    {
        st->old_pitch_buf[i] = (float)L_SUBFR;
    }

    st->upd_cnt = MAX_UPD_CNT;

    mvr2r( GEWB_Ave, st->lsfoldbfi0, M );
    mvr2r( GEWB_Ave, st->lsfoldbfi1, M );
    mvr2r( GEWB_Ave, st->lsf_adaptive_mean, M );

    st->seed_acelp = RANDOM_INITSEED;
    st->seed = RANDOM_INITSEED;
    st->nbLostCmpt = 1;
    st->decision_hyst = 0;

    /* fast recovery */
    set_f( st->old_exc2, 0, L_EXC_MEM );
    set_f( st->old_syn2, 0, L_EXC_MEM );

    /* Stationary noise UV modification */
    st->unv_cnt = 0;
    st->ge_sm = 10;
    st->uv_count = 0;
    st->act_count = 3;
    mvr2r( st->lsp_old, st->lspold_s, M );
    st->noimix_seed = RANDOM_INITSEED;
    st->min_alpha = 1;
    st->exc_pe = 0;

    /*-----------------------------------------------------------------*
     * LD music post-filter
     *-----------------------------------------------------------------*/

    set_f( st->dct_post_old_exc, 0, DCT_L_POST-OFFSET2 );
    st->LDm_enh_min_ns_gain = (float)pow(10.0f, -12/20.0f);
    st->LDm_last_music_flag = 0;
    set_f(st->LDm_lt_diff_etot, 0, MAX_LT);
    st->LDm_thres[0] = TH_0_MIN;
    st->LDm_thres[1] = TH_1_MIN;
    st->LDm_thres[2] = TH_2_MIN;
    st->LDm_thres[3] = TH_3_MIN;
    st->LDm_nb_thr_1 = 0;
    st->LDm_nb_thr_3 = 0;
    st->LDm_mem_etot = 0.0f;

    for (i = 0; i < VOIC_BINS_HR; i++)
    {
        st->LDm_enh_lp_gbin[i] = 1.0f;
        st->LDm_enh_lf_EO [i] = 0.01f;
    }

    for (i = 0; i < MBANDS_GN_LD; i++)
    {
        st->LDm_bckr_noise[i] = E_MIN;
    }

    set_f(st->filt_lfE, 1.0f, DCT_L_POST);
    st->last_nonfull_music = 0;

    /*-----------------------------------------------------------------*
     * CNG and DTX
     *-----------------------------------------------------------------*/
    st->CNG = 0;                                       /* RTXDTX handler CNG=1  nonCNG= 0,*/
    st->prev_ft_speech = 1;                             /* RXDTX handeler  previous frametype flag for  G.192 format AMRWB  SID_FIRST  detection  */
    st->cng_seed = RANDOM_INITSEED;
    st->cng_ener_seed = RANDOM_INITSEED;
    st->cng_ener_seed1 = RANDOM_INITSEED;
    st->old_enr_index = -1;
    st->Enew = 0.0f;
    st->first_CNG = 0;
    mvr2r( st->lsp_old, st->lspCNG, M );
    st->last_allow_cn_step = 0;
    st->shb_cng_ener = -6.02f;
    st->wb_cng_ener = -6.02f;
    st->last_wb_cng_ener = -6.02f;
    st->last_shb_cng_ener = -6.02f;
    st->swb_cng_seed = RANDOM_INITSEED;
    st->ho_hist_ptr = -1;
    st->ho_sid_bw = 0;
    set_f( st->ho_lsp_hist, 0, HO_HIST_SIZE*M );
    set_f( st->ho_ener_hist, 0, HO_HIST_SIZE );
    set_f( st->ho_env_hist, 0, HO_HIST_SIZE*NUM_ENV_CNG );
    st->ho_hist_size = 0;
    st->act_cnt = 0;
    st->ho_circ_ptr = -1;
    set_f( st->ho_lsp_circ, 0, HO_HIST_SIZE*M );
    set_f( st->ho_ener_circ, 0, HO_HIST_SIZE );
    set_f( st->ho_env_circ, 0, HO_HIST_SIZE*NUM_ENV_CNG );
    st->ho_circ_size = 0;

    set_s( st->ho_16k_lsp, 0, HO_HIST_SIZE );
    st->CNG_mode = -1;
    st->last_active_brate = ACELP_7k20;
    st->last_CNG_L_frame = L_FRAME;
    st->act_cnt2 = 0;
    st->num_ho = 0;
    st->cng_type = -1;
    st->last_cng_type = -1;
    set_f( st->lp_env, 0.0f, NUM_ENV_CNG );
    set_f( st->exc_mem, 0.0f, 24 );
    set_f( st->exc_mem1, 0.0f, 30 );
    set_f( st->old_env, 0.0f, NUM_ENV_CNG );

    for ( i=0; i<LPC_SHB_ORDER; i++ )
    {
        st->lsp_shb_prev[i] = 0.5f * ((float) i)/((float) LPC_SHB_ORDER);
        st->lsp_shb_prev_prev[i] = st->lsp_shb_prev[i];
    }

    st->shb_dtx_count = 0;
    st->last_vad = 0;
    st->trans_cnt = 0;
    st->burst_cnt = 0;
    st->last_shb_ener = 0.001f;

    /* HF (6-7kHz) BWE */
    st->seed2 = RANDOM_INITSEED;

    /*-----------------------------------------------------------------*
     * HR SWB BWE parameters
     *-----------------------------------------------------------------*/

    set_f( st->t_audio_prev, 0, 2*END_FREQ_BWE_FULL_FB/50 - NUM_NONTRANS_START_FREQ_COEF );
    st->old_is_transient_hr_bwe = 0;
    st->bwe_highrate_seed = 12345;
    st->mem_EnergyLT = 0.0f;

    /*-----------------------------------------------------------------*
     * HQ core parameters
     *-----------------------------------------------------------------*/

    set_f( st->old_out, 0, L_FRAME48k );
    set_f( st->old_outLB, 0, L_FRAME32k );
    set_f( st->old_coeffs, 0, L_FRAME8k );
    set_s(st->old_is_transient, 0, 3);
    st->old_bfi_cnt = 0;
    set_f( st->old_auOut_2fr, 0, L_FRAME8k*2 );
    set_f( st->old_out_pha[0], 0, N_LEAD_NB );
    set_f( st->old_out_pha[1], 0, N_LEAD_NB );
    st->prev_old_bfi = 0;
    st->phase_mat_flag = 0;
    st->phase_mat_next = 0;
    st->old_Min_ind = 0;
    st->diff_energy = 0.f;
    set_f( st->oldIMDCTout, 0.f, L_FRAME8k/2 );
    set_f( st->prev_oldauOut, 0.f, L_FRAME8k);
    st->stat_mode_out = 0;
    st->stat_mode_old = 0;
    st->oldHqVoicing=0;

    for( i=0; i<MAX_SB_NB; i++ )
    {
        for( j=0; j<MAX_PGF; j++ )
        {
            st->ynrm_values[i][j] = 0.f;
        }
        for( j=0; j<MAX_ROW; j++ )
        {
            st->r_p_values[i][j] = 0.f;
        }
    }
    set_f( st->Norm_gain, 1.f, SFM_N_NB );
    set_f( st->energy_MA_Curr, 100.f, 2 );
    st->HQ_FEC_seed = RANDOM_INITSEED;
    set_f( st->delay_buf_out, 0, HQ_DELTA_MAX*HQ_DELAY_COMP );
    set_f( st->previoussynth, 0, L_FRAME48k);
    set_f( st->old_synth_sw, 0.0f, NS2SA(48000,FRAME_SIZE_NS-ACELP_LOOK_NS-DELAY_BWE_TOTAL_NS) );
    set_f( st->prev_noise_level, 0.0f, 2 );
    st->prev_R = 0;
    set_f( st->prev_coeff_out, 0, L_HQ_WB_BWE );
    set_s( st->prev_SWB_peak_pos, 0, SPT_SHORTEN_SBNUM );

    /* HQ GENERIC */
    st->hq_generic_seed = RANDOM_INITSEED;

    st->mem_norm[0] = 31;
    set_s( st->mem_norm+1, 39, SFM_N_ENV_STAB-1 );
    st->mem_env_delta = 0;
    st->no_att_hangover = 0;
    st->energy_lt = 300.0f;

    st->HqVoicing = 0;
    set_f( st->fer_samples, 0, L_FRAME48k );
    set_f( st->prev_env, 0, SFM_N_WB );
    set_f( st->prev_normq, 0, SFM_N_WB );
    st->prev_hqswb_clas = HQ_NORMAL;

    set_f( st->last_ni_gain, 0, BANDS_MAX );
    set_f( st->last_env, 0, BANDS_MAX );
    st->last_max_pos_pulse = 0;
    st->prev_frm_hfe2 = 0;
    st->prev_stab_hfe2 = 0;
    st->prev_ni_ratio = 0.5f;
    set_f(st->prev_En_sb, 0.0f, NB_SWB_SUBBANDS );

    /* pre-echo reduction */
    reset_preecho_dec( st );

    /*----------------------------------------------------------------------------------*
     * HQ FEC
     *----------------------------------------------------------------------------------*/

    st->old_synthFB = st->synth_history + NS2SA(st->output_Fs, PH_ECU_MEM_NS);
    st->prev_good_synth = st->old_synthFB + NS2SA(st->output_Fs, PH_ECU_LOOKAHEAD_NS);

    set_f( st->X_sav, 0.0f, PH_ECU_SPEC_SIZE );
    st->num_p = 0;
    st->ph_ecu_active = 0;
    st->ni_seed_forfec = 0;
    st->last_fec = 0;
    st->ph_ecu_HqVoicing = 0;
    set_f( st->oldgapsynth, 0.0f, L_FRAME48k );
    st->env_stab = 0.75f;
    st->mem_norm_hqfec[0] = 31;
    set_s( st->mem_norm_hqfec+1, 39, SFM_N_ENV_STAB-1 );
    st->mem_env_delta_hqfec = 0;
    st->env_stab_plc = 0.0f;
    set_f( st->env_stab_state_p, 1.0f/NUM_ENV_STAB_PLC_STATES, NUM_ENV_STAB_PLC_STATES );
    st->envstabplc_hocnt = 0;

    set_f( st->mag_chg_1st, 1.0f, LGW_MAX );
    set_f( st->Xavg, 0.0f, LGW_MAX );
    st->beta_mute = BETA_MUTE_FAC_INI;

    set_s( st->prev_sign_switch, 0, HQ_FEC_SIGN_SFM );
    set_s( st->prev_sign_switch_2, 0, HQ_FEC_SIGN_SFM );

    st->time_offs = 0;
    st->ber_occured_in_pvq = 0;

    /*-----------------------------------------------------------------*
     * SWB BWE parameters
     *-----------------------------------------------------------------*/

    set_f( st->old_wtda_swb, 0, L_FRAME48k);
    set_f( st->old_syn_12k8_16k, 0, NS2SA(16000, DELAY_FD_BWE_ENC_NS) );
    st->prev_mode = NORMAL;
    set_f( st->prev_SWB_fenv, 0, SWB_FENV );
    st->prev_Energy = 0.0f;
    st->prev_L_swb_norm = 8;
    st->Seed = 21211;
    st->prev_frica_flag = 0;
    set_f( st->mem_imdct, 0, L_FRAME48k );
    st->prev_td_energy = 0.0f;
    st->prev_weight = 0.2f;
    st->prev_flag = 0;
    st->prev_coder_type = GENERIC;
    st->last_wb_bwe_ener = 0.0f;
    st->tilt_wb = 0.0f;
    st->prev_Energy_wb = 0.0f;

    /*-----------------------------------------------------------------*
     * TBE parameters
     *-----------------------------------------------------------------*/

    InitSWBdecBuffer( st );
    ResetSHBbuffer_Dec( st );

    if( st->output_Fs == 48000 )
    {
        set_f( st->fbbwe_hpf_mem[0], 0, 4 );
        set_f( st->fbbwe_hpf_mem[1], 0, 4 );
        set_f( st->fbbwe_hpf_mem[2], 0, 4 );
        set_f( st->fbbwe_hpf_mem[3], 0, 4 );
    }

    set_f( st->mem_resamp_HB, 0, INTERP_3_1_MEM_LEN );
    set_f( st->mem_resamp_HB_32k, 0, 2*ALLPASSSECTIONS_STEEP+1 );
    set_f( st->prev_synth_buffer, 0, NS2SA(48000, DELAY_BWE_TOTAL_NS - DELAY_CLDFB_NS) );
    set_f( st->hb_prev_synth_buffer, 0, NS2SA(48000, DELAY_BWE_TOTAL_NS) );
    st->old_bwe_delay = -1;

    st->tilt_mem = 0.0f;
    set_f( st->prev_lsf_diff, 0.5f, LPC_SHB_ORDER-2 );
    st->prev_tilt_para = 0.0f;
    set_f( st->cur_sub_Aq, 0.0f, M+1 );
    set_f( st->int_3_over_2_tbemem_dec, 0.0f, INTERP_3_2_MEM_LEN );
    set_f( st->interpol_3_2_cng_dec, 0.0f, INTERP_3_2_MEM_LEN );

    /* TD BWE post-processing */
    st->ptr_mem_stp_swb = st->mem_stp_swb + LPC_SHB_ORDER - 1;
    set_f( st->mem_zero_swb, 0, LPC_SHB_ORDER );

    for( i=0; i<LPC_SHB_ORDER; i++ )
    {
        st->swb_lsp_prev_interp[i] = (float)cos( (float)i * EVS_PI / (float)10.0f );
    }

    st->prev1_shb_ener_sf = 1.0f;
    st->prev2_shb_ener_sf = 1.0f;
    st->prev3_shb_ener_sf = 1.0f;
    st->prev_res_shb_gshape = 0.125f;
    st->prev_mixFactors = 0.5f;
    st->prev_GainShape = 0.0f;
    set_f(st->fb_state_lpc_syn, 0, LPC_SHB_ORDER);
    st->fb_tbe_demph = 0.0f;

    /*-----------------------------------------------------------------*
     * WB/SWB bandwidth switching parameters
     *-----------------------------------------------------------------*/

    st->tilt_swb = 0.0f;
    st->prev_ener = 0.0f;
    st->prev_ener_shb = 0.0f;
    st->prev_enerLH = 0.0f;
    st->enerLH = 0.0f;
    st->enerLL = 0.0f;
    st->prev_enerLL = 0.0f;
    st->prev_fractive = 0;
    st->prev_bws_cnt = 0;
    st->bws_cnt = N_WS2N_FRAMES;
    st->bws_cnt1 = N_NS2W_FRAMES;
    st->attenu1 = 0.1f;
    st->last_inner_frame = L_FRAME8k;
    st->last_bwidth = 0;

    st->prev_weight1 = 0.5f;
    st->GainFrame_prevfrm = 0.0f;

    /*-----------------------------------------------------------------*
     * channel-aware mode parameters
     *-----------------------------------------------------------------*/

    set_f( st->tilt_code_dec, 0.0f, NB_SUBFR16k );

    st->use_partial_copy = 0;
    st->prev_use_partial_copy = 0;
    st->rf_flag = 0;
    st->rf_flag_last = 0;
    st->prev_rf_frame_type = 0;
    st->next_coder_type = 0;

    st->rf_target_bits = 0;

    st->rf_indx_nelp_fid = 0;
    st->rf_indx_nelp_iG1 = 0;
    st->rf_indx_nelp_iG2[0] = 0;
    st->rf_indx_nelp_iG2[1] = 0;
    st->rf_indx_tbeGainFr = 0;

    /*-----------------------------------------------------------------*
     * Improvement of unvoiced and audio signals in AMR-WB IO mode parameters
     *-----------------------------------------------------------------*/

    st->UV_cnt = 30;
    st->LT_UV_cnt = 60.0f;
    set_f(st->lt_diff_etot, 0, MAX_LT);
    st->Last_ener = 0.0f;
    set_f(st->old_Aq, 0, NB_SUBFR * (M+1) );
    st->old_Aq[0] = 1.0f;
    st->old_Aq[M+1] = 1.0f;
    st->old_Aq[2*(M+1)] = 1.0f;
    st->old_Aq[3*(M+1)] = 1.0f;
    st->lt_voice_fac = 0.0f;

    /*-----------------------------------------------------------------*
     * Bass post-filter parameters
     *-----------------------------------------------------------------*/

    bass_psfilter_init( st->pst_old_syn, &(st->pst_mem_deemp_err), &(st->pst_lp_ener) );
    st->bpf_off = 0;
    set_s( st->Track_on_hist, 0, L_TRACK_HIST );
    set_s( st->vibrato_hist, 0, L_TRACK_HIST );
    set_f( st->mem_mean_pit, 80, L_TRACK_HIST );
    st->psf_att = 1.0f;

    /*-----------------------------------------------------------------*
     * FD BPF & resampling tools parameters
     *-----------------------------------------------------------------*/

    /* open analysis for max. SR 48kHz */
    openCldfb ( &st->cldfbAna, CLDFB_ANALYSIS, 48000);

    /* open analysis BPF for max. SR 16kHz */
    openCldfb ( &st->cldfbBPF, CLDFB_ANALYSIS, 16000);

    /* open synthesis for output SR */
    openCldfb ( &st->cldfbSyn, CLDFB_SYNTHESIS, st->output_Fs);

    st->last_active_bandsToZero_bwdec = 0;
    st->flag_NB_bwddec = 0;
    st->perc_bwddec = 0.0f;
    st->last_flag_filter_NB = 0;
    st->active_frame_cnt_bwddec = 0;
    st->total_frame_cnt_bwddec = 0;
    set_s(st->flag_buffer, 0, 20);
    st->avg_nrg_LT = 0.0f;

    /*-----------------------------------------------------------------*
     * Noise gate parameters
     *-----------------------------------------------------------------*/

    st->ng_ener_ST = -51.0f;

    st->Last_frame_ener = (float)MAX_32;
    st->old_Es_pred = 0;
    set_f(st->old_Aq_12_8 + 1, 0, M );
    st->old_Aq_12_8[0] = 1;

    /*-----------------------------------------------------------------*
     * SC-VBR parameters
     *-----------------------------------------------------------------*/

    st->FadeScale = 1.0f;
    st->last_ppp_mode_dec = 0;
    st->old_ppp_mode = 0;
    st->ppp_mode_dec = 0;
    st->last_nelp_mode_dec = 0;
    st->nelp_mode_dec = 0;
    st->nelp_dec_seed = 0;
    st->firstTime_voiceddec = 1;
    st->prev_gain_pit_dec = 0.0f;
    st->prev_tilt_code_dec = 0.0f;
    st->vbr_hw_BWE_disable_dec = 0;
    st->last_vbr_hw_BWE_disable_dec = 0;
    set_f( st->old_hb_synth, 0, L_FRAME48k );

    /* DTFS variables */
    set_f( st->dtfs_dec_a, 0, MAXLAG_WI );
    set_f( st->dtfs_dec_b, 0, MAXLAG_WI );
    st->dtfs_dec_lag = 0;
    st->dtfs_dec_nH = 0;
    st->dtfs_dec_nH_4kHz = 0;
    st->dtfs_dec_upper_cut_off_freq_of_interest = 0;
    st->dtfs_dec_upper_cut_off_freq = 0;
    st->ph_offset_D = 0;
    st->lastLgainD = 0;
    st->lastHgainD = 0;
    set_f( st->lasterbD, 0, NUM_ERB_WB );

    /* NELP decoder variables */
    set_f( st->bp1_filt_mem_nb_dec, 0, 14 );
    set_f( st->bp1_filt_mem_wb_dec, 0, 8 );
    set_f( st->shape1_filt_mem_dec, 0, 20 );
    set_f( st->shape2_filt_mem_dec, 0, 20 );
    set_f( st->shape3_filt_mem_dec, 0, 20 );

    /*-----------------------------------------------------------------*
     * Mode 2 initialization
     *-----------------------------------------------------------------*/

    /* IGF */
    st->igf = 0;
    memset( &st->hIGFDec, 0, sizeof(st->hIGFDec) );
    st->hIGFDec.igfData.igfInfo.nfSeed = 9733;

    st->enablePlcWaveadjust = 0;

    /* Init Core Decoder */
    open_decoder_LPD( st, st->total_brate, st->bwidth );

    /* PLC mode initialization */
    st->m_decodeMode = DEC_NO_FRAM_LOSS;

    /* Init bandwidth / frame_type */
    st->m_frame_type     = ACTIVE_FRAME;
    st->m_old_frame_type = ACTIVE_FRAME;

    resampleCldfb( st->cldfbAna, st->L_frame*50 );
    resampleCldfb( st->cldfbBPF, st->L_frame*50 );

    /* Create FD_CNG instance */
    createFdCngDec( &st->hFdCngDec );

    /* Init FD-CNG */
    initFdCngDec( st->hFdCngDec, st->cldfbSyn->scale );

    st->cngTDLevel = 0.f;

    st->lp_noise = -20.0f;

    st->force_lpd_reset = 0;


    return;
}


/*----------------------------------------------------------------------*
 * reset_preecho_dec()
 *
 * Initialization of static variables for pre-echo
 *----------------------------------------------------------------------*/

void reset_preecho_dec(
    Decoder_State *st       /* i/o: Decoder static variables structure */
)
{
    st->memfilt_lb = 0;
    st->mean_prev_hb = 0;
    st->smoothmem = 1;
    st->mean_prev = 0;
    st->mean_prev_nc = 0;
    st->wmold_hb = 1;
    st->prevflag = 0;
    st->pastpre = 0;

    return;
}


/*----------------------------------------------------------------------*
 * destroy_decoder()
 *
 * Free memory which was allocated in init_decoder()
 *----------------------------------------------------------------------*/

void destroy_decoder(
    Decoder_State *st       /* o:   Decoder static variables structure */
)
{
    /* CLDFB BPF & resampling tools */
    deleteCldfb( &st->cldfbAna );           /* delete analysis for max. SR 16kHz */
    deleteCldfb( &st->cldfbBPF );           /* delete analysis BPF for max. SR 16kHz */
    deleteCldfb( &st->cldfbSyn );           /* delete synthesis for output SR */

    deleteFdCngDec( &st->hFdCngDec );

    return;
}
