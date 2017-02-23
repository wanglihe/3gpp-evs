/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "rom_enc.h"
#include "prot.h"


/*-----------------------------------------------------------------------*
 * init_encoder()
 *
 * Initialization of state variables
 *-----------------------------------------------------------------------*/

void init_encoder(
    Encoder_State *st        /* i/o: Encoder static variables structure            */
)
{
    short i;

    st->nb_bits_tot = 0;
    /*-----------------------------------------------------------------*
     * ACELP core parameters
     *-----------------------------------------------------------------*/

    if ( st->Opt_AMR_WB )
    {
        st->last_core = AMR_WB_CORE;
    }
    else
    {
        st->last_core = -1;
    }

    st->L_frame = L_FRAME;
    st->last_coder_type = GENERIC;
    st->last_7k2_coder_type = GENERIC;
    st->last_total_brate = st->total_brate;
    st->last_total_brate_cng = -1;
    st->last_core_brate = st->total_brate;
    st->extl = -1;
    st->last_extl = -1;
    st->last_L_frame = L_FRAME;
    st->rate_switching_reset = 0;
    st->rate_switching_reset_16kHz = 0;

    mvr2r( GEWB_Ave, st->mem_AR, M );
    mvr2r( GEWB_Ave, st->lsfoldbfi0, M );
    mvr2r( GEWB_Ave, st->lsfoldbfi1, M );
    mvr2r( GEWB_Ave, st->lsf_adaptive_mean, M );
    init_lvq( st->offset_scale1, st->offset_scale2, st->offset_scale1_p, st->offset_scale2_p, st->no_scales, st->no_scales_p );
    st->next_force_safety_net = 0;

    st->pstreaklen = 0;
    st->streaklimit = 1.0f;
    set_f(st->mem_MA, 0, M );

    init_gp_clip(st->clip_var);
    pitch_ol_init( &st->old_thres, &st->old_pitch, &st->delta_pit, &st->old_corr) ;

    hf_cod_init( st->mem_hp400_enc, st->mem_hf_enc, st->mem_syn_hf_enc, st->mem_hf2_enc, &st->gain_alpha );

    st->LPDmem.tilt_code = 0.0f;
    st->LPDmem.gc_threshold = 0.0f;

    st->clas = UNVOICED_CLAS;
    set_f( st->old_inp_12k8, 0, L_INP_MEM );
    set_f( st->old_wsp, 0, L_WSP_MEM );
    set_f( st->LPDmem.old_exc, 0, L_EXC_MEM );
    set_f( st->old_wsp2, 0, (L_WSP_MEM - L_INTERPOL)/OPL_DECIM );
    set_f( st->old_inp_16k, 0, L_INP_MEM );

    st->mem_deemph = 0.0f;
    st->mem_preemph = 0.0f;
    st->mem_preemph16k = 0.0f;
    st->mem_preemph_enc = 0.0;

    /* AVQ pre-quantizer memory */
    st->mem_preemp_preQ = 0.0f;
    st->mem_deemp_preQ = 0.0f;
    st->last_nq_preQ = 0;
    st->use_acelp_preq = 0;

    /* (Decimated) Weighted Speech Memory */
    st->mem_wsp_enc = 0.0;

    set_f( st->mem_decim16k, 0, 2*L_FILT_MAX );
    st->mem_wsp = 0.0f;
    st->LPDmem.mem_w0 = 0.0f;
    set_f( st->LPDmem.mem_syn, 0, M );
    set_f( st->mem_syn1, 0, M );
    st->mem_deemph_old_syn = 0.0f;
    set_f( st->LPDmem.mem_syn2, 0, M );
    set_f( st->mem_decim, 0, 2*L_FILT_MAX );
    set_f( st->mem_decim2, 0, 3 );
    set_f( st->Bin_E, 0, L_FFT );
    set_f( st->Bin_E_old, 0, L_FFT/2 );
    set_f( st->LPDmem.mem_syn3, 0, M );

    st->ini_frame = 0;
    st->ee_old = 10.0f;
    st->Nb_ACELP_frames = 0;
    st->audio_frame_cnt = AUDIO_COUNTER_INI;	/* Initializatin of the audio frame counter mildly into the audio mode      */

    /* adaptive lag window memory */
    st->old_pitch_la = 0;
    st->old_voicing_la = 0;

    set_f( st->mem_hp20_in, 0.0f, 4 );

    set_f( st->dispMem, 0, 8 );

    /* HF (6-7kHz) BWE */
    st->seed2_enc = RANDOM_INITSEED;

    for (i=0; i<GAIN_PRED_ORDER; i++)
    {
        st->past_qua_en[i] = -14.0f;   /* gain quantization memory (used in AMR-WB IO mode) */
    }

    if( st->input_Fs == 8000 )
    {
        st->min_band = 1;
        st->max_band = 16;
    }
    else
    {
        st->min_band = 0;
        st->max_band = 19;
    }

    for( i=0; i<NB_BANDS; i++ )
    {
        st->fr_bands1[i] = 1e-5f;
        st->fr_bands2[i] = 1e-5f;
        st->ave_enr2[i] = E_MIN;
    }

    if ( st->Opt_AMR_WB )
    {
        mvr2r( mean_isf_amr_wb, st->lsf_old, M );
        isf2isp( st->lsf_old, st->lsp_old1, M, INT_FS_12k8 );
    }
    else
    {
        mvr2r( GEWB_Ave, st->lsf_old, M );
        lsf2lsp( st->lsf_old, st->lsp_old1, M, INT_FS_12k8 );
    }

    mvr2r( st->lsf_old, st->lsf_old1, M );
    mvr2r( st->lsp_old1, st->lsp_old, M );
    mvr2r( st->lsp_old, st->lsp_old16k, M );
    mvr2r( st->lsp_old, st->lspold_enc, M );

    st->stab_fac = 0.0f;

    MDCT_selector_reset( st );

    /* Bass post-filter memories - enceder side of MODE2 */
    st->bpf_off = 0;
    st->pst_mem_deemp_err = 0.0f;
    st->pst_lp_ener = 0.0f;

    /* TC mode */
    st->tc_cnt = 0;
    st->mCb1 = 0;

    /* AC mode */
    st->seed_tcx = 15687;
    st->cor_strong_limit = 1;
    set_f( st->last_exc_dct_in, 0, L_FRAME );
    st->last_ener = 0.0f;
    set_s( st->last_bitallocation_band, 0, 6 );

    st->mem_last_pit_band = BAND1k2+1;

    st->old_dE1 = 0.0f;
    st->old_ind_deltaMax = 0;
    set_f( st->old_enr_ssf, 0.0f, 2*NB_SSF );
    st->spike_hyst = -1;
    st->music_hysteresis = 0;           /* Counter of frames after AUDIO frame to prevent UC */
    st->last_harm_flag_acelp = 0;
    st->GSC_noisy_speech = 0;

    /*-----------------------------------------------------------------*
     * speech/music classifier
     *-----------------------------------------------------------------*/

    st->inact_cnt = 0;
    set_s( st->past_dec, 0, HANG_LEN-1 );
    set_f( st->past_dlp, 0, HANG_LEN-1 );

    for( i=0; i<NB_BANDS_SPMUS; i++ )
    {
        st->past_log_enr[i] = (float)log(E_MIN);
    }

    st->sp_mus_state = -8;
    st->wdrop = 0.0f;
    st->wdlp_0_95_sp = 0.0f;
    set_f( st->last_lsp, 0.0f, M_LSP_SPMUS );
    st->last_cor_map_sum = 0.0f;
    st->last_non_sta = 0.0f;
    set_f( st->past_PS, 0.0f, HIGHEST_FBIN-LOWEST_FBIN );
    st->past_ps_diff = 0;
    st->past_epsP2 = 01;

    st->gsc_thres[0] = TH_0_MIN;
    st->gsc_thres[1] = TH_1_MIN;
    st->gsc_thres[2] = TH_2_MIN;
    st->gsc_thres[3] = TH_3_MIN;
    set_f( st->gsc_lt_diff_etot, 0.0f, MAX_LT );
    st->gsc_mem_etot = 0.0f;
    st->gsc_last_music_flag = 0;
    st->gsc_nb_thr_1 = 0;
    st->gsc_nb_thr_3 = 0;
    st->mold_corr = 0.9f;
    st->lt_gpitch = 0.0f;
    st->mean_avr_dyn  = 0.5f;
    st->last_sw_dyn = 10.0f;
    st->pit_exc_hangover = 0;
    st->Last_pulse_pos = 0;

    /* speech/music classifier improvement */
    for ( i=0; i<BUF_LEN; i++ )
    {
        st->buf_flux[i] = -100;
        st->buf_pkh[i] = 0;
        st->buf_epsP_tilt[i] = 0;
        st->buf_cor_map_sum[i] = 0;
        st->buf_Ntonal[i] = 0;
        st->buf_Ntonal2[i] = 0;
        st->buf_Ntonal_lf[i] = 0;
    }

    set_f( st->lpe_buf, 0, HANG_LEN_INIT );
    set_f( st->voicing_buf, 0, HANG_LEN_INIT );
    st->gsc_hangover = 0;
    set_f( st->sparse_buf, 0, HANG_LEN_INIT );
    set_f( st->hf_spar_buf, 0, HANG_LEN_INIT );
    st->LT_sparse = 0.0f;
    st->gsc_cnt = 0;
    st->last_vad_spa = 0;

    set_f( st->old_Bin_E, 0.0f, 3*N_OLD_BIN_E );
    set_f( st->buf_etot, 0, 4 );
    set_f( st->buf_dlp, 0, 10 );

    st->UV_cnt1 = 300;
    st->LT_UV_cnt1 = 250.0f;
    st->onset_cnt = 0;
    st->attack_hangover = 0;
    st->dec_mov = 0.0f;
    st->dec_mov1 = 0.0f;
    st->mov_log_max_spl = 200.0f;
    st->old_lt_diff[0] = 0.0f;
    st->old_lt_diff[1] = 0.0f;

    st->Etot_h = 0.0f;
    st->Etot_l = 0.0f;
    st->Etot_l_lp = 0.0f;
    st->Etot_last = 0.0f;
    st->Etot_v_h2 = 0.0f;
    st->sign_dyn_lp = 0.0f;

    /*-----------------------------------------------------------------*
     * GSC
     *-----------------------------------------------------------------*/

    /* GSC - pitch excitation parameters */
    st->mem_w0_tmp = 0.0f;
    set_f(st->mem_syn_tmp, 0.0f, M);
    st->high_stable_cor = 0;
    set_f(st->var_cor_t, 0.0f, VAR_COR_LEN);

    st->lps = 0.0f;
    st->lpm = 0.0f;
    st->Last_frame_ener = (float)MAX_32;
    st->lt_dec_thres = 10.0f;
    st->ener_RAT = 0.0f;
    st->mid_dyn = 40.0f;
    st->noise_lev = NOISE_LEVEL_SP0;
    st->past_dyn_dec = 0;

    /*-----------------------------------------------------------------*
     * VAD & noise estimator
     *-----------------------------------------------------------------*/

    wb_vad_init( &st->nb_active_frames, &st->hangover_cnt, &st->lp_speech, &st->nb_active_frames_he, &st->hangover_cnt_he,
                 &st->bcg_flux, &st->soft_hangover, &st->voiced_burst, &st->bcg_flux_init, &st->nb_active_frames_he1, &st->hangover_cnt_he1,
                 &st->vad_flag_reg_H, &st->vad_flag_reg_L, &st->vad_prim_reg, &st->vad_flag_cnt_50, &st->vad_prim_cnt_16,
                 &st->hangover_cnt_dtx, &st->flag_noisy_speech_snr, &st->hangover_cnt_music );

    st->nb_active_frames_HE_SAD = 0;

    /* Noise estimator */
    noise_est_init( &st->totalNoise, &st->first_noise_updt, st->bckr, st->enrO, st->ave_enr, &st->pitO, &st->aEn,
                    &st->harm_cor_cnt, &st->bg_cnt, &st->lt_tn_track, &st->lt_tn_dist,
                    &st->lt_Ellp_dist,&st->lt_haco_ev,&st->low_tn_track_cnt
                    ,&st->Etot_st_est,&st->Etot_sq_st_est
                  );

    st->epsP_0_2_lp = 1.0f;
    st->epsP_0_2_ad_lp = 0.0f;
    st->epsP_2_16_lp = 1.0f;
    st->epsP_2_16_lp2 = 1.0f;
    st->epsP_2_16_dlp_lp = 0.0f;
    st->epsP_2_16_dlp_lp2 = 0.0f;
    st->lt_aEn_zero = 0.0f;
    st->prim_act_quick = 0.0f;
    st->prim_act_slow = 0.0f;
    st->prim_act = 0.0f;
    st->prim_act_quick_he = 0.0f;
    st->prim_act_slow_he = 0.0f;
    st->prim_act_he = 0.0f;
    st->bckr_tilt_lt = 0.f;

    /*-----------------------------------------------------------------*
     * WB, SWB and FB bandwidth detector
     *-----------------------------------------------------------------*/

    st->lt_mean_NB  = 0;
    st->lt_mean_WB  = 0;
    st->lt_mean_SWB = 0;
    st->count_WB  = BWD_COUNT_MAX;
    st->count_SWB = BWD_COUNT_MAX;
    st->count_FB  = BWD_COUNT_MAX;
    st->bwidth = st->max_bwidth;
    st->last_input_bwidth = st->bwidth;
    st->last_bwidth = st->bwidth;
    st->last_bwidth_cng = st->bwidth;

    /*-----------------------------------------------------------------*
     *
     *-----------------------------------------------------------------*/

    /* Tonal detector */
    for ( i=0; i<L_FFT/2; i++)
    {
        st->old_S[i] = 1;
    }
    set_f(st->cor_map, 0, L_FFT/2 );
    st->act_pred = 1;
    st->noise_char = 0;
    st->multi_harm_limit = THR_CORR;
    st->coder_type_raw = VOICED;
    st->last_coder_type_raw = st->coder_type_raw;

    /* Stationary noise UV modification  */
    st->ge_sm = 10;
    st->uv_count = 0;
    st->act_count = 3;
    mvr2r(st->lsp_old, st->lspold_s, M);
    st->noimix_seed = RANDOM_INITSEED;
    st->min_alpha = 1;
    st->exc_pe = 0;

    /*-----------------------------------------------------------------*
     * CNG and DTX
     *-----------------------------------------------------------------*/

    st->lp_noise = 0.0f;
    mvr2r( st->lsp_old1, st->lspCNG, M );
    st->cng_seed = RANDOM_INITSEED;
    st->cng_ener_seed = RANDOM_INITSEED;
    st->cng_ener_seed1 = RANDOM_INITSEED;
    st->lp_ener = 0.0f;
    st->first_CNG = 0;
    st->cnt_SID = 0;
    st->max_SID = 2;
    st->old_enr_index = -1;
    st->Enew = 0.0f;
    st->VarDTX_cnt_voiced = 0;
    st->lt_ener_voiced = 0.0f;
    st->VarDTX_cnt_noise = 0;
    st->lt_ener_noise = 0.0f;
    st->lt_ener_last_SID = 0.0f;
    if( st->var_SID_rate_flag )
    {
        st->interval_SID = 12;
    }
    st->lp_sp_enr = 0.0f;
    st->last_allow_cn_step = 0;

    st->fd_cng_reset_flag = 0;

    if( st->Opt_DTX_ON )
    {
        st->cng_hist_ptr = -1;
        set_f( st->cng_lsp_hist, 0, DTX_HIST_SIZE*M );
        set_f( st->cng_ener_hist, 0, DTX_HIST_SIZE );
        st->cng_cnt = 0;
        st->ho_hist_ptr = -1;
        st->ho_sid_bw = 0;
        set_f( st->ho_lsp_hist, 0, HO_HIST_SIZE*M );
        set_f( st->ho_ener_hist, 0, HO_HIST_SIZE );
        set_f( st->ho_env_hist, 0, HO_HIST_SIZE*NUM_ENV_CNG );
        st->ho_hist_size = 0;
        st->act_cnt = 0;
    }

    st->active_cnt = 0;
    st->cng_type = -1;

    st->CNG_mode = -1;
    st->last_active_brate = ACELP_7k20;
    st->last_CNG_L_frame = L_FRAME;
    set_s( st->ho_16k_lsp, 0, HO_HIST_SIZE );
    st->act_cnt2 = 0;
    st->num_ho = 0;
    st->hangover_terminate_flag = 0;

    st->ho_circ_ptr = -1;
    set_f( st->ho_lsp_circ, 0, HO_HIST_SIZE*M );
    set_f( st->ho_ener_circ, 0, HO_HIST_SIZE );
    set_f( st->ho_env_circ, 0, HO_HIST_SIZE*NUM_ENV_CNG );
    st->ho_circ_size = 0;
    st->burst_ho_cnt = 0;
    st->cng_buf_cnt = 0;

    if( st->var_SID_rate_flag || ((!st->var_SID_rate_flag) && (st->interval_SID >= DTX_HIST_SIZE)) )
    {
        st->cng_hist_size = DTX_HIST_SIZE;
    }
    else
    {
        st->cng_hist_size = st->interval_SID;
    }
    set_f(st->lp_env, 0.0f, 20);
    set_f(st->cng_res_env, 0.0f, 20*8);
    set_f(st->exc_mem, 0.0f, 24);
    set_f(st->exc_mem1, 0.0f, 30);
    set_f(st->exc_mem2, 0.0f, 30);
    set_f(st->old_env, 0.0f, NUM_ENV_CNG);
    /* SWB CNG/DTX */
    st->last_wb_cng_ener = -6.02f;
    st->last_shb_cng_ener = -6.02f;
    st->mov_wb_cng_ener = -6.02f;
    st->mov_shb_cng_ener = -6.02f;
    st->shb_cng_ini_cnt = 1;
    st->shb_NO_DATA_cnt = 0;
    st->last_SID_bwidth = min( st->max_bwidth, SWB );
    st->last_vad = 0;
    /* FEC */
    st->last_clas = UNVOICED_CLAS;

    for (i=0; i<2*NB_SUBFR16k; i++)
    {
        st->old_pitch_buf[i] = L_SUBFR;
    }

    st->old_Es_pred = 0;
    set_f( st->old_Aq_12_8 + 1, 0, M );
    st->old_Aq_12_8[0] = 1;

    /*-----------------------------------------------------------------*
     * CLDFB & resampling tools parameters
     *-----------------------------------------------------------------*/

    openCldfb( &st->cldfbAnaEnc, CLDFB_ANALYSIS, st->input_Fs );

    st->currEnergyLookAhead = 6.1e-5f;

    /*-----------------------------------------------------------------*
     * SC-VBR parameters
     *-----------------------------------------------------------------*/

    st->nelp_enc_seed = 0;
    st->last_nelp_mode = 0;
    st->pppcountE = 0;
    st->last_ppp_mode = 0;
    st->last_last_ppp_mode = 0;
    st->firstTime_voicedenc = 1;
    st->prev_ppp_gain_pit = 0.0;
    st->prev_tilt_code = 0.0;

    st->ppp_mode = 0;
    st->nelp_mode = 0;

    /* stable short pitch detection */
    st->voicing0_sm = 0;
    st->voicing_sm = 0;
    st->LF_EnergyRatio_sm = 1;
    st->predecision_flag = 0;
    st->diff_sm = 0;
    st->energy_sm = 0;

    st->pattern_m = 0;
    st->Last_Resort = 0;
    st->set_ppp_generic = 0;
    st->Q_to_F = 0;

    st->numactive = 0;			        /* keep the count of the frames inside current 600 frame bloack.*/
    st->sum_of_rates = 0.0f;			/* sum of the rates of past 600 active frames*/
    st->global_avr_rate = 0.0f;	        /* global rate upto current time. recorded a (rate in kbps) *6000*/
    st->global_frame_cnt = 0;		    /* 600 active frame block count. Used to update the global rate*/

    st->rate_control = 0;
    st->SNR_THLD = 67.0f;
    st->mode_QQF = 1;
    st->last_Opt_SC_VBR = 0;

    set_f(st->shape1_filt_mem, 0, 20);
    set_f(st->shape2_filt_mem, 0, 20);
    set_f(st->shape3_filt_mem, 0, 20);
    set_f(st->txlpf1_filt1_mem, 0, 20);
    set_f(st->txlpf1_filt2_mem, 0, 20);
    set_f(st->txhpf1_filt1_mem, 0, 20);
    set_f(st->txhpf1_filt2_mem, 0, 20);

    /*-----------------------------------------------------------------*
     * SWB BWE parameters
     *-----------------------------------------------------------------*/

    set_f( st->new_input_hp, 0, NS2SA(16000, ACELP_LOOK_NS + DELAY_FD_BWE_ENC_12k8_NS + DELAY_FIR_RESAMPL_NS - DELAY_CLDFB_NS));
    set_f( st->old_input, 0, NS2SA(48000, DELAY_FD_BWE_ENC_12k8_NS + DELAY_FIR_RESAMPL_NS) );
    set_f( st->old_input_wb, 0, NS2SA(16000, DELAY_FD_BWE_ENC_12k8_NS) );
    set_f( st->old_input_lp, 0, NS2SA(16000, ACELP_LOOK_NS + DELAY_FD_BWE_ENC_16k_NS) );
    set_f( st->old_syn_12k8_16k, 0, NS2SA(16000, DELAY_FD_BWE_ENC_NS) );

    st->prev_mode = NORMAL;
    set_f( st->old_wtda_swb, 0, L_FRAME48k );
    st->prev_L_swb_norm1 = 8;
    st->prev_global_gain = 0.0f;
    st->modeCount = 0;
    st->EnergyLF = 0.0f;


    /*-----------------------------------------------------------------*
     * TBE parameters
     *-----------------------------------------------------------------*/

    InitSWBencBuffer(st);
    ResetSHBbuffer_Enc(st);
    set_f( st->old_speech_shb, 0.0f, L_LOOK_16k + L_SUBFR16k );
    set_f( st->old_speech_wb, 0.0f, (L_LOOK_12k8 + L_SUBFR) * 5/16 );
    set_f( st->old_input_fhb, 0.0f, NS2SA(48000, ACELP_LOOK_NS + DELAY_FD_BWE_ENC_12k8_NS + DELAY_FIR_RESAMPL_NS) - L_FRAME48k/2);

    for( i=0; i < LPC_SHB_ORDER; i++ )
    {
        st->prev_lsp_shb[i] = i/20.0f;
    }

    st->cldfbHBLT = 1.0f;
    st->prev_gainFr_SHB = 0;
    set_f( st->lsp_shb_slow_interpl, 0, LPC_SHB_ORDER );
    set_f( st->lsp_shb_fast_interpl, 0, LPC_SHB_ORDER );
    set_f( st->shb_inv_filt_mem, 0, LPC_SHB_ORDER );
    set_f( st->lsp_shb_spacing, 0.1f, 3);
    st->prev_swb_GainShape = 0;
    st->prev_frGainAtten = 0;
    st->prev_wb_GainShape = 0;
    set_f(st->fb_state_lpc_syn, 0, LPC_SHB_ORDER);
    st->fb_tbe_demph = 0.0f;
    st->tilt_mem = 0.0f;

    openCldfb( &st->cldfbSynTd, CLDFB_SYNTHESIS, 16000);

    st->prev_coder_type = GENERIC;
    set_f( st->prev_lsf_diff, 0.5f, LPC_SHB_ORDER-2 );
    st->prev_tilt_para = 0.0f;
    set_f( st->cur_sub_Aq, 0.0f, M+1 );

    /* TD BWE post-processing */
    st->ptr_mem_stp_swb = st->mem_stp_swb + LPC_SHB_ORDER - 1;
    set_f( st->mem_zero_swb, 0, LPC_SHB_ORDER );

    for( i=0; i<LPC_SHB_ORDER; i++ )
    {
        st->swb_lsp_prev_interp[i] = (float)cos( (float)i * EVS_PI / (float)10.0f );
    }

    set_f( st->dec_2_over_3_mem, 0.0f, 12 );
    set_f( st->dec_2_over_3_mem_lp, 0.0f, 6 );
    set_f( st->old_fdbwe_speech, 0.0f, L_FRAME48k );

    /*-----------------------------------------------------------------*
     * HQ core parameters
     *-----------------------------------------------------------------*/

    st->input = st->input_buff+L_FRAME48k+NS2SA(48000, DELAY_FIR_RESAMPL_NS);
    set_zero( st->input_buff+L_FRAME48k, L_FRAME48k+NS2SA(48000, DELAY_FIR_RESAMPL_NS) );
    st->old_input_signal = st->input - NS2SA(st->input_Fs, DELAY_FIR_RESAMPL_NS) - (short)(st->input_Fs / 50);

    st->old_hpfilt_in = 0.0f;
    st->old_hpfilt_out = 0.0f;
    st->EnergyLT = 0.0f;
    st->Energy_Old = 0;
    st->TransientHangOver = 0;

    set_f( st->old_out, 0, L_FRAME32k );

    st->mode_count = 0;
    st->mode_count1 = 0;

    st->hq_generic_speech_class = 0;
    st->prev_Npeaks = 0;
    set_s( st->prev_peaks, 0, HVQ_MAX_PEAKS );
    st->hvq_hangover = 0;
    st->prev_hqswb_clas = HQ_NORMAL;
    set_s( st->prev_SWB_peak_pos, 0, SPT_SHORTEN_SBNUM );
    st->clas_sec_old = 1.0f;
    st->clas_final_old = 1;
    st->last_gain1 = 0.0f;
    st->last_gain2 = 0.0f;

    /* speech/music classification */
    set_s( st->lt_old_mode, 1, 3 );
    st->lt_voicing = 0.5f;
    st->lt_corr = 0.5f;
    st->lt_tonality = 0;
    set_s( st->lt_corr_pitch, 0, 3 );
    st->lt_hangover = 0;
    st->lowrate_pitchGain = 0;

    st->lt_music_hangover = 0;
    set_f( st->tonality2_buf, 0, HANG_LEN_INIT );
    set_f( st->tonality3_buf, 0, HANG_LEN_INIT );
    set_f( st->LPCErr_buf, 0,HANG_LEN_INIT );
    st->lt_music_state    = 0;
    st->lt_speech_state     = 0;
    st->lt_speech_hangover  = 0;
    st->consec_inactive = 0;
    st->spectral_tilt_reset = 1;
    st->running_avg = 0;
    st->ra_deltasum = 0;
    st->trigger_SID = 0;
    st->snr_sum_vad = 0;

    set_s( st->prev_frm_index, -1, NB_SWB_SUBBANDS_HAR_SEARCH_SB );
    st->prev_frm_hfe2 = 0;
    st->prev_stab_hfe2 = 0;
    st->prev_ni_ratio = 0.5f;
    set_f( st->prev_En_sb, 0.0f, NB_SWB_SUBBANDS );
    set_s( st->last_bitalloc_max_band, 0, 2 );
    set_f( st->last_ni_gain, 0, BANDS_MAX );
    set_f( st->last_env, 0, BANDS_MAX );
    st->last_max_pos_pulse = 0;


    /*-----------------------------------------------------------------*
     * Channel-aware mode
     *-----------------------------------------------------------------*/

    if( !st->Opt_RF_ON || (st->bwidth != WB && st->bwidth != SWB) || st->total_brate != ACELP_13k20 )
    {
        if ( st->Opt_RF_ON )
        {
            printf("\nWarning: Channel-aware mode only available for 13.2 kbps WB/SWB\n");
            printf("         Switched to normal mode!\n");
            st->Opt_RF_ON = 0;
            st->rf_fec_offset = 0;
        }
        st->rf_mode = 0;
    }
    else
    {
        st->rf_mode = st->Opt_RF_ON;
    }

    st->rf_mode_last = st->rf_mode;

    /* initialize RF indice buffers */
    reset_rf_indices( st );

    /*-----------------------------------------------------------------*
     * Mode 2 initialization
     *-----------------------------------------------------------------*/

    st->last_sr_core = st->last_L_frame * 50;

    if( st->codec_mode == MODE2 )
    {
        st->igf = getIgfPresent( st->total_brate, st->bwidth, st->rf_mode );
    }
    else
    {
        st->igf = 0;
    }

    /* FD-CNG encoder */
    createFdCngEnc( &st->hFdCngEnc );
    initFdCngEnc( st->hFdCngEnc, st->input_Fs, st->cldfbAnaEnc->scale );
    configureFdCngEnc( st->hFdCngEnc, st->bwidth, st->rf_mode&&st->total_brate==13200?9600:st->total_brate );

    /*  INIT CORE CODER  */

    st->last_totalNoise = 0.f;
    set_f( st->totalNoise_increase_hist, 0.f, TOTALNOISE_HIST_SIZE );
    st->totalNoise_increase_len = 0;

    init_coder_ace_plus( st );

    InitTransientDetection( (int)(st->input_Fs / 50), NS2SA(st->input_Fs, DELAY_FIR_RESAMPL_NS), &st->transientDetection );

    reset_indices_enc( st );


    st->vbr_generic_ho = 0;

    st->sharpFlag = 0;

    st->Local_VAD = 0;
    set_f( st->nelp_lp_fit_mem, 0, NELP_LP_ORDER*2 );

    return;
}



/*-----------------------------------------------------------------------*
 * destroy_encoder()
 *
 * Free memory which was allocated in init_encoder()
 *-----------------------------------------------------------------------*/

void destroy_encoder(
    Encoder_State *st        /* i/o: Encoder static variables structure   */
)
{
    deleteCldfb( &st->cldfbSynTd );
    deleteCldfb( &st->cldfbAnaEnc );
    deleteFdCngEnc( &st->hFdCngEnc );

    /* Close Core */

    return;
}
