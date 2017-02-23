/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*---------------------------------------------------------------------*
 * core_switching_pre_dec()
 *
 * Preprocessing/preparation for ACELP/HQ core switching
 *---------------------------------------------------------------------*/

void core_switching_pre_dec(
    Decoder_State *st,            /* i/o: decoder state structure     */
    const short output_frame    /* i  : frame length                */
)
{
    short oldLenClasBuff, newLenClasBuff;


    /* Codec mode switching */
    if( st->last_codec_mode == MODE2 )
    {
        mvr2r( st->mem_syn2, st->mem_syn1, M );
        set_f( st->agc_mem2, 0, 2 );
        st->mem_deemph = st->syn[M];
        st->bpf_off = 1;
        set_f( st->pst_old_syn, 0, NBPSF_PIT_MAX);
        st->pst_mem_deemp_err = 0;
        st->psf_lp_noise = st->lp_noise;

        /* reset old HB synthesis buffer */
        if( st->last_L_frame == L_FRAME )
        {
            st->old_bwe_delay = NS2SA( st->output_Fs, MAX_DELAY_TBE_NS - DELAY_SWB_TBE_12k8_NS );
        }
        else
        {
            st->old_bwe_delay = NS2SA( st->output_Fs, MAX_DELAY_TBE_NS - DELAY_SWB_TBE_16k_NS );
        }
        set_f( st->hb_prev_synth_buffer, 0, NS2SA(48000, DELAY_BWE_TOTAL_NS) );

        /* reset upd_cnt */
        st->upd_cnt = MAX_UPD_CNT;
        st->igf = 0;

        if( st->last_core != ACELP_CORE )
        {
            /* reset BWE memories */
            set_f( st->old_bwe_exc, 0, PIT16k_MAX*2 );
            st->bwe_non_lin_prev_scale = 0.0f;
        }

        if( st->output_Fs >= 16000 )
        {
            hf_synth_reset( &st->seed2, st->mem_hf, st->mem_syn_hf, st->mem_hp400, st->mem_hp_interp, st->delay_syn_hf );
        }
        set_f( st->old_syn_12k8_16k, 0, NS2SA(16000, DELAY_FD_BWE_ENC_NS) );

        set_f( st->prev_env, 0, SFM_N_WB );
        set_f( st->prev_normq, 0, SFM_N_WB );

        set_f( st->last_ni_gain, 0, BANDS_MAX );
        set_f( st->last_env, 0, BANDS_MAX );
        st->last_max_pos_pulse = 0;

        if( st->output_Fs > 16000 )
        {
            set_f( st->prev_coeff_out, 0, L_HQ_WB_BWE );
        }

        /* pre-echo */
        st->pastpre = 0;
        /* reset the GSC pre echo energy threshold in case of switching */
        st->Last_frame_ener = (float)MAX_32;

        if( st->last_core == TCX_20_CORE || st->last_core == TCX_10_CORE )
        {
            st->last_core = HQ_CORE;
            mvr2r( st->FBTCXdelayBuf, st->prev_synth_buffer, NS2SA(st->output_Fs, DELAY_BWE_TOTAL_NS - DELAY_CLDFB_NS));
            set_f( st->last_ni_gain, 0, BANDS_MAX );
            set_f( st->last_env, 0, BANDS_MAX );
            st->last_max_pos_pulse = 0;

            set_s( st->prev_SWB_peak_pos, 0, SPT_SHORTEN_SBNUM );
            st->prev_frm_hfe2 = 0;
            st->prev_stab_hfe2 = 0;
        }

        if( st->prev_bfi != 0 )
        {
            short delay_comp;

            /*switch off Hq Voicing as it was not uodated in MODE2*/
            st->oldHqVoicing = 0;
            st->HqVoicing = 0;

            delay_comp = NS2SA(st->output_Fs, DELAY_CLDFB_NS);

            if( !st->last_con_tcx && st->last_core_bfi == ACELP_CORE && st->core == HQ_CORE )
            {
                short i;
                float *realBuffer[CLDFB_NO_COL_MAX], *imagBuffer[CLDFB_NO_COL_MAX];
                float realBufferTmp[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX], imagBufferTmp[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX];

                for( i=0; i<CLDFB_NO_COL_MAX; i++ )
                {
                    set_f( realBufferTmp[i], 0, CLDFB_NO_CHANNELS_MAX );
                    set_f( imagBufferTmp[i], 0, CLDFB_NO_CHANNELS_MAX );
                    realBuffer[i] = realBufferTmp[i];
                    imagBuffer[i] = imagBufferTmp[i];
                }

                /* CLDFB analysis of the synthesis at internal sampling rate */
                cldfb_save_memory( st->cldfbAna );
                cldfbAnalysis( st->syn_Overl, realBuffer, imagBuffer, delay_comp, st->cldfbAna );
                cldfb_restore_memory( st->cldfbAna );

                /* CLDFB synthesis of the combined signal */
                cldfb_save_memory( st->cldfbSyn );
                cldfbSynthesis( realBuffer, imagBuffer, st->fer_samples, delay_comp, st->cldfbSyn );
                cldfb_restore_memory( st->cldfbSyn );
            }

            if( !st->last_con_tcx && st->last_core_bfi == ACELP_CORE && st->core == HQ_CORE )
            {
                lerp(st->syn_Overl, st->fer_samples+delay_comp,output_frame/2, st->last_L_frame/2);
                /*Set to zero the remaining part*/
                set_f( st->fer_samples+delay_comp+output_frame/2, 0, (output_frame/2)-delay_comp);
            }
        }

        st->use_acelp_preq = 0;
        st->reset_mem_AR=0;
    }

    /*FEC*/
    if( st->L_frame <= L_FRAME16k )
    {
        if( st->last_L_frame <= L_FRAME16k && st->core != HQ_CORE )
        {
            if( st->L_frame!=st->last_L_frame )
            {
                if (st->L_frame>st->last_L_frame)
                {
                    oldLenClasBuff = L_SYN_MEM_CLAS_ESTIM * st->last_L_frame/st->L_frame;
                    newLenClasBuff = L_SYN_MEM_CLAS_ESTIM;
                }
                else
                {
                    oldLenClasBuff = L_SYN_MEM_CLAS_ESTIM;
                    newLenClasBuff = L_SYN_MEM_CLAS_ESTIM * st->L_frame/st->last_L_frame;
                }
                lerp( &st->mem_syn_clas_estim[L_SYN_MEM_CLAS_ESTIM-oldLenClasBuff], &st->mem_syn_clas_estim[L_SYN_MEM_CLAS_ESTIM-newLenClasBuff], newLenClasBuff, oldLenClasBuff );
            }
        }
        else
        {
            set_zero( st->mem_syn_clas_estim, L_SYN_MEM_CLAS_ESTIM );
        }
    }

    /* Here we only handle cases where last_ppp and last_nelp not updated when coming from CodecB or other cores
       within ACELP_CORE if switching from another bitarate to vbr, last_ppp and last_nelp is always updated in the previous frame */
    if( st->core == ACELP_CORE && (st->last_core != ACELP_CORE  || st->last_codec_mode == MODE2 ) )
    {
        st->last_ppp_mode_dec = 0;
        st->last_nelp_mode_dec =0;
    }

    /* Handle state reset of stat_noise_uv_mod memory */
    if( st->core == ACELP_CORE && (st->last_core != ACELP_CORE || st->last_codec_mode == MODE2 || st->last_total_brate <= PPP_NELP_2k80 ) )
    {
        st->act_count = 3;
        st->uv_count = 0;
    }

    if( (st->core == ACELP_CORE || st->core == AMR_WB_CORE ) && st->last_core == HQ_CORE )
    {
        if( st->L_frame == L_FRAME16k )
        {
            mvr2r( TRWB2_Ave, st->lsf_old, M ); /* init of LSP */
            mvr2r( TRWB2_Ave, st->lsfoldbfi1, M );
            mvr2r( TRWB2_Ave, st->lsfoldbfi0, M );
            mvr2r( TRWB2_Ave, st->lsf_adaptive_mean, M );
            lsf2lsp( st->lsf_old, st->lsp_old, M, INT_FS_16k );
        }
        else
        {
            mvr2r( TRWB_Ave, st->lsf_old, M ); /* init of LSP */
            mvr2r( TRWB_Ave, st->lsfoldbfi1, M );
            mvr2r( TRWB_Ave, st->lsfoldbfi0, M );
            mvr2r( TRWB_Ave, st->lsf_adaptive_mean, M );
            lsf2lsp( st->lsf_old, st->lsp_old, M, INT_FS_12k8 );
        }

        set_f( st->agc_mem2, 0, 2 );
        st->mem_deemph = 0;
        if( !st->last_con_tcx )
        {
            set_f( st->mem_syn2, 0.0f, M );
        }
        set_f( st->mem_syn1, 0.0f, M );
        st->bwe_non_lin_prev_scale = 0.0f;

        /* Reset ACELP parameters */
        set_zero( st->mem_MA, M );
        if( st->sr_core == 16000 )
        {
            mvr2r( GEWB2_Ave, st->mem_AR, M );
        }
        else
        {
            mvr2r( GEWB_Ave, st->mem_AR, M );
        }
        st->tilt_code = 0.0f;
        st->gc_threshold = 0.0f;
        set_f( st->dispMem, 0, 8 );

        st->last_coder_type = GENERIC;

        fer_energy( output_frame, UNVOICED_CLAS, st->previoussynth, -1, &st->enr_old, 1 );
        st->lp_gainp = 0.0f;
        st->lp_gainc = (float)sqrt( st->lp_ener );


        st->last_voice_factor = 0;
        st->Last_GSC_noisy_speech_flag = 0;

        /* reset CLDFB memories */
        cldfb_reset_memory( st->cldfbAna );
        cldfb_reset_memory( st->cldfbBPF );
        cldfb_reset_memory( st->cldfbSyn );

        /* reset TBE memories */
        if( !st->last_con_tcx )
        {
            set_f( st->old_exc,0, L_EXC_MEM_DEC );
        }
        else if( st->L_frame < L_FRAME16k )
        {
            /* resample from 16kHz to 12.8kHZ */
            synth_mem_updt2( st->L_frame, L_FRAME16k, st->old_exc, st->mem_syn_r, st->mem_syn2, NULL, DEC );
        }

        set_f( st->old_bwe_exc, 0, PIT16k_MAX*2 );

        if( st->output_Fs >= 16000 )
        {
            hf_synth_reset( &st->seed2, st->mem_hf, st->mem_syn_hf, st->mem_hp400, st->mem_hp_interp, st->delay_syn_hf );
        }
        set_f( st->old_syn_12k8_16k, 0, NS2SA(16000, DELAY_FD_BWE_ENC_NS) );
    }

    if( st->core == HQ_CORE && (st->last_core == ACELP_CORE || st->last_core == AMR_WB_CORE) )
    {
        set_f( st->prev_env, 0, SFM_N_WB );
        set_f( st->prev_normq, 0, SFM_N_WB );

        set_f( st->last_ni_gain, 0, BANDS_MAX );
        set_f( st->last_env, 0, BANDS_MAX );
        st->last_max_pos_pulse = 0;

        set_s( st->prev_SWB_peak_pos, 0, SPT_SHORTEN_SBNUM );
        st->prev_frm_hfe2 = 0;
        st->prev_stab_hfe2 = 0;

        if( st->output_Fs > 16000 )
        {
            set_f( st->prev_coeff_out, 0, L_HQ_WB_BWE );
        }

        set_f( st->old_out, 0, output_frame );
    }

    /* handle switching cases where preecho_sb was not called in the last frame (memory not up to date) */
    st->pastpre--;
    if( st->pastpre <= 0 )
    {
        reset_preecho_dec(st);
    }


    if( st->core_brate == FRAME_NO_DATA )
    {
        st->VAD = 0;
        st->m_frame_type = ZERO_FRAME;
    }
    else if( st->core_brate <= SID_2k40 )
    {
        st->VAD = 0;
        st->m_frame_type = SID_FRAME;
    }
    else
    {
        st->VAD = 1;
        st->m_frame_type = ACTIVE_FRAME;
    }

    /*switch on CNA on active frames*/
    if( st->core != AMR_WB_CORE && st->VAD && st->total_brate <= CNA_MAX_BRATE )
    {
        st->flag_cna = 1;
    }
    else if( st->core == AMR_WB_CORE && st->VAD && st->total_brate <= ACELP_8k85 )
    {
        st->flag_cna = 1;
    }
    else if( st->VAD || ((st->cng_type==FD_CNG) && (st->L_frame == L_FRAME16k)))
    {
        st->flag_cna = 0;
    }

    if( st->core == AMR_WB_CORE )
    {
        st->cng_type = LP_CNG;
    }

    /* Reconfigure CNG */
    if ( st->hFdCngDec && ((st->last_L_frame!=st->L_frame ) ||
                           (st->hFdCngDec->hFdCngCom->frameSize!=st->L_frame) ||
                           st->ini_frame == 0 || st->bwidth != st->last_bwidth))
    {
        /* || st->last_core == AMR_WB_CORE || st->last_codec_mode == MODE2)){*/
        if( st->core != AMR_WB_CORE )
        {
            configureFdCngDec( st->hFdCngDec, st->bwidth, st->rf_flag==1&&st->total_brate==13200?9600:st->total_brate, st->L_frame );
        }
        else
        {
            configureFdCngDec( st->hFdCngDec, WB, ACELP_8k00, st->L_frame );

            if( st->VAD )
            {
                st->hFdCngDec->hFdCngCom->CngBitrate = st->total_brate;
            }
        }
        if( st->last_L_frame!=st->L_frame && st->L_frame<=320 && st->last_L_frame<=320 )
        {
            lerp( st->hFdCngDec->hFdCngCom->olapBufferSynth2, st->hFdCngDec->hFdCngCom->olapBufferSynth2, st->L_frame*2, st->last_L_frame*2 );
            if( st->total_brate<=SID_2k40 && st->last_total_brate<=SID_2k40 )
            {
                lerp( st->hFdCngDec->hFdCngCom->olapBufferSynth, st->hFdCngDec->hFdCngCom->olapBufferSynth, st->L_frame*2, st->last_L_frame*2 );
                if( st->L_frame==L_FRAME )
                {
                    short n;
                    for (n=0; n < st->L_frame*2; n++)
                    {
                        st->hFdCngDec->hFdCngCom->olapBufferSynth[n] = st->hFdCngDec->hFdCngCom->olapBufferSynth[n]*0.6250f;
                    }
                }
                else
                {
                    short n;
                    for (n=0; n < st->L_frame*2; n++)
                    {
                        st->hFdCngDec->hFdCngCom->olapBufferSynth[n] =st->hFdCngDec->hFdCngCom->olapBufferSynth[n]*1.6f;
                    }
                }
            }
        }
    }


    return;
}

/*---------------------------------------------------------------------*
 * core_switching_post_dec()
 *
 * Postprocessing for ACELP/HQ core switching
 *---------------------------------------------------------------------*/

void core_switching_post_dec(
    Decoder_State *st,                    /* i/o: decoder state structure     */
    float *synth,                 /* i/o: output synthesis            */
    const short output_frame,           /* i  : frame length                */
    const short core_switching_flag,    /* i  : ACELP->HQ switching flag    */
    const short coder_type              /* i  : ACELP coder type            */
)
{
    short i, delay_comp, delta;
    float tmpF;
    float tmpDelta;
    float synth_subfr_out[SWITCH_MAX_GAP], synth_subfr_bwe[SWITCH_MAX_GAP];
    float mem_synth[NS2SA(16000, DELAY_CLDFB_NS)+2];
    short nZeros;


    if( st->core == ACELP_CORE && st->bfi )
    {
        acelp_core_switch_dec_bfi( st, st->fer_samples, coder_type );
    }

    /* set multiplication factor according to the sampling rate */
    delta = 1;
    if( output_frame == L_FRAME16k )
    {
        delta = 2;
    }
    else if( output_frame == L_FRAME32k )
    {
        delta = 4;
    }
    else if( output_frame == L_FRAME48k )
    {
        delta = 6;
    }

    /* set delay compensation between HQ synthesis and ACELP synthesis */
    delay_comp = delta * HQ_DELAY_COMP;

    if( st->core == HQ_CORE )
    {
        st->use_acelp_preq = 0;
        st->mem_deemph_old_syn = 0.0f;

        if( core_switching_flag && st->last_L_frame == st->last_L_frame_ori && (st->last_core == ACELP_CORE || st->last_core == AMR_WB_CORE))
        {
            acelp_core_switch_dec( st, synth_subfr_out, synth_subfr_bwe, output_frame, core_switching_flag, mem_synth );
        }

        if( core_switching_flag && st->last_core == HQ_CORE && st->prev_bfi )
        {
            mvr2r( st->delay_buf_out, synth_subfr_out, delay_comp );
        }

        /* delay HQ synthesis to synchronize with ACELP synthesis */
        mvr2r( synth, synth + delay_comp, output_frame);
        mvr2r( st->delay_buf_out, synth, delay_comp );
        mvr2r( synth + output_frame, st->delay_buf_out, delay_comp );

        if( core_switching_flag && st->last_L_frame == st->last_L_frame_ori && (st->last_core == ACELP_CORE || st->last_core == AMR_WB_CORE))
        {
            core_switching_OLA( mem_synth, st->last_L_frame, st->output_Fs, synth, synth_subfr_out, synth_subfr_bwe, output_frame, st->bwidth );
        }
        else if( core_switching_flag && st->last_core == HQ_CORE && st->prev_bfi ) /* HQ | ACELP | TRANSITION  with ACELP frame lost */
        {
            nZeros = (short)(NS2SA(st->output_Fs,N_ZERO_MDCT_NS));
            /* Overlapp between old->out[] (stocked in st->fer_samples[]) and good HQ frame on L/2 */
            tmpDelta = 1.0f/(float)(output_frame>>1);
            for( i=0; i<output_frame>>1; i++ )
            {
                tmpF = (float)i * tmpDelta;
                synth[i+delay_comp] = (1-tmpF) * st->fer_samples[i+nZeros] + synth[i+delay_comp] * tmpF;
            }
        }
        else if ( ( !core_switching_flag && st->core == HQ_CORE && (st->last_core == ACELP_CORE || st->last_core == AMR_WB_CORE) ) || /* ACELP | TRANSITION | HQ with TRANSITION lost */
                  ( core_switching_flag && st->prev_bfi && st->last_L_frame != st->last_L_frame_ori ) )                            /* ACELP@12k8 | ACELP@16k | TRANSITION with ACELP@16k lost */
        {
            /* Overlapp between CELP estimation (BFI) and good HQ frame on L/2 */
            tmpDelta = 1.0f/(float)(output_frame>>1);
            for( i=0; i<(output_frame>>1); i++ )
            {
                tmpF = (float)i * tmpDelta;
                synth[i] = synth[i] * tmpF + (1-tmpF) * st->fer_samples[i];
            }
        }

        st->bwe_non_lin_prev_scale = 0.0;
        if ( !(inner_frame_tbl[st->bwidth] == L_FRAME16k && st->core_brate == HQ_32k ) )
        {
            set_f( st->prev_env, 0, SFM_N_WB );
            set_f( st->prev_normq, 0, SFM_N_WB );
        }

        mvr2r( synth, st->previoussynth, output_frame );

        /*Set post-filtering flag to zero*/
        st->pfstat.on = 0;
    }
    else
    {
        /* MDCT to ACELP transition */
        if( st->last_core == HQ_CORE )
        {
            nZeros = (short)(NS2SA(st->output_Fs,N_ZERO_MDCT_NS));
            mvr2r( st->delay_buf_out, synth, delay_comp );  /* copy the HQ/ACELP delay synchroniation buffer at the beginning of ACELP frame */

            if(st->prev_bfi && st->HqVoicing)
            {
                mvr2r( st->fer_samples, st->old_out+nZeros, NS2SA(st->output_Fs,3000000) );
            }

            tmpF = 1.0f/(float)NS2SA(st->output_Fs,3000000);
            for( i=0; i<NS2SA(st->output_Fs,3000000); i++ )
            {
                synth[i+delay_comp] = (1-tmpF*(float)i)*st->old_out[i+nZeros] + tmpF*(float)i*synth[i+delay_comp];
            }
        }

        set_f( st->delay_buf_out, 0, HQ_DELTA_MAX*HQ_DELAY_COMP );
        st->oldHqVoicing = 0;
    }

    /* reset SWB BWE buffers */
    if ( st->bws_cnt == 0 || ( st->bws_cnt > 0 && coder_type != INACTIVE && coder_type != AUDIO ) )
    {
        st->attenu1 = 0.1f;
    }
    if( ( st->last_extl != SWB_BWE && st->extl == SWB_BWE ) || ( st->last_extl != FB_BWE && st->extl == FB_BWE ) ||
            ((st->last_core == HQ_CORE || st->last_extl == SWB_TBE) && st->extl < 0 && st->core != HQ_CORE)
            || (st->last_core == ACELP_CORE && st->core == ACELP_CORE
                && ((st->prev_coder_type != INACTIVE && coder_type != INACTIVE) || (st->prev_coder_type != AUDIO && coder_type == AUDIO))
                && st->bws_cnt >0)
      )
    {
        set_f(st->old_wtda_swb, 0, output_frame);

        if( st->last_extl != WB_BWE )
        {
            st->prev_mode = NORMAL;
        }

        st->prev_Energy = 0.0f;
        st->prev_L_swb_norm = 8;
        st->prev_frica_flag = 0;
        set_f( st->mem_imdct, 0, L_FRAME48k );
        st->prev_td_energy = 0.0f;
        st->prev_weight = 0.2f;
        st->prev_fb_ener_adjust = 0.0f;
    }

    /* reset WB BWE buffers */
    if( st->last_extl != WB_BWE && st->extl == WB_BWE )
    {
        set_f(st->old_wtda_swb, 0, output_frame);
        if ( st->last_extl != SWB_BWE && st->last_extl != FB_BWE )
        {
            st->prev_mode = NORMAL;
        }
        st->prev_Energy_wb = 0.0f;
        st->prev_L_swb_norm = 8;
        set_f( st->mem_imdct, 0, L_FRAME48k );
        st->prev_flag = 0;
    }

    /* reset SWB TBE buffers */
    if( (( st->extl == SWB_TBE || st->extl == FB_TBE || st->extl == SWB_CNG ) &&
            (st->L_frame != st->last_L_frame || ( st->last_extl != SWB_TBE && st->last_extl != FB_TBE ) || st->last_core == HQ_CORE )) ||
            ( st->bwidth < st->last_bwidth && st->last_extl != SWB_TBE ) || st->old_ppp_mode
            || ((st->prev_coder_type == AUDIO || st->prev_coder_type == INACTIVE) && st->bws_cnt > 0 )
            || (st->bws_cnt == 0 && st->prev_bws_cnt == N_WS2N_FRAMES)
      )
    {
        swb_tbe_reset( st->mem_csfilt, st->mem_genSHBexc_filt_down_shb, st->state_lpc_syn,
                       st->syn_overlap, st->state_syn_shbexc, &(st->tbe_demph), &(st->tbe_premph), st->mem_stp_swb, &(st->gain_prec_swb) );

        /* reset GainShape delay for SWB TBE FEC */
        set_f( st->GainShape_Delay, 0, NUM_SHB_SUBFR/2 );

        swb_tbe_reset_synth( st->genSHBsynth_Hilbert_Mem, st->genSHBsynth_state_lsyn_filt_shb_local );

        if( output_frame == L_FRAME16k )
        {
            set_f( st->mem_resamp_HB_32k, 0, 2*ALLPASSSECTIONS_STEEP+1 );        /* reset in case that SWB TBE layer is transmitted, but the output is 16kHz sampled */
        }

        set_f(st->int_3_over_2_tbemem_dec, 0.0f, INTERP_3_2_MEM_LEN);
    }
    else if( (st->extl == SWB_TBE || st->extl == FB_TBE) &&
             ( (st->last_total_brate != st->total_brate) || (st->last_bwidth != st->bwidth) ||
               (st->last_codec_mode != MODE1) || (st->rf_flag != st->rf_flag_last) ) )
    {
        set_f( st->state_lpc_syn, 0.0f, LPC_SHB_ORDER );
        set_f( st->state_syn_shbexc, 0.0f, L_SHB_LAHEAD );
        set_f( st->mem_stp_swb, 0.0f, LPC_SHB_ORDER );
        set_f( st->mem_zero_swb, 0, LPC_SHB_ORDER );
        st->gain_prec_swb = 1.0f;
    }

    /* Interp_3_2 CNG buffers reset */
    if(st->output_Fs == 48000 && ((st->last_core_brate > SID_2k40) && (st->core_brate == FRAME_NO_DATA || st->core_brate == SID_2k40)) )
    {
        set_f(st->interpol_3_2_cng_dec, 0.0f, INTERP_3_2_MEM_LEN);
    }

    /* reset FB TBE buffers */
    if( (st->L_frame != st->last_L_frame || st->last_extl != FB_TBE ) && st->extl == FB_TBE )
    {
        set_f(st->fb_state_lpc_syn, 0, LPC_SHB_ORDER);
        st->fb_tbe_demph = 0;
        fb_tbe_reset_synth( st->fbbwe_hpf_mem,&st->prev_fbbwe_ratio );
    }

    /* reset WB TBE buffers */
    if( st->last_extl != WB_TBE && st->extl == WB_TBE )
    {
        wb_tbe_extras_reset( st->mem_genSHBexc_filt_down_wb2, st->mem_genSHBexc_filt_down_wb3 );
        wb_tbe_extras_reset_synth( st->state_lsyn_filt_shb, st->state_lsyn_filt_dwn_shb, st->mem_resamp_HB );

        set_f( st->state_syn_shbexc, 0, L_SHB_LAHEAD / 4 );
        set_f( st->syn_overlap, 0, L_SHB_LAHEAD );
        set_f( st->mem_csfilt, 0, 2 );
    }


    return;
}


/*---------------------------------------------------------------------*
 * core_switching_hq_prepare_dec()
 *
 * Preprocessing in the first HQ frame after ACELP frame
 * Modify bit allocation for HQ core by removing ACELP subframe budget
*---------------------------------------------------------------------*/

void core_switching_hq_prepare_dec(
    Decoder_State *st,                /* i/o: encoder state structure */
    short *num_bits,          /* i/o: bit budget update       */
    const short output_frame        /* i  : output frame length     */
)
{
    long cbrate;

    if( st->last_core == HQ_CORE && st->prev_bfi )
    {
        mvr2r( st->old_out, st->fer_samples, output_frame );
    }

    /* set switching frame bit-rate */
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

        /* subtract ACELP switching frame bits */
        if( st->core_brate >= ACELP_11k60 )
        {
            (*num_bits)--; /* LP_FLAG bit */
        }

        *num_bits -= ACB_bits_tbl[BIT_ALLOC_IDX(cbrate, GENERIC, 0, 0)];     /* pitch bits*/
        *num_bits -= gain_bits_tbl[BIT_ALLOC_IDX(cbrate, TRANSITION, 0, 0)];   /* gain bits */
        *num_bits -= FCB_bits_tbl[BIT_ALLOC_IDX(cbrate, GENERIC, 0, 0)];    /* FCB bits  */
    }
    else  /* L_frame == L_FRAME16k */
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
            cbrate = min(st->core_brate,ACELP_22k60);
        }

        /* subtract ACELP switching frame bits */
        if( st->core_brate >= ACELP_11k60 )
        {
            /* subtract one bit for LP filtering flag */
            (*num_bits)--;
        }
        *num_bits -= ACB_bits_16kHz_tbl[BIT_ALLOC_IDX_16KHZ(cbrate, GENERIC, 0, 0)];     /* pitch bits*/
        *num_bits -= gain_bits_16kHz_tbl[BIT_ALLOC_IDX_16KHZ(cbrate, GENERIC, 0, 0)];    /* gain bits */
        *num_bits -= FCB_bits_16kHz_tbl[BIT_ALLOC_IDX_16KHZ(cbrate, GENERIC, 0, 0)];     /* FCB bits  */
    }

    /* subtract BWE bits */
    if( !( (inner_frame_tbl[st->bwidth] == L_FRAME16k && st->last_L_frame==L_FRAME16k) || inner_frame_tbl[st->bwidth]==L_FRAME8k ) )
    {
        *num_bits -= (NOOFGAINBITS1 + AUDIODELAYBITS);
    }

    /* reset state of old_out if switching */
    set_f( st->old_out, 0.0f, output_frame );

    return;

}


/*---------------------------------------------------------------------*
 * bandwidth_switching_detect()
 *
 * Classification for band-width switching
 *---------------------------------------------------------------------*/

void bandwidth_switching_detect(
    Decoder_State *st                /* i/o: encoder state structure */
)
{
    /* update band-width switching counter */
    if( st->bws_cnt1 >= N_NS2W_FRAMES )
    {
        st->bws_cnt1 = 0;
    }
    else if( st->total_brate > ACELP_9k60 && st->last_core_brate < ACELP_9k60 && st->bwidth == SWB && st->last_bwidth == WB )
    {
        st->bws_cnt1++;
    }
    else if( st->bws_cnt1 > 0 )
    {
        st->bws_cnt = st->bwidth < st->last_bwidth ? 2*(N_NS2W_FRAMES - st->bws_cnt1)-1 : 0;
        st->bws_cnt1 = st->bwidth < st->last_bwidth ? 0 : ((st->bwidth == SWB) ? st->bws_cnt1+1 : 0); ;
    }

    /* update band-width switching counter */
    if( st->bws_cnt >= N_WS2N_FRAMES )
    {
        st->bws_cnt = 0;
    }
    else if( st->total_brate < ACELP_9k60 && st->last_core_brate > ACELP_9k60 && st->bwidth < st->last_bwidth && st->bwidth == WB )
    {
        st->bws_cnt++;
    }
    else if( st->bws_cnt > 0 )
    {
        st->bws_cnt1 = st->bwidth > st->last_bwidth ? ((N_WS2N_FRAMES - st->bws_cnt) >> 1) : 0;
        st->bws_cnt = st->bwidth > st->last_bwidth ? 0 :(( st->bwidth == WB) ? st->bws_cnt+1 : 0);
    }

    return;
}


/*---------------------------------------------------------------------*
 * bw_switching_pre_proc()
 *
 * Band-width switching pre-processing
 *---------------------------------------------------------------------*/

void bw_switching_pre_proc(
    Decoder_State *st,                  /* i/o: decoder state structure                  */
    const float *old_syn_12k8_16k     /* i  : ACELP core synthesis at 12.8kHz or 16kHz */
)
{
    short i;
    float syn_dct[L_FRAME];

    if( st->core == ACELP_CORE )
    {
        /*----------------------------------------------------------------------*
         * Calculate tilt of the ACELP core synthesis
         *----------------------------------------------------------------------*/
        calc_tilt_bwe( old_syn_12k8_16k, &st->tilt_wb, st->L_frame);
        /*-------------------------------------------------------------------------------*
         * Calculate frequency energy of 0~3.2kHz and 3.2~6.4kHz the ACELP core synthesis
         *-------------------------------------------------------------------------------*/

        edct( old_syn_12k8_16k, syn_dct, L_FRAME );

        st->enerLL = EPSILON;
        for ( i=0; i<L_FRAME/2; i++ )
        {
            st->enerLL += syn_dct[i] * syn_dct[i];
        }
        st->enerLL = (float)sqrt(st->enerLL/128);

        st->enerLH = EPSILON;
        for (; i<L_FRAME; i++ )
        {
            st->enerLH += syn_dct[i] * syn_dct[i];
        }
        st->enerLH = (float)sqrt(st->enerLH/128);
    }
    else
    {
        if( st->old_is_transient[0] )
        {
            st->enerLL = EPSILON;
            for ( i=0; i<32; i++ )
            {
                st->enerLL += st->t_audio_q[i] * st->t_audio_q[i];
            }
            st->enerLL = (float)sqrt(st->enerLL/32);

            st->enerLH = EPSILON;
            for (; i<64; i++ )
            {
                st->enerLH += st->t_audio_q[i] * st->t_audio_q[i];
            }
            st->enerLH = (float)sqrt(st->enerLH/32);
        }
        else
        {
            st->enerLL = EPSILON;
            for ( i=0; i<L_FRAME/2; i++ )
            {
                st->enerLL += st->t_audio_q[i] * st->t_audio_q[i];
            }
            st->enerLL = (float)sqrt(st->enerLL/128);

            st->enerLH = EPSILON;
            for (; i<L_FRAME; i++ )
            {
                st->enerLH += st->t_audio_q[i] * st->t_audio_q[i];
            }
            st->enerLH = (float)sqrt(st->enerLH/128);
        }
    }

    if( st->last_bwidth == 0 && st->extl <= SWB_CNG )
    {
        st->prev_ener_shb = 0.0f;
        set_f(st->prev_SWB_fenv, 0, SWB_FENV);
    }
    else if( ((st->core == ACELP_CORE && st->last_core == HQ_CORE) || (st->core == st->last_core && st->extl != st->last_extl)) && st->last_bwidth >= SWB )
    {
        st->attenu1 = 0.1f;
    }
    if( st->last_core == HQ_CORE || ( st->last_core == ACELP_CORE && !(st->last_extl == WB_TBE || st->last_extl == SWB_TBE || st->last_extl == FB_TBE ) && st->core_brate > ACELP_8k00) )
    {
        st->prev_fractive = 0;
    }

    return;
}
