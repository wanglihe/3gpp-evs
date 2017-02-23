/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * updt_enc()
 *
 * Common updates (all frame types)
 *-------------------------------------------------------------------*/

void updt_enc(
    Encoder_State *st,               /* i/o: state structure                          */
    const short L_frame,           /* i  : length of the frame                      */
    const short coder_type,        /* i  : speech coder type                        */
    const float *old_exc,          /* i  : buffer of excitation                     */
    const float *pitch_buf,        /* i  : floating pitch for each subframe         */
    const float Es_pred,           /* i  : predicited scaled innovation energy      */
    const float *Aq,               /* i  : A(z) quantized for all subframes         */
    const float *lsf_new,          /* i  : current frame LSF vector                 */
    const float *lsp_new,          /* i  : current frame LSP vector                 */
    const float *old_bwe_exc       /* i  : buffer of excitation                     */
)
{
    short i;

    /* update old excitation buffer */
    mvr2r( &old_exc[L_frame], st->LPDmem.old_exc, L_EXC_MEM );
    if( !st->Opt_AMR_WB )
    {
        mvr2r( &old_bwe_exc[L_FRAME32k], st->old_bwe_exc, PIT16k_MAX * 2 );
    }

    /* update old LSP and LSF vector */
    mvr2r( lsp_new, st->lsp_old, M );
    mvr2r( lsf_new, st->lsf_old, M );

    /* update last coder type */
    st->last_coder_type = coder_type;
    if ( coder_type == INACTIVE || (st->bpf_off == 1 && coder_type != AUDIO && coder_type != TRANSITION) )
    {
        st->last_coder_type = UNVOICED;
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
    }

    /* AC mode (GSC) - in speech we can consider that the last pitch band reached the max */
    if ( coder_type != AUDIO && coder_type != INACTIVE )
    {
        st->mem_last_pit_band = 10 + BAND1k2;
        st->past_dyn_dec = NOISE_LEVEL_SP0-1;   /* tends to speech */
        st->noise_lev = NOISE_LEVEL_SP0-1;      /* tends to speech */
        st->mid_dyn = 40.0f * 0.5f + st->mid_dyn * 0.5f;
    }

    /* convert old LSP vector from 12kHz domain to 16kHz domain (needed in case of ACELP@12k8 <-> ACELP@16kHz switching) */
    if( L_frame == L_FRAME )
    {
        mvr2r( st->lsp_old, st->lsp_old16k, M );

        st->rate_switching_reset_16kHz=lsp_convert_poly( st->lsp_old16k, L_FRAME16k, st->Opt_AMR_WB );
    }

    /* update buffer of old subframe pitch values */
    if( st->last_L_frame != L_frame )
    {
        if( L_frame == L_FRAME )
        {
            for( i=0; i<NB_SUBFR; i++ )
            {
                st->old_pitch_buf[NB_SUBFR+i] = 0.8f * st->old_pitch_buf[NB_SUBFR+i+1];
            }
        }
        else
        {
            for( i=NB_SUBFR; i>0; i-- )
            {
                st->old_pitch_buf[NB_SUBFR+i] = 1.25f * st->old_pitch_buf[NB_SUBFR+i-1];
            }
            st->old_pitch_buf[2*NB_SUBFR16k-1] = st->old_pitch_buf[2*NB_SUBFR16k-2];
        }
    }

    mvr2r( &st->old_pitch_buf[L_frame/L_SUBFR], st->old_pitch_buf, L_frame/L_SUBFR );
    mvr2r( pitch_buf, &st->old_pitch_buf[L_frame/L_SUBFR], L_frame/L_SUBFR );

    /* SC-VBR */
    st->last_Opt_SC_VBR = st->Opt_SC_VBR;
    st->last_last_ppp_mode = st->last_ppp_mode;
    st->last_ppp_mode = st->ppp_mode;
    st->last_nelp_mode = st->nelp_mode;

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

void updt_IO_switch_enc(
    Encoder_State *st,            /* i/o: state structure             */
    const short input_frame       /* i  : input frame length          */
)
{
    float xsp_tmp[M];

    if( st->last_core == AMR_WB_CORE )  /* switching to EVS primary mode */
    {
        /* reset onset detection counter */
        st->tc_cnt = -1;

        /* force safety-net LSFQ in the first frames after the switching */
        st->Nb_ACELP_frames = 0;

        /* AMR-WB IO mode uses ISF(ISP), but EVS primary mode LSF(LSP) */
        mvr2r( stable_LSP, xsp_tmp, M );
        isf2lsf( st->lsf_old, st->lsf_old, xsp_tmp, M, INT_FS_12k8 );
        mvr2r( stable_LSP, xsp_tmp, M );
        isp2lsp( st->lsp_old, st->lsp_old, xsp_tmp, M );
        isp2lsp( st->lsp_old1, st->lsp_old1, xsp_tmp, M );

        mvr2r( stable_LSP, xsp_tmp, M );
        isp2lsp( st->lspCNG, st->lspCNG, xsp_tmp, M );
        if( st->old_enr_index >= 0 )
        {
            st->old_enr_index = min( (short)((float)st->old_enr_index / STEP_AMR_WB_SID * STEP_SID), 127 );
        }
        /* Perform preemphasis of the old input signal @16kHz */
        st->mem_preemph16k = 0;
        preemph( st->old_inp_16k, PREEMPH_FAC_16k, L_INP_MEM, &(st->mem_preemph16k) );

        /* reset TD BWE buffers */
        set_f( st->old_speech_shb, 0.0f, L_LOOK_16k + L_SUBFR16k );
        set_f( st->old_speech_wb, 0.0f, (L_LOOK_12k8 + L_SUBFR) * 5/16 );
        set_f( st->old_bwe_exc, 0.0f, PIT16k_MAX * 2 );
        set_f( st->old_bwe_exc_extended, 0.0f, NL_BUFF_OFFSET );

        st->bwe_non_lin_prev_scale = 0.0;
        set_f( st->decim_state1, 0.0f, (2*ALLPASSSECTIONS_STEEP+1) );
        set_f( st->decim_state2, 0.0f, (2*ALLPASSSECTIONS_STEEP+1) );
        set_f( st->old_wtda_swb, 0, L_FRAME16k );
        set_f( st->old_input_wb, 0, NS2SA(16000, DELAY_FD_BWE_ENC_NS) );

        wb_tbe_extras_reset( st->mem_genSHBexc_filt_down_wb2, st->mem_genSHBexc_filt_down_wb3 );
        if(  input_frame >= L_FRAME32k )
        {
            swb_tbe_reset( st->mem_csfilt, st->mem_genSHBexc_filt_down_shb, st->state_lpc_syn,
                           st->syn_overlap, st->state_syn_shbexc, &(st->tbe_demph),&(st->tbe_premph), st->mem_stp_swb, &(st->gain_prec_swb) );
        }

        if( input_frame == L_FRAME48k )
        {
            set_f( st->fb_state_lpc_syn, 0, LPC_SHB_ORDER );
            st->fb_tbe_demph = 0;
            fb_tbe_reset_enc( st->elliptic_bpf_2_48k_mem, &st->prev_fb_energy );
        }

        /* reset FD BWE buffers */
        st->prev_mode = NORMAL;

        /* reset the unvoiced/audio signal improvement  memories */
        st->seed_tcx = 15687;

        st->use_acelp_preq = 0;

    }
    else            /* switching to AMR-WB IO mode */
    {
        set_f( st->mem_MA, 0, M );

        /* AMR-WB IO mode uses ISF(ISP), but EVS primary mode LSF(LSP) */
        mvr2r( stable_ISP, xsp_tmp, M );
        lsf2isf( st->lsf_old, st->lsf_old, xsp_tmp, M, INT_FS_12k8 );
        mvr2r( stable_ISP, xsp_tmp, M );
        lsp2isp( st->lsp_old, st->lsp_old, xsp_tmp, M );
        mvr2r( st->lsp_old, st->lsp_old1, M );
        lsp2isp( st->lsp_old1, st->lsp_old1, xsp_tmp, M );
        mvr2r( stable_ISP, xsp_tmp, M );
        lsp2isp( st->lspCNG, st->lspCNG, xsp_tmp, M );
        if( st->old_enr_index >= 0 )
        {
            st->old_enr_index = min( (short)((float)st->old_enr_index / STEP_SID * STEP_AMR_WB_SID), 63 );
        }

        /* gain quantization memory */
        set_f(st->past_qua_en, -14.0f, GAIN_PRED_ORDER );

        /* reset VBR signalling */
        st->ppp_mode = 0;
        st->nelp_mode = 0;

        /* reset the unvoiced/audio signal improvement  memories */
        st->seed_tcx = 15687;
    }

    /* Force SID in case of AMR-WB IO mode/EVS primary mode switching */
    st->cnt_SID = 0;

    /* CNG - reset */
    st->cng_cnt = 0;
    st->ho_hist_size = 0;
    st->burst_ho_cnt = 0;

    /* LP memories */
    mvr2r( UVWB_Ave, st->mem_AR, M );

    /* FEC - update adaptive LSF mean vector */
    mvr2r( st->lsf_old, st->lsfoldbfi0, M );
    mvr2r( st->lsf_old, st->lsfoldbfi1, M );
    mvr2r( st->lsf_old, st->lsf_adaptive_mean, M );

    return;
}
