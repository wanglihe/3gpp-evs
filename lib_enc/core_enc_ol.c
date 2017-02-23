/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"
#include "basop_proto_func.h"

/*-------------------------------------------------------------------*
 * Local functions
 *-------------------------------------------------------------------*/

static void BITS_ALLOC_ACELP_config_rf( const short coder_type, float *tilt_code, short *rf_frame_type,
                                        short *rf_target_bits, short nb_subfr, short rf_fec_indicator, float *pitch_buf );

static void BITS_ALLOC_TCX_config_rf( short *rf_frame_type, short *rf_target_bits, short PLC_Mode,
                                      short coder_type, short last_core, int   TD_Mode );

static void closest_centroid_rf( const float *data, const float *weights, const float *quantizer,
                                 const short centroids, const short  length, short *ind_vec );

/*-------------------------------------------------------------------*
 * core_encode_openloop()
 *
 * Open-loop core encoder
 *-------------------------------------------------------------------*/

void core_encode_openloop(
    Encoder_State *st,                    /* i/o: encoder state structure             */
    const short coder_type,             /* i  : coding type                         */
    const short pitch[3],               /* i  : open-loop pitch values for quantiz. */
    const float voicing[3],             /* i  : open-loop pitch gains               */
    const float Aw[NB_SUBFR16k*(M+1)],  /* i  : weighted A(z) unquant. for subframes*/
    const float *lsp_new,               /* i  : LSPs at the end of the frame        */
    const float *lsp_mid,               /* i  : LSPs at the middle of the frame     */
    float *pitch_buf,             /* i/o: floating pitch values for each subfr*/
    float *voice_factors,         /* o  : voicing factors                     */
    float *ptr_bwe_exc,           /* o  : excitation for SWB TBE              */
    const short vad_hover_flag
)
{
    float lsf_q[M], lsp_q[M], lspmid_q[M];
    Word16 lspq_ind[M];
    Word16 A_q_ind[M+1];
    float A_q_ace[NB_SUBFR16k*(M+1)];
    float A_q_tcx[NB_SUBFR16k*(M+1)];
    int param_lpc[NPRM_LPC_NEW];
    int nbits_lpc;
    int param_core[2*NPRM_DIV];
    int target_bits;
    float stab_fac;
    int indexBuffer[N_MAX+1];
    CONTEXT_HM_CONFIG hm_cfg;
    float lsp_tcx_q[M], lsf_tcx_q[M];
    int tcx_lpc_cdk;
    float A_w[M+1];
    float gain_pitch_buf[NB_SUBFR16k];
    float gain_code_buf[NB_SUBFR16k];
    short bits_param_lpc[10], no_param_lpc;

    /* lsf quant parameters */
    float lsp_q_rf[M];
    float Aq_rf[NB_SUBFR*(M+1)];
    float stab_fac_rf;
    float *exc_rf;
    float *syn_rf;
    short rf_PLC_Mode;
    short TD_Mode;
    short rf_tcx_lpc_cdk;
    float lsp[M], lsf[M];
    float rf_mem_MA[M];
    float exc_buf_rf[L_EXC_MEM + L_FRAME + 1];
    float syn_buf_rf[M+L_FRAME16k+L_FRAME16k/2];

    float w_rf[M], lsf_uq_rf[M];
    float lsf_q_1st_rf[M], lsf_q_d_rf[M], lsf_q_rf[M];
    float lsp_old_q_rf[M], lsf_old_q_rf[M];

    /*--------------------------------------------------------------*
     * back up parameters for RF
     *---------------------------------------------------------------*/

    /* back up the old LSPs and LSFs */
    mvr2r( st->lsp_old, lsp_old_q_rf, M );
    mvr2r( st->lsf_old, lsf_old_q_rf, M );

    /* back up old exc before primary encoding */
    set_f( exc_buf_rf, 0, (L_EXC_MEM+L_FRAME+1) );
    exc_rf = exc_buf_rf + L_EXC_MEM;
    mvr2r(st->LPDmem.old_exc, exc_buf_rf, L_EXC_MEM);

    /* back up old synthesis before primary encoding */
    set_f( syn_buf_rf, 0, (M+L_FRAME16k+L_FRAME16k/2) );
    syn_rf = syn_buf_rf + M;
    mvr2r(st->LPDmem.mem_syn, syn_buf_rf, M);

    /* back up syn2 mem */
    mvr2r(st->LPDmem.mem_syn2, st->rf_mem_syn2, M);

    /* back up LPD mem_w0 target generation memory */
    st->rf_mem_w0 = st->LPDmem.mem_w0;

    /* back up clip gain memory */
    mvr2r( st->clip_var, st->rf_clip_var, 6 );

    /* back up tilt code */
    st->rf_tilt_code = st->LPDmem.tilt_code;

    /* back up dispMem */
    mvr2r( st->dispMem, st->rf_dispMem, 8 );

    /* back up gc_threshold for noise addition */
    st->rf_gc_threshold = st->LPDmem.gc_threshold;



    /*--------------------------------------------------------------*
     * Initializations
     *---------------------------------------------------------------*/

    tcx_lpc_cdk = 0;
    set_i( param_lpc, 0, NPRM_LPC_NEW );
    set_i( param_core, 0, 2*NPRM_DIV );
    mvi2i( st->tcxltp_param, &param_core[1+NOISE_FILL_RANGES], LTPSIZE );

    no_param_lpc = 0;     /* avoid MSVC warnings */
    nbits_lpc = 0;        /* avoid MSVC warnings */
    stab_fac = 0.0f;      /* avoid MSVC warnings */

    hm_cfg.indexBuffer = indexBuffer;

    /*--------------------------------------------------------------*
     * LPC Quantization
     *---------------------------------------------------------------*/

    if( st->lpcQuantization == 1 && coder_type == VOICED )
    {
        (&(st->acelp_cfg))->midLpc = 0;
    }
    else
    {
        (&(st->acelp_cfg))->midLpc=st->acelp_cfg.midLpc_enable;
    }

    if( st->core == ACELP_CORE || !st->enableTcxLpc )
    {
        if( st->envWeighted )
        {
            /* Unweight the envelope */
            E_LPC_lsp_unweight( st->lsp_old, st->lsp_old, st->lsf_old, 1.0f/st->gamma );
            st->envWeighted = 0;
        }

        if( st->core == TCX_20_CORE )
        {
            lpc_quantization( st, st->core, st->lpcQuantization, st->lsf_old, lsp_new, lsp_mid, lsp_q, lsf_q,
                              lspmid_q, st->mem_MA, st->mem_AR, st->narrowBand, AUDIO, st->acelp_cfg.midLpc, param_lpc, &nbits_lpc,
                              &(st->seed_acelp), st->sr_core, st->Bin_E, st->Bin_E_old, bits_param_lpc, &no_param_lpc );
        }
        else
        {
            lpc_quantization( st, st->core, st->lpcQuantization, st->lsf_old, lsp_new, lsp_mid, lsp_q, lsf_q,
                              lspmid_q, st->mem_MA, st->mem_AR, st->narrowBand, coder_type, st->acelp_cfg.midLpc, param_lpc, &nbits_lpc,
                              &(st->seed_acelp), st->sr_core, st->Bin_E, st->Bin_E_old, bits_param_lpc, &no_param_lpc );
        }

        /*-------------------------------------------------------------*
         * Rate switching: reset
         *-------------------------------------------------------------*/

        if( st->rate_switching_reset )
        {
            mvr2r( lsp_q, st->lsp_old, M );
            mvr2r( lsf_q, st->lsf_old, M );
            mvr2r( lsp_q, lspmid_q, M );
        }

        /*--------------------------------------------------------------*
        * LPC Interpolation
        *---------------------------------------------------------------*/

        stab_fac = lsf_stab( lsf_q, st->lsf_old, 0, st->L_frame );
    }




    /*--------------------------------------------------------------*
    * Run ACELP
    *---------------------------------------------------------------*/

    if( st->core == ACELP_CORE )
    {
        if( st->acelp_cfg.midLpc )
        {
            int_lsp4( st->L_frame, st->lsp_old, lspmid_q, lsp_q, A_q_ace, M, 0 );
        }
        else
        {
            int_lsp( st->L_frame, st->lsp_old, lsp_q, A_q_ace, M, interpol_frac_12k8, 0 );
        }

        /* Calculate target bits */
        target_bits = st->bits_frame_core - nbits_lpc - st->nb_bits_header_ace;

        if( st->rf_mode )
        {
            /* joint bit allocation for redundant frame and TBE */
            /* calculate target bits for core coding */
            target_bits -= st->rf_target_bits_write;
        }

        if( st->igf )
        {
            target_bits -= get_tbe_bits( st->total_brate, st->bwidth, st->rf_mode );
        }

        if( st->acelp_cfg.midLpc )
        {
            target_bits -= MIDLSF_NBITS;
        }

        if( st->plcExt.enableGplc )
        {
            target_bits -= st->plcExt.nBits;
        }

        /* reset TBE buffers previous frame frame wasn't ACELP*/
        if( st->last_core != ACELP_CORE )
        {
            TBEreset_enc( st, st->bwidth );
        }

        /* Run ACELP encoder */
        coder_acelp( &(st->acelp_cfg), coder_type, Aw, A_q_ace, st->speech_enc_pe,
                     st->synth, &(st->LPDmem), voicing, pitch, param_core, stab_fac, st, &st->plcExt, target_bits,
                     gain_pitch_buf, gain_code_buf, pitch_buf, voice_factors, ptr_bwe_exc
                   );

        st->glr_idx[0] = encSideSpecPowDiffuseDetector( st->plcExt.last_lsf_ref, st->plcExt.last_lsf_con,
                         st->last_sr_core, &(st->prev_lsf4_mean), st->glr, coder_type );

        mvr2r( lsf_q, st->plcExt.last_lsf_ref, M );
        mvr2r( st->plcExt.lsf_con, st->plcExt.last_lsf_con, M );

        updateSpecPowDiffuseIdx( st, gain_pitch_buf, gain_code_buf );

        if( st->last_stab_fac > 0.02 )
        {
            st->glr_idx[0] = 0;
        }
        st->last_stab_fac = stab_fac;

        st->plcExt.LPDmem = &(st->LPDmem);

        encoderSideLossSimulation( st, &st->plcExt, lsf_q, stab_fac, st->plcExt.calcOnlylsf, st->L_frame );

        st->tcxltp_norm_corr_past = voicing[1];
        st->tcx_cfg.tcx_curr_overlap_mode = st->tcx_cfg.tcx_last_overlap_mode = ALDO_WINDOW;

    }



    /*--------------------------------------------------------------*
    * Run TCX20
    *---------------------------------------------------------------*/

    if( st->core == TCX_20_CORE )
    {
        if( st->enableTcxLpc )
        {
            if( st->rf_mode )
            {
                mvr2r( st->mem_MA, rf_mem_MA, M );
            }

            tcx_lpc_cdk = tcxlpc_get_cdk( st->tcx_cfg.coder_type );

            /* Get the envelope corresponding to the current frame */
            E_LPC_int_lpc_tcx( st->lspold_enc, lsp_new, A_q_tcx );

            /* Weight the envelope */
            weight_a( A_q_tcx, A_q_tcx, st->gamma, M );

            /* Save the weighted envelope */
            mvr2r( A_q_tcx, A_w, M+1 );

            /* Convert to lsp and lsf */
            a2lsp_stab( A_q_tcx, lsp, lsp_new );
            lsp2lsf( lsp, lsf, M, INT_FS_12k8 );

            /* Quantize */
            Q_lsf_tcxlpc( lsf, lsf_tcx_q, lspq_ind, param_lpc, st->narrowBand, tcx_lpc_cdk, st->mem_MA, st->tcx_cfg.coder_type, st->Bin_E );

            /* Account for consumed bits */
            nbits_lpc = TCXLPC_NUMBITS;
            if( param_lpc[0] )
            {
                nbits_lpc += TCXLPC_IND_NUMBITS;
            }

            /* Convert quantized lsf to lsp and A */
            lsf2lsp( lsf_tcx_q, lsp_tcx_q, M, INT_FS_12k8 );
            lsp2a_stab( lsp_tcx_q, A_q_tcx, M );
        }
        else
        {
            E_LPC_int_lpc_tcx( st->lsp_old, lsp_q, A_q_tcx );
        }

        if( st->tcx_lpc_shaped_ari )
        {
            basop_E_LPC_f_lsp_a_conversion(lspq_ind, A_q_ind, M);
        }

        /* Calculate target bits */
        target_bits = st->bits_frame_core - nbits_lpc - st->nb_bits_header_tcx;
        if( st->rf_mode )
        {
            /* joint bit allocation for redundant frame and TBE */
            /* calculate target bits for core coding */
            target_bits -= st->rf_target_bits_write;
        }

        if( st->mdct_sw == MODE1 )
        {
            /* Account for core mode signalling bits difference: bandwidth and ACELP/TCX signaling bit are replaced */
            target_bits += (FrameSizeConfig[st->frame_size_index].bandwidth_bits + 1) - signalling_mode1_tcx20_enc(st, 0);
        }
        else if( st->mdct_sw_enable == MODE2 )
        {
            --target_bits;
        }

        if( st->plcExt.enableGplc )
        {
            target_bits -= st->plcExt.nBits;
        }

        /* subtract bits for TCX overlap mode (1 bit: full, 2 bits: half or no overlap) */
        target_bits -= (st->tcx_cfg.tcx_curr_overlap_mode == HALF_OVERLAP || st->tcx_cfg.tcx_curr_overlap_mode == MIN_OVERLAP) ? 2 : 1;

        target_bits -= st->tcxltp_bits;

        /* Run TCX20 encoder */

        coder_tcx( 0, &(st->tcx_cfg),
                   A_q_tcx, A_q_ind,
                   st->synth, st->L_frame, st->L_frameTCX, st->tcx_cfg.tcx_coded_lines, target_bits,
                   st->tcxonly, st->spectrum_long, &(st->LPDmem), param_core, st, &hm_cfg);

        coder_tcx_post( st, &(st->LPDmem), &(st->tcx_cfg), st->synth, A_q_tcx, Aw, st->wspeech_enc );


        st->plcExt.LPDmem = &(st->LPDmem);

        GplcTcxEncSetup( st, &st->plcExt );

        if( st->enableTcxLpc )
        {
            E_LPC_lsp_unweight(lsp_tcx_q, lsp_q, lsf_q, 1.0f/st->gamma ); /* Update lsf_q for encoderSideLossSimulation() */
        }
        encoderSideLossSimulation( st, &st->plcExt, lsf_q, stab_fac, 1, st->L_frame );

    }



    /* Update lsp/lsf memory */
    mvr2r( lsp_new, st->lspold_enc, M );

    if( st->enableTcxLpc && st->core != ACELP_CORE )
    {
        /* Update lsf / lsp memory */
        mvr2r(lsf_tcx_q, st->lsf_old, M);
        mvr2r(lsp_tcx_q, st->lsp_old, M);
        st->envWeighted = 1;

        /* Update ACELP quantizer state */
        lsf_update_memory( st->narrowBand, st->lsf_old, st->mem_MA, st->mem_MA );
        st->pstreaklen = 0;
        st->streaklimit = 1.0f;
        /* check resonance for pitch clipping algorithm */
        gp_clip_test_lsf( st->core_brate, st->lsf_old, st->clip_var, 0 );
        mvr2r( st->lsf_old, st->mem_AR, M );
    }
    else
    {
        mvr2r( lsf_q, st->lsf_old, M );
        mvr2r( lsp_q, st->lsp_old, M );
    }

    /* Update MODE1 CNG parameters */
    if( st->Opt_DTX_ON && vad_hover_flag )
    {
        st->burst_ho_cnt++;
        if( st->burst_ho_cnt > HO_HIST_SIZE )
        {
            st->burst_ho_cnt = HO_HIST_SIZE;
        }
    }
    else
    {
        st->burst_ho_cnt = 0;
    }

    if( st->Opt_DTX_ON )
    {
        /* update CNG parameters in active frames */
        if( st->bwidth == NB && st->enableTcxLpc && st->core != ACELP_CORE )
        {
            float buf[L_LP], res[L_FRAME], A[M+1], r[M+1], tmp, lsptmp[M];

            assert( st->L_frame==L_FRAME );

            mvr2r( st->synth+L_FRAME-L_LP, buf, L_LP );
            tmp = st->synth[L_FRAME-L_LP-1];
            preemph( buf, st->preemph_fac, L_LP, &tmp );
            autocorr( buf, r, M, L_LP, LP_assym_window, 0, 0, 0 );
            lag_wind( r, M, INT_FS_12k8, LAGW_WEAK );
            lev_dur( A, r, M, NULL );
            a2lsp_stab( A, lsptmp, lsp_new );

            residu( A, M, buf+L_LP-L_FRAME, res, L_FRAME );

            cng_params_upd( lsptmp, res, st->L_frame, &st->ho_circ_ptr, st->ho_ener_circ,
                            &st->ho_circ_size, st->ho_lsp_circ, ENC, st->ho_env_circ,
                            &st->cng_buf_cnt, st->cng_exc2_buf, st->cng_brate_buf, st->last_active_brate );
        }
        else
        {
            cng_params_upd( lsp_new, st->LPDmem.old_exc+L_EXC_MEM-st->L_frame, st->L_frame,
                            &st->ho_circ_ptr, st->ho_ener_circ, &st->ho_circ_size,
                            st->ho_lsp_circ, ENC, st->ho_env_circ, &st->cng_buf_cnt,
                            st->cng_exc2_buf, st->cng_brate_buf, st->last_active_brate );
        }

        if( st->L_frame == L_FRAME )
        {
            /* store LSPs@16k, potentially to be used in CNG@16k */
            mvr2r( st->lsp_old16k, &(st->ho_lsp_circ2[(st->ho_circ_ptr)*M]), M );
        }

        /* Set 16k LSP flag for CNG buffer */
        st->ho_16k_lsp[st->ho_circ_ptr] = (st->L_frame == L_FRAME ? 0 : 1 );

        /* efficient DTX hangover control */
        if ( st->burst_ho_cnt > 1 )
        {
            dtx_hangover_control( st, lsp_new );
        }
    }

    /*--------------------------------------------------------------*
    * Adaptive Bass Post-filter
    *---------------------------------------------------------------*/


    if( st->core > ACELP_CORE || st->rate_switching_reset )
    {
        /*TCX mode: copy values*/
        set_zero( st->mem_bpf, 2*L_FILT16k );           /*TCX->no gain*/
        set_zero( st->mem_error_bpf, 2*L_FILT16k );     /*TCX->no gain*/
        st->bpf_gain_param = 0;
    }
    else if( st->acelp_cfg.bpf_mode >= 1 )
    {
        /*ACELP: estimate bpf parameter with delay=0*/

        /*Estimate bpf parameter*/
        bass_pf_enc( st->speech_enc, st->synth, pitch_buf, gain_pitch_buf, st->L_frame, L_SUBFR, st->mem_bpf, st->mem_error_bpf,
                     &(st->bpf_gain_param), st->acelp_cfg.bpf_mode ,&(st->pst_lp_ener), &(st->pst_mem_deemp_err) );

    }


    /*--------------------------------------------------------------*
      * Analysis Print Out
      *---------------------------------------------------------------*/



    /*--------------------------------------------------------------*
    * Generate Bitstream
    *---------------------------------------------------------------*/

    enc_prm( coder_type, param_core, param_lpc, st, st->L_frame, &hm_cfg, bits_param_lpc, no_param_lpc );

    /* Channel-aware mode - encode partial copy */
    if( st->rf_mode )
    {
        set_f( lsf_q_1st_rf, 0.0f, M);

        if( st->core == ACELP_CORE )
        {
            /* convert lsp to lsf */
            lsp2lsf( lsp_new, lsf_uq_rf, M, st->sr_core );

            /* first stage VQ, 8 bits; reuse TCX high rate codebook */
            st->rf_indx_lsf[0][0] = vlpc_1st_cod( lsf_uq_rf, lsf_q_1st_rf, st->sr_core, w_rf );
            v_sub( lsf_uq_rf, lsf_q_1st_rf, lsf_q_d_rf, M );

            /* second stage vq */
            closest_centroid_rf( lsf_q_d_rf, w_rf, lsf_q_diff_cb_8b_rf, (1<<8), M, &st->rf_indx_lsf[0][1] );

            /* quantized lsf from two stages */
            v_add( lsf_q_1st_rf, lsf_q_diff_cb_8b_rf + M * st->rf_indx_lsf[0][1], lsf_q_rf, M );

            v_sort( lsf_q_rf, 0, M-1 );
            reorder_lsf( lsf_q_rf, LSF_GAP, M, st->sr_core );
        }
        else
        {
            rf_tcx_lpc_cdk = tcxlpc_get_cdk( GENERIC );

            /* Quantize */
            Q_lsf_tcxlpc( lsf, lsf_tcx_q, lspq_ind, param_lpc, st->narrowBand, rf_tcx_lpc_cdk, rf_mem_MA, GENERIC, st->Bin_E );

            /* VQ, 5+4+4 bits; reuse TCX low rate codebook */
            st->rf_indx_lsf[0][0] = param_lpc[1];
            st->rf_indx_lsf[0][1] = param_lpc[2];
            st->rf_indx_lsf[0][2] = param_lpc[3];
        }

        if (st->core == ACELP_CORE)
        {
            /* current n-th ACELP frame and its corresponding partial copy */
            lsf2lsp( lsf_q_rf, lsp_q_rf, M, st->sr_core );

            /* Interpolate LSPs and convert to LPC */
            int_lsp( st->L_frame, lsp_old_q_rf, lsp_q_rf, Aq_rf, M, interpol_frac_12k8, 0 );

            /* stability estimation */
            stab_fac_rf = lsf_stab( lsf_q_rf, lsf_old_q_rf, 0, st->L_frame );

            /* Configure partial copy estimation of the current n-th frame to be packed in future with n+fec_offset frame */
            /* o: rf_frame_type, o: rf_target_bits */
            BITS_ALLOC_ACELP_config_rf( coder_type, st->rf_tilt_buf, &st->rf_frame_type, &st->rf_target_bits, st->nb_subfr, st->rf_fec_indicator, pitch_buf );

            /* RF frame type in the buffer */
            st->rf_indx_frametype[0] = st->rf_frame_type;
            st->rf_targetbits_buff[0] = st->rf_target_bits;

            if( st->rf_frame_type != RF_NO_DATA )
            {
                /* coder_acelp_rf does the partial copy encoding based on the rf frame type chosen for the RF encoding */
                coder_acelp_rf( st->rf_target_bits, st->speech_enc_pe, coder_type, st->rf_frame_type, Aw, Aq_rf,
                                voicing, pitch, stab_fac_rf, st, &(st->acelp_cfg_rf), exc_rf, syn_rf );
            }
        }
        else
        {

            TD_Mode = 1;
            st->rf_clas[0] = st->clas;
            st->rf_gain_tcx[0] = param_core[0];

            /* attenuate somewhat the gain for onset when the correlation with previous frame is too low: avoid preecho */
            if( st->rf_gain_tcx[1]!= 0 && st->rf_gain_tcx[0] > 1.6*st->rf_gain_tcx[1] && st->tcxltp_gain <= 0.2 )
            {
                st->rf_gain_tcx[0] = 1.6*st->rf_gain_tcx[1];

                if( st->rf_gain_tcx[0] > 127 )
                {
                    st->rf_gain_tcx[0] = 127;
                }
            }

            /* get concealment decision*/
            rf_PLC_Mode = 0;
            if ((st->core == TCX_20_CORE)
                    &&(st->last_core == TCX_20_CORE)
                    && (st->rf_second_last_core == TCX_20_CORE)
                    && ((st->tcxltp_pitch_int <= 0.5f*st->L_frame) || (st->tcxltp_gain <= 0.4f))
                    && (st->tcxltp_pitch_int == st->rf_tcxltp_pitch_int_past)
                    && !st->rf_last_tns_active
                    && !st->rf_second_last_tns_active
                    && !(st->tcx_cfg.fIsTNSAllowed & st->fUseTns[0])
               )
            {
                rf_PLC_Mode = 1;
            }
            else if (st->last_core != 0)
            {
                if( (st->clas <= UNVOICED_TRANSITION || st->last_clas <= UNVOICED_TRANSITION || st->tcxltp_gain <= 0.4f) && st->last_core!=-1 )
                {
                    rf_PLC_Mode = st->last_core;
                }
            }

            /* call TD1 when the gain drop compare to previous frame*/
            if (rf_PLC_Mode == 0 && st->rf_gain_tcx[1]!= 0 &&
                    ((st->transientDetection.transientDetector.bIsAttackPresent && st->rf_gain_tcx[0] < 0.97*st->rf_gain_tcx[1]) ||
                     st->rf_gain_tcx[0] < 0.90*st->rf_gain_tcx[1]))
            {
                TD_Mode = 0;
            }
            else
            {
                TD_Mode = 1;
            }

            /* updates */
            st->rf_tcxltp_pitch_int_past  = st->tcxltp_pitch_int;
            st->rf_second_last_tns_active = st->rf_last_tns_active;
            st->rf_last_tns_active        = st->tcx_cfg.fIsTNSAllowed & st->fUseTns[0];
            st->rf_second_last_core       = st->last_core;

            st->rf_tcxltp_param[0] = st->tcxltp_param[1];

            /* Configure partial copy estimation of the current n-th frame to be packed in future with n+fec_offset frame */
            /* o: rf_frame_type, o: rf_target_bits */
            BITS_ALLOC_TCX_config_rf( &st->rf_frame_type, &st->rf_target_bits, rf_PLC_Mode, coder_type, st->last_core, TD_Mode );

            /* RF frame type in the buffer */
            st->rf_indx_frametype[0] = st->rf_frame_type;
            st->rf_targetbits_buff[0] = st->rf_target_bits;
        }
    }


    return;
}

/*-------------------------------------------------------------------*
* closest_centroid_rf()
*
* Determine a set of closest VQ centroids for a given input
*-------------------------------------------------------------------*/

static void closest_centroid_rf(
    const float *data,            /* i  : input data                        */
    const float *weights,         /* i  : weights                           */
    const float *quantizer,       /* i  : quantizer table                   */
    const short  centroids,       /* i  : number of centroids               */
    const short  length,          /* i  : dimension of quantiser            */
    short *ind_vec          /* o  : list of best match indice vectors */
)
{
    short i,j;
    float tmp, werr, best_werr;

    ind_vec[0] = 0;
    best_werr = 1.0E20f;

    for( i=0; i<centroids; i++ )
    {
        werr = 0.0f;
        for( j=0; j<length; j++ )
        {
            tmp = (float) *(data + j) - quantizer[i * length + j];
            werr += (float) (*(weights + j) * tmp * tmp);
        }

        if( werr < best_werr )
        {
            ind_vec[0] = i;
            best_werr = werr;
        }
    }
    return;
}


/*-------------------------------------------------------------------*
 * core_acelp_tcx20_switching()
 *
 * Open-loop ACELP/TCX20 core decision
 *-------------------------------------------------------------------*/

void core_acelp_tcx20_switching(
    Encoder_State *st,            /* i/o: encoder state structure             */
    const short vad_flag,
    short sp_aud_decision0,
    float non_staX,
    short *pitch,         /* i  : open-loop pitch values for quantiz. */
    float *pitch_fr,      /* i/o: fraction pitch values               */
    float *voicing_fr,    /* i/o: fractional voicing values           */
    const float currFlatness,   /* i  : flatness                            */
    const float lsp_mid[M],     /* i  : LSPs at the middle of the frame     */
    const float stab_fac        /* i  : LP filter stability                 */
)
{
    int i, j;
    float A_q_tcx[NB_SUBFR16k*(M+1)];
    float dsnr, snr_tcx, snr_acelp;
    int   iter;
    float xn_buf[L_MDCT_OVLP_MAX+L_FRAME_PLUS+L_MDCT_OVLP_MAX];
    float Ap[M+1];
    float gainlpc[FDNS_NPTS];
    float en[N_MAX/4];
    float sqGain, ener, tmp, fac, offset;
    int L_frame = st->L_frame;
    int overlap;
    int tcx_offset = st->tcx_cfg.tcx_offset;
    float *x = st->spectrum_long;
    float target;
    int T0;
    float gain, noise, scale;
    float *pt_ener_sfr, ener_sfr[NB_SUBFR16k];

    /* Check minimum pitch for quantization */
    for( i = 0; i < 3; i++ )
    {
        /* check minimum pitch for quantization */
        if( pitch[i] < PIT_MIN_SHORTER )
        {
            pitch[i] *= 2;
        }

        /* convert pitch values to 16kHz domain */
        if ( st->L_frame == L_FRAME16k )
        {
            pitch[i] = (short)(pitch[i] * 1.25f + 0.5f);
        }
    }
    if( st->narrowBand == 1 )
    {
        pitchDoubling_det( st->wspeech_enc, pitch, pitch_fr, voicing_fr );
    }

    lsp2a_stab(lsp_mid, A_q_tcx, M);

    tcx_ltp_encode( st->tcxltp, st->tcxonly, TCX_20, st->L_frame, L_SUBFR, st->speech_enc+st->encoderLookahead_enc, st->speech_ltp+st->encoderLookahead_enc,
                    st->wspeech_enc+st->encoderLookahead_enc, pitch[1], st->tcxltp_param, &st->tcxltp_bits, &st->tcxltp_pitch_int, &st->tcxltp_pitch_fr,
                    &st->tcxltp_gain, &st->tcxltp_pitch_int_past, &st->tcxltp_pitch_fr_past, &st->tcxltp_gain_past, &st->tcxltp_norm_corr_past, st->last_core,
                    st->pit_min, st->pit_fr1, st->pit_fr2, st->pit_max, st->pit_res_max, &st->transientDetection, 0, A_q_tcx, M );

    /* Force TCX when TCX20 in MODE1 is selected */
    if( st->mdct_sw == MODE1 )
    {
        st->core = TCX_20_CORE;
    }
    else
    {
        /*--------------------------------------------------------------*
         * Estimate TCX SNR
         *---------------------------------------------------------------*/

        target = 1000.f;
        if ( st->sr_core == 16000.f ) target = 850.f;
        if ( st->sr_core == 12800.f ) target = 850.f;
        if ( st->narrowBand == 1 ) target = 500.f;

        if( st->last_core == ACELP_CORE )
        {
            L_frame += tcx_offset;

            if( st->tcx_cfg.lfacNext < 0 )
            {
                L_frame -= st->tcx_cfg.lfacNext;
                tcx_offset = st->tcx_cfg.lfacNext;
            }
            else
            {
                tcx_offset = 0;
            }
        }

        overlap = st->tcx_cfg.tcx_mdct_window_delay;

        mvr2r(st->speech_ltp-(overlap>>1)+tcx_offset, xn_buf, L_frame+overlap);

        if( st->last_core == ACELP_CORE )
        {
            if( tcx_offset < 0 )
            {
                set_f( xn_buf, 0.0f, overlap>>1 );
            }
        }
        else
        {
            for (i = 0; i < overlap; i++)
            {
                xn_buf[i] *= st->tcx_cfg.tcx_mdct_window[i];
            }
        }

        for (i = 0; i < overlap; i++)
        {
            xn_buf[L_frame+i] *= st->tcx_cfg.tcx_mdct_window[overlap-1-i];
        }

        TCX_MDCT( xn_buf, x, overlap, L_frame-overlap, overlap );

        for( i = 0; i < L_frame; i++ )
        {
            x[i] *= (float)L_frame / sqrt(2*NORM_MDCT_FACTOR);
        }

        weight_a( A_q_tcx, Ap, st->gamma, M );

        lpc2mdct( Ap, M, gainlpc );

        mdct_preShaping( x, L_frame, gainlpc );

        if( st->narrowBand == 1 )
        {
            j = (int)( (float)L_frame*0.625f );

            set_f( x + j, 0.0f, L_frame - j );
        }

        for( i = 0; i < L_frame; i+=4 )
        {
            ener = 0.01f + x[i]*x[i] + x[i+1]*x[i+1] + x[i+2]*x[i+2] + x[i+3]*x[i+3];
            en[i/4] = 9.0f + 10.0f*(float)log10(ener);
        }

        fac = 128.0f;
        offset = fac;

        for( iter = 0; iter < 10; iter++ )
        {
            fac *= 0.5f;
            offset -= fac;
            ener = 0.0f;

            for( i=0; i<L_frame/4; i++ )
            {
                tmp = en[i] - offset;

                if( tmp > 3.0f )
                {
                    ener += tmp;
                }

                if( ener > target )
                {
                    offset += fac;
                    break;
                }
            }
        }

        if( offset <= 32.f )
        {
            offset = -128.f;
        }

        sqGain = (float)pow(10.0f, offset/20.0f);
        ener = sqGain*sqGain/12.f*(float)sqrt(2.f)/(float)L_frame;

        snr_tcx = 0.0f;
        pt_ener_sfr = ener_sfr;

        for( i = 0; i < st->L_frame; i += L_SUBFR )
        {
            *pt_ener_sfr = sum2_f( st->wspeech_enc + i, L_SUBFR ) + 1e-6f;

            snr_tcx += (float)log10( *pt_ener_sfr/(ener*L_SUBFR) );
            pt_ener_sfr++;
        }
        snr_tcx *= ((float)(10*L_SUBFR))/(float)st->L_frame;


        /*--------------------------------------------------------------*
        * Estimate ACELP SNR
        *---------------------------------------------------------------*/

        scale = 0.055f;
        if ( st->sr_core == 16000 ) scale = 0.092f;
        if ( st->sr_core == 12800 ) scale = 0.059f;
        if ( st->narrowBand ) scale = 0.15f;

        snr_acelp = 0.0f;
        fac = (float)st->sr_core/12800.f;
        pt_ener_sfr = ener_sfr;

        for( i = 0; i < st->L_frame; i += L_SUBFR )
        {
            T0 = (int)( (fac*pitch_fr[(int)((float)(i/L_SUBFR)/fac+0.5f)]) + 0.5f );
            gain = get_gain( st->wspeech_enc+i, st->wspeech_enc+i-T0, L_SUBFR, NULL );

            noise = 1e-6f;
            for( j = 0; j < L_SUBFR; j++ )
            {
                tmp = st->wspeech_enc[i+j] - gain * st->wspeech_enc[i+j-T0];
                noise += tmp * tmp;
            }

            noise *= scale;
            snr_acelp += (float)log10( *pt_ener_sfr/noise );
            pt_ener_sfr++;
        }

        snr_acelp *= ((float)(10*L_SUBFR))/(float)st->L_frame;



        /*--------------------------------------------------------------*
        * Switching Decision
        *---------------------------------------------------------------*/

        dsnr = 0.0f;
        /* hysteresis for very small SNR differences between ACELP and TCX */

        /* try to use TCX instead of ACELP on temporally stationary frames */
        if( (snr_acelp > snr_tcx) &&
                (snr_acelp < snr_tcx + 2.0f) &&
                (st->prevTempFlatness + currFlatness < 3.25f || stab_fac == 1.0f || (st->sr_core == 12800 && sp_aud_decision0 == 1 && st->prevTempFlatness + currFlatness < 20.f)) &&
                (st->acelpFramesCount <= 6))
        {
            dsnr = -2.0f;
        }

        /* try to use ACELP instead of TCX on transient and "buzzy" frames */
        if( (snr_acelp < snr_tcx) &&
                (snr_acelp > snr_tcx - 2.0f) &&
                (st->prevTempFlatness + currFlatness > 3.25f) &&
                (st->acelpFramesCount >= 6))
        {
            dsnr =  2.0f;
        }

        if( (st->sr_core == 12800) && (offset < 74.0f) && (non_staX > 5.0f) && (snr_acelp >= snr_tcx - 4) && st->acelpFramesCount >= 1 && (((st->lps > st->lpm) && mean(voicing_fr, 4) >= 0.3f) || (st->acelpFramesCount >= 6 && (st->lps > st->lpm - 1.5f))) && (sp_aud_decision0 == 0) && vad_flag )
        {
            /* Fine tuned across various databases based on various metrics to detect TCX frames in speech.*/
            dsnr = 4.0f;
        }

        if( st->flag_noisy_speech_snr )
        {
            if( vad_flag || st->Opt_DTX_ON )
            {
                dsnr += 2.f;
            }
            else
            {
                dsnr -= 2.f;
            }
        }

        if (st->sr_core == 12800 && (non_staX < 2.f || (st->flag_noisy_speech_snr==0&&vad_flag==1&&offset==-128.f&&st->acelpFramesCount>=6)) && (st->last_core==ACELP_CORE||st->last_core==TCX_20_CORE))
        {
            st->core = st->last_core;
        }
        else if( snr_acelp + dsnr > snr_tcx )
        {
            st->core = ACELP_CORE;
            st->acelpFramesCount = min(32767-1, st->acelpFramesCount + 1);
        }
        else
        {
            st->core = TCX_20_CORE;
            st->acelpFramesCount = 0;
        }
    }

    /* Fixed Decision (using -C) */
    if( st->acelpEnabled == 1 && st->tcx20Enabled == 0 )
    {
        st->core = ACELP_CORE;
    }

    if( st->acelpEnabled == 0 && st->tcx20Enabled == 1 )
    {
        st->core = TCX_20_CORE;
    }


    st->prevTempFlatness = currFlatness;

    return;
}

/*-------------------------------------------------------------------*
 * BITS_ALLOC_ACELP_config_rf()
 *
 * configure channel aware mode
 *-------------------------------------------------------------------*/
static void BITS_ALLOC_ACELP_config_rf(
    const short coder_type,
    float *tilt_code,
    short *rf_frame_type,
    short *rf_target_bits,
    short nb_subfr,
    short rf_fec_indicator,
    float *pitch_buf
)
{
    float mean_tc, min_tilt_code, max_tilt_code;
    short nrgMode, ltfMode, ltpMode, gainsMode;

    short en_partial_red = 1;
    float dpit1, dpit2, dpit3;

    /* Init */
    *rf_target_bits = 0;

    /* ----------------------------------------*
     * RF frame type selection                 *
     *-----------------------------------------*/

    /* Mean tilt code estimation */
    mean_tc = 0;
    mean_tc = mean(tilt_code, nb_subfr);

    /* Maximum tilt code estimation */
    max_tilt_code = tilt_code[0];
    maximum(tilt_code, nb_subfr, &max_tilt_code);

    /* Minimum tilt code estimation */
    min_tilt_code=tilt_code[0];
    minimum(tilt_code, nb_subfr, &min_tilt_code);

    /* ----------------------------------------*/
    /* Decide Criticality                      */
    /*-----------------------------------------*/
    dpit1 = (float)fabs( pitch_buf[0] - pitch_buf[1] );
    dpit2 = (float)fabs( pitch_buf[1] - pitch_buf[2] );
    dpit3 = (float)fabs( pitch_buf[2] - pitch_buf[3] );

    if ( rf_fec_indicator == 1)
    {
        if  ( max_tilt_code > 0.48f && dpit1 <= 0.0f && dpit2 <= 0.0f && dpit3 <= 0.0f && coder_type == VOICED )
        {
            en_partial_red = 0;
        }
        if  ( max_tilt_code > 0.47f  && dpit1 <= 1.0f && dpit2 <= 1.0f && dpit3 <= 1.0f && coder_type == GENERIC)
        {
            en_partial_red = 0;
        }
    }
    else
    {
        if  ( max_tilt_code > 0.47 && dpit1 <= 0.25f && dpit2 <= 0.25f && dpit3 <= 0.25f && coder_type == VOICED )
        {
            en_partial_red = 0;
        }
        if  ( max_tilt_code > 0.45  && dpit1 <= 1.25f && dpit2 <= 1.25f && dpit3 <= 1.25f && coder_type == GENERIC)
        {
            en_partial_red = 0;
        }
    }


    /* ---------------------------------------------------------*
     * Identify number of bits required as per rf frame type    *
     * ---------------------------------------------------------*/

    /* rf_mode, 1 bit */
    *rf_target_bits += 1;

    /* rf_fec_offset 2 bits */
    *rf_target_bits += 2;

    /* rf_frame_type, 3 bits */
    *rf_target_bits += 3;

    /* LSF bits 8 + 8 bits */
    *rf_target_bits += 16;

    /* Intialize the RF mode frame type to all-pred */
    *rf_frame_type = RF_ALLPRED;

    if ( coder_type == INACTIVE || en_partial_red == 0 )
    {
        *rf_frame_type = RF_NO_DATA;
    }
    else if ( coder_type == UNVOICED || coder_type == INACTIVE )
    {
        *rf_frame_type = RF_NELP;
    }
    else if( ( coder_type == GENERIC ) && max_tilt_code < 0.05f)
    {
        *rf_frame_type = RF_NOPRED;
    }
    else if( ( coder_type == GENERIC ) && mean_tc < 0.3f)
    {
        *rf_frame_type = RF_GENPRED;
    }

    nrgMode = ACELP_NRG_MODE[1][1][*rf_frame_type];
    ltfMode = ACELP_LTF_MODE[1][1][*rf_frame_type];
    ltpMode = ACELP_LTP_MODE[1][1][*rf_frame_type];
    gainsMode = ACELP_GAINS_MODE[1][1][*rf_frame_type];

    /* Number of RF bits for different RF coder types */
    switch (*rf_frame_type)
    {
    case RF_ALLPRED:
        /* Es_pred bits 3 bits, LTF: 1, pitch: 8,5,5,5, FCB: 0, gain: 7,0,7,0, Diff GFr: 4*/
        *rf_target_bits += (  ACELP_NRG_BITS[nrgMode]
                              + ACELP_LTF_BITS[ltfMode]
                              + ACELP_LTP_BITS_SFR[ltpMode][0] + ACELP_LTP_BITS_SFR[ltpMode][1] + ACELP_LTP_BITS_SFR[ltpMode][2] + ACELP_LTP_BITS_SFR[ltpMode][3]
                              + ACELP_GAINS_BITS[gainsMode] + ACELP_GAINS_BITS[gainsMode]
                              + 2 /*2 bits for PartialCopy GainFrame*/ );
        break;

    case RF_NOPRED:
        /* Es_pred bits 3 bits, LTF: 0, pitch: 0, FCB: 7,7,7,7, gain: 6,0,6,0, Diff GFr: 2*/
        /*bits += (3 + 0 + 0 + 28 + 12 + 2); */ /* 64 rf bits */
        *rf_target_bits += (  ACELP_NRG_BITS[nrgMode]
                              + ACELP_LTF_BITS[ltfMode]
                              + 28
                              + ACELP_GAINS_BITS[gainsMode] + ACELP_GAINS_BITS[gainsMode]
                              + 2 /*2 bits for PartialCopy GainFrame*/ );
        break;

    case RF_GENPRED:
        /* Es_pred bits 3 bits, LTF: 0, pitch: 8,0,8,0, FCB: 6,7,5,5, gain: 5,0,5,0, Diff GFr: 0*/
        /*bits += (3 + 0 + 16 + 23 + 10 + 0);  */   /* 72 rf bits */
        *rf_target_bits += (  ACELP_NRG_BITS[nrgMode]
                              + ACELP_LTF_BITS[ltfMode]
                              + ACELP_LTP_BITS_SFR[ltpMode][0] + ACELP_LTP_BITS_SFR[ltpMode][1] + ACELP_LTP_BITS_SFR[ltpMode][2] + ACELP_LTP_BITS_SFR[ltpMode][3]
                              + 14
                              + ACELP_GAINS_BITS[gainsMode] + ACELP_GAINS_BITS[gainsMode]
                              + 2 /*2 bits for PartialCopy GainFrame*/ );
        break;

    case RF_NELP:
        /* gain: 19, Diff GFr: 5 */
        /*bits += (19 + 5);    */
        *rf_target_bits +=  (19 + NUM_BITS_SHB_FRAMEGAIN);
        break;

    case RF_NO_DATA:
        *rf_target_bits  = 6;
        break;
    default:
        assert(!"RF_Frame_type does not belong to ACELP Partial copy frame types possible!");
        break;
    }

    return;

}

/*-------------------------------------------------------------------*
  * BITS_ALLOC_TCX_config_rf()
  *
  * configure channel aware mode
  *-------------------------------------------------------------------*/
static void BITS_ALLOC_TCX_config_rf(
    short *rf_frame_type,
    short *rf_target_bits,
    short PLC_Mode,
    short coder_type,
    short last_core,
    int   TD_Mode
)
{
    /* Init: rf_mode + rf_fec_offset + rf_frame_type */
    *rf_target_bits = 1 + 2 + 3;

    if( coder_type == INACTIVE || last_core == ACELP_CORE )
    {
        *rf_frame_type = RF_NO_DATA;
    }
    else
    {
        /* classification */
        *rf_target_bits += 2;

        if( PLC_Mode )
        {
            /* TCX global gain  = 7 bits */
            *rf_target_bits += 7;
            *rf_frame_type = RF_TCXFD;
        }
        else
        {
            /* pitch and gain */
            /* LTP data */
            if( TD_Mode )
            {
                *rf_target_bits += 9;
                *rf_frame_type = RF_TCXTD2;
            }
            else
            {
                *rf_target_bits += 9;
                *rf_frame_type = RF_TCXTD1;
            }
        }

        if( *rf_frame_type == RF_TCXFD )
        {
            /* TCXFD: LSF bits 5 + 4 + 4 bits     */
            /* only embed LSF for FD concealment */
            *rf_target_bits += TCXLPC_NUMBITS;
        }
    }
    return;
}
