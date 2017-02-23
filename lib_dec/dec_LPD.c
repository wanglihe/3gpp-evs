/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"
#include "cnst.h"
#include "basop_proto_func.h"
#include "stat_com.h"

#if defined(_MSC_VER) && (_MSC_VER <= 1200)  /* disable global optimizations to overcome an internal compiler error */
#pragma optimize("g", off)
#endif

/*-------------------------------------------------------------------*
* decoder_LPD()
*
* Core decoder
*--------------------------------------------------------------------*/

void decoder_LPD(
    float signal_out[],     /* output: signal with LPD delay (7 subfrs) */
    float signal_outFB[],
    short *total_nbbits,    /* i/o:    number of bits / decoded bits    */
    Decoder_State *st,              /* i/o:    decoder memory state pointer     */
    float *bpf_noise_buf,   /* i/o: BPF noise buffer                    */
    short bfi,              /* i  : BFI flag                            */
    short *bitsRead,        /* o  : number of read bits                 */
    short *coder_type,      /* o  : coder type                          */
    int   param[],          /* o  : buffer of parameters                */
    float *pitch_buf,       /* i/o: floating pitch values for each subfr*/
    float *voice_factors,   /* o  : voicing factors                     */
    float *ptr_bwe_exc      /* o  : excitation for SWB TBE              */
)
{
    int   *prm, param_lpc[NPRM_LPC_NEW];
    float synth_buf[OLD_SYNTH_SIZE_DEC+L_FRAME_PLUS+M];
    float *synth;
    float synth_bufFB[OLD_SYNTH_SIZE_DEC+L_FRAME_PLUS+M];
    float *synthFB;
    float lsf[(NB_DIV+1)*M], lsp[(NB_DIV+1)*M], lspmid[M], lsfmid[M];
    float Aq[(NB_SUBFR16k+1)*(M+1)];
    int   pitch[NB_SUBFR16k];
    float pit_gain[NB_SUBFR16k];
    short i, k;
    short  past_core_mode;
    short L_frame, nb_subfr, L_frameTCX;
    Word16 Aind[M+1], lspind[M];
    float tmp_old[M+1], tmp_new[M+1], enr_old, enr_new;
    float lspnew_uw[NB_DIV*M], lsfnew_uw[NB_DIV*M];
    float lsf_q_1st_rf[M], lsf_q_rf[M], lsp_q_rf[M];
    float lsp_diff;
    short LSF_Q_prediction;  /* o  : LSF prediction mode                 */
    short tcx_last_overlap_mode, tcx_current_overlap_mode;


    /*--------------------------------------------------------------------------------*
     * Initializations
     *--------------------------------------------------------------------------------*/

    prm = param;  /* to avoid compilation warnings */
    LSF_Q_prediction = -1;  /* to avoid compilation warnings */
    enr_old = 0.0f;
    enr_new = 0.0f;

    if ( st->use_partial_copy && st->rf_frame_type >= RF_TCXFD && st->rf_frame_type <= RF_TCXTD2 )
    {
        bfi = st->bfi;
    }

    past_core_mode = st->last_core_bfi;

    /*Adjust bit per frame*/
    if( !bfi )
    {
        st->bits_frame_core = st->bits_frame - (*bitsRead);
    }

    /* Framing parameters */
    L_frame  = st->L_frame;
    L_frameTCX = st->L_frameTCX;
    nb_subfr = st->nb_subfr;

    /* Initialize pointers */
    synth = synth_buf + st->old_synth_len;
    synthFB = synth_bufFB + st->old_synth_lenFB;
    mvr2r( st->old_synth, synth_buf, st->old_synth_len );
    mvr2r( st->old_synthFB, synth_bufFB, st->old_synth_lenFB );
    set_zero( synth, L_FRAME_PLUS+M );
    set_zero( synthFB, L_FRAME_PLUS+M );


    /*For post-processing (post-filtering+blind BWE)*/
    if( st->tcxonly == 0 )
    {
        /* for bass postfilter */
        set_i( pitch, L_SUBFR, nb_subfr );
        set_zero( pit_gain, nb_subfr );
    }

    /* PLC: [Common: Memory update]
     * PLC: Update the number of lost frames */
    if( bfi )
    {
        st->nbLostCmpt++;
    }


    /*--------------------------------------------------------------------------------*
     * BITSTREAM DECODING
     *--------------------------------------------------------------------------------*/


    if( !bfi )
    {
        st->second_last_core = st->last_core;
        tcx_last_overlap_mode = st->tcx_cfg.tcx_last_overlap_mode;
        tcx_current_overlap_mode = st->tcx_cfg.tcx_curr_overlap_mode;

        dec_prm( &(st->core), &(st->last_core), coder_type, param, param_lpc, total_nbbits, st, L_frame, bitsRead );

        if(!st->rate_switching_init && (st->last_codec_mode) == MODE2 && st->BER_detect)
        {
            *coder_type= st->last_coder_type;
            st->last_core = st->second_last_core;
            st->tcx_cfg.tcx_last_overlap_mode = tcx_last_overlap_mode;
            st->tcx_cfg.tcx_curr_overlap_mode = tcx_current_overlap_mode;
            st->bfi = 1;
            bfi = 1;
            st->flagGuidedAcelp = 0;
            st->nbLostCmpt++;
            st->core_brate = st->last_core_brate;
            st->core = GetPLCModeDecision( st );
        }
    }
    else
    {
        if( st->use_partial_copy && st->rf_frame_type >= RF_TCXFD && st->rf_frame_type <= RF_TCXTD2 )
        {
            dec_prm( &(st->core), &(st->last_core), coder_type, param, param_lpc, total_nbbits, st, L_frame, bitsRead );
        }

        if( st->nbLostCmpt > 1 )
        {
            st->flagGuidedAcelp = 0;
        }
        /* PLC: [Common: mode decision]
         * PLC: Decide which Concealment to use. Update pitch lags if needed */
        st->core = GetPLCModeDecision( st );
    }

    /* PLC: [Common: Memory update]
     * PLC: Update the number of lost frames */
    if( !bfi )
    {
        if( st->prev_bfi == 1 )
        {
            st->prev_nbLostCmpt = st->nbLostCmpt;
        }
        else
        {
            st->prev_nbLostCmpt = 0;
        }

        st->nbLostCmpt = 0;
    }


    /*--------------------------------------------------------------------------------*
     * LPC PARAMETERS
     *--------------------------------------------------------------------------------*/


    if( (bfi == 0 ) || ( bfi == 1 && st->use_partial_copy && st->rf_frame_type == RF_TCXFD) )
    {
        if(st->use_partial_copy && ( st->rf_frame_type < RF_TCXFD || st->rf_frame_type > RF_TCXTD2))
        {
            if( st->envWeighted )
            {
                mvr2r( st->lspold_uw, st->lsp_old, M );
                mvr2r( st->lsfold_uw, st->lsf_old, M );
                st->envWeighted = 0;
            }

            /* first stage VQ, 8 bits; reuse TCX high rate codebook */
            set_f(lsf_q_1st_rf, 0.0f, M);
            vlpc_1st_dec(param_lpc[0], lsf_q_1st_rf, st->sr_core );

            /* second stage vq */
            /* quantized lsf from two stages */
            v_add( lsf_q_1st_rf, lsf_q_diff_cb_8b_rf + M * param_lpc[1], lsf_q_rf, M );

            v_sort( lsf_q_rf, 0, M-1 );
            reorder_lsf( lsf_q_rf, LSF_GAP, M, st->sr_core );

            /* current n-th ACELP frame and its corresponding partial copy */
            lsf2lsp( lsf_q_rf, lsp_q_rf, M, st->sr_core );

            /* copy the old and current lsfs and lsps into the lsf[] and lsp[] buffer for interpolation */
            mvr2r(st->lsf_old, &lsf[0], M);
            mvr2r(st->lsp_old, &lsp[0], M);
            mvr2r(lsf_q_rf, &lsf[M], M);
            mvr2r(lsp_q_rf, &lsp[M], M);
            lsp_diff = 0.0f;

            for( i = 0; i < M; i++ )
            {
                lsp_diff += (float)fabs(lsp[i+M] - lsp[i]);
            }

            if( st->core == ACELP_CORE && st->last_core == ACELP_CORE
                    && lsp_diff < 1.6f && lsp_diff > 0.12f && st->next_coder_type == GENERIC
                    && !st->prev_use_partial_copy && st->last_coder_type == UNVOICED && st->rf_frame_type >= RF_GENPRED)
            {
                mvr2r( &lsp[0], &lsp[M], M );
            }

            /* update mem_MA and mem_AR memories */
            lsf_update_memory( st->narrowBand, &lsf[M], st->mem_MA, st->mem_MA );
            mvr2r(&lsf[M], st->mem_AR, M);

            for( k=0; k<st->numlpc; ++k )
            {
                mvr2r( &lsp[(k+1)*M], &lspnew_uw[k*M], M );
                mvr2r( &lsf[(k+1)*M], &lsfnew_uw[k*M], M );
            }
        }
        else if( (st->enableTcxLpc && st->core != ACELP_CORE) || (bfi && st->use_partial_copy && st->rf_frame_type == RF_TCXFD) )
        {
            int tcx_lpc_cdk;

            if( bfi && st->use_partial_copy && st->rf_frame_type == RF_TCXFD )
            {
                tcx_lpc_cdk = tcxlpc_get_cdk(GENERIC);
            }
            else
            {
                tcx_lpc_cdk = tcxlpc_get_cdk(*coder_type);
            }

            mvr2r( st->lsf_old, &lsf[0], M );
            mvr2r( st->lsp_old, &lsp[0], M );

            D_lsf_tcxlpc( param_lpc, &lsf[M], lspind, st->narrowBand, tcx_lpc_cdk, st->mem_MA );

            lsf2lsp( &lsf[M], &lsp[M], M, st->sr_core );

            lsf_update_memory( st->narrowBand, &lsf[M], st->mem_MA, st->mem_MA );
            mvr2r(&lsf[M], st->mem_AR, M);
            st->envWeighted = 1;

            E_LPC_lsp_unweight( &lsp[M], lspnew_uw, lsfnew_uw, 1.0f/st->gamma );
        }
        else
        {
            if( st->envWeighted )
            {
                mvr2r( st->lspold_uw, st->lsp_old, M );
                mvr2r( st->lsfold_uw, st->lsf_old, M );
                st->envWeighted = 0;
            }

            /* Unquantize LPC */
            if( st->core == TCX_20_CORE )
            {
                lpc_unquantize( st, st->lsf_old, st->lsp_old, lsf, lsp, st->lpcQuantization, param_lpc, st->numlpc, st->core, st->mem_MA,
                                lspmid, lsfmid, AUDIO, st->acelp_cfg.midLpc, st->narrowBand, &(st->seed_acelp), st->sr_core,
                                &st->mid_lsf_int, st->prev_bfi, &LSF_Q_prediction, &st->safety_net );
            }
            else
            {
                lpc_unquantize( st, st->lsf_old, st->lsp_old, lsf, lsp, st->lpcQuantization, param_lpc, st->numlpc, st->core, st->mem_MA,
                                lspmid, lsfmid, *coder_type, st->acelp_cfg.midLpc, st->narrowBand, &(st->seed_acelp), st->sr_core,
                                &st->mid_lsf_int, st->prev_bfi, &LSF_Q_prediction, &st->safety_net );

                if( st->prev_use_partial_copy && st->last_core == ACELP_CORE && st->core == ACELP_CORE && st->prev_rf_frame_type >= RF_GENPRED && *coder_type == UNVOICED )
                {
                    if( st->lpcQuantization && st->acelp_cfg.midLpc )
                    {
                        mvr2r( lspmid, &lsp[0], M );
                        mvr2r( &lsp[M], lspmid, M );
                    }
                }
            }

            for( k=0; k<st->numlpc; ++k )
            {
                mvr2r( &lsp[(k+1)*M], &lspnew_uw[k*M], M );
                mvr2r( &lsf[(k+1)*M], &lsfnew_uw[k*M], M );
            }
        }

        /* PLC: [LPD: LPC concealment] built the moving average for the LPC concealment */
        for( k=0; k<st->numlpc; k++ )
        {
            for( i=0; i<M; i++ )
            {
                st->lsf_adaptive_mean[i] = (st->lsfoldbfi1[i]+st->lsfoldbfi0[i]+lsfnew_uw[k*M+i])/3;
                st->lsfoldbfi1[i] = st->lsfoldbfi0[i];
                st->lsfoldbfi0[i] = lsfnew_uw[k*M+i];
            }
        }
    }
    else
    {
        /* PLC: [LPD: LPC concealment] Conceal the LPC from the lost frame */
        float const* lsfBase;                      /* base for differential lsf coding */

        if( st->tcxonly == 0 || st->core < TCX_10_CORE )
        {
            st->numlpc = 1;
        }
        else
        {
            st->numlpc = 2;
        }

        if( st->nbLostCmpt == 1)
        {
            mvr2r( st->lsf_old, st->old_lsf_q_cng, M );
            mvr2r( st->lsp_old, st->old_lsp_q_cng, M );
        }


        lsfBase = PlcGetlsfBase( st->lpcQuantization, st->narrowBand, st->sr_core );

        dlpc_bfi( st->L_frame, lsfnew_uw, st->lsfold_uw, st->last_good, st->nbLostCmpt, st->mem_MA, st->mem_AR, &(st->stab_fac),
                  st->lsf_adaptive_mean, st->numlpc, st->lsf_cng, st->plcBackgroundNoiseUpdated, st->lsf_q_cng, st->old_lsf_q_cng, lsfBase);

        st->envWeighted = 0;

        mvr2r( st->lspold_uw, lsp, M );
        mvr2r( st->lsfold_uw, lsf, M );

        for( k = 0; k < st->numlpc; k++ )
        {
            mvr2r( &lsfnew_uw[k*M], &lsf[(k+1)*M], M );

            lsf2lsp( &lsf[(k+1)*M], &lsp[(k+1)*M], M, st->sr_core );
            lsf2lsp( st->lsf_q_cng, st->lsp_q_cng, M, st->sr_core );

            mvr2r( &lsp[(k+1)*M], &lspnew_uw[k*M], M );
        }
    }


    /*--------------------------------------------------------------*
      * Rate switching
      *---------------------------------------------------------------*/
    if( st->rate_switching_reset )
    {
        mvr2r( &(lsf[M]),&(lsf[0]), M );
        mvr2r( &(lsp[M]),&(lsp[0]), M );
        mvr2r( &(lsf[M]),st->lsf_old, M );
        mvr2r( &(lsp[M]),st->lsp_old, M );
        mvr2r( &(lsf[M]),lsfmid, M );
        mvr2r( &(lsp[M]),lspmid, M );
        lsp2a_stab( st->lsp_old, st->old_Aq_12_8, M );
    }

    if(st->enablePlcWaveadjust)
    {
        if (st->core == ACELP_CORE)
        {
            st->tonality_flag = 0;
        }
        if (bfi)
        {
            st->plcInfo.nbLostCmpt++;
        }
    }

    /*--------------------------------------------------------------------------------*
     * ACELP
     *--------------------------------------------------------------------------------*/


    if( (st->prev_bfi!=0) && (bfi==0) && (*coder_type==VOICED) && (st->prev_nbLostCmpt>4) )
    {
        st->dec_glr_idx = 1;
        st->reset_mem_AR = 1;
    }

    if( st->core == ACELP_CORE )
    {
        if( !st->tcxonly )
        {
            /* Set pointer to parameters */
            prm = param;

            /* Stability Factor */
            if( !bfi )
            {
                st->stab_fac = lsf_stab( &lsf[M], &lsf[0], 0, st->L_frame );
            }

            if( !bfi && st->prev_bfi )
            {
                /* check if LSP interpolation can be relaxed or if LPC power can be diffused*/
                lsp2a_stab( &lsp[0], tmp_old, M);
                enr_old = enr_1_Az( tmp_old, 2*L_SUBFR );

                lsp2a_stab( &lsp[M], tmp_new, M);
                enr_new = enr_1_Az( tmp_new, 2*L_SUBFR );
            }

            if( !bfi && (st->dec_glr_idx == 1 || (!(st->safety_net) && enr_new>=256.f && enr_new > 2*enr_old)) && st->prev_bfi  )
            {
                RecLpcSpecPowDiffuseLc( &lsp[M], &lsp[0], &lsf[M], st, st->dec_glr_idx == 1 ? 1 : 0 );
                int_lsp( L_frame, &lsp[0], &lsp[M], Aq, M, interpol_frac_12k8, 0 );
                mvr2r( &lsf[M], lsfnew_uw, M );
            }
            else
            {
                /* LPC Interpolation for ACELP */
                if( !bfi && st->acelp_cfg.midLpc )
                {
                    st->relax_prev_lsf_interp = 0;

                    if( st->prev_bfi )
                    {
                        /* check if LSP interpolation can be relaxed */
                        if ( enr_new < (0.25f * enr_old) )
                        {
                            /* don't use safety_net as this is getting impacted with lsf_restruct */
                            if ( (st->clas_dec == UNVOICED_CLAS) || (st->clas_dec == SIN_ONSET) || (st->clas_dec == INACTIVE_CLAS) || (*coder_type == GENERIC) || (*coder_type == TRANSITION) )
                            {
                                st->relax_prev_lsf_interp = 1;
                            }
                            else
                            {
                                st->relax_prev_lsf_interp = -1;
                            }
                        }
                    }

                    if( st->stab_fac == 0 && st->old_bfi_cnt > 0 && st->clas_dec != VOICED_CLAS && st->clas_dec != ONSET && st->relax_prev_lsf_interp == 0 )
                    {
                        st->relax_prev_lsf_interp = 2;
                    }
                    int_lsp4( L_frame, &lsp[0], lspmid, &lsp[M], Aq, M, st->relax_prev_lsf_interp );
                }
                else
                {
                    int_lsp( L_frame, &lsp[0], &lsp[M], Aq, M, interpol_frac_12k8, 0 );
                    int_lsp( L_frame, st->old_lsp_q_cng, st->lsp_q_cng, st->Aq_cng, M, interpol_frac_12k8, 0 );
                }
            }
        }

        if( bfi && st->last_core != ACELP_CORE )
        {
            /* PLC: [TCX: TD PLC] */
            con_tcx( st, &synthFB[0] );
            lerp(synthFB, synth, st->L_frame, st->L_frameTCX);
            st->con_tcx = 1;
            set_f( &st->mem_pitch_gain[2], st->lp_gainp, st->nb_subfr);
        }
        else
        {
            /* ACELP decoder */
            if (st->L_frame == L_FRAME)
            {
                mvr2r(Aq+2*(M+1), st->cur_sub_Aq, (M+1));
            }
            else
            {
                mvr2r(Aq+3*(M+1), st->cur_sub_Aq, (M+1));
            }

            if( bfi )
            {
                /* PLC: [ACELP: general]
                 * PLC: Use the ACELP like concealment */
                con_acelp( Aq, st->core_ext_mode, &synth[0], pitch, pit_gain, st->stab_fac, st, pitch_buf, voice_factors, ptr_bwe_exc );
                mvr2r( &st->mem_pitch_gain[2], &st->mem_pitch_gain[st->nb_subfr+2], st->nb_subfr );
                set_zero( &st->mem_pitch_gain[2],st->nb_subfr );
            }
            else
            {
                decoder_acelp( st, *coder_type, prm, Aq, st->acelp_cfg, &synth[0], pitch, pit_gain, st->stab_fac, pitch_buf, voice_factors, LSF_Q_prediction, ptr_bwe_exc );

                if( st->flagGuidedAcelp > 0 )
                {
                    st->guidedT0 = min(st->T0_4th + st->guidedT0, NBPSF_PIT_MAX);
                }

                for (i = 0; i < st->nb_subfr; i++)
                {
                    st->mem_pitch_gain[2+(2*st->nb_subfr-1)-i] = st->mem_pitch_gain[2+ (st->nb_subfr-1) -i];
                    st->mem_pitch_gain[2+(st->nb_subfr-1)-i] = pit_gain[i];
                }
            }
        }

        /* LPC for ACELP/BBWE */
        if( st->narrowBand || st->sr_core==12800 || st->sr_core==16000 )
        {
            mvr2r( Aq, st->mem_Aq, nb_subfr*(M+1) );
        }

        /* PLC: [TCX: Tonal Concealment] */
        /* Signal that this frame is not TCX */
        TonalMDCTConceal_UpdateState( &st->tonalMDCTconceal, 0, 0, 0, 0 );

        if( !bfi )
        {
            st->second_last_tns_active = st->last_tns_active;
            st->last_tns_active = 0;
            st->tcxltp_last_gain_unmodified = 0.0f;
        }

    }



    /*--------------------------------------------------------------------------------*
     * TCX20
     *--------------------------------------------------------------------------------*/


    if( st->core == TCX_20_CORE )
    {
        /* Set pointer to parameters */
        prm = param;

        /* Stability Factor */
        if( !bfi )
        {
            st->stab_fac = lsf_stab( &lsf[M], &lsf[0], 0, st->L_frame );
        }

        if (st->enableTcxLpc)
        {
            /* Convert quantized lsp to A */
            lsp2a_stab( &lsp[M], Aq, M );
        }
        else
        {
            if (!st->tcxonly)
            {
                if( !bfi && st->prev_bfi && !(st->safety_net) && st->rate_switching_reset )
                {
                    /* diffuse LPC power on rate switching*/
                    RecLpcSpecPowDiffuseLc( &lsp[M], &lsp[0], &lsf[M], st, 0 );
                    int_lsp( L_frame, &lsp[0], &lsp[M], Aq, M, interpol_frac_12k8, 0 );
                    mvr2r( &lsf[M], lsfnew_uw, M );
                }
                else
                {
                    /* LPC Interpolation for TCX */
                    E_LPC_int_lpc_tcx( &lsp[0], &lsp[M], Aq );
                }
            }
            else
            {
                lsp2a_stab( &lsp[M], Aq, M );
            }
        }

        if( !bfi && st->tcx_lpc_shaped_ari )
        {
            basop_E_LPC_f_lsp_a_conversion(lspind, Aind, M);
        }

        /* TCX decoder */
        decoder_tcx( &st->tcx_cfg, prm, Aq, Aind, L_frame, L_frameTCX, st->tcx_cfg.tcx_coded_lines,
                     &synth[0], &synthFB[0], st, *coder_type, bfi, 0, st->stab_fac );

    }

    /*--------------------------------------------------------------------------------*
     * TCX10
     *--------------------------------------------------------------------------------*/

    if( st->core == TCX_10_CORE )
    {
        prm = NULL;                   /* just to avoid MSVC warnings */

        for (k=0; k<2; k++)
        {
            /* Set pointer to parameters */
            prm = param + (k*DEC_NPRM_DIV);

            /* Stability Factor */
            if( !bfi )
            {
                st->stab_fac = lsf_stab( &lsf[(k+1)*M], &lsf[k*M], 0, st->L_frame );
            }

            lsp2a_stab( &lsp[(k+1)*M], Aq, M );
            IGFDecRestoreTCX10SubFrameData( &st->hIGFDec, k );

            /* TCX decoder */
            decoder_tcx( &st->tcx_cfg, prm, Aq, Aind, L_frame/2, L_frameTCX/2, st->tcx_cfg.tcx_coded_lines/2,
                         &synth[k*L_frame/2], &synthFB[k*L_frameTCX/2], st, *coder_type, bfi, k, st->stab_fac );
        }

    }

    if (st->core == TCX_10_CORE || st->core == TCX_20_CORE)
    {
        if( st->enablePlcWaveadjust ||           /* bfi      */
                ( st->last_total_brate >= HQ_48k &&  /* recovery */
                  st->last_codec_mode == MODE2 ) )
        {
            /* waveform adjustment */
            concealment_signal_tuning( bfi, st->core, synthFB, &st->plcInfo, st->nbLostCmpt, st->prev_bfi,
                                       st->tonalMDCTconceal.secondLastPcmOut, past_core_mode,st->tonalMDCTconceal.lastPcmOut, st );

            if( (bfi || st->prev_bfi) && st->plcInfo.Pitch && st->plcInfo.concealment_method == TCX_NONTONAL )
            {
                lerp( synthFB, synth, L_frame, L_frameTCX );

                if( !bfi && st->prev_bfi )
                {
                    st->plcInfo.Pitch = 0;
                }
            }
        }

        if (!bfi)
        {
            TonalMDCTConceal_SaveTimeSignal( &st->tonalMDCTconceal, synthFB, L_frameTCX );
        }

        decoder_tcx_post( st, synth, synthFB, Aq, bfi );

        if( st->core == TCX_20_CORE )
        {
            /* LPC Interpolation for BBWE/post-processing */
            if( st->narrowBand || st->sr_core==12800 || st->sr_core==16000 )
            {
                int_lsp( L_frame, st->lspold_uw, lspnew_uw, Aq, M, interpol_frac_12k8, 0 );
                mvr2r( Aq, st->mem_Aq, nb_subfr*(M+1) );
            }
        }
    }


    /* PLC: [Common: Classification] */
    /* the classifier buffer is always updated if the sr is at
       16000 or below - the classification itself is just performed if(!st->tcxonly ) */
    if( st->sr_core <= 16000 )
    {
        if( st->core == TCX_20_CORE || st->core == TCX_10_CORE || (st->tcxonly && st->bfi) )
        {
            float pitch_C[4];

            /* note: the classifier needs the pitch only for tcx_only == 0 , i.e. not for TCX10 */
            pitch_C[0] = pitch_C[1] = pitch_C[2] = pitch_C[3] = (float)floor(st->old_fpitch+0.5f);

            FEC_clas_estim( synth, pitch_C, st->L_frame, st->tcxonly? GENERIC : st->core_ext_mode, st->codec_mode, st->mem_syn_clas_estim,
                            &st->clas_dec, &st->lp_ener_bfi, st->core_brate, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                            st->tcxltp ? st->tcxltp_last_gain_unmodified : -1.0f, st->narrowBand, 1,
                            bfi, st->preemph_fac, st->tcxonly, st->last_core_brate );
        }
    }




    /*--------------------------------------------------------------------------------*
     * Updates
     *--------------------------------------------------------------------------------*/


    if( bfi && st->last_core != ACELP_CORE )
    {
        /* Update FEC_scale_syn parameters */
        if (st->tcxltp_gain == 0)
        {
            fer_energy( L_frame, UNVOICED, synth, L_frame/2, &st->enr_old, L_frame );
        }
        else
        {
            fer_energy( L_frame, st->clas_dec, synth, st->old_fpitch, &st->enr_old, L_frame );
        }
    }

    if( !bfi && st->clas_dec >= VOICED_TRANSITION && st->clas_dec < INACTIVE_CLAS)
    {
        short offset;

        if( st->core == ACELP_CORE )
        {
            offset = (st->nb_subfr-1)*(M+1);
        }
        else
        {
            /* TCX */
            offset = 0;
        }

        /* use latest LPC set */
        st->old_enr_LP = enr_1_Az( Aq+offset, L_SUBFR );
    }


    /* Update */
    mvr2r( synth_buf+L_frame, st->old_synth, st->old_synth_len );
    mvr2r( st->old_synthFB + L_frameTCX - NS2SA(st->output_Fs, PH_ECU_MEM_NS), st->synth_history, NS2SA(st->output_Fs, PH_ECU_MEM_NS) );
    mvr2r( synth_bufFB+L_frameTCX, st->old_synthFB, st->old_synth_lenFB );
    mvr2r( st->old_out+NS2SA(st->output_Fs,N_ZERO_MDCT_NS), st->old_synthFB+st->old_synth_lenFB, NS2SA(st->output_Fs, PH_ECU_LOOKAHEAD_NS));

    mvr2r( &lspnew_uw[(st->numlpc-1)*M], st->lspold_uw, M );
    mvr2r( &lsfnew_uw[(st->numlpc-1)*M], st->lsfold_uw, M );

    if( bfi == 1 )
    {
        mvr2r( st->lspold_uw, st->lsp_old, M ); /* for recovery */
        mvr2r( st->lsfold_uw, st->lsf_old, M ); /* for recovery */
    }
    else
    {
        mvr2r( &lsp[st->numlpc*M], st->lsp_old, M );
        mvr2r( &lsf[st->numlpc*M], st->lsf_old, M );
    }
    mvr2r( st->lsp_q_cng, st->old_lsp_q_cng, M );
    mvr2r( st->lsf_q_cng, st->old_lsf_q_cng, M );

    /* Update MODE1 CNG parameters */
    if( !st->tcxonly )
    {
        /* update CNG parameters in active frames */
        if ( st->bwidth == NB && st->enableTcxLpc && st->core != ACELP_CORE )
        {
            float buf[L_LP], res[L_FRAME], A[M+1], r[M+1], tmp, lsptmp[M];
            assert(st->L_frame==L_FRAME);
            mvr2r(synth+L_FRAME-L_LP, buf, L_LP);
            tmp = synth[L_FRAME-L_LP-1];
            preemph(buf, st->preemph_fac, L_LP, &tmp);
            autocorr( buf, r, M, L_LP, LP_assym_window, 0, 0, 0 );
            lag_wind(r, M, INT_FS_12k8, LAGW_WEAK);
            lev_dur( A, r, M, NULL );
            a2lsp_stab( A, lsptmp, &lspnew_uw[0] );
            residu(A, M, buf+L_LP-L_FRAME, res, L_FRAME );

            cng_params_upd( lsptmp, res, st->L_frame, &st->ho_circ_ptr, st->ho_ener_circ, &st->ho_circ_size, st->ho_lsp_circ,
                            DEC, st->ho_env_circ, NULL, NULL, NULL, st->last_active_brate );
        }
        else
        {
            cng_params_upd( &lsp[M], st->old_exc+L_EXC_MEM_DEC-st->L_frame, st->L_frame, &st->ho_circ_ptr, st->ho_ener_circ,
                            &st->ho_circ_size, st->ho_lsp_circ, DEC, st->ho_env_circ, NULL, NULL, NULL, st->last_active_brate );
        }

        /* Set 16k LSP flag for CNG buffer */
        st->ho_16k_lsp[st->ho_circ_ptr] = (st->L_frame == L_FRAME ? 0 : 1 );
    }

    st->last_is_cng = 0;

    /* Postfiltering */
    post_decoder( st, *coder_type, synth_buf, pit_gain, pitch, signal_out, bpf_noise_buf );

    if( signal_outFB )
    {
        mvr2r( synthFB, signal_outFB, L_frameTCX );
    }

    if( st->enablePlcWaveadjust)
    {
        if(!bfi)
        {
            st->plcInfo.nbLostCmpt = 0;
        }

        if( st->core == ACELP_CORE ) /* may happen only if bfi==1 */
        {
            set_state( st->plcInfo.Transient, st->core, MAX_POST_LEN );
        }
    }


    return;
}
