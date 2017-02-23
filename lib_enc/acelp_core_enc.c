/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_enc.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * acelp_core_enc()
 *
 * ACELP core encoder
 *--------------------------------------------------------------------*/

void acelp_core_enc(
    Encoder_State *st,                      /* i/o: encoder state structure             */
    LPD_state *mem,                     /* i/o: acelp memories                      */
    const float inp[],                    /* i  : input signal of the current frame   */
    const short vad_flag,
    const float ener,                     /* i  : residual energy from Levinson-Durbin*/
    const short pitch[3],                 /* i  : open-loop pitch values for quantiz. */
    const float voicing[3],               /* i  : Open-loop pitch gains               */
    float A[NB_SUBFR16k*(M+1)],     /* i  : A(z) unquantized for the 4 subframes*/
    float Aw[NB_SUBFR16k*(M+1)],    /* i  : weighted A(z) unquant. for subframes*/
    const float epsP[M+1],                /* i  : LP prediction errors                */
    float lsp_new[M],               /* i  : LSPs at the end of the frame        */
    float lsp_mid[M],               /* i  : LSPs in the middle of the frame     */
    short coder_type,               /* i  : coding type                         */
    const short sharpFlag,                /* i  : formant sharpening flag             */
    short vad_hover_flag,
    const short attack_flag,              /* i  : flag signalling attack encoded by AC mode (GSC) */
    float bwe_exc_extended[],       /* i/o: bandwidth extended excitation       */
    float *voice_factors,           /* o  : voicing factors                     */
    float old_syn_12k8_16k[],       /* o  : intermediate ACELP synthesis at 12.8kHz or 16kHz to be used by SWB BWE */
    float pitch_buf[NB_SUBFR16k],   /* o  : floating pitch for each subframe    */
    short *unbits                   /* o  : number of unused bits               */
)
{
    short nBits;                                              /* reserved bits                        */
    short i;
    float old_exc[L_EXC], *exc;                               /* excitation signal buffer             */
    float lsf_new[M];                                         /* ISFs at the end of the frame         */
    float Aq[NB_SUBFR16k*(M+1)];                              /* A(z)   quantized for the 4 subframes */
    float syn[L_FRAME16k];                                    /* synthesis signal buffer              */
    float res[L_FRAME16k];                                    /* Residual signal for FER protection   */
    float exc2[L_FRAME16k];                                   /* enhanced excitation                  */
    float Es_pred;                                            /* predicited scaled innovation energy  */
    float tmp_noise;                                          /* NB post-filter long-term noise energy*/
    short tc_subfr;                                           /* TC sub-frame indication              */
    float old_bwe_exc[(PIT16k_MAX + (L_FRAME16k + 1) + L_SUBFR16k) * 2]; /* excitation buffer */
    float *bwe_exc;                                           /* excitation for SWB TBE */
    short allow_cn_step;
    float int_fs;
    float att;
    float lim;
    short T_op[3];
    short nb_bits;      /* parameters handling */
    int indice;

    /* SC-VBR - back-up memories for LSF quantizer and synthesis filter */
    short mCb1, pstreaklen;
    float mem_MA[M], mem_AR[M], Bin_E[L_FFT], Bin_E_old[L_FFT/2], lsp_new_bck[M], lsp_mid_bck[M], mem_syn_bck[M];
    float clip_var, mem_w0_bck, streaklimit;

    float q_env[NUM_ENV_CNG];
    short sid_bw = -1;
    float exc3[L_FRAME16k];
    float syn1[L_FRAME16k];
    float enr;
    float enr_index;

    float tilt_code_bck;
    float gc_threshold_bck;
    float clip_var_bck[6];
    short next_force_sf_bck;



    /*------------------------------------------------------------------*
     * Initialization
     *------------------------------------------------------------------*/

    Es_pred = 0;

    mvs2s( pitch, T_op, 3 );

    /* convert pitch values to 16kHz domain */
    if ( st->L_frame == L_FRAME16k )
    {
        T_op[0] = (short)(T_op[0] * 1.25f + 0.5f);
        T_op[1] = (short)(T_op[1] * 1.25f + 0.5f);
        T_op[2] = T_op[1];
    }

    exc = old_exc + L_EXC_MEM;                                  /* pointer to excitation signal in the current frame */
    mvr2r( mem->old_exc, old_exc, L_EXC_MEM );

    bwe_exc = old_bwe_exc + PIT16k_MAX * 2;                     /* pointer to BWE excitation signal in the current frame */
    mvr2r( st->old_bwe_exc, old_bwe_exc, PIT16k_MAX * 2 );

    st->bpf_off = 0;
    if( st->last_core == HQ_CORE || st->last_codec_mode == MODE2 )
    {
        /* in case of HQ->ACELP switching, do not apply BPF */
        st->bpf_off = 1;
        /* reset the GSC pre echo energy threshold in case of switching */
        st->Last_frame_ener = (float)MAX_32;
    }

    /* force safety-net LSFQ in the first frames after CNG segment */
    if( st->last_core_brate <= SID_2k40 )
    {
        st->Nb_ACELP_frames = 0;
    }
    st->Nb_ACELP_frames++;

    if( st->L_frame == L_FRAME )
    {
        int_fs = INT_FS_12k8;
    }
    else
    {
        int_fs = INT_FS_16k;
    }

    tmp_noise = 0;
    tc_subfr = 0;

    /* SC-VBR temporary variables */
    mCb1 = 0;
    pstreaklen = 0;
    clip_var = 0;
    mem_w0_bck = 0;
    streaklimit = 0;

    /* channel-aware mode */
    reset_rf_indices(st);

    /*-----------------------------------------------------------------*
     * ACELP@12k8 / ACELP@16k switching
     *-----------------------------------------------------------------*/

    if( st->last_L_frame != st->L_frame && st->last_core != HQ_CORE )
    {
        /* in case of switching, do not apply BPF (flag employed also in updt_enc()) */
        st->bpf_off = 1;

        /* force safety-net LSFQ in the first frames after ACELP@12k8/ACELP@16k switching */
        st->Nb_ACELP_frames = 1;

        /* convert old quantized LSP vector */
        if( st->L_frame == L_FRAME )
        {
            st->rate_switching_reset = lsp_convert_poly( st->lsp_old, st->L_frame, 0 );
        }
        else
        {
            st->rate_switching_reset = st->rate_switching_reset_16kHz;
            mvr2r( st->lsp_old16k, st->lsp_old, M );
        }

        /* convert old quantized LSF vector */
        lsp2lsf( st->lsp_old, st->lsf_old, M, int_fs );

        /* interpolation of unquantized ISPs */
        if( st->rate_switching_reset )
        {
            /*extrapolation in case of unstable LSP*/
            int_lsp4( st->L_frame, lsp_mid, lsp_mid, lsp_new, A, M, 0 );
        }
        else
        {
            int_lsp4( st->L_frame, st->lsp_old, lsp_mid, lsp_new, A, M, 0 );
        }

        /* Reset LPC mem */
        mvr2r( GEWB_Ave, st->mem_AR, M );
        set_zero( st->mem_MA, M );

        /* update synthesis filter memories */
        synth_mem_updt2( st->L_frame, st->last_L_frame, mem->old_exc, mem->mem_syn_r, mem->mem_syn2, mem->mem_syn, ENC );
        mvr2r( mem->old_exc, old_exc, L_EXC_MEM );
        mvr2r( mem->mem_syn2, st->mem_syn1, M );
        mvr2r( mem->mem_syn2, mem->mem_syn3, M );

        /* update Aw[] coefficients */
        weight_a_subfr( st->L_frame/L_SUBFR, A, Aw, st->gamma, M );
    }

    if( st->last_bwidth == NB && st->bwidth != NB && st->ini_frame != 0 )
    {
        st->rate_switching_reset = 1;
    }

    /*----------------------------------------------------------------*
     * Encoding of CNG frames
     *----------------------------------------------------------------*/

    if( st->core_brate == SID_2k40 || st->core_brate == FRAME_NO_DATA )
    {
        if( st->cng_type == LP_CNG )
        {
            /* Run CNG post parameter update */
            cng_params_postupd( st->ho_circ_ptr, &st->cng_buf_cnt, st->cng_exc2_buf,
                                st->cng_brate_buf, st->ho_env_circ);
            /* encode CNG parameters */
            CNG_enc( st, st->L_frame, Aq, inp, ener, lsp_new, lsf_new , &allow_cn_step, st->burst_ho_cnt, q_env, &sid_bw, st->exc_mem2 );
            /* comfort noise generation */
            CNG_exc( st->core_brate, st->L_frame, &st->Enew, &st->cng_seed, exc, exc2, &st->lp_ener,
                     st->last_core_brate, &st->first_CNG, &st->cng_ener_seed, bwe_exc, allow_cn_step, &st->last_allow_cn_step, st->num_ho,
                     q_env, st->lp_env, st->old_env, st->exc_mem, st->exc_mem1, &sid_bw, &st->cng_ener_seed1, exc3, st->Opt_AMR_WB );
        }
        else
        {
            if( st->core_brate == SID_2k40 )
            {
                FdCng_encodeSID( st->hFdCngEnc, st, st->preemph_fac );
                st->last_CNG_L_frame = st->L_frame;
            }

            generate_comfort_noise_enc( st );

            FdCng_exc( st->hFdCngEnc->hFdCngCom, &st->CNG_mode, st->L_frame, st->lsp_old, st->first_CNG, st->lspCNG, Aq, lsp_new,lsf_new, exc, exc2, bwe_exc);
            mvr2r( exc2, exc3, st->L_frame );
            if( st->core_brate == SID_2k40 )
            {
                enr = dotp( exc, exc, st->L_frame ) / st->L_frame;
                enr = (float)log10( enr + 0.1f ) / (float)log10( 2.0f );

                /* decrease the energy in case of WB input */
                if( st->bwidth != NB )
                {
                    if( st->bwidth == WB )
                    {
                        if( st->CNG_mode >= 0 )
                        {
                            /* Bitrate adapted attenuation */
                            att = ENR_ATT[st->CNG_mode];
                        }
                        else
                        {
                            /* Use least attenuation for higher bitrates */
                            att = ENR_ATT[4];
                        }
                    }
                    else
                    {
                        att = 1.5f;
                    }

                    enr -= att;
                }

                enr_index = (short)( (enr + 2.0f) * STEP_SID );
                if( enr_index > 127 )
                {
                    enr_index = 127;
                }

                if( enr_index < 0 )
                {
                    enr_index = 0;
                }
                st->old_enr_index = enr_index;
            }
        }


        /* synthesis at 12.8kHz sampling rate */
        syn_12k8( st->L_frame, Aq, exc3, syn1, mem->mem_syn3, 1 );

        /* reset the encoder */
        CNG_reset_enc( st, mem, pitch_buf, voice_factors, 0 );

        /* update st->mem_syn1 for ACELP core switching */
        mvr2r( mem->mem_syn3, st->mem_syn1, M );

        /* update ACELP core synthesis filter memory */
        mvr2r( mem->mem_syn3, mem->mem_syn, M );

        /* update old synthesis buffer - needed for ACELP internal sampling rate switching */
        mvr2r( syn1 + st->L_frame - L_SYN_MEM, mem->mem_syn_r, L_SYN_MEM );


        /* save and delay synthesis to be used by SWB BWE */
        save_old_syn( st->L_frame, syn1, old_syn_12k8_16k, st->old_syn_12k8_16k, st->preemph_fac, &st->mem_deemph_old_syn );

        /*Update MODE2 core switching memory*/
        deemph( syn1, st->preemph_fac, st->L_frame,  &(st->LPDmem.syn[M]) );
        mvr2r( syn1+st->L_frame-M-1, st->LPDmem.syn, M+1 );
    }

    /*----------------------------------------------------------------*
     * Encoding of all other frames
     *----------------------------------------------------------------*/

    else
    {
        /*-----------------------------------------------------------------*
         * After inactive period, use the most up-to-date ISPs
         *-----------------------------------------------------------------*/

        if ( st->last_core_brate == FRAME_NO_DATA || st->last_core_brate == SID_2k40 )
        {
            mvr2r( st->lspCNG, st->lsp_old, M );

            lsp2lsf( st->lspCNG, st->lsf_old, M, int_fs );
        }

        /*-----------------------------------------------------------------*
         * Reset higher ACELP pre-quantizer in case of switching
         *-----------------------------------------------------------------*/

        if( !st->use_acelp_preq )
        {
            st->mem_deemp_preQ = 0.0f;
            st->mem_preemp_preQ = 0.0f;
            st->last_nq_preQ = 0;
        }
        st->use_acelp_preq = 0;

        /*-----------------------------------------------------------------*
         * LSF Quantization
         * A[z] calculation
         *-----------------------------------------------------------------*/

        /* SC-VBR & channel-aware mode - back-up memories for LSF quantizer and synthesis filter */
        lsf_syn_mem_backup( st, mem, &tilt_code_bck, &gc_threshold_bck, clip_var_bck, &next_force_sf_bck, lsp_new, lsp_mid, &clip_var,
                            mem_AR, mem_MA, lsp_new_bck, lsp_mid_bck, &mCb1, Bin_E, Bin_E_old, mem_syn_bck, &mem_w0_bck, &streaklimit, &pstreaklen );

        lsf_enc( st, st->L_frame, coder_type, lsf_new, lsp_new, lsp_mid, Aq, &st->stab_fac, st->Nb_ACELP_frames );

        /*---------------------------------------------------------------*
         * Calculation of LP residual (filtering through A[z] filter)
         *---------------------------------------------------------------*/

        calc_residu( inp, res, Aq, st->L_frame );

        /* smoothing in case of CNG */
        if( st->Opt_DTX_ON && vad_hover_flag )
        {
            st->burst_ho_cnt++;
            if(st->burst_ho_cnt > HO_HIST_SIZE)
            {
                st->burst_ho_cnt = HO_HIST_SIZE;
            }
            if( st->bwidth != NB )
            {
                if( st->bwidth == WB && st->CNG_mode >= 0 )
                {
                    lim = HO_ATT[st->CNG_mode];
                }
                else
                {
                    lim = 0.6f;
                }

                att = lim/6.0f;
                att = 1.0f/(1 + att * st->burst_ho_cnt);

                if ( att < lim )
                {
                    att = lim;
                }

                for( i = 0; i < st->L_frame; i++ )
                {
                    res[i] *= att;
                }
            }
        }
        else
        {
            st->burst_ho_cnt = 0;
        }

        /*---------------------------------------------------------------*
         * Calculation of prediction for scaled innovation energy
         * (for memory-less gain quantizer)
         *---------------------------------------------------------------*/

        if( ( coder_type != UNVOICED && coder_type != AUDIO && coder_type != INACTIVE && !(st->core_brate <= ACELP_8k00 && coder_type != TRANSITION) )
                || (coder_type == INACTIVE && st->total_brate >= ACELP_32k) )
        {
            nb_bits = Es_pred_bits_tbl[BIT_ALLOC_IDX(st->core_brate, coder_type, -1, -1)];
            Es_pred_enc( &Es_pred, &indice, st->L_frame, L_SUBFR, res, voicing, nb_bits, 0 );

            push_indice( st, IND_ES_PRED, indice, nb_bits );
        }

        /*------------------------------------------------------------*
         * Encode excitation according to coding type
         *------------------------------------------------------------*/

        if( st->nelp_mode )
        {
            /* SC-VBR - NELP frames */
            encod_nelp( st, mem, inp, Aw, Aq, res, syn, &tmp_noise, exc, exc2, pitch_buf, voice_factors, bwe_exc );
        }
        else if( coder_type == UNVOICED )
        {
            /* UNVOICED frames (Gauss. excitation) */
            encod_unvoiced( st, mem, inp, Aw, Aq, vad_flag, res, syn, &tmp_noise, exc, pitch_buf, voice_factors, bwe_exc );
        }
        else if( coder_type == TRANSITION )
        {
            tc_subfr = encod_tran( st, mem, st->L_frame, inp, Aw, Aq, coder_type, Es_pred, T_op, voicing, res, syn,
                                   exc, exc2, pitch_buf, voice_factors, bwe_exc, attack_flag, unbits, sharpFlag );
        }
        else if( st->ppp_mode )
        {
            /* SC-VBR - PPP frames */
            encod_ppp( st, mem, inp, Aw, Aq, &coder_type, sharpFlag, T_op, voicing, res, syn, exc, exc2, pitch_buf, voice_factors, bwe_exc );

            if( st->bump_up )   /* PPP failed, bump up */
            {
                /* restore memories of LSF quantizer and synthesis filter */
                lsf_syn_mem_restore( st,  mem, tilt_code_bck,
                                     gc_threshold_bck,
                                     clip_var_bck, next_force_sf_bck,
                                     lsp_new, lsp_mid, clip_var, mem_AR, mem_MA, lsp_new_bck, lsp_mid_bck, mCb1, Bin_E,Bin_E_old,mem_syn_bck, mem_w0_bck, streaklimit, pstreaklen );

                /* redo LSF quantization */
                lsf_enc( st, st->L_frame, coder_type, lsf_new, lsp_new, lsp_mid, Aq, &st->stab_fac, st->Nb_ACELP_frames );

                /* recalculation of LP residual (filtering through A[z] filter) */
                calc_residu( inp, res, Aq, st->L_frame );
                st->burst_ho_cnt = 0;

                /* VOICED frames in SC-VBR */
                encod_gen_voic( st, mem, st->L_frame, sharpFlag, inp, Aw, Aq, coder_type, Es_pred, T_op, voicing, res, syn,
                                exc, exc2, pitch_buf, voice_factors, bwe_exc, unbits );
            }
        }
        else if( coder_type == AUDIO || ( coder_type == INACTIVE && st->core_brate <= ACELP_24k40 ) )
        {
            /* AUDIO and INACTIVE frames (coded by GSC technology) */
            encod_audio( st, mem, inp, Aw, Aq, T_op, voicing, res, syn, exc, pitch_buf, voice_factors, bwe_exc, attack_flag, coder_type, lsf_new, &tmp_noise );
        }
        else
        {
            /* GENERIC, VOICED and INACTIVE frames (coded by AVQ technology) */
            encod_gen_voic( st, mem, st->L_frame, sharpFlag, inp, Aw, Aq, coder_type, Es_pred, T_op, voicing, res, syn,
                            exc, exc2, pitch_buf, voice_factors, bwe_exc, unbits );
        }


        /* update st->mem_syn1 for ACELP core switching */
        mvr2r( mem->mem_syn, st->mem_syn1, M );

        /* update old synthesis buffer - needed for ACELP internal sampling rate switching */
        mvr2r( syn + st->L_frame - L_SYN_MEM, mem->mem_syn_r, L_SYN_MEM );

        /* save and delay synthesis to be used by SWB BWE */
        save_old_syn( st->L_frame, syn, old_syn_12k8_16k, st->old_syn_12k8_16k, st->preemph_fac, &st->mem_deemph_old_syn );

        /*Update MODE2 core switching memory*/
        mvr2r( syn, syn1, st->L_frame );
        deemph( syn1, st->preemph_fac, st->L_frame,  &(st->LPDmem.syn[M]) );
        mvr2r( syn1+st->L_frame-M-1, st->LPDmem.syn, M+1 );


        /*--------------------------------------------------------------------------------------*
         * Modify the excitation signal when the noise is stationary
         *--------------------------------------------------------------------------------------*/

        if ( st->nelp_mode != 1 )
        {
            /* exc2 buffer is needed only for updating of Aq[] which is needed for core switching */
            mvr2r( exc, exc2, st->L_frame);
            stat_noise_uv_enc( st, coder_type, epsP, lsp_new, lsp_mid, Aq, exc2 );
        }

        /*-----------------------------------------------------------------*
         * Encode supplementary information for Frame Error Concealment
         *-----------------------------------------------------------------*/

        FEC_encode( st, syn, coder_type, st->clas, pitch_buf, res, &st->Last_pulse_pos, st->L_frame, st->total_brate, st->core_brate );

        if( st->L_frame == L_FRAME )
        {
            mvr2r( Aq+2*(M+1), st->cur_sub_Aq, (M+1) );
        }
        else
        {
            mvr2r( Aq+3*(M+1), st->cur_sub_Aq, (M+1) );
        }


    }   /* end of active inp coding */


    /*-----------------------------------------------------------------*
     * Write ACELP unused bits
     *-----------------------------------------------------------------*/

    if ( st->core_brate != SID_2k40 && st->core_brate != FRAME_NO_DATA && st->core_brate != PPP_NELP_2k80 )
    {
        /* unused bits */
        if ( coder_type == AUDIO || ( coder_type == INACTIVE && st->core_brate <= ACELP_24k40 ) )
        {
            nBits = 0;
        }
        else if(st->L_frame == L_FRAME )
        {
            nBits = reserved_bits_tbl[BIT_ALLOC_IDX(st->core_brate, coder_type, -1, TC_SUBFR2IDX(tc_subfr))];
        }
        else
        {
            nBits = 0;
        }

        while( nBits > 0 )
        {
            i = min(nBits, 16);
            push_indice( st, IND_UNUSED, 0, i );
            nBits -= i;
        }
    }

    /*-----------------------------------------------------------------*
     * Apply non linearity in case of SWB TBE
     *-----------------------------------------------------------------*/

    if( (st->last_Opt_SC_VBR == 1 && st->Opt_SC_VBR == 0) || ((st->extl == SWB_TBE || st->extl == WB_TBE || st->extl == FB_TBE) && st->last_extl != SWB_TBE && st->last_extl != WB_TBE && st->last_extl != FB_TBE) )
    {
        st->bwe_non_lin_prev_scale = 0.0f;
        set_f( st->old_bwe_exc_extended, 0.0f, NL_BUFF_OFFSET );
    }

    if( !st->Opt_SC_VBR )
    {
        /* Apply a non linearity to the SHB excitation */
        non_linearity( bwe_exc, bwe_exc_extended, st->old_bwe_exc_extended, L_FRAME32k, &st->bwe_non_lin_prev_scale, coder_type,  voice_factors, st->L_frame );
    }

    if ( st->core_brate == SID_2k40 || st->core_brate == FRAME_NO_DATA )
    {
        st->bwe_non_lin_prev_scale = 0.0f;
    }

    /*-----------------------------------------------------------------*
     * Updates
     *-----------------------------------------------------------------*/

    updt_enc( st, st->L_frame, coder_type, old_exc, pitch_buf, Es_pred, Aq, lsf_new, lsp_new, old_bwe_exc );

    if( st->Opt_DTX_ON && st->core_brate > SID_2k40 )
    {
        /* update CNG parameters in active frames */
        cng_params_upd( lsp_new, exc, st->L_frame, &st->ho_circ_ptr, st->ho_ener_circ,
                        &st->ho_circ_size, st->ho_lsp_circ, ENC, st->ho_env_circ,
                        &st->cng_buf_cnt, st->cng_exc2_buf, st->cng_brate_buf, st->last_active_brate );

        if( st->L_frame == L_FRAME )
        {
            /* store LSPs@16k, potentially to be used in CNG@16k */
            mvr2r( st->lsp_old16k, &(st->ho_lsp_circ2[(st->ho_circ_ptr)*M]), M );
        }

        /* set LSP@16k flag for the first buffer */
        st->ho_16k_lsp[st->ho_circ_ptr] = (st->L_frame == L_FRAME ? 0 : 1 );

        /* efficient DTX hangover control */
        if ( st->burst_ho_cnt > 1 )
        {
            dtx_hangover_control( st, lsp_new );
        }
    }

    /* SC-VBR update of average data rate */
    if ( vad_flag == 1 )
    {
        update_average_rate( st );

    }



    return;
}
