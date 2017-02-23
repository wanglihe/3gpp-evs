/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "prot.h"
#include "options.h"
#include "rom_com.h"


/*-------------------------------------------------------------------*
 * decoder_acelp()
 *
 * Decode ACELP frame
 *-------------------------------------------------------------------*/

void decoder_acelp(
    Decoder_State *st,            /* i/o: coder memory state          */
    const short coder_type,     /* i  : coder type                  */
    int   prm[],          /* i  : parameters                  */
    const float A[],            /* i  : coefficients NxAz[M+1]      */
    ACELP_config acelp_cfg,      /* i  : ACELP config                */
    float synth[],        /* i/o: synthesis                   */
    int   *pT,            /* o  : pitch for all subframe      */
    float *pgainT,        /* o  : pitch gain for all subfr    */
    const float stab_fac,       /* i  : stability of isf            */
    float *pitch_buffer,  /* i/o: pitch values for each subfr.*/
    float *voice_factors, /* o  : voicing factors             */
    const short LSF_Q_prediction,/* i  : LSF prediction mode        */
    float *bwe_exc        /* o  : excitation for SWB TBE      */
)
{
    short i, i_subfr, L_frame;
    int   T0, T0_frac, T0_min, T0_min_frac, T0_max, T0_max_frac, T0_res;
    float tmp, gain_pit, gain_code, Es_pred;
    float code[L_SUBFR];
    float mem_syn[M], *syn, syn_buf[M+L_FRAME16k+L_FRAME16k/2];
    float *exc, exc_buf[L_EXC_MEM_DEC+L_FRAME16k+1];
    float exc2[L_FRAME16k];
    const float *p_A;
    float *pt_pitch, pitch_buf[NB_SUBFR16k];
    float gain_inov;
    float mem_back[M];
    float h1[L_FRAME16k/4+1];
    float mem[M];
    const float *pA;
    float gain_code2;
    float code2[L_SUBFR];
    short lp_flag;
    int offset;
    float error = 0.0f;
    float gain_preQ = 0;                /* Gain of prequantizer excitation   */
    float code_preQ[L_SUBFR];           /* Prequantizer excitation           */
    PulseConfig config;
    float weights[5];

    float prev_gain_pit;
    float tmp_noise;   /* Long term temporary noise energy */

    float gain_code_tmp;
    float gain_pit_tmp;
    float gain_code_pre;

    gain_inov = 0;    /* to avoid compilation warnings */
    T0 = 0;           /* to avoid compilation warnings */
    T0_frac = 0;      /* to avoid compilation warnings */
    T0_res = 0;       /* to avoid compilation warnings */

    /*------------------------------------------------------------------------*
     * Initializations                                                        *
     *------------------------------------------------------------------------*/

    gain_code_pre = 0;
    set_f( code_preQ, 0.f, L_SUBFR );

    gain_pit = 0;
    gain_code = 0;
    gain_code2 = 0.f;

    prev_gain_pit = 0;
    tmp_noise = 0;

    if( st->nb_subfr == NB_SUBFR )
    {
        weights[0] = 0.1f;
        weights[1] = 0.2f;
        weights[2] = 0.3f;
        weights[3] = 0.4f;
    }
    else  /*nb_subfr == NB_SUBFR16k */
    {
        weights[0] = (float)1/15;
        weights[1] = (float)2/15;
        weights[2] = (float)3/15;
        weights[3] = (float)4/15;
        weights[4] = (float)5/15;
    }

    st->lp_gainp = 0;
    st->lp_gainc = 0;


    /* Framing parameters */
    L_frame = st->L_frame;

    /*------------------------------------------------------------------------*
     * Previous frame is TCX                                                  *
     *------------------------------------------------------------------------*/
    /* Reset phase dispersion */

    if( st->last_core_bfi > ACELP_CORE )
    {
        set_zero( st->dispMem, 8 );
    }

    /* Update of synthesis filter memories in case of 12k8 core */
    if( st->prev_bfi && st->last_con_tcx && st->L_frame < L_FRAME16k )
    {
        synth_mem_updt2( st->L_frame, L_FRAME16k, st->old_exc, st->mem_syn_r, st->mem_syn2, NULL, DEC );
    }

    if( st->last_con_tcx && st->old_enr_LP )
    {
        float enr_LP, ratio;

        /* rescale excitation buffer if LPC energies differs too much */
        enr_LP = enr_1_Az( A, L_SUBFR );

        ratio = st->old_enr_LP/enr_LP;
        if (ratio < 0.8)
        {
            v_multc( st->old_exc, ratio, st->old_exc, L_EXC_MEM_DEC);
        }
    }

    /*------------------------------------------------------------------------*
     * Initialize buffers                                                     *
     *------------------------------------------------------------------------*/

    mvr2r( st->mem_syn2, mem_back, M );

    /* set ACELP synthesis memory */
    mvr2r( st->mem_syn2, mem_syn, M );

    /* set excitation memory*/
    exc=exc_buf+L_EXC_MEM_DEC;
    mvr2r( st->old_exc, exc_buf, L_EXC_MEM_DEC );
    *(exc+L_frame) = 0.f;

    /* Init syn buffer */
    syn = syn_buf + M;
    mvr2r( st->mem_syn2, syn_buf, M );

    /*------------------------------------------------------------------------*
     * Fast recovery flag
     *------------------------------------------------------------------------*/

    if( st->prev_bfi && coder_type == VOICED )
    {
        /*Force BPF to be applied fully*/
        st->bpf_gain_param = 3;
    }

    /*------------------------------------------------------------------------*
     * - decode mean_ener_code for gain decoder (d_gain2.c)                   *
     *------------------------------------------------------------------------*/

    if ( acelp_cfg.nrg_mode > 0 )
    {
        Es_pred_dec( &Es_pred, prm[0], acelp_cfg.nrg_bits,acelp_cfg.nrg_mode>1 );
        prm++;
    }
    else
    {
        Es_pred = 0.f;
    }

    /*------------------------------------------------------------------------*
     *          Loop for every subframe in the analysis frame                 *
     *------------------------------------------------------------------------*
     *  To find the pitch and innovation parameters. The subframe size is     *
     *  L_SUBFR and the loop is repeated L_frame/L_SUBFR times.               *
     *     - compute impulse response of weighted synthesis filter (h1[])     *
     *     - compute the target signal for pitch search                       *
     *     - find the closed-loop pitch parameters                            *
     *     - encode the pitch delay                                           *
     *     - update the impulse response h1[] by including fixed-gain pitch   *
     *     - find target vector for codebook search                           *
     *     - correlation between target vector and impulse response           *
     *     - codebook search                                                  *
     *     - encode codebook address                                          *
     *     - VQ of pitch and codebook gains                                   *
     *     - find synthesis speech                                            *
     *     - update states of weighting filter                                *
     *------------------------------------------------------------------------*/

    p_A = A;
    pt_pitch = pitch_buf;

    for( i_subfr = 0; i_subfr < L_frame; i_subfr += L_SUBFR )
    {
        if( st->use_partial_copy && st->rf_frame_type == RF_NELP )
        {
            if( i_subfr == 0 )
            {
                decod_nelp( st, UNVOICED, &tmp_noise, pt_pitch, exc, exc2, voice_factors, bwe_exc, 0 /*st->bfi*/, pgainT );
                set_f(pitch_buffer, L_SUBFR, NB_SUBFR);
            }
        }
        else
        {
            /*-------------------------------------------------------*
             * - Decode adaptive codebook.                           *
             *-------------------------------------------------------*/

            if( st->use_partial_copy && st->acelp_cfg.gains_mode[i_subfr/L_SUBFR] == 0 )
            {
                gain_pit = prev_gain_pit;
            }

            if( acelp_cfg.ltp_bits != 0 )
            {
                /*if( st->use_partial_copy
                    && st->rf_frame_type == RF_GENPRED
                    && ( i_subfr == L_SUBFR || i_subfr == 3*L_SUBFR ) )
                {
                    *pt_pitch = (float)T0 + (float)T0_frac/(float)T0_res;
                }
                else*/
                {
                    /* pitch lag decoding */
                    *pt_pitch = Mode2_pit_decode( acelp_cfg.ltp_mode, i_subfr, L_SUBFR, &prm, &T0, &T0_frac, &T0_res,
                                                  &T0_min, &T0_min_frac, &T0_max, &T0_max_frac, st->pit_min,
                                                  st->pit_fr1, st->pit_fr1b, st->pit_fr2, st->pit_max, st->pit_res_max );
                }

                /* find pitch excitation */
                if( st->pit_res_max == 6 && !(st->use_partial_copy) )
                {
                    if ( T0_res == (st->pit_res_max>>1) )
                    {
                        pred_lt4( &exc[i_subfr], &exc[i_subfr], T0, T0_frac<<1, L_SUBFR+1, inter6_2, PIT_L_INTERPOL6_2, PIT_UP_SAMP6 );
                    }
                    else
                    {
                        pred_lt4( &exc[i_subfr], &exc[i_subfr], T0, T0_frac, L_SUBFR+1, inter6_2, PIT_L_INTERPOL6_2, PIT_UP_SAMP6 );
                    }
                }
                else
                {
                    if( T0_res == (st->pit_res_max>>1) )
                    {
                        pred_lt4( &exc[i_subfr], &exc[i_subfr], T0, T0_frac<<1, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP );
                    }
                    else
                    {
                        pred_lt4( &exc[i_subfr], &exc[i_subfr], T0, T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP );
                    }
                }

                /* LP filtering of the adaptive excitation*/
                lp_flag = acelp_cfg.ltf_mode;

                if( acelp_cfg.ltf_mode == NORMAL_OPERATION )
                {
                    lp_flag = *prm;
                    prm++;
                }

                lp_filt_exc_dec( st, MODE2, st->core_brate, 0, coder_type, i_subfr, L_SUBFR, L_frame, lp_flag, exc );
            }
            else
            {
                /* No adaptive codebook (UC) */
                set_zero( exc + i_subfr, L_SUBFR );

                T0 = L_SUBFR;
                T0_frac = 0;
                T0_res = 1;
                pitch_buf[i_subfr/L_SUBFR] = (float)L_SUBFR;
            }

            if( st->igf )
            {
                if( st->sr_core == 12800 )
                {
                    offset = T0 * HIBND_ACB_L_FAC + (int) ((float) T0_frac * 0.25f * HIBND_ACB_L_FAC + 2 * HIBND_ACB_L_FAC + 0.5f) - 2 * HIBND_ACB_L_FAC;
                    for (i = 0; i < L_SUBFR * HIBND_ACB_L_FAC; i++)
                    {
                        bwe_exc[i + i_subfr * HIBND_ACB_L_FAC] = bwe_exc[i + i_subfr * HIBND_ACB_L_FAC - offset + (int) error];
                    }
                    error += (float) offset - (float) T0 * HIBND_ACB_L_FAC  - 0.25f * HIBND_ACB_L_FAC * (float) T0_frac;
                }
                else
                {
                    offset = T0 * 2 + (int) ((float) T0_frac * 0.5f + 4 + 0.5f) - 4;
                    for (i=0; i<L_SUBFR * 2; i++)
                    {
                        bwe_exc[i + i_subfr * 2] = bwe_exc[i + i_subfr * 2 - offset + (int) error];
                    }
                    error += (float) offset - (float) T0 * 2 - 0.5f * (float) T0_frac;
                }
            }

            pitch_buffer[i_subfr/L_SUBFR] = (float)T0 + (float)T0_frac/(float)T0_res;

            /*-------------------------------------------------------*
             * - Decode innovative codebook.                         *
             *-------------------------------------------------------*/

            if( st->use_partial_copy && ( st->rf_frame_type == RF_ALLPRED || ( st->rf_frame_type == RF_GENPRED && ( i_subfr == L_SUBFR || i_subfr == 3*L_SUBFR )) ) )
            {
                set_f(code, 0.0f, L_SUBFR);
            }
            else
            {
                config = PulseConfTable[acelp_cfg.fixed_cdk_index[i_subfr/L_SUBFR]];
                D_ACELP_indexing( code, config, NB_TRACK_FCB_4T, prm, &st->BER_detect );
                (prm) += 8;

                /*-------------------------------------------------------*
                * - Add the fixed-gain pitch contribution to code[].    *
                *-------------------------------------------------------*/

                cb_shape( acelp_cfg.pre_emphasis, acelp_cfg.pitch_sharpening, acelp_cfg.phase_scrambling, acelp_cfg.formant_enh,
                          acelp_cfg.formant_tilt, acelp_cfg.formant_enh_num, acelp_cfg.formant_enh_den, p_A, code, st->tilt_code, pitch_buf[i_subfr/L_SUBFR] );
            }

            /*-------------------------------------------------------*
             * - Generate Gaussian excitation                        *
             *-------------------------------------------------------*/

            if( acelp_cfg.gains_mode[i_subfr/L_SUBFR] == 7 && !st->use_partial_copy )
            {
                gaus_L2_dec( code2, st->tilt_code, p_A, acelp_cfg.formant_enh_num, &(st->seed_acelp) );
            }
            else
            {
                gain_code2 = 0.f;
                set_zero( code2, L_SUBFR );
            }

            /*-------------------------------------------------*
             * - Decode codebooks gains.                       *
             *-------------------------------------------------*/

            if( st->acelp_cfg.gains_mode[i_subfr/L_SUBFR] != 0 )
            {
                decode_acelp_gains( code, acelp_cfg.gains_mode[i_subfr/L_SUBFR], Es_pred, &gain_pit, &gain_code,
                                    &prm, &(st->past_gpit), &(st->past_gcode), &gain_inov, L_SUBFR, code2, &gain_code2 );
            }

            if( st->use_partial_copy && st->rf_frame_type == RF_ALLPRED )
            {
                st->past_gcode = 0.0f;
            }

            if( st->use_partial_copy && st->rf_frame_type == RF_NOPRED )
            {
                st->past_gpit = 0.004089f;
            }

            /*----------------------------------------------------------*
             * Update parameters for the next subframe.                 *
             * - tilt of code: 0.0 (unvoiced) to 0.5 (voiced)           *
             *----------------------------------------------------------*/

            st->tilt_code = est_tilt( exc+i_subfr, gain_pit, code, gain_code, &(st->voice_fac), L_SUBFR, acelp_cfg.voice_tilt );

            pgainT[i_subfr/L_SUBFR] = gain_pit;

            /*-------------------------------------------------------*
             * - Find the total excitation.                          *
             *-------------------------------------------------------*/

            gain_code_tmp = gain_code;
            gain_pit_tmp = gain_pit;
            if ( st->core == ACELP_CORE && st->last_core == ACELP_CORE && ( st->use_partial_copy || st->prev_use_partial_copy))
            {
                if ( i_subfr > 0 && gain_pit > 1.23f && st->prev_tilt_code_dec > 0.2f && st->next_coder_type == VOICED && (st->use_partial_copy || st->prev_use_partial_copy ) )
                {
                    gain_pit *= (0.8f - i_subfr/640.0f);
                }

                else if( !st->prev_use_partial_copy && st->last_coder_type == UNVOICED && st->next_coder_type != UNVOICED && gain_code < gain_code_pre)
                {
                    gain_code = 0.0f;
                }
            }

            gain_code_pre = gain_code;
            st->tilt_code_dec[i_subfr/L_SUBFR] = st->tilt_code;

            for( i = 0; i < L_SUBFR; i++ )
            {
                exc2[i+i_subfr] = gain_pit*exc[i+i_subfr];
                exc2[i+i_subfr] +=  gain_code2*code2[i];
                exc[i+i_subfr] = exc2[i+i_subfr] + gain_code*code[i];
            }

            /*-----------------------------------------------------------------*
             * Prepare TBE excitation
             *-----------------------------------------------------------------*/

            gain_code = gain_code_tmp;
            gain_pit = gain_pit_tmp;

            if( st->igf != 0 )
            {
                prep_tbe_exc( L_frame, i_subfr, gain_pit, gain_code, code, st->voice_fac, &voice_factors[i_subfr/L_SUBFR],
                              bwe_exc, gain_preQ, code_preQ, T0, coder_type, st->core_brate );
            }

            /*---------------------------------------------------------*
             * Enhance the excitation                                  *
             *---------------------------------------------------------*/

            enhancer( MODE2, -1, acelp_cfg.fixed_cdk_index[i_subfr/L_SUBFR], 0, coder_type, L_frame, st->voice_fac, stab_fac,
                      st->past_gcode, gain_inov, &(st->gc_threshold), code, &exc2[i_subfr], gain_pit, st->dispMem );

        } /* !RF_NELP frame partial copy */

        /*----------------------------------------------------------*
         * - compute the synthesis speech                           *
         *----------------------------------------------------------*/

        syn_filt( p_A, M,&exc2[i_subfr], &syn[i_subfr], L_SUBFR, mem_syn, 1 );

        /*-----------------------------------------------------------------*
         * update lp_filtered gains for the case of frame erasure
         *-----------------------------------------------------------------*/

        st->lp_gainp += weights[i_subfr/L_SUBFR] * st->past_gpit;
        st->lp_gainc += weights[i_subfr/L_SUBFR] * st->past_gcode;

        /*----------------------------------------------------------*
         * - update pitch lag for guided ACELP                      *
         *----------------------------------------------------------*/

        if( st->enableGplc && (i_subfr/L_SUBFR) == (L_frame/L_SUBFR)-1 )
        {
            st->T0_4th = T0;
        }

        /*----------------------------------------------------------*
         * - Update LPC coeffs                                      *
         *----------------------------------------------------------*/

        p_A += (M+1);
        pt_pitch++;
        /* copy current gain for next subframe use, in case there is no explicit encoding */
        prev_gain_pit = gain_pit;

    } /* end of subframe loop */

    if(st->BER_detect)
    {
        for (i=0; i<L_frame; i++)
        {
            exc[i] = 0.0f;
        }
        int_lsp( L_frame, st->old_lsp_q_cng, st->lsp_q_cng, st->Aq_cng, M, interpol_frac_12k8, 0 );

        p_A =st->Aq_cng;
        if(st->last_good < UNVOICED_TRANSITION )
        {
            mvr2r(st->mem_syn2, mem_syn, M );
        }
        else
        {
            set_zero( mem_syn, M );
        }

        for (i_subfr = 0; i_subfr < L_frame; i_subfr += L_SUBFR)
        {
            syn_filt(p_A, M,&exc[i_subfr], &syn[i_subfr], L_SUBFR, mem_syn, 1);
            p_A += (M+1);
        }
    }

    tmp = 0;
    pA = A + (st->nb_subfr-1)*(M+1);
    set_zero( h1, L_SUBFR+1 );
    set_zero( mem, M );
    h1[0] = 1.0f;
    syn_filt( pA, M,h1, h1, L_SUBFR, mem, 0);     /* impulse response of LPC     */
    deemph( h1, st->preemph_fac, L_SUBFR, &tmp);    /* impulse response of deemph  */
    /* impulse response level = gain introduced by synthesis+deemphasis */
    st->last_gain_syn_deemph = (float)sqrt(dotp( h1, h1, L_SUBFR) );

    /*-----------------------------------------------------------*
     * PLC: [ACELP: Fade-out]
     * PLC: update the background level
     *-----------------------------------------------------------*/

    /* Do the classification */
    FEC_clas_estim( syn, pitch_buf, st->L_frame, st->core_ext_mode, st->codec_mode, st->mem_syn_clas_estim, &st->clas_dec,
                    &st->lp_ener_bfi, st->core_brate, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -1.0f,
                    st->narrowBand, 0,
                    0, st->preemph_fac, st->tcxonly, st->last_core_brate );

    /* Update Pitch Lag memory */
    mvr2r( &st->old_pitch_buf[L_frame/L_SUBFR], st->old_pitch_buf, L_frame/L_SUBFR );
    mvr2r( pitch_buf, &st->old_pitch_buf[L_frame/L_SUBFR], L_frame/L_SUBFR );

    FEC_scale_syn( st->L_frame, st->clas_dec, st->last_good, syn, pitch_buf, st->enr_old, 0, coder_type,
                   LSF_Q_prediction, &st->scaling_flag, &st->lp_ener_FEC_av, &st->lp_ener_FEC_max, st->bfi,
                   st->total_brate, st->prev_bfi, st->last_core_brate, exc, exc2, A, &st->old_enr_LP, mem_back, mem_syn
                   , st->last_con_tcx && (st->L_frameTCX_past != st->L_frame) && (st->last_core != 0)
                   , (st->clas_dec == ONSET || (st->last_good >= VOICED_TRANSITION && st->last_good < INACTIVE_CLAS) ) );

    /* update ACELP synthesis memory */
    mvr2r( mem_syn, st->mem_syn2, M );
    mvr2r( syn+L_frame-L_SYN_MEM, st->mem_syn_r, L_SYN_MEM );

    tmp = st->syn[M];
    deemph( syn, st->preemph_fac, L_frame, &tmp );
    mvr2r( syn+L_frame-(L_frame/2), st->old_syn_Overl, L_frame/2 );
    mvr2r( syn+L_frame-M-1, st->syn, M+1 );
    mvr2r( syn, synth, L_frame );

    mvr2r( syn, st->old_core_synth, L_frame );


    /* update old_exc */
    mvr2r( exc_buf+L_frame, st->old_exc,  L_EXC_MEM_DEC );

    /* Output pitch parameters for bass post-filter */
    for (i_subfr = 0; i_subfr < L_frame; i_subfr += L_SUBFR )
    {
        *pT++ = (int)( pitch_buf[i_subfr/L_SUBFR] + 0.5f);
    }

    /* Update TCX-LTP */
    st->tcxltp_last_gain_unmodified = 0.f;

    /*Update MODE1*/
    mvr2r( p_A-(M+1), st->old_Aq_12_8, M+1 );
    st->old_Es_pred = Es_pred;

    st->tcxltp_third_last_pitch = st->tcxltp_second_last_pitch;
    st->tcxltp_second_last_pitch = st->old_fpitch;
    st->old_fpitch = pitch_buf[(L_frame/L_SUBFR) - 1];

    return;
}
