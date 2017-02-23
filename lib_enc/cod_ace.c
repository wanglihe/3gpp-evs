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
 * coder_acelp()
 *
 * Encode ACELP frame
 *-------------------------------------------------------------------*/

void coder_acelp(
    ACELP_config *acelp_cfg,      /* i/o: configuration of the ACELP  */
    const short coder_type,     /* i  : coding type                 */
    const float A[],            /* i  : coefficients 4xAz[M+1]      */
    const float Aq[],           /* i  : coefficients 4xAz_q[M+1]    */
    const float speech[],       /* i  : speech[-M..lg]              */
    float synth[],        /* o  : synthesis                   */
    LPD_state *LPDmem,         /* i/o: ACELP memories              */
    const float voicing[],      /* input: open-loop LTP gain        */
    const short T_op[],         /* input: open-loop LTP lag         */
    int *prm,           /* output: acelp parameters         */
    const float stab_fac,
    Encoder_State *st,            /* i/o : coder memory state         */
    HANDLE_PLC_ENC_EVS hPlc_Ext,
    const short target_bits,
    float *gain_pitch_buf,/* o  : gain pitch values           */
    float *gain_code_buf, /* o  : gain code values            */
    float *pitch_buf,     /* o  : pitch values for each subfr.*/
    float *voice_factors, /* o  : voicing factors             */
    float *bwe_exc        /* o  : excitation for SWB TBE      */
)
{
    short i, i_subfr;
    int T0, T0_min, T0_min_frac, T0_max, T0_max_frac, T0_res, T0_frac;
    float tmp, Es_pred;
    float gain_pit, gain_code, voice_fac;
    ACELP_CbkCorr g_corr;
    float g_corr2[6];
    const float *p_A, *p_Aq;
    float h1[L_SUBFR];        /* weighted impulse response of LP */
    float code[L_SUBFR];
    float cn[L_SUBFR];
    float xn[L_SUBFR];
    float xn2[L_SUBFR];
    float y1[L_SUBFR];        /* Filtered adaptive excitation       */
    float y2[L_SUBFR];        /* Filtered adaptive excitation       */
    float res_save;
    float exc_buf[L_EXC_MEM+L_FRAME16k+1], *exc;
    float exc2[L_SUBFR];
    float *syn,syn_buf[M+L_FRAME16k+L_FRAME16k/2];  /*128 for the memory, L_FRAME for the current synth and 128 for the ZIR for next TCX*/
    float syn2[L_FRAME16k];
    float norm_gain_code, gain_inov;
    short clip_gain;
    float gain_code2;
    float code2[L_SUBFR];
    float y22[L_SUBFR];        /* Filtered adaptive excitation       */
    float lp_select;
    float *pt_pitch, *pt_gain_pitch, *pt_gain_code;
    int offset;
    float error;
    float gain_preQ;                    /* Gain of prequantizer excitation   */
    float code_preQ[L_SUBFR];           /* Prequantizer excitation           */

    /* Configure ACELP */
    BITS_ALLOC_config_acelp( target_bits, coder_type, &(st->acelp_cfg), st->narrowBand, st->nb_subfr );

    /*------------------------------------------------------------------------*
     * Initialize buffers                                                     *
     *------------------------------------------------------------------------*/

    set_f( code_preQ, 0.f, L_SUBFR );
    gain_preQ = 0.0f;

    /* Reset phase dispersion */
    if( st->last_core > ACELP_CORE )
    {
        set_zero( st->dispMem, 8 );
    }

    /* set excitation memory*/
    exc = exc_buf + L_EXC_MEM;
    mvr2r( LPDmem->old_exc, exc_buf, L_EXC_MEM );
    *(exc+st->L_frame) = 0.f; /*to solve a warning*/

    /* Init syn buffer */
    syn = syn_buf + M;
    mvr2r( LPDmem->mem_syn, syn_buf, M );

    pt_pitch = pitch_buf;
    pt_gain_pitch = gain_pitch_buf;
    pt_gain_code = gain_code_buf;

    T0 = 0;
    T0_res = 0;
    T0_frac = 0;

    error = 0.0f;

    /*---------------------------------------------------------------*
     * Calculation of LP residual (filtering through A[z] filter)
     *---------------------------------------------------------------*/

    calc_residu( speech, exc, Aq, st->L_frame );

    /*------------------------------------------------------------------------*
     * Find and quantize mean_ener_code for gain quantizer                    *
     *------------------------------------------------------------------------*/

    if( acelp_cfg->nrg_mode > 0 )
    {
        Es_pred_enc( &Es_pred, prm, st->L_frame, L_SUBFR, exc, voicing, acelp_cfg->nrg_bits, acelp_cfg->nrg_mode>1 );
        prm++;
    }
    else
    {
        Es_pred = 0.f;
    }

    if( st->L_frame == L_FRAME )
    {
        mvr2r( Aq+2*(M+1), st->cur_sub_Aq, (M+1) );
    }
    else
    {
        mvr2r( Aq+3*(M+1), st->cur_sub_Aq, (M+1) );
    }

    /*------------------------------------------------------------------------*
     *          Loop for every subframe in the analysis frame                 *
     *------------------------------------------------------------------------*
     *  To find the pitch and innovation parameters. The subframe size is     *
     *  L_SUBFR and the loop is repeated L_FRAME_PLUS/L_SUBFR times.          *
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
    p_Aq = Aq;

    res_save = exc[0];

    for( i_subfr = 0; i_subfr < st->L_frame; i_subfr += L_SUBFR )
    {
        /* Restore exc[i_subfr] and save next exc[L_SUBFR+i_subfr] */
        exc[i_subfr] = res_save;
        res_save = exc[L_SUBFR+i_subfr];

        /*--------------------------------------------------------------------------*
         * Find target for pitch search (xn[]), target for innovation search (cn[]) *
         * and impulse response of the weighted synthesis filter (h1[]).            *
         *--------------------------------------------------------------------------*/

        find_targets( speech, &syn[i_subfr-M], i_subfr, &LPDmem->mem_w0, p_Aq,exc, L_SUBFR, p_A, st->preemph_fac, xn, cn, h1 );

        /*-----------------------------------------------------------------*
         * Gain clipping test to avoid unstable synthesis on frame erasure
         * or in case of floating point encoder & fixed p. decoder
         *-----------------------------------------------------------------*/

        clip_gain = gp_clip( st->core_brate, voicing, i_subfr, coder_type, xn, st->clip_var );

        /*-----------------------------------------------------------------*
         * - find unity gain pitch excitation (adaptive codebook entry)    *
         *   with fractional interpolation.                                *
         * - find filtered pitch exc. y1[]=exc[] convolved with h1[])      *
         * - compute pitch gain1                                           *
         *-----------------------------------------------------------------*/

        if ( acelp_cfg->ltp_bits != 0 )
        {
            /* pitch lag coding */
            Mode2_pit_encode( acelp_cfg->ltp_mode, i_subfr, &prm, &exc[i_subfr],
                              T_op, &T0_min, &T0_min_frac, &T0_max, &T0_max_frac, &T0, &T0_frac, &T0_res,
                              h1, xn, st->pit_min, st->pit_fr1, st->pit_fr1b, st->pit_fr2, st->pit_max, st->pit_res_max );

            /* find pitch excitation */
            if( st->pit_res_max == 6 )
            {
                if ( T0_res == (st->pit_res_max>>1) )
                {
                    pred_lt4( &exc[i_subfr], &exc[i_subfr], T0, T0_frac<<1, L_SUBFR+1, inter6_2, PIT_L_INTERPOL6_2, PIT_UP_SAMP6);
                }
                else
                {
                    pred_lt4( &exc[i_subfr], &exc[i_subfr], T0, T0_frac, L_SUBFR+1, inter6_2, PIT_L_INTERPOL6_2, PIT_UP_SAMP6);
                }
            }
            else
            {
                if( T0_res == (st->pit_res_max>>1) )
                {
                    pred_lt4( &exc[i_subfr], &exc[i_subfr], T0, T0_frac<<1, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP);
                }
                else
                {
                    pred_lt4( &exc[i_subfr], &exc[i_subfr], T0, T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP);
                }
            }

            /* filter adaptive codebook */
            lp_select =  lp_filt_exc_enc( MODE2, st->core_brate, 0, coder_type, i_subfr, exc, h1, xn, y1, xn2, L_SUBFR,
                                          st->L_frame, g_corr2, clip_gain, &(gain_pit), &(acelp_cfg->ltf_mode) );

            if( acelp_cfg->ltf_mode == NORMAL_OPERATION )
            {
                *prm = lp_select;
                prm++;
            }

            g_corr.y1y1 = g_corr2[0];
            g_corr.xy1 = -0.5f*(g_corr2[1]-0.01f)+0.01f;
        }
        else
        {
            /* No adaptive codebook (UC) */
            gain_pit=0.f;
            g_corr.xy1=0.f;
            g_corr.y1y1=0.f;
            set_zero( y1, L_SUBFR );
            set_zero( exc+i_subfr, L_SUBFR );
            T0 = L_SUBFR;
            T0_frac = 0;
            T0_res = 1;
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
                error += (float) offset - (float) T0 * 2  - 0.5f * (float) T0_frac;
            }
        }

        pitch_buf[i_subfr/L_SUBFR] = (float)T0 + (float)T0_frac/(float)T0_res;

        /*----------------------------------------------------------------------*
         *                 Encode the algebraic innovation                      *
         *----------------------------------------------------------------------*/

        E_ACELP_innovative_codebook(  exc, T0, T0_frac, T0_res, gain_pit, LPDmem->tilt_code,
                                      acelp_cfg->fixed_cdk_index[i_subfr/L_SUBFR],
                                      acelp_cfg->pre_emphasis,
                                      acelp_cfg->pitch_sharpening,
                                      acelp_cfg->phase_scrambling,
                                      acelp_cfg->formant_enh,
                                      acelp_cfg->formant_tilt,
                                      acelp_cfg->formant_enh_num,
                                      acelp_cfg->formant_enh_den,
                                      i_subfr, p_Aq, h1, xn, cn, y1, y2,
                                      st->acelp_autocorr, &prm, code
                                      ,st->L_frame, st->last_L_frame, st->total_brate
                                   );

        E_corr_xy2( xn, y1, y2, g_corr2, L_SUBFR );
        g_corr.y2y2 = 0.01F + g_corr2[2];
        g_corr.xy2 = 0.01F + -0.5f*g_corr2[3];
        g_corr.y1y2 = 0.01F + 0.5f*g_corr2[4];

        g_corr.xx = 0.01F + dotp( xn, xn, L_SUBFR );

        /*----------------------------------------------------------------------*
         *                 Add Gaussian excitation                              *
         *----------------------------------------------------------------------*/

        if(acelp_cfg->gains_mode[i_subfr/L_SUBFR]==7)
        {
            gauss_L2( h1, code2, y2, y22, &gain_code2, g_corr2, gain_pit,
                      LPDmem->tilt_code, p_Aq, acelp_cfg->formant_enh_num, &(st->seed_acelp) );

            g_corr.y1y1 = g_corr2[0];
            g_corr.y1y2 = g_corr2[4];
        }
        else
        {
            gain_code2 = 0.f;
            set_zero( code2, L_SUBFR );
            set_zero( y22, L_SUBFR );
        }

        /*----------------------------------------------------------*
         *  - Compute the fixed codebook gain                       *
         *  - quantize fixed codebook gain                          *
         *----------------------------------------------------------*/

        encode_acelp_gains( code, acelp_cfg->gains_mode[i_subfr/L_SUBFR], Es_pred, clip_gain, &g_corr, &gain_pit, &gain_code,
                            &prm, &norm_gain_code, &gain_inov, L_SUBFR, code2, &gain_code2, st->flag_noisy_speech_snr );

        gp_clip_test_gain_pit( st->core_brate, gain_pit, st->clip_var);

        /*----------------------------------------------------------*
         * - voice factor (for codebook tilt sharpening)            *
         *----------------------------------------------------------*/

        LPDmem->tilt_code = est_tilt( exc+i_subfr, gain_pit, code, gain_code, &voice_fac, L_SUBFR, acelp_cfg->voice_tilt );

        st->rf_tilt_buf[i_subfr/L_SUBFR] = LPDmem->tilt_code;


        /*-----------------------------------------------------------------*
         * Update memory of the weighting filter
         *-----------------------------------------------------------------*/

        LPDmem->mem_w0 = xn[L_SUBFR-1] - gain_pit * y1[L_SUBFR-1] - gain_code * y2[L_SUBFR-1] - gain_code2*y22[L_SUBFR-1];

        /*-------------------------------------------------------*
         * - Find the total excitation.                          *
         *-------------------------------------------------------*/

        for( i = 0; i < L_SUBFR; i++ )
        {
            exc2[i] = gain_pit*exc[i+i_subfr];
            exc2[i] +=  gain_code2*code2[i];
            exc[i+i_subfr] = exc2[i] + gain_code*code[i];
        }

        /*-----------------------------------------------------------------*
        * Prepare TBE excitation
        *-----------------------------------------------------------------*/

        prep_tbe_exc( st->L_frame, i_subfr, gain_pit, gain_code, code, voice_fac, &voice_factors[i_subfr/L_SUBFR],
                      bwe_exc, gain_preQ, code_preQ, T0, coder_type, st->core_brate );

        /*---------------------------------------------------------*
         * Enhance the excitation                                  *
         *---------------------------------------------------------*/

        enhancer( MODE2, -1, acelp_cfg->fixed_cdk_index[i_subfr/L_SUBFR], 0, coder_type, st->L_frame, voice_fac, stab_fac,
                  norm_gain_code, gain_inov, &(LPDmem->gc_threshold), code, exc2, gain_pit, st->dispMem );

        /*----------------------------------------------------------*
         * - compute the synthesis speech                           *
         *----------------------------------------------------------*/

        syn_filt( p_Aq, M, exc2, &syn2[i_subfr], L_SUBFR, LPDmem->mem_syn2, 1 );

        syn_filt( p_Aq, M, &exc[i_subfr], &syn[i_subfr], L_SUBFR, &syn[i_subfr-M], 0 );

        /*----------------------------------------------------------*
         * Save buffers for BPF                                     *
         *----------------------------------------------------------*/

        *pt_pitch = ((float)T0+(float)T0_frac/(float)T0_res+0.5f);
        *pt_gain_pitch = gain_pit;
        *pt_gain_code = gain_code;

        /*----------------------------------------------------------*
         * Update                                                   *
         *----------------------------------------------------------*/

        p_A += (M+1);
        p_Aq += (M+1);
        pt_pitch++;
        pt_gain_pitch++;
        pt_gain_code++;

        if( hPlc_Ext != NULL )
        {
            hPlc_Ext->T0_4th = T0;
        }

    } /* end of subframe loop */

    p_A -= (M+1);
    p_Aq -= (M+1);


    /*----------------------------------------------------------*
     * Update LPD memory                                        *
     *----------------------------------------------------------*/

    mvr2r( exc + st->L_frame - L_EXC_MEM, LPDmem->old_exc, L_EXC_MEM );
    mvr2r( syn + st->L_frame - M, LPDmem->mem_syn, M );
    mvr2r( syn + st->L_frame - L_SYN_MEM, LPDmem->mem_syn_r, L_SYN_MEM );

    if( hPlc_Ext != NULL )
    {
        mvr2r( exc + st->L_frame - L_EXC_MEM - 8, hPlc_Ext->old_exc, 8 );
    }

    /*----------------------------------------------------------*
     * ZIR at the end of the ACELP frame (for TCX)              *
     *----------------------------------------------------------*/

    mvr2r( syn2, syn, st->L_frame );
    tmp = LPDmem->syn[M];
    deemph( syn, st->preemph_fac, st->L_frame, &tmp );
    mvr2r( syn + st->L_frame/2, LPDmem->Txnq, st->L_frame/2 );
    mvr2r( syn + st->L_frame - (M+1), LPDmem->syn, M+1 );
    mvr2r( syn, synth, st->L_frame );


    /*Update MODE1*/
    mvr2r( p_Aq, st->old_Aq_12_8, M+1 );
    st->old_Es_pred = Es_pred;

    return;
}
