/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*-------------------------------------------------------------------*
 * encod_gen_voic()
 *
 * Encode excitation signal
 *-------------------------------------------------------------------*/

void encod_gen_voic(
    Encoder_State *st,                    /* i/o: state structure                         */
    LPD_state *mem,                   /* i/o: encoder memories                        */
    const short L_frame,                /* i  : length of the frame                     */
    const short sharpFlag,              /* i  : formant sharpening flag                 */
    const float speech[],               /* i  : input speech                            */
    const float Aw[],                   /* i  : weighted A(z) unquantized for subframes */
    const float Aq[],                   /* i  : LP coefficients                         */
    const short coder_type,             /* i  : coding type                             */
    const float Es_pred,                /* i  : predicted scaled innov. energy          */
    const short T_op[],                 /* i  : open loop pitch                         */
    const float voicing[],              /* i  : voicing                                 */
    const float *res,                   /* i  : residual signal                         */
    float *syn,                   /* i/o: core synthesis                          */
    float *exc,                   /* i/o: current non-enhanced excitation         */
    float *exc2,                  /* i/o: current enhanced excitation             */
    float *pitch_buf,             /* i/o: floating pitch values for each subframe */
    float *voice_factors,         /* o  : voicing factors                         */
    float *bwe_exc,               /* o  : excitation for SWB TBE                  */
    short *unbits                 /* i/o: number of unused bits                   */
)
{
    float xn[L_SUBFR];            /* Target vector for pitch search    */
    float xn2[L_SUBFR];           /* Target vector for codebook search */
    float cn[L_SUBFR];            /* Target vector in residual domain  */
    float h1[L_SUBFR+(M+1)];      /* Impulse response vector           */
    float code[L_SUBFR];          /* Fixed codebook excitation         */
    float y1[L_SUBFR];            /* Filtered adaptive excitation      */
    float y2[L_SUBFR];            /* Filtered algebraic excitation     */
    float gain_pit;               /* Pitch gain                        */
    float voice_fac;              /* Voicing factor                    */
    float gain_code;              /* Gain of code                      */
    float gain_inov;              /* inovation gain                    */
    short i, i_subfr;             /* tmp variables                     */
    short T0, T0_frac;            /* close loop integer pitch and fractional part */
    short T0_min, T0_max;         /* pitch variables                   */
    float *pt_pitch;              /* pointer to floating pitch buffer  */
    float g_corr[6];              /* ACELP correl, values + gain pitch */
    float gains_mem[2*(NB_SUBFR-1)]; /* pitch gain and code gain from previous subframes */
    short clip_gain;              /* LSF clip gain                     */
    const float *p_Aw, *p_Aq;     /* pointer to LP filter coeff. vector*/
    float error;
    int offset;
    float gain_preQ;              /* Gain of prequantizer excitation   */
    float code_preQ[L_SUBFR];     /* Prequantizer excitation           */
    short unbits_PI;              /* number of unused bits for EVS_PI  */
    float norm_gain_code;
    short pitch_limit_flag;
    short harm_flag_acelp;
    short lp_select, lp_flag;

    /*------------------------------------------------------------------*
     * Initializations
     *------------------------------------------------------------------*/

    gain_pit = 0;
    gain_code = 0;
    gain_preQ = 0;
    unbits_PI = 0;
    error = 0.0f;

    if( L_frame == L_FRAME )
    {
        T0_max = PIT_MAX;
        T0_min = PIT_MIN;
    }
    else /* L_frame == L_FRAME16k */
    {
        T0_max = PIT16k_MAX;
        T0_min = PIT16k_MIN;
    }

    p_Aw = Aw;
    p_Aq = Aq;
    pt_pitch = pitch_buf;
    gain_preQ = 0;
    set_f( code_preQ, 0, L_SUBFR );

    /* set and write harmonicity flag */
    harm_flag_acelp = 0;

    if( st->core_brate > ACELP_24k40 && st->core_brate <= ACELP_32k && coder_type == GENERIC )
    {
        if( st->last_harm_flag_acelp > 2 )
        {
            harm_flag_acelp = 1;
        }

        push_indice( st, IND_HARM_FLAG_ACELP, harm_flag_acelp, 1 );
    }

    /*------------------------------------------------------------------*
     * ACELP subframe loop
     *------------------------------------------------------------------*/

    for( i_subfr=0; i_subfr<L_frame; i_subfr+=L_SUBFR )
    {
        /*----------------------------------------------------------------*
         * Find the the excitation search target "xn" and innovation
         *   target in residual domain "cn"
         * Compute impulse response, h1[], of weighted synthesis filter
         *----------------------------------------------------------------*/

        mvr2r( &res[i_subfr], &exc[i_subfr], L_SUBFR );

        find_targets( speech, mem->mem_syn, i_subfr, &mem->mem_w0, p_Aq, res, L_SUBFR, p_Aw, st->preemph_fac, xn, cn, h1 );

        /*----------------------------------------------------------------*
         * Close-loop pitch search and quantization
         * Adaptive exc. construction
         *----------------------------------------------------------------*/

        *pt_pitch = pit_encode( st, st->core_brate, 0, L_frame, coder_type, &pitch_limit_flag, i_subfr,
                                exc, L_SUBFR, T_op, &T0_min, &T0_max, &T0, &T0_frac, h1, xn );

        if( L_frame == L_FRAME )
        {
            offset = tbe_celp_exc_offset( T0, T0_frac );
            for (i=0; i<L_SUBFR * HIBND_ACB_L_FAC; i++)
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

        /*-----------------------------------------------------------------*
         * Find adaptive exitation
         *-----------------------------------------------------------------*/

        pred_lt4( &exc[i_subfr], &exc[i_subfr], T0, T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP);

        /*-----------------------------------------------------------------*
         * Gain clipping test to avoid unstable synthesis on frame erasure
         *-----------------------------------------------------------------*/

        clip_gain = gp_clip( st->core_brate, voicing, i_subfr, coder_type, xn, st->clip_var );

        if( coder_type == INACTIVE )
        {
            /* in case of AVQ inactive, limit the gain to 0.65 */
            clip_gain = 2;
        }

        /*-----------------------------------------------------------------*
         * LP filtering of the adaptive excitation, codebook target computation
         *-----------------------------------------------------------------*/

        lp_select = lp_filt_exc_enc( MODE1, st->core_brate, 0, coder_type, i_subfr, exc, h1, xn, y1, xn2, L_SUBFR,
                                     L_frame, g_corr, clip_gain, &gain_pit, &lp_flag );

        if( lp_flag == NORMAL_OPERATION )
        {
            push_indice( st, IND_LP_FILT_SELECT, lp_select, 1 );
        }

        /* update long-term pitch gain for speech/music classifier */
        st->lowrate_pitchGain = 0.9f * st->lowrate_pitchGain + 0.1f * gain_pit;

        /*-----------------------------------------------------------------*
         * Transform-domain contribution (active frames)
         *-----------------------------------------------------------------*/

        if( st->core_brate > ACELP_24k40 && coder_type != INACTIVE )
        {
            transf_cdbk_enc( st, st->core_brate, st->extl, coder_type, harm_flag_acelp, i_subfr, -1, cn, exc,
                             p_Aq, p_Aw, h1, xn, xn2, y1, y2, Es_pred, &gain_pit, gain_code, g_corr, clip_gain,
                             &(st->mem_deemp_preQ), &(st->mem_preemp_preQ), &gain_preQ, code_preQ, unbits );
        }

        /*-----------------------------------------------------------------*
         * Innovation encoding
         *-----------------------------------------------------------------*/

        inov_encode( st, st->core_brate, 0, L_frame, st->last_L_frame, coder_type, st->bwidth, sharpFlag, i_subfr, -1, p_Aq,
                     gain_pit, cn, exc, h1, mem->tilt_code, *pt_pitch, xn2, code, y2, &unbits_PI );

        /*-----------------------------------------------------------------*
         * Gain encoding
         *-----------------------------------------------------------------*/

        if ( st->core_brate <= ACELP_8k00 )
        {
            gain_enc_lbr( st, st->core_brate, coder_type, i_subfr, xn, y1, y2, code,
                          &gain_pit, &gain_code, &gain_inov, &norm_gain_code, g_corr, gains_mem, clip_gain );
        }
        else if ( st->core_brate > ACELP_32k )
        {
            gain_enc_SQ( st, st->core_brate, coder_type, i_subfr, -1, xn, y1, y2, code, Es_pred,
                         &gain_pit, &gain_code, &gain_inov, &norm_gain_code, g_corr, clip_gain );
        }
        else
        {
            gain_enc_mless( st, st->core_brate, L_frame, coder_type, i_subfr, -1, xn, y1, y2, code, Es_pred,
                            &gain_pit, &gain_code, &gain_inov, &norm_gain_code, g_corr, clip_gain );
        }

        if ( st->last_ppp_mode == 1 )
        {
            /* SC-VBR - all other st->clip_var values will be updated even in a PPP frame */
            st->clip_var[1] = gain_pit;
        }

        /*-----------------------------------------------------------------*
         * update LP-filtered gains for the case of frame erasures
         *-----------------------------------------------------------------*/

        gp_clip_test_gain_pit( st->core_brate, gain_pit, st->clip_var);
        mem->tilt_code = est_tilt( exc+i_subfr, gain_pit, code, gain_code, &voice_fac, L_SUBFR, 0 );

        /*-----------------------------------------------------------------*
         * Transform-domain contribution (inactive frames)
         *-----------------------------------------------------------------*/

        if ( st->core_brate > ACELP_24k40 && coder_type == INACTIVE )
        {
            transf_cdbk_enc( st, st->core_brate, st->extl, coder_type, 0, i_subfr, -1, cn, exc,
                             p_Aq, p_Aw, h1, xn, xn2, y1, y2, Es_pred, &gain_pit, gain_code, g_corr, clip_gain,
                             &(st->mem_deemp_preQ), &(st->mem_preemp_preQ), &gain_preQ, code_preQ, unbits );
        }

        /*-----------------------------------------------------------------*
         * Update memory of the weighting filter
         *-----------------------------------------------------------------*/

        mem->mem_w0 = xn[L_SUBFR-1] - (gain_pit*y1[L_SUBFR-1]) - (gain_code*y2[L_SUBFR-1]);

        /*-----------------------------------------------------------------*
         * Construct adaptive part of the excitation
         * Save the non-enhanced excitation for FEC_exc
         *-----------------------------------------------------------------*/

        for ( i = 0; i < L_SUBFR;  i++ )
        {
            exc2[i+i_subfr] = gain_pit * exc[i+i_subfr];
            exc[i+i_subfr] = exc2[i+i_subfr] + gain_code * code[i];
        }

        /*-----------------------------------------------------------------*
         * Add the ACELP pre-quantizer contribution
         *-----------------------------------------------------------------*/
        if( st->core_brate > ACELP_24k40 )
        {
            for ( i = 0; i < L_SUBFR; i++ )
            {
                exc2[i+i_subfr] += gain_preQ * code_preQ[i];
                exc[i+i_subfr] += gain_preQ * code_preQ[i];
            }
        }

        /*-----------------------------------------------------------------*
         * Prepare TBE excitation
         *-----------------------------------------------------------------*/

        prep_tbe_exc( L_frame, i_subfr, gain_pit, gain_code, code, voice_fac, &voice_factors[i_subfr/L_SUBFR],
                      bwe_exc, gain_preQ, code_preQ, T0, coder_type, st->core_brate );


        /*-----------------------------------------------------------------*
         * Synthesize speech to update mem_syn[].
         * Update A(z) filters
         *-----------------------------------------------------------------*/

        syn_filt( p_Aq, M, &exc[i_subfr], &syn[i_subfr], L_SUBFR, mem->mem_syn, 1 );

        p_Aw += (M+1);
        p_Aq += (M+1);
        pt_pitch++;
    }

    /* write reserved bits */
    while( unbits_PI > 0 )
    {
        i = min(unbits_PI, 16);
        push_indice( st, IND_UNUSED, 0, i );
        unbits_PI -= i;
    }

    /* SC-VBR */
    st->prev_ppp_gain_pit = gain_pit;
    st->prev_tilt_code = mem->tilt_code;

    return;
}
