/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * enc_pit_exc()
 *
 * Encode pitch only contribution
 *-------------------------------------------------------------------*/

void enc_pit_exc(
    Encoder_State *st,                   /* i/o: State structure                                  */
    LPD_state *mem,                  /* i/o: encoder memories                                 */
    const float *speech,               /* i  : Input speech                                     */
    const float Aw[],                  /* i  : weighted A(z) unquantized for subframes          */
    const float *Aq,                   /* i  : 12k8 Lp coefficient                              */
    const float Es_pred,               /* i  : predicted scaled innov. energy                   */
    const short *T_op,                 /* i  : open loop pitch                                  */
    const float *voicing,              /* i  : voicing                                          */
    const float *res,                  /* i  : residual signal                                  */
    float *synth,                /* i/o: core synthesis                                   */
    float *exc,                  /* i/o: current non-enhanced excitation                  */
    short *T0,                   /* i/o: close loop integer pitch                         */
    short *T0_frac,              /* i/o: close-loop pitch period - fractional part        */
    float *pitch_buf,            /* i/o: Fractionnal per subframe pitch                   */
    const short nb_subfr,              /* i  : Number of subframe considered                    */
    float *gpit                  /* o  : pitch mean gpit                                  */
)
{
    float xn[PIT_EXC_L_SUBFR];          /* Target vector for pitch search    */
    float xn2[PIT_EXC_L_SUBFR];         /* Target vector for codebook search */
    float cn[PIT_EXC_L_SUBFR];          /* Target vector in residual domain  */
    float h1[PIT_EXC_L_SUBFR+(M+1)];    /* Impulse response vector           */
    float y1[PIT_EXC_L_SUBFR];          /* Filtered adaptive excitation      */
    float code[L_SUBFR];                /* Fixed codebook excitation         */
    float y2[L_SUBFR];                  /* Filtered algebraic excitation     */
    float voice_fac;                    /* Voicing factor                    */
    float gain_code;                    /* Gain of code                      */
    float gain_inov;                    /* inovation gain                    */
    float gain_pit;                     /* Pitch gain                        */
    short pit_idx, i_subfr;             /* tmp variables                     */
    short T0_min, T0_max;               /* pitch variables                   */
    float g_corr[10];                   /* ACELP correlation values + gain pitch */
    short clip_gain, i;                 /* LSF clip gain and LP flag         */
    const float *p_Aw, *p_Aq;           /* pointer to LP filter coefficient vector */
    float *pt_pitch;                    /* pointer to floating pitch         */
    short L_subfr;
    float cum_gpit, gpit_tmp;
    short Local_BR, Pitch_BR, Pitch_CT;
    short unbits_PI = 0;                /* saved bits for EVS_PI             */
    float norm_gain_code;
    short pitch_limit_flag;
    short lp_select, lp_flag;

    /*------------------------------------------------------------------*
     * Initialization
     *------------------------------------------------------------------*/

    pitch_limit_flag = 1;  /* always extended pitch Q range */

    if( st->GSC_noisy_speech )
    {
        Local_BR = ACELP_7k20;
        Pitch_CT = GENERIC;
        Pitch_BR = ACELP_7k20;
    }
    else
    {
        Local_BR = ACELP_7k20;
        Pitch_CT = AUDIO;
        Pitch_BR = st->core_brate;
    }

    gain_code = 0;

    T0_max = PIT_MAX;
    T0_min = PIT_MIN;

    cum_gpit = 0.0f;
    L_subfr = L_FRAME/nb_subfr;


    /*------------------------------------------------------------------*
     * ACELP subframe loop
     *------------------------------------------------------------------*/

    p_Aw = Aw;
    p_Aq = Aq;
    pt_pitch = pitch_buf;       /* pointer to the pitch buffer */
    for( i_subfr = 0; i_subfr < L_FRAME; i_subfr += L_subfr )
    {

        /*----------------------------------------------------------------*
         * Find the the excitation search target "xn" and innovation
         *   target in residual domain "cn"
         * Compute impulse response, h1[], of weighted synthesis filter
         *----------------------------------------------------------------*/

        mvr2r( &res[i_subfr], &exc[i_subfr], L_subfr );

        if( L_subfr == L_SUBFR )
        {
            find_targets( speech, st->mem_syn_tmp, i_subfr, &st->mem_w0_tmp, p_Aq, res, L_subfr, p_Aw, st->preemph_fac, xn, cn, h1 );
        }
        else
        {
            find_targets( speech, st->mem_syn_tmp, i_subfr, &st->mem_w0_tmp, p_Aq, res, L_subfr, p_Aw, st->preemph_fac, xn, NULL, h1 );
        }

        /*----------------------------------------------------------------*
         * Close-loop pitch search and quantization
         * Adaptive exc. construction
         *----------------------------------------------------------------*/

        *pt_pitch = pit_encode( st, Pitch_BR, 0, L_FRAME, Pitch_CT, &pitch_limit_flag, i_subfr, exc,
                                L_subfr, T_op, &T0_min, &T0_max, T0, T0_frac, h1, xn );

        /*-----------------------------------------------------------------*
         * Find adaptive exitation
         *-----------------------------------------------------------------*/

        pred_lt4( &exc[i_subfr], &exc[i_subfr], *T0, *T0_frac, L_subfr+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP);

        /*-----------------------------------------------------------------*
         * Gain clipping test to avoid unstable synthesis on frame erasure
         * or in case of floating point encoder & fixed p. decoder
         *-----------------------------------------------------------------*/

        clip_gain = gp_clip( st->core_brate, voicing, i_subfr, AUDIO, xn, st->clip_var );

        /*-----------------------------------------------------------------*
         * Codebook target computation
         * (No LP filtering of the adaptive excitation)
         *-----------------------------------------------------------------*/

        lp_select = lp_filt_exc_enc( MODE1, st->core_brate, 0, AUDIO, i_subfr, exc, h1, xn, y1, xn2, L_subfr,
                                     L_FRAME, g_corr, clip_gain, &gain_pit, &lp_flag );

        if( lp_flag == NORMAL_OPERATION )
        {
            push_indice( st, IND_LP_FILT_SELECT, lp_select, 1 );
        }

        /* update long-term pitc hgain for speech/music classifier */
        st->lowrate_pitchGain = 0.9f * st->lowrate_pitchGain + 0.1f * gain_pit;

        gpit_tmp =  gain_pit;

        if( st->GSC_noisy_speech == 0 || L_subfr != L_SUBFR )
        {
            pit_idx = (short) vquant(&gain_pit, mean_gp, &gain_pit, dic_gp, 1, 16);
            push_indice( st, IND_PIT_IDX, pit_idx, 4 );
        }

        if( st->GSC_noisy_speech && L_subfr == L_SUBFR )
        {
            /*-----------------------------------------------------------------*
             * Innovation & gain encoding
             *-----------------------------------------------------------------*/

            inov_encode( st, Local_BR, 0, L_FRAME, st->last_L_frame, LOCAL_CT, WB, 1, i_subfr, -1,
                         p_Aq, gain_pit, cn, exc, h1, mem->tilt_code, *pt_pitch, xn2, code, y2, &unbits_PI );

            gain_enc_mless( st, Local_BR, L_FRAME, LOCAL_CT, i_subfr, -1, xn, y1, y2, code, Es_pred,
                            &gain_pit, &gain_code, &gain_inov, &norm_gain_code, g_corr, clip_gain );
        }

        gp_clip_test_gain_pit( st->core_brate, gain_pit, st->clip_var );

        if( st->GSC_noisy_speech )
        {
            mem->tilt_code = est_tilt( exc+i_subfr, gain_pit, code, gain_code, &voice_fac, L_SUBFR, 0 );
        }
        else
        {
            mem->tilt_code  = 0.0f;
        }

        /*-----------------------------------------------------------------*
         * Update memory of the weighting filter
         *-----------------------------------------------------------------*/

        if( st->GSC_noisy_speech )
        {
            st->mem_w0_tmp = xn[L_subfr-1] - (gain_pit*y1[L_subfr-1]) - (gain_code*y2[L_subfr-1]);
        }
        else
        {
            st->mem_w0_tmp = xn[L_subfr-1] - (gain_pit*y1[L_subfr-1]);
        }

        /*-----------------------------------------------------------------*
         * Construct adaptive part of the excitation
         * Save the non-enhanced excitation for FEC_exc
         *-----------------------------------------------------------------*/

        if( st->GSC_noisy_speech )
        {
            for ( i = 0; i < L_subfr; i++ )
            {
                exc[i+i_subfr] = gain_pit * exc[i+i_subfr] + gain_code * code[i];
            }
        }
        else
        {
            for ( i = 0; i < L_subfr; i++ )
            {
                exc[i+i_subfr] = gain_pit * exc[i+i_subfr];
            }
        }

        /*-----------------------------------------------------------------*
         * Synthesize speech to update mem_syn[].
         * Update A(z) filters
         *-----------------------------------------------------------------*/

        syn_filt( p_Aq, M, &exc[i_subfr], &synth[i_subfr], L_subfr, st->mem_syn_tmp, 1 );

        if( L_subfr == 2*L_SUBFR )
        {
            if( i_subfr == 0 )
            {
                cum_gpit = gpit_tmp*.5f;
            }
            else
            {
                cum_gpit += gpit_tmp*.5f;
            }

            p_Aw += 2*(M+1);
            p_Aq += 2*(M+1);
            pt_pitch++;
            *pt_pitch = *(pt_pitch-1);
            pt_pitch++;
        }
        else if( L_subfr == 4*L_SUBFR )
        {
            cum_gpit = gpit_tmp;

            pt_pitch++;
            *pt_pitch = *(pt_pitch-1);
            pt_pitch++;
            *pt_pitch = *(pt_pitch-1);
            pt_pitch++;
            *pt_pitch = *(pt_pitch-1);
            pt_pitch++;

            p_Aw += 4*(M+1);
            p_Aq += 4*(M+1);
        }
        else
        {
            if( i_subfr == 0 )
            {
                cum_gpit = gpit_tmp*.25f;
            }
            else
            {
                cum_gpit += gpit_tmp*.25f;
            }

            pt_pitch++;
            p_Aw += (M+1);
            p_Aq += (M+1);
        }
    }

    *gpit = 0.1f * *gpit + 0.9f * cum_gpit;

    return;
}
