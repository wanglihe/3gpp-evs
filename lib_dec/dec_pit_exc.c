/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*-------------------------------------------------------------------*
 * dec_pit_exc()
 *
 * Decode pitch-only contribution (used by the GSC technology)
 *-------------------------------------------------------------------*/

void dec_pit_exc(
    Decoder_State *st,                  /* i/o: decoder static memory                     */
    const short L_frame,              /* i  : length of the frame                       */
    const float *Aq,                  /* i  : LP filter coefficient                     */
    const float Es_pred,              /* i  : predicted scaled innov. energy            */
    float *pitch_buf,           /* o  : floating pitch values for each subframe   */
    float *code,                /* o  : innovation                                */
    float *exc,                 /* i/o: adapt. excitation exc                     */
    const short nb_subfr,             /* i  : Number of subframe considered             */
    float *gain_buf
)
{
    short T0, T0_frac, T0_min, T0_max;/* integer pitch variables                        */
    float gain_pit;            /* pitch gain                                            */
    float gain_code;           /* gain/normalized gain of the algebraic excitation      */
    float norm_gain_code;      /* normalized gain of the algebraic excitation           */
    float gain_inov;           /* Innovation gain                                       */
    float voice_fac;           /* voicing factor                                        */
    short L_subfr, pit_idx;
    const float *p_Aq;         /* Pointer to frame LP coefficient                       */
    float *pt_pitch;           /* pointer to floating pitch                             */
    short i_subfr, i;          /* tmp variables                                         */
    short Local_BR, Pitch_BR, Pitch_CT;
    short pitch_limit_flag;
    short nbits;
    float *pt_gain;            /* Pointer to floating gain values for each subframe     */

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
    pitch_limit_flag = 1;  /* always extended pitch Q range */

    /*------------------------------------------------------------------*
     * ACELP subframe loop
     *------------------------------------------------------------------*/

    L_subfr = st->L_frame/nb_subfr;
    p_Aq = Aq;                  /* pointer to interpolated LPC parameters */
    pt_pitch = pitch_buf;       /* pointer to the pitch buffer */
    pt_gain = gain_buf;         /* pointer to the gain buffer  */

    for ( i_subfr = 0; i_subfr < L_FRAME; i_subfr += L_subfr )
    {
        /*----------------------------------------------------------------------*
         *  Decode pitch lag
         *----------------------------------------------------------------------*/

        *pt_pitch = pit_decode( st, Pitch_BR, 0, L_frame, i_subfr, Pitch_CT, &pitch_limit_flag, &T0, &T0_frac, &T0_min, &T0_max, L_subfr );

        /*--------------------------------------------------------------*
         * Find the adaptive codebook vector.
         *--------------------------------------------------------------*/

        pred_lt4( &exc[i_subfr], &exc[i_subfr], T0, T0_frac, L_subfr+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP );

        /*--------------------------------------------------------------*
         * Innovation decoding
         *--------------------------------------------------------------*/

        if( st->GSC_noisy_speech )
        {
            inov_decode( st, Local_BR, 0, L_frame, LOCAL_CT, 1, i_subfr, -1, p_Aq, st->tilt_code, *pt_pitch, code );

            /*--------------------------------------------------------------*
             * Gain decoding
             * Estimate spectrum tilt and voicing
             *--------------------------------------------------------------*/

            gain_dec_mless( st, Local_BR, L_frame, LOCAL_CT, i_subfr, -1, code, Es_pred, &gain_pit, &gain_code, &gain_inov, &norm_gain_code );

            st->tilt_code = est_tilt( exc+i_subfr, gain_pit, code, gain_code, &voice_fac,L_SUBFR,0 );
        }
        else
        {
            nbits = 4;

            set_f(code, 0, L_SUBFR);
            gain_code = 0.0f;
            st->tilt_code = 0.0f;

            pit_idx = (short)get_next_indice( st, nbits );

            gain_pit = 0.5853f + dic_gp[pit_idx];

            if( st->BER_detect )  /* Bitstream is corrupted, use the past pitch gain */
            {
                gain_pit = st->lp_gainp;
            }
        }

        /*----------------------------------------------------------------------*
         * Find the total excitation
         *----------------------------------------------------------------------*/

        if( st->GSC_noisy_speech )
        {
            for (i = 0; i < L_subfr; i++)
            {
                exc[i+i_subfr] = gain_pit * exc[i+i_subfr] + gain_code * code[i];
            }
        }
        else
        {
            for (i = 0; i < L_subfr; i++)
            {
                exc[i+i_subfr] = gain_pit * exc[i+i_subfr];
            }
        }

        if( L_subfr == 2*L_SUBFR )
        {
            p_Aq += 2*(M+1);
            pt_pitch++;
            *pt_pitch = *(pt_pitch-1);
            pt_pitch++;
            *pt_gain = gain_pit;
            pt_gain++;
            *pt_gain = *(pt_gain-1);
            pt_gain++;

            if( i_subfr == 0 )
            {
                /* update gains for FEC - equivalent to lp_gain_updt() */
                st->lp_gainp = (3.0f/10.0f) * gain_pit;
                st->lp_gainc = 0;
            }
            else
            {
                /* update gains for FEC - equivalent to lp_gain_updt() */
                st->lp_gainp = (7.0f/10.0f) * gain_pit;
                st->lp_gainc = 0;
            }
        }
        else if( L_subfr == 4*L_SUBFR )
        {
            pt_pitch++;
            *pt_pitch = *(pt_pitch-1);
            pt_pitch++;
            *pt_pitch = *(pt_pitch-1);
            pt_pitch++;
            *pt_pitch = *(pt_pitch-1);
            pt_pitch++;
            *pt_gain = gain_pit;
            pt_gain++;
            *pt_gain = *(pt_gain-1);
            pt_gain++;
            *pt_gain = *(pt_gain-1);
            pt_gain++;
            *pt_gain = *(pt_gain-1);
            pt_gain++;
            p_Aq += 4*(M+1);

            /* update gains for FEC - equivalent to lp_gain_updt() */
            st->lp_gainp = gain_pit;
            st->lp_gainc = 0;
        }
        else
        {
            p_Aq += (M+1);
            pt_pitch++;
            *pt_gain = gain_pit;
            pt_gain++;

            lp_gain_updt( i_subfr, gain_pit, 0, &st->lp_gainp, &st->lp_gainc, L_frame );
        }
    }

    return;
}
