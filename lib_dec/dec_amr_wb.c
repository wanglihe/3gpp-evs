/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*---------------------------------------------------------------------*
 * decod_amr_wb()
 *
 * Decode excitation signal in AMR-WB IO mode
 *---------------------------------------------------------------------*/

void decod_amr_wb(
    Decoder_State *st,                  /* i/o: decoder static memory                     */
    const float *Aq,                  /* i  : LP filter coefficients                    */
    float *pitch_buf,           /* o  : floating pitch values for each subframe   */
    float *exc,                 /* i/o: adapt. excitation exc                     */
    float *exc2,                /* i/o: adapt. excitation/total exc               */
    short hf_gain[NB_SUBFR],    /* o  : decoded HF gain                           */
    float *voice_factors        /* o  : voicing factors                           */
    , float *gain_buf
)
{
    short T0, T0_frac, T0_min, T0_max;/* integer pitch variables                               */
    float gain_pit;                   /* pitch gain                                            */
    float gain_code;                  /* gain/normalized gain of the algebraic excitation      */
    float norm_gain_code;             /* normalized gain of the algebraic excitation           */
    float gain_inov;                  /* Innovation gain                                       */
    float voice_fac;                  /* voicing factor                                        */
    float code[L_SUBFR];              /* algebraic codevector                                  */
    const float *p_Aq;                /* Pointer to frame LP coefficient                       */
    float *pt_pitch;                  /* pointer to floating pitch                             */
    short i_subfr, i;                 /* tmp variables                                         */
    short pitch_limit_flag;

    /*------------------------------------------------------------------*
     * ACELP subframe loop
     *------------------------------------------------------------------*/

    p_Aq = Aq;                  /* pointer to interpolated LPC parameters */
    pt_pitch = pitch_buf;       /* pointer to the pitch buffer */
    st->lt_voice_fac = 0.0f;
    pitch_limit_flag = 0;  /* always restrained pitch Q range in IO mode */

    for( i_subfr = 0; i_subfr < L_FRAME; i_subfr += L_SUBFR )
    {
        /*----------------------------------------------------------------------*
         * Decode pitch lag
         *----------------------------------------------------------------------*/

        *pt_pitch = pit_decode( st, st->core_brate, 1, L_FRAME, i_subfr, -1, &pitch_limit_flag, &T0, &T0_frac, &T0_min, &T0_max, L_SUBFR );

        /*--------------------------------------------------------------*
         * Find the adaptive codebook vector
         *--------------------------------------------------------------*/

        pred_lt4( &exc[i_subfr], &exc[i_subfr], T0, T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP );

        /*--------------------------------------------------------------*
         * LP filtering of the adaptive excitation
         *--------------------------------------------------------------*/

        lp_filt_exc_dec( st, MODE1, st->core_brate, 1, -1, i_subfr, L_SUBFR, L_FRAME, 0,exc );

        /*--------------------------------------------------------------*
         * Innovation decoding
         *--------------------------------------------------------------*/

        inov_decode( st, st->core_brate, 1, L_FRAME, -1, 0, i_subfr, -1, p_Aq, st->tilt_code, *pt_pitch, code );

        /*--------------------------------------------------------------*
         * Gain decoding
         * Estimate spectrum tilt and voicing
         *--------------------------------------------------------------*/

        gain_dec_amr_wb( st, st->core_brate, &gain_pit, &gain_code, st->past_qua_en, &gain_inov, code, &norm_gain_code );

        /* update LP filtered gains for the case of frame erasures */
        lp_gain_updt( i_subfr, gain_pit, norm_gain_code, &st->lp_gainp, &st->lp_gainc, L_FRAME );

        st->tilt_code = est_tilt( exc+i_subfr, gain_pit, code, gain_code, &voice_fac, L_SUBFR, 0 );

        /*----------------------------------------------------------------------*
         * Find the total excitation
         *----------------------------------------------------------------------*/

        for (i = 0; i < L_SUBFR;  i++)
        {
            exc2[i+i_subfr] = gain_pit*exc[i+i_subfr];
            exc[i+i_subfr] = exc2[i+i_subfr] + gain_code*code[i];
        }

        /*----------------------------------------------------------------*
         * Excitation enhancements
         *----------------------------------------------------------------*/

        enhancer( MODE1, st->core_brate, -1, 1, -1, L_FRAME, voice_fac, st->stab_fac,
                  norm_gain_code, gain_inov, &st->gc_threshold, code, exc2 + i_subfr, gain_pit, st->dispMem );

        /*-----------------------------------------------------------------*
         * HF gain modification factors at 23.85 kbps
         *-----------------------------------------------------------------*/

        if ( st->core_brate == ACELP_23k85 )
        {
            hf_gain[i_subfr/L_SUBFR] = (short)get_next_indice( st,4);
        }

        voice_fac = VF_0th_PARAM + VF_1st_PARAM * voice_fac + VF_2nd_PARAM * voice_fac * voice_fac;
        voice_factors[i_subfr/L_SUBFR] = min( max(0.0f, voice_fac), 1.0f);

        p_Aq += (M+1);
        pt_pitch++;

        st->lt_voice_fac += 0.25f*voice_fac;
        gain_buf[i_subfr/L_SUBFR] = gain_pit;
    }

    return;
}
