/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prot.h"

/*-------------------------------------------------------------------*
 * decod_unvoiced()
 *
 * Decode unvoiced (UC) frames
 *-------------------------------------------------------------------*/

void decod_unvoiced(
    Decoder_State *st,            /* i/o: decoder static memory                   */
    const float *Aq,            /* i  : LP filter coefficient                   */
    const short coder_type,     /* i  : coding type                             */
    float *tmp_noise,     /* o  : long term temporary noise energy        */
    float *pitch_buf,     /* o  : floating pitch values for each subframe */
    float *voice_factors, /* o  : voicing factors                         */
    float *exc,           /* o  : adapt. excitation exc                   */
    float *exc2,          /* o  : adapt. excitation/total exc             */
    float *bwe_exc        /* i/o: excitation for SWB TBE                  */
    , float *gain_buf
)
{
    float gain_pit = 0;         /* Quantized pitch gain                         */
    float gain_code;            /* Quantized algebraic codeebook gain           */
    float gain_inov;            /* inovation gain                               */
    float norm_gain_code;       /* normalized algebraic codeebook gain          */
    float voice_fac;            /* Voicing factor                               */
    float code[L_SUBFR];        /* algebraic codevector                         */
    short i_subfr;
    float *pt_pitch;
    const float *p_Aq;

    if( st->last_ppp_mode_dec == 1 || st->last_nelp_mode_dec == 1 )
    {
        /* SC_VBR - reset the decoder, to avoid memory not updated issue for this unrealistic case */
        CNG_reset_dec( st, pitch_buf, voice_factors );
    }

    p_Aq = Aq;                  /* pointer to interpolated LPC parameters       */
    pt_pitch = pitch_buf;       /* pointer to the pitch buffer                  */

    for( i_subfr=0; i_subfr<L_FRAME; i_subfr+=L_SUBFR )
    {
        /*----------------------------------------------------------------*
         * Unvoiced subframe processing
         *----------------------------------------------------------------*/

        gaus_dec( st, st->core_brate, i_subfr, code, &norm_gain_code,
                  &st->lp_gainp, &st->lp_gainc, &gain_inov, &st->tilt_code,
                  &voice_fac, &gain_pit, pt_pitch, exc, &gain_code, exc2 );

        *tmp_noise = norm_gain_code;

        /*----------------------------------------------------------------*
         * Excitation enhancements (update of total excitation signal)
         *----------------------------------------------------------------*/

        enhancer( MODE1, st->core_brate, -1, 0, coder_type, L_FRAME, voice_fac, st->stab_fac,
                  norm_gain_code, gain_inov, &st->gc_threshold, code, exc2 + i_subfr, gain_pit, st->dispMem );

        voice_factors[i_subfr/L_SUBFR] = 0.0f;

        interp_code_5over2( &exc[i_subfr], &bwe_exc[i_subfr * HIBND_ACB_L_FAC], L_SUBFR );

        p_Aq += (M+1);
        pt_pitch++;
        st->tilt_code_dec[i_subfr/L_SUBFR] = st->tilt_code;
    }

    /* SC-VBR */
    st->prev_gain_pit_dec = gain_pit;

    set_f( gain_buf, 0, NB_SUBFR );

    return;
}
