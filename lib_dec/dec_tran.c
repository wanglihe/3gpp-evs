/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * decod_tran()
 *
 * Decode transition (TC) frames
 *-------------------------------------------------------------------*/

void decod_tran(
    Decoder_State *st,                 /* i/o: decoder static memory                  */
    const short L_frame,             /* i  : length of the frame                    */
    const short tc_subfr,            /* i  : TC subframe index                      */
    const float *Aq,                 /* i  : LP filter coefficient                  */
    const short coder_type,          /* i  : coding type                            */
    const float Es_pred,             /* i  : predicted scaled innov. energy         */
    float *pitch_buf,          /* o  : floating pitch values for each subframe*/
    float *voice_factors,      /* o  : voicing factors                        */
    float *exc,                /* i/o: adapt. excitation exc                  */
    float *exc2,               /* i/o: adapt. excitation/total exc            */
    float *bwe_exc,            /* i/o: excitation for SWB TBE                 */
    short *unbits,             /* i/o: number of unused bits                  */
    const short sharpFlag,           /* i  : formant sharpening flag                */
    float *gain_buf
)
{
    short T0, T0_frac, T0_min, T0_max; /* integer pitch variables               */
    float gain_code;             /* Quantized algebraic codeebook gain          */
    float norm_gain_code;        /* normalized algebraic codeebook gain         */
    float gain_pit = 0;          /* Quantized pitch gain                        */
    float voice_fac;             /* Voicing factor                              */
    float gain_inov;             /* inovation gain                              */
    float code[L_SUBFR];         /* algebraic codevector                        */
    const float *p_Aq;           /* pointer to lp filter coefficient            */
    float *pt_pitch;             /* pointer to floating pitch                   */
    short i_subfr, i;            /* tmp variables                               */
    short position;              /* TC related flag                             */
    float gain_preQ = 0;         /* Gain of prequantizer excitation             */
    float code_preQ[L_SUBFR];    /* Prequantizer excitation                     */
    short Jopt_flag;             /* flag indicating zero adaptive contribtuion  */
    float norm_gain_preQ;

    gain_preQ = 0;
    set_f( code_preQ, 0, L_SUBFR );
    /*----------------------------------------------------------------*
     * ACELP subframe loop
     *----------------------------------------------------------------*/

    p_Aq = Aq;
    pt_pitch = pitch_buf;
    Jopt_flag = 0;
    norm_gain_preQ = 0.0f;

    for (i_subfr = 0; i_subfr < L_frame; i_subfr += L_SUBFR)
    {
        /*------------------------------------------------------------*
         * TC : subframe determination &
         * adaptive/glottal part of excitation construction
         *------------------------------------------------------------*/

        transition_dec( st, st->core_brate, 0, L_frame, i_subfr, coder_type, tc_subfr, &Jopt_flag, exc,
                        &T0, &T0_frac, &T0_min, &T0_max, &pt_pitch, &position, bwe_exc );

        /*-----------------------------------------------------------------*
         * Transform domain contribution decoding - active frames
         *-----------------------------------------------------------------*/

        if( st->core_brate > ACELP_24k40 )
        {
            transf_cdbk_dec( st, st->core_brate, coder_type, 0, i_subfr, tc_subfr, Es_pred, 0,
                             &st->mem_preemp_preQ, &gain_preQ, &norm_gain_preQ, code_preQ, unbits );
        }

        /*-----------------------------------------------------------------*
         * ACELP codebook search + pitch sharpening
         *-----------------------------------------------------------------*/

        inov_decode( st, st->core_brate, 0, L_frame, coder_type, sharpFlag, i_subfr, tc_subfr, p_Aq, st->tilt_code, *pt_pitch, code );

        /*-----------------------------------------------------------------*
         * De-quantize the gains
         * Update tilt of code: 0.0 (unvoiced) to 0.5 (voiced)
         *-----------------------------------------------------------------*/

        if( Jopt_flag == 0 )
        {
            /* 2/3-bit decoding */
            gain_dec_tc( st, st->core_brate, L_frame, i_subfr, tc_subfr, Es_pred, code, &gain_pit, &gain_code, &gain_inov, &norm_gain_code );
        }
        else
        {
            /* 5-bit decoding */
            if ( st->core_brate > ACELP_32k )
            {
                gain_dec_SQ( st, st->core_brate, coder_type, i_subfr, tc_subfr, code, Es_pred, &gain_pit, &gain_code, &gain_inov, &norm_gain_code );
            }
            else
            {
                gain_dec_mless( st, st->core_brate, L_frame, coder_type, i_subfr, tc_subfr, code, Es_pred, &gain_pit, &gain_code, &gain_inov, &norm_gain_code );
            }
        }

        /* update LP filtered gains for the case of frame erasures */
        lp_gain_updt( i_subfr, gain_pit, norm_gain_code + norm_gain_preQ, &st->lp_gainp, &st->lp_gainc, L_frame );

        st->tilt_code = est_tilt( exc+i_subfr, gain_pit, code, gain_code, &voice_fac,L_SUBFR,0 );

        /*----------------------------------------------------------------------*
         * Find the total excitation
         *----------------------------------------------------------------------*/

        for (i = 0; i < L_SUBFR;  i++)
        {
            exc2[i+i_subfr] = gain_pit*exc[i+i_subfr];
            exc[i+i_subfr] = exc2[i+i_subfr] + gain_code*code[i];
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

        /*----------------------------------------------------------------*
         * Excitation enhancements (update of total excitation signal)
         *----------------------------------------------------------------*/
        if( st->core_brate > ACELP_32k )
        {
            mvr2r( exc+i_subfr, exc2+i_subfr, L_SUBFR );
        }
        else
        {
            enhancer( MODE1, st->core_brate, -1, 0, coder_type, L_frame, voice_fac, st->stab_fac,
                      norm_gain_code, gain_inov, &st->gc_threshold, code, exc2 + i_subfr, gain_pit, st->dispMem );
        }

        p_Aq += (M+1);
        pt_pitch++;
        st->tilt_code_dec[i_subfr/L_SUBFR] = st->tilt_code;
        gain_buf[i_subfr/L_SUBFR] = gain_pit;
    }

    /* SC-VBR */
    st->prev_gain_pit_dec = gain_pit;

    return;
}
