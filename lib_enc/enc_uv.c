/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * encod_unvoiced()
 *
 * Encode unvoiced (UC) frames
 *-------------------------------------------------------------------*/

void encod_unvoiced(
    Encoder_State *st,                /* i/o: state structure                         */
    LPD_state *mem,               /* i/o: encoder memories                        */
    const float *speech,            /* i  : input speech                            */
    const float Aw[],               /* i  : weighted A(z) unquantized for subframes */
    const float *Aq,                /* i  : LP coefficients                         */
    const short vad_flag,
    const float *res,               /* i  : residual signal                         */
    float *syn,               /* o  : core synthesis                          */
    float *tmp_noise,         /* o  : long-term noise energy                  */
    float *exc,               /* i/o: current non-enhanced excitation         */
    float *pitch_buf,         /* o  : floating pitch values for each subframe */
    float *voice_factors,     /* o  : voicing factors                         */
    float *bwe_exc            /* i/o: excitation for SWB TBE                  */
)
{
    float xn[L_SUBFR];              /* Target vector for pitch search               */
    float h1[L_SUBFR];              /* Impulse response vector                      */
    float code[L_SUBFR];            /* Fixed codebook excitation                    */
    float y2[L_SUBFR];              /* Filtered algebraic excitation                */
    float *pt_pitch;                /* pointer to floating pitch buffer             */
    float gain_pit;                 /* Pitch gain                                   */
    float voice_fac;                /* Voicing factor                               */
    float gain_code;                /* gain of code                                 */
    float gain_inov;                /* inovative gain                               */
    const float *p_Aw, *p_Aq;       /* pointer to LP filter coeff. vector           */
    short i_subfr;
    float norm_gain_code;

    /*------------------------------------------------------------------*
     * Initializations
     *------------------------------------------------------------------*/

    gain_pit = 0;

    if( st->Opt_SC_VBR && vad_flag == 0 && (st->last_ppp_mode == 1 || st->last_nelp_mode == 1) )
    {
        /* SC_VBR - reset the encoder, to avoid memory not updated issue for the
           case when UNVOICED mode is used to code inactive speech */
        CNG_reset_enc( st, mem, pitch_buf, voice_factors, 1 );

    }

    p_Aw = Aw;
    p_Aq = Aq;
    pt_pitch = pitch_buf;

    /*----------------------------------------------------------------*
     * subframe loop
     *----------------------------------------------------------------*/

    for( i_subfr=0; i_subfr<L_FRAME; i_subfr+=L_SUBFR )
    {
        /*----------------------------------------------------------------*
         * Bandwidth expansion of A(z) filter coefficients
         * Find the excitation search target "xn" and innovation target in residual domain "cn"
         * Compute impulse response, h1[], of weighted synthesis filter
         *----------------------------------------------------------------*/

        mvr2r( &res[i_subfr], &exc[i_subfr], L_SUBFR );

        find_targets( speech, mem->mem_syn, i_subfr, &mem->mem_w0, p_Aq,res, L_SUBFR, p_Aw, st->preemph_fac, xn, NULL, h1 );

        /*----------------------------------------------------------------*
         * Unvoiced subframe processing
         *----------------------------------------------------------------*/

        *pt_pitch = gaus_encode( st, i_subfr, h1, xn, exc, &mem->mem_w0, st->clip_var, &mem->tilt_code, code, &gain_code,
                                 y2, &gain_inov, &voice_fac, &gain_pit, &norm_gain_code, st->core_brate );

        *tmp_noise = norm_gain_code;

        voice_factors[i_subfr/L_SUBFR] =  0.0f;

        interp_code_5over2( &exc[i_subfr], &bwe_exc[i_subfr * HIBND_ACB_L_FAC], L_SUBFR );

        /*-----------------------------------------------------------------*
         * Synthesize speech to update mem_syn[].
         * Update A(z) filters
         *-----------------------------------------------------------------*/

        syn_filt( p_Aq, M, &exc[i_subfr], &syn[i_subfr], L_SUBFR, mem->mem_syn, 1 );

        p_Aw += (M+1);
        p_Aq += (M+1);
        pt_pitch++;
    }

    /* SC-VBR */
    st->prev_ppp_gain_pit = gain_pit;
    st->prev_tilt_code = mem->tilt_code;

    return;
}
