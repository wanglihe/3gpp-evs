/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*-------------------------------------------------------------------*
 * encod_amr_wb()
 *
 * Encode excitation signal in AMR-WB IO mode
 *-------------------------------------------------------------------*/

void encod_amr_wb(
    Encoder_State *st,               /* i/o: state structure                                  */
    LPD_state *mem,              /* i/o: encoder state structure                          */
    const float speech[],          /* i  : input speech                                     */
    const float Aw[],              /* i  : weighted A(z) unquantized for subframes          */
    const float Aq[],              /* i  : 12k8 Lp coefficient                              */
    const short pitch[3],          /* i  : open-loop pitch values for quantiz.              */
    const float voicing[],         /* i  : voicing                                          */
    const float *res,              /* i  : residual signal                                  */
    float *syn,              /* i/o: core synthesis                                   */
    float *exc,              /* i/o: current non-enhanced excitation                  */
    float *exc2,             /* i/o: current enhanced excitation                      */
    float *pitch_buf,        /* i/o: floating pitch values for each subframe          */
    short hf_gain[NB_SUBFR], /* o  : decoded HF gain                                  */
    const float *speech16k         /* i  : input speech @16kHz                              */
)
{
    float xn[L_SUBFR];                  /* Target vector for pitch search    */
    float xn2[L_SUBFR];                 /* Target vector for codebook search */
    float cn[L_SUBFR];                  /* Target vector in residual domain  */
    float h1[L_SUBFR+(M+1)];            /* Impulse response vector           */
    float code[L_SUBFR];                /* Fixed codebook excitation         */
    float y1[L_SUBFR];                  /* Filtered adaptive excitation      */
    float y2[L_SUBFR];                  /* Filtered algebraic excitation     */
    float gain_pit ;                    /* Pitch gain                        */
    float voice_fac;                    /* Voicing factor                    */
    float gain_code;                    /* Gain of code                      */
    float gain_inov;                    /* inovation gain                    */
    short i, i_subfr;                   /* tmp variables                     */
    short T_op[3];                      /* pitch period for quantization     */
    short T0, T0_frac;                  /* close loop integer pitch and fractional part */
    short T0_min, T0_max;               /* pitch variables                   */
    float *pt_pitch;                    /* pointer to floating pitch buffer  */
    float g_corr[6];                    /* ACELP correl, values + gain pitch */
    short clip_gain;                    /* LSF clip gain                     */
    const float *p_Aw, *p_Aq;           /* pointer to LP filter coeff. vector*/
    short unbits = 0;
    float norm_gain_code;
    short pitch_limit_flag;
    short lp_select, lp_flag;


    /*------------------------------------------------------------------*
     * Initializations
     *------------------------------------------------------------------*/

    pitch_limit_flag = 0;  /* always restrained pitch Q range in IO mode */
    T0_max = PIT_MAX;
    T0_min = PIT_MIN;

    p_Aw = Aw;
    p_Aq = Aq;
    pt_pitch = pitch_buf;

    mvs2s( pitch, T_op, 2 );
    if( T_op[0] <= PIT_MIN )
    {
        T_op[0] *= 2;
    }

    if( T_op[1] <= PIT_MIN )
    {
        T_op[1] *= 2;
    }

    /*------------------------------------------------------------------*
     * ACELP subframe loop
     *------------------------------------------------------------------*/

    for ( i_subfr = 0; i_subfr < L_FRAME; i_subfr += L_SUBFR )
    {
        /*----------------------------------------------------------------*
         * Bandwidth expansion of A(z) filter coefficients
         * Find the the excitation search target "xn" and innovation
         *   target in residual domain "cn"
         * Compute impulse response, h1[], of weighted synthesis filter
         *----------------------------------------------------------------*/

        mvr2r( &res[i_subfr], &exc[i_subfr], L_SUBFR );

        find_targets( speech, mem->mem_syn, i_subfr, &mem->mem_w0, p_Aq, res, L_SUBFR, p_Aw, TILT_FAC, xn, cn, h1 );

        /*----------------------------------------------------------------*
         * Close-loop pitch search and quantization
         * Adaptive exc. construction
         *----------------------------------------------------------------*/

        *pt_pitch = pit_encode( st, st->core_brate, 1, L_FRAME, -1, &pitch_limit_flag, i_subfr, exc,
                                L_SUBFR, T_op, &T0_min, &T0_max, &T0, &T0_frac, h1, xn );

        /*-----------------------------------------------------------------*
         * Find adaptive exitation
         *-----------------------------------------------------------------*/

        pred_lt4( &exc[i_subfr], &exc[i_subfr], T0, T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP);

        /*-----------------------------------------------------------------*
         * Gain clipping test to avoid unstable synthesis on frame erasure
         *   or in case of floating point encoder & fixed p. decoder
         *-----------------------------------------------------------------*/

        clip_gain = gp_clip( st->core_brate, voicing, i_subfr, 0, xn, st->clip_var );

        /*-----------------------------------------------------------------*
         * LP filtering of the adaptive excitation, codebook target computation
         *-----------------------------------------------------------------*/

        lp_select = lp_filt_exc_enc( MODE1, st->core_brate, 1, -1, i_subfr, exc, h1, xn, y1, xn2,
                                     L_SUBFR, L_FRAME, g_corr, clip_gain, &gain_pit, &lp_flag );

        if( lp_flag == NORMAL_OPERATION )
        {
            push_indice( st, IND_LP_FILT_SELECT, lp_select, 1 );
        }

        /*-----------------------------------------------------------------*
         * Innovation encoding
         *-----------------------------------------------------------------*/

        inov_encode( st, st->core_brate, 1, L_FRAME, st->last_L_frame, -1, -1, 0, i_subfr, -1, p_Aq,
                     gain_pit, cn, exc, h1, mem->tilt_code, *pt_pitch, xn2, code, y2, &unbits);

        /*-----------------------------------------------------------------*
         * Gain encoding
         * Pitch gain clipping test
         * Estimate spectrum tilt and voicing
         *-----------------------------------------------------------------*/

        gain_enc_amr_wb( st, xn, y1, y2, code, st->core_brate, &gain_pit, &gain_code,
                         &gain_inov, &norm_gain_code, g_corr, clip_gain, st->past_qua_en );

        gp_clip_test_gain_pit( st->core_brate, gain_pit, st->clip_var);

        mem->tilt_code = est_tilt( exc+i_subfr, gain_pit, code, gain_code, &voice_fac, L_SUBFR, 0 );

        /*-----------------------------------------------------------------*
         * Update memory of the weighting filter
         *-----------------------------------------------------------------*/

        mem->mem_w0 = xn[L_SUBFR-1] - gain_pit * y1[L_SUBFR-1] - gain_code * y2[L_SUBFR-1];

        /*-----------------------------------------------------------------*
         * Find the total excitation
         *-----------------------------------------------------------------*/

        for ( i = 0; i < L_SUBFR;  i++ )
        {
            exc2[i+i_subfr] = gain_pit * exc[i+i_subfr];
            exc[i+i_subfr] = exc2[i+i_subfr] + gain_code * code[i];
        }


        /*-----------------------------------------------------------------*
         * Synthesize speech to update mem_syn[]
         * Update A(z) filters
         *-----------------------------------------------------------------*/

        syn_filt( p_Aq, M, &exc[i_subfr], &syn[i_subfr], L_SUBFR, mem->mem_syn, 1 );

        /*-----------------------------------------------------------------*
         * HF gain modification factors at 23.85 kbps
         *-----------------------------------------------------------------*/

        if ( st->core_brate == ACELP_23k85 )
        {
            if( st->input_Fs >= 16000 )
            {
                hf_cod( st->core_brate, &speech16k[i_subfr * L_SUBFR16k/L_SUBFR], p_Aq, &exc[i_subfr], &syn[i_subfr],
                        &st->seed2_enc, st->mem_hp400_enc, st->mem_syn_hf_enc, st->mem_hf_enc, st->mem_hf2_enc, &st->hangover_cnt,
                        &st->gain_alpha, &hf_gain[i_subfr/L_SUBFR] );
            }

            push_indice( st, IND_HF_GAIN_MODIFICATION, hf_gain[i_subfr/L_SUBFR], 4 );
        }

        p_Aw += (M+1);
        p_Aq += (M+1);
        pt_pitch++;
    }

    return;
}
