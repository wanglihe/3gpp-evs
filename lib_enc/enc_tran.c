/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "prot.h"
#include "rom_com.h"

/*-------------------------------------------------------------------*
 * encod_tran()
 *
 * Encode transition (TC) frames
 *-------------------------------------------------------------------*/

short encod_tran(
    Encoder_State *st,                  /* i/o: state structure                                   */
    LPD_state *mem,                 /* i/o: encoder memories                                  */
    const short L_frame,              /* i  : length of the frame                               */
    const float speech[],             /* i  : input speech                                      */
    const float Aw[],                 /* i  : weighted A(z) unquantized for subframes           */
    const float Aq[],                 /* i  : LP coefficients                                   */
    const short coder_type,           /* i  : coding type                                       */
    const float Es_pred,              /* i  : predicted scaled innov. energy                    */
    const short T_op[],               /* i  : open loop pitch                                   */
    const float voicing[],            /* i  : voicing                                           */
    const float *res,                 /* i  : residual signal                                   */
    float *syn,                 /* i/o: core synthesis                                    */
    float *exc,                 /* i/o: current non-enhanced excitation                   */
    float *exc2,                /* i/o: current enhanced excitation                       */
    float *pitch_buf,           /* i/o: floating pitch values for each subframe           */
    float *voice_factors,       /* o  : voicing factors                                   */
    float *bwe_exc,             /* i/o: excitation for SWB TBE                            */
    const short attack_flag,          /* i  : Flag to indicate when an audio attack is dealt with TM */
    short *unbits,              /* i/o: number of unused bits                             */
    const short sharpFlag             /* i  : formant sharpening flag                           */
)
{
    float xn[L_SUBFR];                /* Target vector for pitch search                         */
    float xn2[L_SUBFR];               /* Target vector for codebook search                      */
    float cn[L_SUBFR];                /* Target vector in residual domain                       */
    float h1[L_SUBFR+(M+1)];          /* Impulse response vector                                */
    float code[L_SUBFR];              /* Fixed codebook excitation                              */
    float y1[L_SUBFR];                /* Filtered adaptive excitation                           */
    float y2[L_SUBFR];                /* Filtered algebraic excitation                          */
    float gain_pit;                   /* Pitch gain                                             */
    float voice_fac;                  /* Voicing factor                                         */
    float gain_code;                  /* Gain of code                                           */
    float gain_inov;                  /* inovation gain                                         */
    short i, i_subfr, tc_subfr;       /* tmp variables                                          */
    short position, T0_min, T0_max;   /* pitch and TC variables                                 */
    short T0, T0_frac;                /* close loop integer pitch and fractional part           */
    float *pt_pitch;                  /* pointer to floating pitch buffer                       */
    float g_corr[6];                  /* ACELP correlation values  and gain pitch               */
    short clip_gain;                  /* LSF clip gain                                          */
    const float *p_Aw, *p_Aq;         /* pointer to LP filter coefficient vector                */
    float gain_preQ;                  /* Gain of prequantizer excitation                        */
    float code_preQ[L_SUBFR];         /* Prequantizer excitation                                */
    short Jopt_flag;                  /* joint optimization flag                                */
    short unbits_PI;                  /* saved bits for EVS_PI                                  */
    float norm_gain_code;

    /*------------------------------------------------------------------*
     * Initializations
     *------------------------------------------------------------------*/

    gain_pit = 0;
    gain_code = 0;
    gain_preQ = 0;
    unbits_PI = 0;

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

    Jopt_flag = 0;
    tc_subfr = -1;
    if( attack_flag )
    {
        tc_subfr = 3*L_SUBFR;
    }

    p_Aw = Aw;
    p_Aq = Aq;
    pt_pitch = pitch_buf;
    gain_preQ = 0;
    set_f( code_preQ, 0, L_SUBFR );

    /*----------------------------------------------------------------*
     * ACELP subframe loop
     *----------------------------------------------------------------*/

    for ( i_subfr=0; i_subfr<L_frame; i_subfr+=L_SUBFR )
    {
        /*----------------------------------------------------------------*
         * Find the the excitation search target "xn" and innovation
         *   target in residual domain "cn"
         * Compute impulse response, h1[], of weighted synthesis filter
         *----------------------------------------------------------------*/

        mvr2r( &res[i_subfr], &exc[i_subfr], L_SUBFR );

        find_targets( speech, mem->mem_syn, i_subfr, &mem->mem_w0, p_Aq, res, L_SUBFR, p_Aw, st->preemph_fac, xn, cn, h1 );

        /*-----------------------------------------------------------------*
         * TC: subframe determination &
         * adaptive/glottal part of excitation construction
         *-----------------------------------------------------------------*/

        transition_enc( st, st->core_brate, L_frame, coder_type, i_subfr, &tc_subfr, &Jopt_flag, &position, voicing, T_op, &T0,
                        &T0_frac, &T0_min, &T0_max, exc, y1, res, h1, xn, xn2, st->clip_var, &gain_pit, g_corr, &clip_gain, &pt_pitch, bwe_exc );

        /*-----------------------------------------------------------------*
         * Transform domain contribution encoding - active frames
         *-----------------------------------------------------------------*/

        if( st->core_brate > ACELP_24k40 )
        {
            transf_cdbk_enc( st, st->core_brate, st->extl, coder_type, 0, i_subfr, tc_subfr, cn, exc,
                             p_Aq, p_Aw, h1, xn, xn2, y1, y2, Es_pred, &gain_pit, gain_code, g_corr, clip_gain,
                             &(st->mem_deemp_preQ), &(st->mem_preemp_preQ), &gain_preQ, code_preQ, unbits );
        }

        /*-----------------------------------------------------------------*
         * ACELP codebook search + pitch sharpening
         *-----------------------------------------------------------------*/

        inov_encode( st, st->core_brate, 0, L_frame, st->last_L_frame, coder_type, st->bwidth, sharpFlag, i_subfr, tc_subfr, p_Aq,
                     gain_pit, cn, exc, h1, mem->tilt_code, *pt_pitch, xn2, code, y2, &unbits_PI );

        if( (st->L_frame == L_FRAME16k) && (tc_subfr == 0) && (i_subfr == L_SUBFR) && (T0 == 2*L_SUBFR) )
        {
            Jopt_flag = 1;
        }

        /*-----------------------------------------------------------------*
         * Quantize the gains
         * Test quantized gain of pitch for pitch clipping algorithm
         * Update tilt of code: 0.0 (unvoiced) to 0.5 (voiced)
         *-----------------------------------------------------------------*/

        if( Jopt_flag == 0 )
        {
            /* SQ gain_code */
            gain_enc_tc( st, st->core_brate, L_frame, i_subfr, tc_subfr, xn, y2, code, Es_pred,
                         &gain_pit, &gain_code, &gain_inov, &norm_gain_code );
        }
        else
        {
            if( st->core_brate > ACELP_32k )
            {
                /* SQ gain_pit and gain_code */
                gain_enc_SQ( st, st->core_brate, coder_type, i_subfr, tc_subfr, xn, y1, y2, code, Es_pred,
                             &gain_pit, &gain_code, &gain_inov, &norm_gain_code, g_corr, clip_gain );
            }
            else
            {
                /* VQ gain_pit and gain_code */
                gain_enc_mless( st, st->core_brate, L_frame, coder_type, i_subfr, tc_subfr, xn, y1, y2, code, Es_pred,
                                &gain_pit, &gain_code, &gain_inov, &norm_gain_code, g_corr, clip_gain );
            }
        }

        /*-----------------------------------------------------------------*
         * update LP-filtered gains for the case of frame erasures
         *-----------------------------------------------------------------*/

        gp_clip_test_gain_pit( st->core_brate, gain_pit, st->clip_var);
        mem->tilt_code = est_tilt( exc+i_subfr, gain_pit, code, gain_code, &voice_fac, L_SUBFR, 0 );

        /*-----------------------------------------------------------------*
         * Update memory of the weighting filter
         *-----------------------------------------------------------------*/

        mem->mem_w0 = xn[L_SUBFR-1] - (gain_pit*y1[L_SUBFR-1]) - (gain_code*y2[L_SUBFR-1]);

        /*-----------------------------------------------------------------*
         * Construct adaptive part of the excitation
         * Save the non-enhanced excitation for FEC_exc
         *-----------------------------------------------------------------*/

        for( i = 0; i < L_SUBFR;  i++ )
        {
            exc2[i+i_subfr] = gain_pit * exc[i+i_subfr];
            exc[i+i_subfr] = exc2[i+i_subfr] + gain_code * code[i];
        }

        /*-----------------------------------------------------------------*
         * Add the ACELP pre-quantizer contribution
         *-----------------------------------------------------------------*/

        if( st->core_brate > ACELP_24k40 )
        {
            for( i = 0; i < L_SUBFR; i++ )
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
        i = min( unbits_PI, 16 );
        push_indice( st, IND_UNUSED, 0, i );
        unbits_PI -= i;
    }

    /* write TC configuration */
    if( L_frame == L_FRAME )
    {
        if( tc_subfr == TC_0_0 )
        {
            push_indice( st, IND_TC_SUBFR, 1, 1 );
        }
        else if( tc_subfr == TC_0_64 )
        {
            push_indice( st, IND_TC_SUBFR, 0, 1 );
            push_indice( st, IND_TC_SUBFR, 1, 1 );
            push_indice( st, IND_TC_SUBFR, 0, 1 );
            push_indice( st, IND_TC_SUBFR, 1, 1 );
        }
        else if( tc_subfr == TC_0_128 )
        {
            push_indice( st, IND_TC_SUBFR, 0, 1 );
            push_indice( st, IND_TC_SUBFR, 1, 1 );
            push_indice( st, IND_TC_SUBFR, 0, 1 );
            push_indice( st, IND_TC_SUBFR, 0, 1 );
        }
        else if( tc_subfr == TC_0_192 )
        {
            push_indice( st, IND_TC_SUBFR, 0, 1 );
            push_indice( st, IND_TC_SUBFR, 1, 1 );
            push_indice( st, IND_TC_SUBFR, 1, 1 );
        }
        else if( tc_subfr == L_SUBFR )
        {
            push_indice( st, IND_TC_SUBFR, 0, 1 );
            push_indice( st, IND_TC_SUBFR, 0, 1 );
            push_indice( st, IND_TC_SUBFR, 1, 1 );
        }
        else if( tc_subfr == 2*L_SUBFR )
        {
            push_indice( st, IND_TC_SUBFR, 0, 1 );
            push_indice( st, IND_TC_SUBFR, 0, 1 );
            push_indice( st, IND_TC_SUBFR, 0, 1 );
            push_indice( st, IND_TC_SUBFR, 1, 1 );
        }
        else if( tc_subfr == 3*L_SUBFR )
        {
            push_indice( st, IND_TC_SUBFR, 0, 1 );
            push_indice( st, IND_TC_SUBFR, 0, 1 );
            push_indice( st, IND_TC_SUBFR, 0, 1 );
            push_indice( st, IND_TC_SUBFR, 0, 1 );
        }
    }
    else  /* L_frame == L_FRAME16k */
    {
        if( tc_subfr == 0 )
        {
            push_indice( st, IND_TC_SUBFR, 0, 2 );
        }
        else if( tc_subfr == L_SUBFR )
        {
            push_indice( st, IND_TC_SUBFR, 1, 2 );
        }
        else if( tc_subfr == 2*L_SUBFR )
        {
            push_indice( st, IND_TC_SUBFR, 2, 2 );
        }
        else if( tc_subfr == 3*L_SUBFR )
        {
            push_indice( st, IND_TC_SUBFR, 3, 2 );
            push_indice( st, IND_TC_SUBFR, 0, 1 );
        }
        else if( tc_subfr == 4*L_SUBFR )
        {
            push_indice( st, IND_TC_SUBFR, 3, 2 );
            push_indice( st, IND_TC_SUBFR, 1, 1 );
        }
    }

    /* SC-VBR */
    st->prev_ppp_gain_pit = gain_pit;
    st->prev_tilt_code = mem->tilt_code;

    return tc_subfr;
}
