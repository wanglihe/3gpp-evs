/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"

/*-------------------------------------------------------------------
 * encod_ppp()
 *
 * Encode PPP frames in SC-VBR
 *-------------------------------------------------------------------*/

void encod_ppp(
    Encoder_State *st,              /* i/o: state structure                         */
    LPD_state *mem,             /* i/o: encoder memories                        */
    const float speech[],         /* i  : input speech                            */
    const float Aw[],             /* i  : weighted A(z) unquantized for subframes */
    const float Aq[],             /* i  : 12k8 Lp coefficient                     */
    short *coder_type,      /* i/o: coding type                             */
    const short sharpFlag,        /* i  : formant sharpening flag                 */
    const short T_op[],           /* i  : open loop pitch                         */
    const float voicing[],        /* i  : voicing                                 */
    float *res,             /* i/o: residual signal                         */
    float *synth,           /* i/o: core synthesis                          */
    float *exc,             /* i/o: current non-enhanced excitation         */
    float *exc2,            /* i/o: current enhanced excitation             */
    float *pitch_buf,       /* i/o: floating pitch values for each subframe */
    float *voice_factors,   /* o  : voicing factors                         */
    float *bwe_exc          /* o  : excitation for SWB TBE                  */
)
{
    float xn[L_SUBFR];      /* Target vector for pitch search               */
    float h1[L_SUBFR+(M+1)];/* Impulse response vector                      */
    short i_subfr;          /* tmp variables                                */
    const float *p_Aw,*p_Aq;/* pointer to LP filter coeff. vector           */
    short k;
    float p_Aq_old[M+1], excQ_ppp[L_FRAME], p_Aq_curr[M], pitch[NB_SUBFR];
    float LPC_de_old[M+1], LPC_de_curr[M+1];
    short rate_ctrl;

    rate_ctrl = st->rate_control;

    /*------------------------------------------------------------------*
    * ACELP subframe loop
    *------------------------------------------------------------------*/

    p_Aw = Aw;
    p_Aq = Aq;

    for( i_subfr=0; i_subfr<L_FRAME; i_subfr+=L_SUBFR )
    {
        /*----------------------------------------------------------------*
         * Bandwidth expansion of A(z) filter coefficients
         * Find the the excitation search target "xn" and innovation
         * target in residual domain "cn"
         * Compute impulse response, h1[], of weighted synthesis filter
         *----------------------------------------------------------------*/

        mvr2r( &res[i_subfr], &exc[i_subfr], L_SUBFR );

        find_targets( speech, mem->mem_syn, i_subfr, &mem->mem_w0, p_Aq,res, L_SUBFR, p_Aw, TILT_FAC,xn, NULL, h1 );

        /* call voiced encoder at this point */
        if( i_subfr == 0 ) /* generate the L_FRAME exc */
        {
            for( k=0; k<M; k++ )
            {
                p_Aq_curr[k] = p_Aq[k+(3*(M+1))+1];
            }

            lsp2a_stab( st->lsp_old, p_Aq_old, M );

            deemph_lpc( p_Aq_curr, p_Aq_old, LPC_de_curr, LPC_de_old, 1 );

            /* last frame-end lpc and curr frame-end lpc */
            ppp_voiced_encoder( st, res, excQ_ppp, T_op[1], LPC_de_old, LPC_de_curr, exc, pitch, st->vadsnr );

            if( st->bump_up )
            {
                i_subfr = L_FRAME;
            }
        }

        if( st->bump_up != 1 )
        {
            /*-----------------------------------------------------------------*
             * Gain clipping test to avoid unstable synthesis on frame erasure
             * or in case of floating point encoder & fixed p. decoder
             *-----------------------------------------------------------------*/

            gp_clip( st->core_brate, voicing, i_subfr, *coder_type, xn, st->clip_var );

            /* run the above to maintain gain clipping memories */
            gp_clip_test_gain_pit( st->core_brate, st->prev_ppp_gain_pit, st->clip_var );

            /*-----------------------------------------------------------------*
             * Synthesize speech to update mem_syn[].
             * Update A(z) filters
             *-----------------------------------------------------------------*/

            syn_filt( p_Aq, M, &excQ_ppp[i_subfr], &synth[i_subfr], L_SUBFR, mem->mem_syn, 1 );

            p_Aw += (M+1);
            p_Aq += (M+1);
        }

    }   /* end of subframe loop */


    if( st->bump_up )
    {
        /* PPP failed, bump up */
        st->ppp_mode = 0;
        st->core_brate = ACELP_7k20;
        st->pppcountE = 0;

        if ( st->set_ppp_generic )
        {
            *coder_type = GENERIC;
        }
        else
        {
            *coder_type = VOICED;
        }

        /* delete previous indices */
        reset_indices_enc( st );

        /* signalling matrix (writing of signalling bits) */
        signalling_enc( st, *coder_type, sharpFlag );
    }
    else
    {
        mvr2r( excQ_ppp, exc, L_FRAME );

        /*-----------------------------------------------------------------*
         * Updates: last value of new target is stored in mem_w0
         *-----------------------------------------------------------------*/

        mem->mem_w0 = xn[L_SUBFR-1] - (exc[L_FRAME-1]);

        mvr2r( exc, exc2, L_FRAME );

        st->dispMem[0] = 2;
        st->dispMem[2] = st->prev_ppp_gain_pit;

        for( k=3; k<7; k++ )
        {
            st->dispMem[k] = st->dispMem[k-1];
        }
        mem->tilt_code = st->prev_tilt_code;
        mvr2r( pitch, pitch_buf, NB_SUBFR );
        pitch_buf[NB_SUBFR16k-1] = pitch[NB_SUBFR-1];
        interp_code_5over2( exc2, bwe_exc, L_FRAME );
        set_f( voice_factors, 0.0f, NB_SUBFR16k );
    }

    st->rate_control = rate_ctrl;

    return;
}
