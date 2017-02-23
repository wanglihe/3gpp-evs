/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*---------------------------------------------------------------------*
 * decod_gen_voic()
 *
 * Decode generic (GC), voiced (VC) and AMR-WB IO frames
 *---------------------------------------------------------------------*/

void decod_gen_voic(
    Decoder_State *st,                  /* i/o: decoder static memory                     */
    const short L_frame,              /* i  : length of the frame                       */
    const short sharpFlag,            /* i  : formant sharpening flag                   */
    const float *Aq,                  /* i  : LP filter coefficient                     */
    const short coder_type,           /* i  : coding type                               */
    const float Es_pred,              /* i  : predicted scaled innov. energy            */
    const short do_WI,                /* i  : FEC fast recovery flag                    */
    float *pitch_buf,           /* o  : floating pitch values for each subframe   */
    float *voice_factors,       /* o  : voicing factors                           */
    float *exc,                 /* i/o: adapt. excitation exc                     */
    float *exc2,                /* i/o: adapt. excitation/total exc               */
    float *bwe_exc,             /* o  : excitation for SWB TBE                    */
    short *unbits,              /* number of unused bits                          */
    float *gain_buf
)
{
    short T0 = PIT_MIN, T0_frac = 0, T0_min, T0_max;/* integer pitch variables               */
    float gain_pit = 0.0f;          /* pitch gain                                            */
    float gain_code = 0.0f;         /* gain/normalized gain of the algebraic excitation      */
    float norm_gain_code = 0.0f;    /* normalized gain of the algebraic excitation           */
    float gain_inov = 0;            /* Innovation gain                                       */
    float gains_mem[2*(NB_SUBFR-1)];/* pitch gain and code gain from previous subframes      */
    float voice_fac;                /* voicing factor                                        */
    float code[L_SUBFR];            /* algebraic codevector                                  */

    const float *p_Aq;              /* Pointer to frame LP coefficient                       */
    float *pt_pitch;                /* pointer to floating pitch                             */
    short i_subfr, i;               /* tmp variables                                         */
    int   offset;
    float error = 0.0f;
    float gain_preQ = 0;            /* Gain of prequantizer excitation                       */
    float code_preQ[L_SUBFR];       /* Prequantizer excitation                               */
    float norm_gain_preQ;
    short pitch_limit_flag;

    DTFS_STRUCTURE *PREVP, *CURRP;
    short shft_prev = 0, shft_curr = 0;
    float ph_offset, dummy2[2], out[L_FRAME16k], enratio = 0.0f;
    float sp_enratio, curr_spch_nrg, prev_spch_nrg, curr_res_nrg, prev_res_nrg, syn_tmp[L_FRAME16k], mem_tmp[M];
    short harm_flag_acelp;

    /* read harmonicity flag */
    harm_flag_acelp = 0;
    if( st->core_brate > ACELP_24k40 && st->core_brate <= ACELP_32k && coder_type == GENERIC )
    {
        harm_flag_acelp = (short)get_next_indice( st, 1 );
    }

    /*------------------------------------------------------------------*
     * ACELP subframe loop
     *------------------------------------------------------------------*/

    p_Aq = Aq;                  /* pointer to interpolated LPC parameters */
    pt_pitch = pitch_buf;       /* pointer to the pitch buffer */
    norm_gain_preQ = 0.0f;
    gain_preQ = 0;
    set_f( code_preQ, 0, L_SUBFR );

    for( i_subfr = 0; i_subfr < L_frame; i_subfr += L_SUBFR )
    {

        /*----------------------------------------------------------------------*
         * Decode pitch lag
         *----------------------------------------------------------------------*/

        *pt_pitch = pit_decode( st, st->core_brate, 0, L_frame, i_subfr, coder_type, &pitch_limit_flag, &T0, &T0_frac, &T0_min, &T0_max, L_SUBFR );

        /*--------------------------------------------------------------*
         * Find the adaptive codebook vector
         *--------------------------------------------------------------*/

        pred_lt4( &exc[i_subfr], &exc[i_subfr], T0, T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP );

        if( L_frame == L_FRAME )
        {
            offset = tbe_celp_exc_offset(T0, T0_frac);
            for (i=0; i<L_SUBFR * HIBND_ACB_L_FAC; i++)
            {
                bwe_exc[i + i_subfr * HIBND_ACB_L_FAC] = bwe_exc[i + i_subfr * HIBND_ACB_L_FAC - offset + (int) error];
            }
            error += (float) offset - (float) T0 * HIBND_ACB_L_FAC - 0.25f * HIBND_ACB_L_FAC * (float) T0_frac;
        }
        else
        {
            offset = T0 * 2 + (int) ((float) T0_frac * 0.5f + 4 + 0.5f) - 4;
            for (i=0; i<L_SUBFR * 2; i++)
            {
                bwe_exc[i + i_subfr * 2] = bwe_exc[i + i_subfr * 2 - offset + (int) error];
            }
            error += (float) offset - (float) T0 * 2 - 0.5f * (float) T0_frac;
        }

        /*--------------------------------------------------------------*
         * LP filtering of the adaptive excitation
         *--------------------------------------------------------------*/

        lp_filt_exc_dec( st, MODE1, st->core_brate, 0, coder_type, i_subfr, L_SUBFR, L_frame, 0, exc );

        /*-----------------------------------------------------------------*
         * Transform-domain contribution decoding (active frames)
         *-----------------------------------------------------------------*/

        if( st->core_brate > ACELP_24k40 && coder_type != INACTIVE )
        {
            gain_code = 0.0f;
            transf_cdbk_dec( st, st->core_brate, coder_type, harm_flag_acelp, i_subfr, -1, Es_pred,
                             gain_code, &st->mem_preemp_preQ, &gain_preQ, &norm_gain_preQ, code_preQ, unbits );
        }

        /*--------------------------------------------------------------*
         * Innovation decoding
         *--------------------------------------------------------------*/

        inov_decode( st, st->core_brate, 0, L_frame, coder_type, sharpFlag, i_subfr, -1, p_Aq, st->tilt_code,*pt_pitch, code );

        /*--------------------------------------------------------------*
         * Gain decoding
         * Estimate spectrum tilt and voicing
         *--------------------------------------------------------------*/

        if ( st->core_brate <= ACELP_8k00 )
        {
            gain_dec_lbr( st, st->core_brate, coder_type, i_subfr, code, &gain_pit, &gain_code, &gain_inov, &norm_gain_code, gains_mem );
        }
        else if ( st->core_brate > ACELP_32k )
        {
            gain_dec_SQ( st, st->core_brate, coder_type, i_subfr, -1, code, Es_pred, &gain_pit, &gain_code, &gain_inov, &norm_gain_code );
        }
        else
        {
            gain_dec_mless( st, st->core_brate, L_frame, coder_type, i_subfr, -1, code, Es_pred, &gain_pit, &gain_code, &gain_inov, &norm_gain_code );
        }

        st->tilt_code = est_tilt( exc+i_subfr, gain_pit, code, gain_code, &voice_fac,L_SUBFR,0 );

        /*-----------------------------------------------------------------*
         * Transform-domain contribution decoding (inactive frames)
         *-----------------------------------------------------------------*/

        if( st->core_brate > ACELP_24k40 && coder_type == INACTIVE )
        {
            transf_cdbk_dec( st, st->core_brate, coder_type, 0, i_subfr, -1, Es_pred, gain_code,
                             &st->mem_preemp_preQ, &gain_preQ, &norm_gain_preQ, code_preQ, unbits );
        }

        /* update LP filtered gains for the case of frame erasures */
        lp_gain_updt( i_subfr, gain_pit, norm_gain_code + norm_gain_preQ, &st->lp_gainp, &st->lp_gainc, L_frame );

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

        prep_tbe_exc( L_frame, i_subfr, gain_pit, gain_code, code, voice_fac, &voice_factors[i_subfr/L_SUBFR], bwe_exc,
                      gain_preQ, code_preQ, T0, coder_type, st->core_brate );

        /*----------------------------------------------------------------*
         * Excitation enhancements (update of total excitation signal)
         *----------------------------------------------------------------*/
        if( st->core_brate > ACELP_32k || coder_type == INACTIVE )
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
        gain_buf[i_subfr/L_SUBFR] = gain_pit;
        st->tilt_code_dec[i_subfr/L_SUBFR] = st->tilt_code;
    }

    /* FEC fast recovery */
    if ( do_WI )
    {
        shft_prev = L_EXC_MEM - (short) rint_new(st->bfi_pitch);
        prev_res_nrg = sum2_f( st->old_exc2+shft_prev, (short) rint_new(st->bfi_pitch) ) + 1e-6f;
        prev_spch_nrg = sum2_f( st->old_syn2+shft_prev, (short) rint_new(st->bfi_pitch) ) + 1e-6f;

        mvr2r( st->mem_syn2, mem_tmp, M );
        syn_12k8( st->L_frame, Aq, exc2, syn_tmp, mem_tmp, 1 );

        shft_curr = st->L_frame - (short) rint_new(pitch_buf[NB_SUBFR16k-1]);
        curr_res_nrg = sum2_f(exc2+shft_curr, (short) rint_new(pitch_buf[NB_SUBFR16k-1]));
        curr_spch_nrg = sum2_f(syn_tmp+shft_curr,(short) rint_new(pitch_buf[NB_SUBFR16k-1]));

        enratio = curr_res_nrg/prev_res_nrg;
        sp_enratio = curr_spch_nrg/prev_spch_nrg;

        if ( enratio > 0.25f &&
                enratio < 15.0f &&
                sp_enratio > 0.15f &&
                st->bfi_pitch < 150 &&
                pitch_buf[NB_SUBFR16k-1] < 150 )
        {
            PREVP = DTFS_new();
            CURRP = DTFS_new();

            DTFS_to_fs( st->old_exc2+shft_prev, (short)rint_new( st->bfi_pitch ), PREVP, (short)st->output_Fs, do_WI );
            DTFS_to_fs( exc2+shft_curr, (short)rint_new( pitch_buf[NB_SUBFR16k-1] ), CURRP, (short)st->output_Fs, do_WI );

            ph_offset = 0.0f;
            WIsyn( *PREVP, CURRP, dummy2, &ph_offset, out, st->L_frame, 1 );

            mvr2r( out, exc2, st->L_frame);
            mvr2r( exc2, exc, st->L_frame);

            /* update bwe_exc for SWB-TBE */
            for (i_subfr = 0; i_subfr < L_frame; i_subfr += L_SUBFR)
            {
                interp_code_4over2( exc + i_subfr, bwe_exc + (i_subfr*2), L_SUBFR );
            }

            free(PREVP);
            free(CURRP);
        }
    }

    /* SC-VBR */
    st->prev_gain_pit_dec = gain_pit;

    return;
}
