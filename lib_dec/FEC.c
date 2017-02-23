/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "rom_dec.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * Local functions
 *-------------------------------------------------------------------*/

static void gain_dec_bfi( float *past_qua_en );
static void pulseRes_preCalc( Word16* cond1, Word16* cond2, Word32* cond3 ,Word16 new_pit, Word16 Tc, Word16 L_frame );


/*-------------------------------------------------------------------*
 * FEC_exc_estim()
 *
 * Calculation of excitation signal
 *-------------------------------------------------------------------*/

void FEC_exc_estim(
    Decoder_State *st,               /* i/o: Decoder static memory                      */
    const short L_frame,           /* i  : length of the frame                        */
    float *exc,              /* o  : pointer to excitation buffer (with past)   */
    float *exc2,             /* o  : total excitation (for synthesis)           */
    float *exc_dct_in,       /* o  : GSC excitation in DCT domain               */
    float *pitch_buf,        /* o  : pitch buffer for each subframe             */
    float *voice_factors,    /* o  : voicing factors                            */
    float *tmp_tc,           /* o  : FEC pitch                                  */
    float *bwe_exc,          /* o  : excitation for SWB TBE                     */
    float *lsf_new,          /* i  : ISFs at the end of the frame               */
    float *tmp_noise         /* o  : long-term noise energy                     */
)
{
    short i, Tc, new_pit;
    float exc2_buf[L_FRAME16k + MODE1_L_FIR_FER-1];
    float *pt_exc, *pt1_exc, alpha, step, gain_inov, gain, gainCNG, hp_filt[5];
    float fT0, delta, ftmp;
    short Diff_len;
    short Len, max_len;
    short last_bin;
    short extrapolationFailed = 1;
    float predPitchLag;
    Word16 cond1, cond2;
    Word32 cond3;

    /*-----------------------------------------------------------------*
     * Initializations
     *-----------------------------------------------------------------*/

    gainCNG = (float)sqrt( st->lp_ener );

    ftmp = 2.0f * st->lp_gainc;
    if ( gainCNG > ftmp )
    {
        gainCNG = ftmp;
    }

    Diff_len = 0;
    set_f( exc_dct_in, 0.0f, L_FRAME16k );

    /*-----------------------------------------------------------------*
     * Pitch extrapolation
     *-----------------------------------------------------------------*/

    /* pitch extrapolation */
    pitch_pred_linear_fit( st->nbLostCmpt, st->last_good, st->old_pitch_buf,
                           L_frame == L_FRAME ? &st->old_pitch_buf[2*NB_SUBFR-1] : &st->old_pitch_buf[2*NB_SUBFR16k-1],
                           &predPitchLag, L_frame == L_FRAME?PIT_MIN_DOUBLEEXTEND:PIT16k_MIN_EXTEND,
                           L_frame == L_FRAME?PIT_MAX:PIT16k_MAX, st->mem_pitch_gain, 0, 0, &extrapolationFailed, L_frame/L_SUBFR );

    new_pit = (int)(predPitchLag+0.5f);

    /* initialize FEC pitch to the long-term pitch */
    *tmp_tc = st->bfi_pitch;
    if( L_frame == L_FRAME )
    {
        if( ( ( st->old_pitch_buf[2*NB_SUBFR-1] < 1.8f * st->bfi_pitch ) &&
                ( st->old_pitch_buf[2*NB_SUBFR-1] > 0.6f * st->bfi_pitch ) ) ||     /* last pitch coherent with the past  */
                ( st->upd_cnt >= MAX_UPD_CNT ) )  /* or last update too far in the past */
        {
            /* take the pitch value of last subframe of the previous frame */
            *tmp_tc = st->old_pitch_buf[2*NB_SUBFR-1];
        }
    }
    else  /* L_frame == L_FRAME16k */
    {
        if( ( ( st->old_pitch_buf[2*NB_SUBFR16k-1] < 1.8f * st->bfi_pitch ) &&
                ( st->old_pitch_buf[2*NB_SUBFR16k-1] > 0.6f * st->bfi_pitch ) ) ||     /* last pitch coherent with the past  */
                ( st->upd_cnt >= MAX_UPD_CNT ) )  /* or last update too far in the past */
        {
            /* take the pitch value of last subframe of the previous frame */
            *tmp_tc = st->old_pitch_buf[2*NB_SUBFR16k-1];
        }
    }
    /* convert pitch period */
    Tc = (short)(*tmp_tc + 0.5f);

    pulseRes_preCalc( &cond1, &cond2, &cond3 , (Word16)new_pit, (Word16)Tc, (Word16)L_frame );

    if( (cond1 < 0 ) && (new_pit > 0) && (cond2 != 0) && (cond3 > 0) && extrapolationFailed == 0 )
    {
        fT0 = *tmp_tc;
        delta = (float)(new_pit - fT0)/(L_frame/L_SUBFR);   /* # of subframes */
        for ( i=0; i<L_frame/L_SUBFR; i++ )                 /* subframe pitch values */
        {
            fT0 += delta;
            pitch_buf[i] = (short)(fT0 + 0.5f);
        }
    }
    else
    {
        for ( i=0; i<L_frame/L_SUBFR; i++ )
        {
            pitch_buf[i] = *tmp_tc;
        }
    }

    /*-----------------------------------------------------------------*
     * Estimate gain damping factor
     *-----------------------------------------------------------------*/

    alpha = ALPHA_VT; /* rapid convergence to 0 */
    if( st->last_coder_type == UNVOICED && st->nbLostCmpt <= 3 )
    {
        /* last good frame was clearly unvoiced */
        alpha = ALPHA_UU;
    }
    else if( st->last_coder_type == AUDIO || st->last_good == INACTIVE_CLAS )
    {
        if( st->Last_GSC_pit_band_idx > 0 && st->nbLostCmpt > 1 )
        {
            alpha = 0.8f;
        }
        else if( st->nbLostCmpt <= 5 )
        {
            alpha = 0.995f;
        }
        else
        {
            alpha = 0.95f;
        }
    }
    else if(st->last_good == UNVOICED_CLAS )
    {
        if( st->nbLostCmpt <= 1 )
        {
            /* If stable, do not decrease the energy, pitch gain = 0 */
            alpha = st->stab_fac * (1.0f - 2.0f*ALPHA_U) + 2.0f*ALPHA_U;  /* [0.8, 1.0] */
        }
        else if ( st->nbLostCmpt == 2 )
        {
            alpha = ALPHA_U * 1.5f;   /* 0.6 */
        }
        else
        {
            alpha = ALPHA_U;   /* 0.4 go rapidly to CNG gain, pitch gain = 0 */
        }
    }
    else if(st->last_good == UNVOICED_TRANSITION )
    {
        alpha = ALPHA_UT;
    }
    else if(st->last_good == ONSET && st->nbLostCmpt <= 3 && (st->last_coder_type == GENERIC || st->last_coder_type == TRANSITION) )
    {
        alpha = 0.8f;         /* mild convergence to 0 for the first 3 erased frames */
    }
    else if( (st->last_good == VOICED_CLAS || st->last_good == ONSET ) && st->nbLostCmpt <= 3 )
    {
        alpha = ALPHA_V;      /* constant for the first 3 erased frames */
    }
    else if(st->last_good == SIN_ONSET )
    {
        alpha = ALPHA_S;
    }

    if(st->last_good >= VOICED_CLAS && st->last_good < INACTIVE_CLAS && st->last_coder_type != AUDIO )
    {
        if( st->nbLostCmpt == 1 )
        {
            /* if this is the first erased frame, move pitch gain towards 1 for voiced to remove energy fluctuations */
            gain = (float)sqrt(st->lp_gainp);
            if( gain > 0.98f )
            {
                gain = 0.98f;
            }
            else if( gain < 0.85f )
            {
                gain = 0.85f;
            }

            alpha *= gain;
        }
        else
        {
            alpha = st->lp_gainp;
        }
    }

    /*-----------------------------------------------------------------*
     * Extrapolate past excitation signal using estimated pitch period
     *-----------------------------------------------------------------*/

    if( (st->last_good >= UNVOICED_TRANSITION && st->last_good < INACTIVE_CLAS) ||
            ( (st->last_coder_type == AUDIO || st->last_good == INACTIVE_CLAS) && st->Last_GSC_pit_band_idx > 0) )
    {
        /* Do pitch-based extrapolation only for stable voiced and audio signals */

        pt_exc = exc;
        pt1_exc = pt_exc - Tc;

        if ( st->nbLostCmpt == 1 )
        {
            /* first pitch cycle is low-pass filtered */
            for( i=0 ; i< Tc; i++ )
            {
                *pt_exc++ = ( 0.18f * pt1_exc[-1] + 0.64f * pt1_exc[0] + 0.18f * pt1_exc[1] );
                pt1_exc++;
            }
        }

        /* last pitch cycle of the previous frame is repeatedly copied up to an extra subframe */
        while( pt_exc < exc + L_frame + L_SUBFR )
        {
            *pt_exc++ = *pt1_exc++;
        }

        /* resynchronization of excitation based on glottal pulse locations in the last good frame */
        if ( new_pit > 0 )
        {
            if ((cond1 < 0 ) && (new_pit > 0) && (cond2 != 0) && (cond3 > 0) && extrapolationFailed == 0 )
            {
                mvr2r( exc, exc-L_frame-L_SUBFR, L_frame+L_SUBFR );
                PulseResynchronization( exc-L_frame-L_SUBFR, exc, L_frame, L_frame/L_SUBFR, Tc, new_pit );
            }
        }

        if(st->last_good == UNVOICED_TRANSITION && ( st->last_coder_type == GENERIC || st->last_coder_type == TRANSITION ) )
        {
            /* start of the frame gain */
            gain = 0.0f;

            /* end of the frame gain */
            st->lp_gainp = 0.0f;

            step = 0.0f;
        }
        else
        {
            /* start of the frame gain */
            gain = 1.0f;

            /* end of the frame gain */
            st->lp_gainp = alpha;

            /* linearly attenuate the gain throughout the frame */
            step = (1.0f/L_frame) * (gain - st->lp_gainp);
        }

        /* scaling of the harmonic part of excitation */
        for( i=0 ; i< L_frame; i++ )
        {
            exc[i] *= gain;
            gain -= step;
        }

        if( (st->last_coder_type == AUDIO || st->last_good == INACTIVE_CLAS) && st->Last_GSC_pit_band_idx > 0 )
        {
            /* in case the last frame was coded by GSC, convert the excitation signal to frequency domain */
            edct( exc, exc_dct_in, st->L_frame);

            /* Reset unvaluable part of the adaptive (pitch) excitation contribution */
            Diff_len = (short)(mfreq_loc[st->Last_GSC_pit_band_idx]/BIN_SIZE);
            max_len = st->L_frame - Diff_len;
            Len = min( max_len, 80 );

            for( i=0; i<Len; i++ )
            {
                exc_dct_in[i + Diff_len] *= sm_table[i];
            }

            for( ; i<max_len; i++ )
            {
                exc_dct_in[i + Diff_len] = 0.0f;
            }

            Diff_len++;
        }
    }

    /*-----------------------------------------------------------------*
     * Replicate the last spectrum in case the last good frame was coded by GSC
     *-----------------------------------------------------------------*/

    if( (st->last_coder_type == AUDIO || st->last_good == INACTIVE_CLAS) && st->total_brate <= ACELP_24k40 && !st->Opt_AMR_WB )
    {
        /* Replication of the last spectrum, with a slight downscaling of its dynamic */
        st->GSC_noisy_speech = st->Last_GSC_noisy_speech_flag;
        gsc_dec( st, exc_dct_in, st->Last_GSC_pit_band_idx, Diff_len, 0, st->L_frame/L_SUBFR, st->last_coder_type, &last_bin, lsf_new, NULL, tmp_noise );

        /* Transform back to time domain */
        edct( exc_dct_in, exc, st->L_frame);
    }

    /*-----------------------------------------------------------------*
     * Construct the random part of excitation
     *-----------------------------------------------------------------*/

    else
    {
        /* generate the random part of the excitation */
        for (i=0; i<L_frame+MODE1_L_FIR_FER-1; i++)
        {
            exc2_buf[i] = (float)own_random( &st->seed );
        }

        /* start of the frame gain */
        gain = st->lp_gainc;

        /* end of the frame gain */
        st->lp_gainc = alpha * st->lp_gainc + (1.0f - alpha) * gainCNG;

        if( st->last_good == UNVOICED_TRANSITION && ( st->last_coder_type == GENERIC || st->last_coder_type == TRANSITION ) && gainCNG > 0 )
        {
            st->lp_gainc = gainCNG;
        }

        /* linearly attenuate the gain throughout the frame */
        step = (1.0f/L_frame) * (gain - st->lp_gainc);

        /* calculate gain to normalize energy */
        pt_exc = exc2_buf + MODE1_L_FIR_FER/2;

        gain_inov = 1.0f / (float)sqrt( dotp( pt_exc, pt_exc, L_frame ) / L_frame + 0.01f );

        /* attenuate somewhat on unstable unvoiced */
        if( (st->last_good == UNVOICED_CLAS || st->last_good == INACTIVE_CLAS) && st->last_coder_type != UNVOICED )
        {
            gain_inov *= 0.8f;
        }

        /* scaling of the random part of excitation */
        pt_exc = exc2_buf;
        for( i=0; i< MODE1_L_FIR_FER/2; i++ )
        {
            /* non-causal ringing of the FIR filter */
            *pt_exc++ *= (gain_inov * gain);
        }

        for( i=0; i< L_frame; i++ )
        {
            /* the inner part of the FIR filter */
            *pt_exc++ *= (gain_inov * gain);
            gain -= step;
        }

        for( i=0; i< MODE1_L_FIR_FER/2; i++ )
        {
            /* causal ringing of the FIR filter */
            *pt_exc++ *= (gain_inov * gain);
        }
    }

    /*-----------------------------------------------------------------*
     * Total excitation
     *-----------------------------------------------------------------*/

    if( (st->last_coder_type == AUDIO || st->last_good == INACTIVE_CLAS) && st->total_brate <= ACELP_24k40 && !st->Opt_AMR_WB )
    {
        /* For GSC - the excitation is already computed */
        mvr2r( exc, exc2, st->L_frame);
    }
    else if(st->last_good >= UNVOICED_TRANSITION && st->last_good < INACTIVE_CLAS )
    {
        /* For voiced and generic signals - prepare a HP filter for the random part of excitation */
        for( i=0; i<MODE1_L_FIR_FER; i++ )
        {
            hp_filt[i] = (1.0f - st->tilt_code) * h_high[i];
        }

        /* HP filter the random part of the excitation and add the adaptive part */
        for( i=0; i< L_frame; i++ )
        {
            exc2[i] = exc[i] + dotp( &exc2_buf[i], hp_filt, MODE1_L_FIR_FER );
        }
    }
    else
    {
        /* For purely unvoiced signals - just copy the unfiltered random part of the excitation */
        mvr2r( exc2_buf + MODE1_L_FIR_FER/2, exc, L_frame );
        mvr2r( exc2_buf + MODE1_L_FIR_FER/2, exc2, L_frame );
    }

    if( L_frame == L_FRAME )
    {
        interp_code_5over2( exc, bwe_exc, L_frame );
    }
    else
    {
        interp_code_4over2( exc, bwe_exc, L_frame );
    }

    /*-----------------------------------------------------------------*
     * Updates
     *-----------------------------------------------------------------*/

    /* Update voicing factors of TBE */
    if( st->last_coder_type == AUDIO || st->last_good == INACTIVE_CLAS )
    {
        if(st->L_frame == L_FRAME )
        {
            set_f( voice_factors, 1.0f, NB_SUBFR );
        }
        else
        {
            set_f( voice_factors, 1.0f, NB_SUBFR16k );
        }
    }
    else
    {
        if(st->L_frame == L_FRAME )
        {
            set_f( voice_factors, st->last_voice_factor, NB_SUBFR);
        }
        else
        {
            set_f( voice_factors, st->last_voice_factor, NB_SUBFR16k );
        }
    }

    if( st->Opt_AMR_WB )
    {
        /* update buffer of gains for the next frame */
        gain_dec_bfi(st->past_qua_en);
    }

    st->bfi_pitch = pitch_buf[L_frame/L_SUBFR - 1];
    st->bfi_pitch_frame = st->L_frame;

    return;
}


/*-------------------------------------------------------------------*
 * gain_dec_bfi()
 *
 * Estimate past quantized gain prediction residual to be used in
 * next frame
 *-------------------------------------------------------------------*/

static void gain_dec_bfi(
    float *past_qua_en    /* i/o: gain quantization memory (4 words)  */
)
{
    short   i;
    float av_pred_en;

    av_pred_en = 0.0f;
    for (i = 0; i < GAIN_PRED_ORDER; i++)
    {
        av_pred_en += past_qua_en[i];
    }

    av_pred_en = (float)(av_pred_en*(1.0f/(float)GAIN_PRED_ORDER)-3.0f);

    if (av_pred_en < -14.0f)
    {
        av_pred_en = -14.0f;
    }

    for (i=GAIN_PRED_ORDER-1; i>0; i--)
    {
        past_qua_en[i] = past_qua_en[i-1];
    }

    past_qua_en[0] = av_pred_en;

    return;
}


/*-------------------------------------------------------------------*
 * pulseRes_preCalc()
 *
 * calculates some conditions for Pulse resynchronization to take place
 *-------------------------------------------------------------------*/

static void pulseRes_preCalc(
    Word16* cond1,
    Word16* cond2,
    Word32* cond3,
    Word16 new_pit,
    Word16 Tc,
    Word16 L_frame
)
{
    Word16 tmp_pit, tmp_pit_e, tmp_frame, tmp_frame_e;
    Word32 tmp_pit2;

    tmp_pit = BASOP_Util_Divide1616_Scale(new_pit/*Q0*/,Tc/*Q0*/,&tmp_pit_e)/*Q15*/;
    tmp_frame = add(  extract_l(L_mult0(L_frame ,  64/*1.f/L_SUBFR Q12*/)/*Q12*/) , 4096/*1.f Q12*/ );/*Q12*/
    tmp_frame = BASOP_Util_Divide1616_Scale(4096/*1.f Q12*/,tmp_frame, &tmp_frame_e);/*Q15*/
    tmp_frame = shl(tmp_frame,add(tmp_frame_e,1));
    tmp_frame = sub(32767/*1.f Q15*/, tmp_frame);/*Q15*/

    BASOP_SATURATE_WARNING_OFF
    /*To calc Q15 threshold, overflow may happen - do negation and compare with negated value to check also highest possible value*/
    tmp_pit = shl(negate(tmp_pit),tmp_pit_e);
    BASOP_SATURATE_WARNING_ON

    *cond1 = sub(tmp_pit, negate(tmp_frame));

    *cond2 = sub(Tc, new_pit);

    tmp_pit_e = BASOP_Util_Add_MantExp(new_pit,15-0,negate(Tc),15-0,&tmp_pit);/*Q15*/
    tmp_pit = abs_s(tmp_pit);
    tmp_pit2 = L_mult(Tc,4915/*0.15f Q15*/);/*Q16*/

    BASOP_SATURATE_WARNING_OFF
    /*To calc Q15 threshold, overflow may happen - do negation and compare with negated value to check also highest possible value*/
    tmp_pit2 = L_shl(L_negate(tmp_pit2),sub(15-16,tmp_pit_e));
    BASOP_SATURATE_WARNING_ON

    *cond3 = L_sub(L_mult0(-1, tmp_pit),tmp_pit2);

    return;
}
