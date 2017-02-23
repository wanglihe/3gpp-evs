/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "prot.h"
#include "options.h"


/*-------------------------------------------------------------------*
* con_acelp()
*
* Concealment function for ACELP and TD-TCX
*--------------------------------------------------------------------*/

void con_acelp(
    float A[],                   /* input: coefficients NxAz[M+1]    */
    int coder_type,              /* input: ACELP coder type          */
    float synth[],               /* i/o:   synthesis                 */
    int *pT,                     /* out:   pitch for all subframe    */
    float *pgainT,               /* out:   pitch gain for all subfr  */
    float stab_fac,              /* input: stability of isf          */
    Decoder_State *st,           /* i/o :  coder memory state        */
    float pitch_buffer[],        /* i/o: floating pitch values for each subframe          */
    float *voice_factors,        /* o  : voicing factors                                  */
    float *bwe_exc               /* o  : excitation for SWB TBE                           */
)
{
    int   i, i_subfr, L_frame, T0;
    float tmp_deemph;
    float mem_syn[M], mem_syn2[M], mem[M], *syn;
    float *noise_buf;
    float *exc, *harmonic_exc_buf, buf[L_EXC_MEM_DEC+M+L_FRAME16k+L_FRAME16k/2], *p_A;
    float pitch_buf[NB_SUBFR16k];
    float alpha = 0.0f;
    float step, gain, gainCNG, gain_inov, ftmp;
    float *pt_exc,tmp_tc,predPitchLag;
    float pc = 0.f;
    short extrapolationFailed,tmpSeed,Tc;
    float *pt1_exc;
    float *w = st->tcx_cfg.tcx_mdct_window;
    int   W1, W2, j;
    int   l_fir_fer = L_FIR_FER;
    float lpFiltAdapt[3];
    float hp_filt[3];
    int   nSubframes;
    float exc_unv[L_FRAME16k+L_FRAME16k/2];
    float syn_unv[L_FRAME16k+L_FRAME16k/2];
    float mem_syn_unv[M];
    float gain_lpc[NB_SUBFR16k];
    float h1[L_SUBFR+1];
    float gainSynthDeemph;
    float tmp = 0.f;
    int   fUseExtrapolatedPitch;

    /* Framing parameters */

    L_frame = st->L_frame;

    fUseExtrapolatedPitch = 0;
    extrapolationFailed = 1;

    /*------------------------------------------------------------------------*
     * Initialize buffers                                                     *
     *------------------------------------------------------------------------*/

    /* set ACELP synthesis memory */
    mvr2r( st->mem_syn2, mem_syn, M);

    /* set excitation memory*/
    harmonic_exc_buf = buf+M;
    exc = harmonic_exc_buf+L_EXC_MEM_DEC;
    mvr2r(st->old_exc, harmonic_exc_buf, L_EXC_MEM_DEC);
    exc[L_frame] = 0.0f;

    /*------------------------------------------------------------------------*
     * PLC: [ACELP:Extrapolate Pitch Lag]
     *------------------------------------------------------------------------*/

    if (st->flagGuidedAcelp == 1)
    {
        T0 = st->guidedT0;
    }

    pitch_pred_linear_fit( st->nbLostCmpt, st->last_good, st->old_pitch_buf, &st->old_fpitch,
                           &predPitchLag, st->pit_min, st->pit_max, st->mem_pitch_gain, 0,
                           st->plc_use_future_lag, &extrapolationFailed, st->nb_subfr );
    T0 = (int)(predPitchLag+0.5f);

    if( extrapolationFailed )
    {
        /*------------------------------------------------------------------------*
         * - Construct adaptive codebook from side information                    *
         *------------------------------------------------------------------------*/

        if( st->flagGuidedAcelp == 0)
        {
            nSubframes = 0;
        }
        else
        {
            nSubframes = 2;
            /* Construct adaptive codebook with T0, T0_frac, T0_res, gain_pit for 2 sub-frames */
            for( i=0; i<2*L_SUBFR; i++ )
            {
                exc[i] = exc[i-st->guidedT0];
            }
        }
    }
    else
    {
        nSubframes = 0;
    }


    if( nSubframes > 0 )
    {
        tmp_tc = (float)st->guidedT0;    /* take the transmit pitch*/
    }
    else
    {
        tmp_tc = st->old_fpitch;         /* take the previous frame last pitch*/
    }

    /* PLC: [ACELP: Fade-out]
     * PLC: calculate damping factor */
    alpha = Damping_fact(coder_type, st->nbLostCmpt, st->last_good, stab_fac, &(st->lp_gainp), 0);

    if (st->nbLostCmpt==1)
    {
        st->cummulative_damping = 1.0f;
    }
    else
    {
        st->cummulative_damping *= alpha;
    }

    /*-----------------------------------------------------------------*
     * PLC: [ACELP: adaptive codebook]
     * PLC: Construct the harmonic part of excitation
     *-----------------------------------------------------------------*/

    if( st->last_good >= UNVOICED_TRANSITION )
    {

        /*---------------------------------------------------------------*
         *  Last pitch cycle of the previous frame is repeatedly copied. *
         *---------------------------------------------------------------*/

        Tc = (short)(tmp_tc + 0.5f);

        if ((T0 > 0) && (T0 != Tc) && (abs(T0-Tc) < 0.15f*Tc) && extrapolationFailed == 0 )
        {
            fUseExtrapolatedPitch = 1;
        }

        if(st->enableGplc)
        {
            pt_exc = &exc[nSubframes*L_SUBFR];
        }
        else
        {
            pt_exc = exc;
        }
        pt1_exc = pt_exc - Tc;

        if (fUseExtrapolatedPitch != 0)
        {
            /* Required because later pt1_exc[1] used in filtering points to exc[0]. To make it safe also for GPL pt_exc is used instead of exc */
            pt_exc[0] = 0;
            pt_exc = harmonic_exc_buf;
            assert(pt_exc < pt1_exc-1);
        }

        if( st->nbLostCmpt == 1 )
        {
            /* pitch cycle is first low-pass filtered */
            /*get filter coefficients*/
            genPlcFiltBWAdap(st->sr_core, &lpFiltAdapt[0], 0, st->cummulative_damping);
            for( i=0 ; i< Tc; i++ )
            {
                *pt_exc++ = ( lpFiltAdapt[0] * pt1_exc[-1] + lpFiltAdapt[1] * pt1_exc[0] + lpFiltAdapt[2] * pt1_exc[1]);
                pt1_exc++;
            }
        }
        else
        {
            /* copy the first pitch cycle without low-pass filtering */
            for( i=0 ; i< Tc; i++ )
            {
                *pt_exc++ = *pt1_exc++;
            }
        }

        if (fUseExtrapolatedPitch != 0)
        {
            pt1_exc = harmonic_exc_buf;
        }

        for (i = 0; i < L_frame+(1-nSubframes)*L_SUBFR-Tc; i++)
        {
            *pt_exc++ = *pt1_exc++;
        }

        /*-------------------------------------------------------*
         *  PLC: [ACELP: adaptive codebook]
         *  PLC: Resync pulse positions.
         *-------------------------------------------------------*/

        if( nSubframes > 0 )
        {
            pitch_buf[0] = (float)st->guidedT0;
            pitch_buf[1] = (float)st->guidedT0;
        }

        if (nSubframes>0)
        {
            pitch_buf[3] = pitch_buf[2] = pitch_buf[1]; /* do not resync on second half of frame */
            if(st->nb_subfr == 5)
            {
                /* for guidedacelp cases and nSubframes=2, set pitch_buf[4] to avoid memory_access issues in post_decoder() */
                pitch_buf[4] = pitch_buf[3];
            }
        }
        else
        {
            if (fUseExtrapolatedPitch != 0)
            {
                get_subframe_pitch(st->nb_subfr, st->old_fpitch, predPitchLag, pitch_buf);

                PulseResynchronization(harmonic_exc_buf, exc, L_frame, st->nb_subfr, st->old_fpitch, predPitchLag);
            }
            else
            {
                set_f( pitch_buf, st->old_fpitch, st->nb_subfr);
            }
        }

        /*------------------------------------------------------------*
         *  PLC: [ACELP: adaptive codebook]
         *  PLC: Create the harmonic part needed for the overlap-add.
         *------------------------------------------------------------*/

        pt_exc = exc+L_frame;
        pt1_exc = pt_exc - ((T0 == 0) ? Tc : T0);

        for (i = 0; i < L_frame/2; i++)
        {
            *pt_exc++ = *pt1_exc++;
        }

        /*-------------------------------------------------------*
         *  PLC: [ACELP: adaptive codebook]
         *  PLC: update the floating point pitch for consecutive loss
         *-------------------------------------------------------*/

        if (fUseExtrapolatedPitch != 0)
        {

            if (st->flagGuidedAcelp == 1)
            {
                st->old_fpitch = (float) T0;
            }
            else
            {
                st->old_fpitch = predPitchLag;
            }
        }

        /*-------------------------------------------------------*
         *  PLC: [ACELP: adaptive BPF]
         *  PLC: Accommodate the BPF
         *-------------------------------------------------------*/

        st->bpf_gain_param = 3 ; /*full BPF*/

        /*-------------------------------------------------------*
         *  PLC: [ACELP: adaptive codebook]
         *  PLC: Calculate the initial gain and fade out step.
         *-------------------------------------------------------*/

        pc = (float) fabs( pitch_buf[3] + pitch_buf[2] - pitch_buf[1] - pitch_buf[0] ) * 256.0f / (float)L_frame;

        if ((st->last_good <= UNVOICED_TRANSITION) && (coder_type == GENERIC) && (pc > 6))
        {
            gain = 0.0f;
            st->lp_gainp = 0.0f;
        }
        else
        {
            gain = 1.0f; /* start-of-the-frame gain */
            st->lp_gainp = alpha;
        }

        step = (1.0f/L_frame) * (gain - st->lp_gainp);

        /*-------------------------------------------------------*
         *  PLC: [ACELP: Fade-out]
         *  Apply fade out
         *-------------------------------------------------------*/

        for (i_subfr = 0; i_subfr < st->nb_subfr; i_subfr++)
        {
            pgainT[i_subfr] = gain;

            for (i = i_subfr*L_SUBFR; i < (i_subfr+1)*L_SUBFR; i++)
            {
                exc[i] *= gain;
                gain -= step;
            }
        }

        for ( ; i < L_frame+(L_frame/2); i++ )
        {
            exc[i] *= gain;
            gain -= step;
        }

        for (i = 0; i < st->nb_subfr; i ++)
        {
            pT[i] = (int)( pitch_buf[i] + 0.5f);
            pitch_buffer[i] = pitch_buf[i];
        }

        /* update old exc without random part*/
        mvr2r(harmonic_exc_buf+L_frame, st->old_exc,  L_EXC_MEM_DEC);

    }
    else
    {
        /* No harmonic part */
        set_zero(&exc[0], L_frame+L_frame/2);

        for (i = 0; i < st->nb_subfr; i ++)
        {
            pitch_buf[i] = st->pit_max;
            pgainT[i] = 0.f;
            pT[i] = L_SUBFR;
            pitch_buffer[i] = L_SUBFR;
        }

        st->bpf_gain_param = 0; /*no BPF*/
    }

    /*-----------------------------------------------------------------*
     * Construct the random part of excitation
     *-----------------------------------------------------------------*/
    noise_buf = buf;
    tmpSeed = st->seed_acelp;

    for (i = 0; i < L_frame+l_fir_fer-1; i++)
    {
        noise_buf[i] = (float)own_random(&tmpSeed);
    }
    st->seed_acelp = tmpSeed;

    for ( ; i < L_frame+(L_frame/2)+l_fir_fer-1; i++)
    {
        noise_buf[i] = (float)own_random(&tmpSeed);
    }

    /*get filter coefficients*/
    genPlcFiltBWAdap(st->sr_core, &hp_filt[0], 1, st->cummulative_damping);

    /* PLC: [ACELP: Fade-out]
     * PLC: retrieve background level */
    tmp = 1.0f;
    gainSynthDeemph=getLevelSynDeemph( &(tmp), A, L_SUBFR, st->preemph_fac, L_frame/L_SUBFR);
    gainCNG = st->cngTDLevel/gainSynthDeemph;

    gain = st->lp_gainc; /* start-of-the-frame gain */

    ftmp = 2.0f * gain;

    if (gainCNG > ftmp)
    {
        gainCNG = ftmp;
    }

    st->lp_gainc = alpha * st->lp_gainc + (1.0f - alpha) * gainCNG;  /* end-of-the-frame gain */

    if( (st->last_good == UNVOICED_TRANSITION) && (coder_type == GENERIC))
    {
        st->lp_gainc = gainCNG;
    }

    highPassFiltering(st->last_good, L_frame+ l_fir_fer/2, noise_buf, hp_filt, l_fir_fer);

    /* Find energy normalization factor */
    pt_exc = noise_buf + l_fir_fer/2;
    gain_inov = 1.0f / (float)sqrt( dotp( pt_exc, pt_exc, L_frame ) / L_frame + 0.01f );

    /* PLC: [ACELP: Fade-out]
     * PLC: Linearly attenuate the gain through the frame */
    step = (1.0f/L_frame) * (gain - st->lp_gainc);

    if ((st->last_good == UNVOICED_CLAS) && (coder_type!=UNVOICED)) /* Attenuate somewhat on unstable unvoiced */
    {
        gain_inov *= 0.8f;
    }

    if (st->last_good >= UNVOICED_TRANSITION)
    {
        float tilt_code;

        tilt_code = (float)(0.10f*(1.0f + st->voice_fac));
        gain_inov *= (1.f - tilt_code);
    }

    pt_exc = noise_buf;

    /* non-causal ringing of the FIR filter   */

    for( i=0 ; i< l_fir_fer/2; i++ )
    {
        *pt_exc++ *= (gain_inov * gain);
    }

    /* Actual filtered random part of excitation */

    for( i=0 ; i< L_frame; i++ )
    {
        *pt_exc++ *= (gain_inov * gain);
        gain -= step;
    }

    for( i=0 ; i< (L_frame/2)+l_fir_fer/2; i++ )
    {
        *pt_exc++ *= (gain_inov * gain);
    }

    st->past_gcode = gain;

    if( st->last_good < UNVOICED_TRANSITION )
    {
        mvr2r(noise_buf+l_fir_fer/2, exc, L_frame+L_frame/2);
        mvr2r(harmonic_exc_buf+L_frame, st->old_exc, L_EXC_MEM_DEC);
        mvr2r( exc,  exc_unv, L_frame+ (L_frame/2));   /* Update exc_unv */
    }
    else
    {
        mvr2r( noise_buf+l_fir_fer/2, exc_unv, L_frame+ (L_frame/2));   /* Update exc_unv */
    }

    /* Compute total excitation in noisebuffer to save memories */
    if(st->last_good >= UNVOICED_TRANSITION)
    {
        for( i=0 ; i< L_frame+1; i++ )
        {
            noise_buf[i] = exc[i]+exc_unv[i];
        }
    }
    else
    {
        noise_buf = exc_unv;
    }

    if(L_frame == L_FRAME)
    {
        interp_code_5over2(noise_buf, bwe_exc, L_frame);
        set_f(voice_factors, st->last_voice_factor, NB_SUBFR);
    }
    else
    {
        interp_code_4over2(noise_buf, bwe_exc, L_frame);
        set_f(voice_factors, st->last_voice_factor, NB_SUBFR16k);
    }

    /*----------------------------------------------------------*
     * - compute the synthesis speech                           *
     *----------------------------------------------------------*/

    /* Init syn buffer */
    syn = buf + M;
    mvr2r(st->mem_syn2, buf, M );

    if (st->nbLostCmpt == 1)
    {

        if(st->last_good < UNVOICED_TRANSITION )
        {
            mvr2r(st->mem_syn2, mem_syn_unv, M );
        }
        else
        {
            set_zero( mem_syn_unv, M );
        }
    }
    else
    {
        mvr2r( st->mem_syn_unv_back, mem_syn_unv, M );
    }

    /* voiced synth */

    if(st->last_good >= UNVOICED_TRANSITION)
    {
        p_A = A;

        for (i_subfr = 0; i_subfr < L_frame; i_subfr += L_SUBFR)
        {
            tmp = 0;
            set_zero( h1, L_SUBFR+1 );
            set_zero( mem, M );
            h1[0] = 1.0f;
            syn_filt(p_A, M,h1, h1, L_SUBFR, mem, 0); /* impulse response of LPC     */
            deemph(h1, st->preemph_fac, L_SUBFR, &tmp); /* impulse response of deemph  */
            /* impulse response level = gain introduced by synthesis+deemphasis */
            gain_lpc[i_subfr/L_SUBFR] = 1.f/(float)sqrt(dotp( h1, h1, L_SUBFR));
            p_A += (M+1);
        }

        j=0;
        for (i_subfr = 0; i_subfr < L_frame; i_subfr += L_SUBFR)
        {

            for (i=0; i<L_SUBFR; i++)
            {
                exc[i_subfr + i] *= st->last_gain_syn_deemph*gain_lpc[j];
            }
            j++;
        }

        for (i=L_frame; i<L_frame+L_frame/2; i++)
        {
            exc[i] *= st->last_gain_syn_deemph*gain_lpc[3];
        }
        p_A = A;

        for (i_subfr = 0; i_subfr < L_frame; i_subfr += L_SUBFR)
        {
            syn_filt(p_A, M,&exc[i_subfr], &syn[i_subfr], L_SUBFR, mem_syn, 1);
            p_A += (M+1);
        }

        mvr2r( mem_syn, mem_syn2, M );

        /* synthesize ola*/
        syn_filt(p_A-(M+1), M, &exc[L_frame], &syn[L_frame], (L_frame/2), mem_syn2, 0);
    }

    /* unvoiced synth */
    tmp = 0;
    p_A = st->Aq_cng;

    for (i_subfr = 0; i_subfr < L_frame; i_subfr += L_SUBFR)
    {
        set_zero(h1, L_SUBFR+1);
        set_zero(mem, M);
        h1[0] = 1.0f;
        syn_filt(p_A, M,h1, h1, L_SUBFR, mem, 0); /* impulse response of LPC     */
        deemph(h1, st->preemph_fac, L_SUBFR, &tmp); /* impulse response of deemph  */
        /* impulse response level = gain introduced by synthesis+deemphasis */
        gain_lpc[i_subfr/L_SUBFR] = 1.f/(float)sqrt(dotp( h1, h1, L_SUBFR));
        p_A += (M+1);
    }

    j=0;
    for (i_subfr = 0; i_subfr < L_frame; i_subfr += L_SUBFR)
    {

        for (i=0; i<L_SUBFR; i++)
        {
            exc_unv[i_subfr + i] *= st->last_gain_syn_deemph*gain_lpc[j];
        }
        j++;
    }

    for (i=L_frame; i<L_frame+L_frame/2; i++)
    {
        exc_unv[i] *= st->last_gain_syn_deemph*gain_lpc[j-1];
    }
    p_A = st->Aq_cng;

    for (i_subfr = 0; i_subfr < L_frame; i_subfr += L_SUBFR)
    {
        syn_filt(p_A, M,&exc_unv[i_subfr], &syn_unv[i_subfr], L_SUBFR, mem_syn_unv, 1);
        p_A += (M+1);
    }

    mvr2r(mem_syn_unv,st->mem_syn_unv_back,M);

    if(st->last_good < UNVOICED_TRANSITION)
    {
        mvr2r(mem_syn_unv,mem_syn,M);
    }

    /* unvoiced for ola */
    syn_filt(p_A-(M+1), M,&exc_unv[i_subfr], &syn_unv[i_subfr], (L_frame/2), mem_syn_unv, 0);

    /* add separate synthesis buffers */
    if(st->last_good >= UNVOICED_TRANSITION)
    {
        for( i=0 ; i < L_frame; i++ )
        {
            syn[i] += syn_unv[i];
        }
    }
    else
    {
        mvr2r(syn_unv,syn,L_frame+L_frame/2);
    }


    /* update buffer for the classification */
    FEC_clas_estim( syn, pitch_buf, st->L_frame, coder_type, st->codec_mode, st->mem_syn_clas_estim, &(st->clas_dec),
                    &st->lp_ener_bfi, st->core_brate, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -1.0f, st->narrowBand, 0,
                    1, st->preemph_fac, st->tcxonly, st->last_core_brate );

    /* Update Pitch Lag memory */
    mvr2r( &st->old_pitch_buf[L_frame/L_SUBFR], st->old_pitch_buf, L_frame/L_SUBFR );
    mvr2r( pitch_buf, &st->old_pitch_buf[L_frame/L_SUBFR], L_frame/L_SUBFR );

    /*updating enr_old parameters*/
    fer_energy(L_frame, st->last_good, syn, tmp_tc, &(st->enr_old), 1);

    /* update ACELP synthesis memory */
    mvr2r( mem_syn, st->mem_syn2, M );
    mvr2r( syn+L_frame-L_SYN_MEM, st->mem_syn_r, L_SYN_MEM );

    /* Deemphasis and output synth */
    tmp_deemph = st->syn[M];
    deemph(syn, st->preemph_fac, L_frame+L_frame/2, &tmp_deemph);

    mvr2r(syn, synth, L_frame);
    mvr2r(syn+L_frame-L_frame/2, st->old_syn_Overl, L_frame/2);


    /* save last half frame if next frame is TCX */
    mvr2r(syn+L_frame, st->syn_Overl_TDAC, L_frame/2);
    mvr2r(syn+L_frame-M-1, st->syn, 1+M);

    /* update old_Aq */
    mvr2r(p_A-(M+1), st->old_Aq_12_8, M+1);


    mvr2r(syn+L_frame, st->syn_Overl, L_frame/2);

    W1 = st->tcx_cfg.tcx_mdct_window_length;
    W2 = st->tcx_cfg.tcx_mdct_window_length/2;

    st->tcx_cfg.tcx_curr_overlap_mode = FULL_OVERLAP;
    {
        short n = (short)((float)L_frame*N_ZERO_MDCT_NS/FRAME_SIZE_NS);

        mvr2r(syn+L_frame-n, st->old_outLB, L_frame-n);

        for (i=0; i<W1; i++)
        {
            st->old_outLB[i+n] *= w[W1-1-i]*w[W1-1-i];
        }
        set_zero(&st->old_outLB[W1+n], n);
    }

    /* create aliasing and windowing */

    for (i=0; i<W1; i++)
    {
        buf[i] = st->syn_Overl_TDAC[i]*w[W1-1-i];
    }

    for (i=0; i<W2; i++)
    {
        st->syn_Overl_TDAC[i] = buf[i]+buf[W1-1-i]; /* A-D */
    }

    for (i=0; i<W2; i++)
    {
        st->syn_Overl_TDAC[W2+i] = buf[W2+i]+buf[W1-1-W2-i]; /* B-C */
    }

    for (i=0; i<W1; i++)
    {
        st->syn_Overl_TDAC[i] *= w[W1-1-i];
    }

    /* update memory for full band */
    lerp(st->syn_Overl_TDAC, st->syn_Overl_TDACFB, st->L_frameTCX/2, L_frame/2);
    lerp(st->syn_Overl, st->syn_OverlFB, st->L_frameTCX/2, L_frame/2);
    lerp(st->old_outLB, st->old_out, st->L_frameTCX, L_frame);

    /* copy total excitation exc2 as 16kHz for acelp mode1 decoding */
    lerp(exc, st->old_exc2, L_EXC_MEM, L_frame);
    lerp(syn, st->old_syn2, L_EXC_MEM, L_frame);
    st->bfi_pitch = pitch_buf[st->nb_subfr-1];
    st->bfi_pitch_frame = L_frame;

    return;
}
