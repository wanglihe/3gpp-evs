/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "options.h"
#include "prot.h"
#include "rom_dec.h"


/*-----------------------------------------------------------------*
 * con_tcx()
 *
 *
 *-----------------------------------------------------------------*/

void con_tcx(
    Decoder_State* st,                    /* i/o: coder memory state        */
    float          synth[]                /* i/o: synth[]                   */
)
{
    short i, n, L_frame, L_subfr, fLowPassFilter, T0;
    int offset = 0;
    int mem_syn_r_size_old, mem_syn_r_size_new;
    float *noise;
    float mem_syn[M], *syn;
    float *exc, buf[OLD_EXC_SIZE_DEC+L_FRAME_MAX+L_FRAME_MAX/NB_SUBFR+1+L_FRAME_MAX/2];
    float pre_emph_buf;
    float pitch_buf[NB_SUBFR16k], hp_filt[L_FIR_FER2];
    float predPitchLag, alpha = 0.0f;
    float tmp_deemph, step, gain, gainCNG, gain_inov, ftmp;
    float *pt_exc, *pt1_exc;
    short Tc, tmpSeed;
    short fUseExtrapolatedPitch=0;
    float *ana_window;
    float r[M+1], A_local[M+1], mem;
    float *w;
    short W1, W2;
    short extrapolationFailed=1;
    float gainSynthDeemph;
    float tmp=0.f;
    float old_pitch_buf[2*NB_SUBFR16k+2];

    /* Framing parameters */
    L_frame = st->L_frameTCX;
    L_subfr = st->L_frameTCX/st->nb_subfr;
    w = st->tcx_cfg.tcx_mdct_windowFB;
    W1 = st->tcx_cfg.tcx_mdct_window_lengthFB;
    W2 = st->tcx_cfg.tcx_mdct_window_lengthFB/2;

    /* take the previous frame last pitch */
    Tc = (short)(st->old_fpitchFB + 0.5f);

    set_zero(buf,sizeof(buf)/sizeof(buf[0]));

    v_multc(st->old_pitch_buf, (float)L_frame / st->L_frame, old_pitch_buf, 2*NB_SUBFR16k+2);

    /* set excitation memory*/
    exc = buf+OLD_EXC_SIZE_DEC;

    tmp_deemph = synth[-1];
    pre_emph_buf = synth[-1];

    if( st->nbLostCmpt == 1 || st->tcxConceal_recalc_exc )
    {
        /* apply pre-emphasis to the signal */
        mem = synth[-(L_frame/2+st->pit_max_TCX+2*M)-1];

        preemph(&synth[-(L_frame/2+st->pit_max_TCX+2*M)], st->preemph_fac, L_frame/2+st->pit_max_TCX+2*M, &mem);
        st->lp_gainc = 0.0f;

        st->lp_gainp = get_gain( synth-2*L_subfr, synth-2*L_subfr-Tc, 2*L_subfr, NULL );

        if(st->lp_gainp < 0.0f)
        {
            st->lp_gainp = 0.0f;
        }

        if(st->lp_gainp > 1.0f)
        {
            st->lp_gainp = 1.0f;
        }
        ana_window = buf;
        ham_cos_window(ana_window, 3*L_frame/4, L_frame/4);

        /* Autocorrelation */
        autocorr( &(synth[-L_frame-1]), r, M, L_frame, ana_window, 0, 0, 0 );

        /* Lag windowing */
        lag_wind( r, M, st->output_Fs, LAGW_STRONG );

        /* Levinson Durbin */
        lev_dur(A_local, r, M, NULL);

        /* copy for multiple frame loss */
        mvr2r( A_local, st->old_Aq_12_8, M+1 );

        /* Residu */
        assert((2*L_subfr+Tc+1+M) <= st->old_synth_lenFB);
        residu(A_local, M, &(synth[-(2*L_subfr+Tc+1+M)]), &(exc[-(2*L_subfr+Tc+1+M)]), 2*L_subfr+Tc+1+M);
    }
    else
    {
        /* apply pre-emphasis to the signal */
        mem = synth[-L_frame-1];

        preemph(&synth[-L_frame], st->preemph_fac, L_frame, &mem);
        mvr2r( st->old_Aq_12_8, A_local, M+1 );
        offset = L_frame/2;

        if(st->last_good >= UNVOICED_TRANSITION )
        {
            i = max(Tc - L_frame/2, 0);
            mvr2r(st->old_excFB, &(exc[-i]), offset+i);
        }
        else
        {
            mvr2r(st->old_excFB, &(exc[-2*L_subfr]), 2*L_subfr+offset);
        }
    }

    /*-----------------------------------------------------------------*
     * PLC: Construct the harmonic part of excitation
     *-----------------------------------------------------------------*/

    if(st->last_good != UNVOICED_CLAS && !((st->last_good == UNVOICED_TRANSITION) && (st->core_ext_mode == GENERIC)) )
    {
        if ( st->nbLostCmpt == 1 || st->tcxConceal_recalc_exc )
        {
            st->lp_gainc = 0.0f;

            for ( i=0; i<2*L_subfr; i++ )
            {
                st->lp_gainc += ( exc[i-2*L_subfr] - st->lp_gainp * exc[i-2*L_subfr-Tc] ) * ( exc[i-2*L_subfr] - st->lp_gainp * exc[i-2*L_subfr-Tc] );
            }
            st->lp_gainc = (float)sqrt(st->lp_gainc / (2.0f*L_subfr) );
        }
        if( ( st->nbLostCmpt == 1 ) && st->rf_frame_type >= RF_TCXFD && st->rf_frame_type <= RF_TCXTD2 && st->use_partial_copy )
        {
            predPitchLag = (st->tcxltp_pitch_int + st->tcxltp_pitch_fr/(float)st->pit_res_max) * (float)st->L_frameTCX / (float)st->L_frame;
            T0 = (int)(predPitchLag+0.5f);

            if ((T0 > 0) && (T0 != Tc) && ((float)abs(T0-Tc) < 0.15f*Tc) )
            {
                fUseExtrapolatedPitch = 1;
            }
        }
        else
        {
            pitch_pred_linear_fit( st->nbLostCmpt, st->last_good, old_pitch_buf, &(st->old_fpitchFB),
                                   &predPitchLag, st->pit_min_TCX, st->pit_max_TCX, st->mem_pitch_gain,
                                   st->output_Fs > 25600, st->plc_use_future_lag, &extrapolationFailed, st->nb_subfr);

            T0 = (int)(predPitchLag+0.5f);

            if ((T0 > 0) && (T0 != Tc) && ((float)abs(T0-Tc) < 0.15f*Tc) && extrapolationFailed == 0 )
            {
                fUseExtrapolatedPitch = 1;
            }
        }

        fLowPassFilter = 0;
        pt_exc = exc + offset;

        pt1_exc = pt_exc - Tc;

        if (fUseExtrapolatedPitch != 0)
        {
            pt_exc = buf;
        }

        if(st->stab_fac < 1 && st->nbLostCmpt == 1 )
        {
            /* pitch cycle is first low-pass filtered */
            for( i=0 ; i< Tc; i++ )
            {
                if (st->output_Fs <= 16000)
                {
                    *pt_exc++ = ( 0.0053f * pt1_exc[-5] +
                                  0.0000f * pt1_exc[-4] +
                                  -0.0440f * pt1_exc[-3] +
                                  0.0000f * pt1_exc[-2] +
                                  0.2637f * pt1_exc[-1] +
                                  0.5500f * pt1_exc[0] +
                                  0.2637f * pt1_exc[1] +
                                  0.0000f * pt1_exc[2] +
                                  -0.0440f * pt1_exc[3] +
                                  0.0000f * pt1_exc[4] +
                                  0.0053f * pt1_exc[5]);
                }
                else /*(st->output_Fs >= 32000)*/
                {
                    *pt_exc++ = (-0.0053f * pt1_exc[-5] +
                                 -0.0037f * pt1_exc[-4] +
                                 -0.0140f * pt1_exc[-3] +
                                 0.0180f * pt1_exc[-2] +
                                 0.2668f * pt1_exc[-1] +
                                 0.4991f * pt1_exc[0] +
                                 0.2668f * pt1_exc[1] +
                                 0.0180f * pt1_exc[2] +
                                 -0.0140f * pt1_exc[3] +
                                 -0.0037f * pt1_exc[4] +
                                 -0.0053f * pt1_exc[5]);
                }
                pt1_exc++;
            }
            fLowPassFilter = 1;
        }
        else
        {
            /* copy the first pitch cycle without low-pass filtering */
            for( i=0 ; i< Tc; i++ )
            {
                *pt_exc++ = *pt1_exc++;
            }
            fLowPassFilter = 1;
        }

        if (fUseExtrapolatedPitch != 0)
        {
            pt1_exc = buf;
        }

        for (i = 0; i < L_frame-fLowPassFilter*Tc+L_subfr; i++)
        {
            *pt_exc++ = *pt1_exc++;
        }

        if (fUseExtrapolatedPitch != 0)
        {
            get_subframe_pitch(st->nb_subfr,  st->old_fpitch, predPitchLag * st->L_frame/L_frame, pitch_buf);

            PulseResynchronization(buf, exc, L_frame, st->nb_subfr, st->old_fpitchFB, predPitchLag);
        }
        else
        {
            set_f(pitch_buf, st->old_fpitch, st->nb_subfr);
        }

        if ( st->nbLostCmpt == 1 )
        {
            pt_exc = exc+L_frame;
            pt1_exc = pt_exc - ((T0 == 0) ? Tc : T0);

            for (i = 0; i < L_frame/2; i++)
            {
                *pt_exc++ = *pt1_exc++;
            }
        }

        if (fUseExtrapolatedPitch != 0)
        {
            st->old_fpitchFB = predPitchLag;
        }
        st->bpf_gain_param = 0;

        /* PLC: calculate damping factor */
        alpha = Damping_fact(st->core_ext_mode, st->nbLostCmpt, st->last_good, st->stab_fac, &(st->lp_gainp), 0 );

        if( st->nbLostCmpt == 1 )
        {
            st->cummulative_damping = 1;
        }
        else
        {
            st->cummulative_damping *= alpha;
        }

        gain = 1.0f;
        if( st->rf_frame_type == RF_TCXTD1 && st->use_partial_copy == 1 )
        {
            gain = 0.5f;
        }

        step = (1.0f/(L_frame+(L_frame/2))) * (gain - alpha);

        /* PLC: Apply fade out */
        for( i=offset; i < L_frame+(L_frame/2); i++ )
        {
            exc[i] *= gain;
            gain -= step;
        }

        offset = max(((short)(st->old_fpitchFB + 0.5f)) - L_frame/2, 0);
        mvr2r(exc+L_frame-offset, st->old_excFB, L_frame/2+offset);

        /* copy old_exc as 16kHz for acelp decoding */
        if(st->nbLostCmpt == 1)
        {
            lerp(exc-L_frame/2, st->old_exc, L_EXC_MEM_DEC, L_frame+L_frame/2);
        }
        else
        {
            mvr2r(st->old_exc+L_FRAME16k, st->old_exc, L_FRAME16k/2);
            lerp(exc, st->old_exc+L_FRAME16k/2, L_FRAME16k, L_frame);
        }
    }
    else
    {
        /* No harmonic part */
        set_zero(&exc[0], L_frame+L_frame/2);

        if ( st->nbLostCmpt == 1 )
        {
            st->lp_gainc = 0.0f;

            for ( i=0; i<2*L_subfr; i++ )
            {
                st->lp_gainc += ( exc[i-2*L_subfr] ) * ( exc[i-2*L_subfr]);
            }
            st->lp_gainc = (float)sqrt(st->lp_gainc / (2.0f*L_subfr) );
        }

        set_f( pitch_buf, (float)L_SUBFR, st->nb_subfr);

        /* PLC: calculate damping factor */
        alpha = Damping_fact(st->core_ext_mode, st->nbLostCmpt, st->last_good, st->stab_fac, &(st->lp_gainp), 0);
    }

    /*-----------------------------------------------------------------*
     * Construct the random part of excitation
     *-----------------------------------------------------------------*/

    tmpSeed = st->seed_acelp;
    noise = buf;

    for (i = 0; i < L_frame+L_FIR_FER2-1; i++)
    {
        noise[i] = (float)own_random(&tmpSeed);
    }
    st->seed_acelp = tmpSeed;

    for ( ; i < L_frame+(L_frame/2)+2*L_FIR_FER2; i++)
    {
        noise[i] = (float)own_random(&tmpSeed);
    }

    if (st->last_good == VOICED_CLAS || st->last_good == ONSET)
    {
        mem = noise[0];
        preemph(&noise[1], st->output_Fs<=16000 ? 0.2f : 0.6f, L_frame+(L_frame/2)+L_FIR_FER2, &mem);
    }

    /* high rate filter tuning */
    if(st->output_Fs<=16000)
    {
        for( i=0; i< L_FIR_FER2; i++ )
        {
            hp_filt[i] = h_high3_16[i];
        }
    }
    else /*(st->output_Fs==32000)*/
    {
        for( i=0; i< L_FIR_FER2; i++ )
        {
            hp_filt[i] = h_high3_32[i];
        }
    }

    if ( st->nbLostCmpt == 1 )
    {
        highPassFiltering(st->last_good, L_frame+L_frame/2+L_FIR_FER2, noise, hp_filt, L_FIR_FER2);
    }
    else
    {
        if(st->last_good > UNVOICED_TRANSITION )
        {
            for( i=0 ; i< L_frame+L_frame/2+L_FIR_FER2; i++ )
            {
                noise[i] = (1-st->cummulative_damping)*noise[i]+ st->cummulative_damping*dotp(&noise[i], hp_filt, L_FIR_FER2 );
            }
        }
    }

    /* PLC: [TCX: Fade-out] retrieve background level */
    tmp = 1.0f;
    gainSynthDeemph = getLevelSynDeemph(&(tmp), A_local, L_frame/4, st->preemph_fac, 1);
    if (st->tcxonly)
    {
        gainCNG = st->CngLevelBackgroundTrace_bfi/gainSynthDeemph;
    }
    else
    {
        gainCNG = st->cngTDLevel/gainSynthDeemph;
    }
    gain = st->lp_gainc; /* start-of-the-frame gain */

    if( st->rf_frame_type == RF_TCXTD1 && st->use_partial_copy == 1 )
    {
        gain *= 0.7f;
    }

    ftmp = 2.0f * gain;

    if (gainCNG > ftmp)
    {
        gainCNG = ftmp;
    }

    st->lp_gainc = alpha * st->lp_gainc + (1.0f - alpha) * gainCNG;  /* end-of-the-frame gain */

    /* PLC: [TCX: Fade-out] Linearly attenuate the gain through the frame */
    step = (1.0f/L_frame) * (gain - st->lp_gainc);
    pt_exc = noise + L_FIR_FER2/2;

    gain_inov = 1.0f / (float)sqrt( dotp( pt_exc, pt_exc, L_frame ) / L_frame + 0.01f );/* normalize energy */

    if ((st->last_good == UNVOICED_CLAS) && (st->core_ext_mode != UNVOICED))
    {
        gain_inov *= 0.8f;
    }
    else if (!((st->last_good == UNVOICED_CLAS) || (st->last_good == UNVOICED_TRANSITION)))
    {
        gain_inov *= (1.1f- 0.75*st->lp_gainp);
    }
    st->lp_gainp = alpha;
    pt_exc = noise;            /* non-causal ringing of the FIR filter   */

    for( i=0 ; i< L_FIR_FER2/2; i++ )
    {
        *pt_exc++ *= (gain_inov * gain);
    }

    for( i=0 ; i< L_frame+L_FIR_FER2/2; i++ )     /* Actual filtered random part of excitation */
    {
        *pt_exc++ *= (gain_inov * gain);
        gain -= step;
    }

    for( i=0 ; i< (L_frame/2); i++ )   /* causal ringing of the FIR filter */
    {
        *pt_exc++ *= (gain_inov * gain);
    }

    /*-----------------------------------------------------------------*
     * Construct the total excitation
     *-----------------------------------------------------------------*/

    if(st->last_good >= UNVOICED_TRANSITION )
    {
        for( i=0 ; i< (L_frame+L_frame/2); i++ )
        {
            exc[i] += noise[i+(L_FIR_FER2/2)];
        }
    }
    else
    {
        mvr2r(noise+L_FIR_FER2/2, exc, L_frame+L_frame/2);
        mvr2r(exc+L_frame-2*L_subfr, st->old_excFB, 2*L_subfr+L_frame/2);

        /* copy old_exc as 16kHz for acelp decoding */
        if(st->nbLostCmpt == 1)
        {
            lerp(exc, st->old_exc, L_EXC_MEM_DEC, L_frame+L_frame/2);
        }
        else
        {
            mvr2r(st->old_exc+L_FRAME16k, st->old_exc, L_FRAME16k/2);
            lerp(exc, st->old_exc+L_FRAME16k/2, L_FRAME16k, L_frame);
        }
    }

    /* Update Pitch Lag memory */
    mvr2r( &st->old_pitch_buf[st->nb_subfr], st->old_pitch_buf, st->nb_subfr );
    mvr2r( pitch_buf, &st->old_pitch_buf[st->nb_subfr], st->nb_subfr );

    /*----------------------------------------------------------*
     * - compute the synthesis speech                           *
     *----------------------------------------------------------*/

    syn = buf + M;
    mvr2r(synth-M, buf, M);
    mvr2r( buf, mem_syn, M);

    syn_filt(A_local, M, &exc[0], &syn[0], L_frame+ (L_frame/2), mem_syn, 1);

    n = (short)((float)L_frame*N_ZERO_MDCT_NS/FRAME_SIZE_NS);

    /* update ACELP synthesis memory */
    mem_syn_r_size_old=(int)(1.25*L_frame/20.f);

    /* copy mem_syn as 16kHz */
    mem_syn_r_size_new=(int)(1.25*L_FRAME16k/20.f);
    mvr2r(syn+L_frame-L_SYN_MEM, st->mem_syn_r, L_SYN_MEM);
    lerp( st->mem_syn_r+L_SYN_MEM-mem_syn_r_size_old, st->mem_syn_r+L_SYN_MEM-mem_syn_r_size_new, mem_syn_r_size_new, mem_syn_r_size_old );
    mvr2r( st->mem_syn_r+L_SYN_MEM-M, st->mem_syn2, M );

    /* Deemphasis and output synth and ZIR */
    deemph(syn, st->preemph_fac, L_frame+L_frame/2, &tmp_deemph);
    mvr2r(syn+L_frame-M-1, st->syn, 1+M);


    lerp( syn+L_frame-L_frame/2, st->old_syn_Overl, st->L_frame/2, L_frame/2 );
    mvr2r(syn+L_frame-n, st->old_out, L_frame-n);

    for (i = 0; i < W1; i++)
    {
        st->old_out[i+n] *= w[W1-1-i]*w[W1-1-i];
    }
    set_zero(&st->old_out[W1+n], n);

    mvr2r(syn, synth, L_frame);

    mvr2r(syn+L_frame, st->syn_OverlFB, L_frame/2);

    /* copy total excitation exc2 as 16kHz for ACELP MODE1 decoding */
    lerp( exc, st->old_exc2, L_EXC_MEM, L_frame );
    lerp( syn, st->old_syn2, L_EXC_MEM, L_frame );
    st->bfi_pitch = pitch_buf[st->nb_subfr-1];
    st->bfi_pitch_frame = st->L_frame;

    /* create aliasing and windowing need for transition to TCX10/5 */
    mvr2r(syn+L_frame, st->syn_Overl_TDACFB, L_frame/2);

    for (i=0; i<W1; i++)
    {
        buf[i] = st->syn_Overl_TDACFB[i]*w[W1-1-i];
    }

    for (i=0; i<W2; i++)
    {
        st->syn_Overl_TDACFB[i] = buf[i]+buf[W1-1-i]; /* A-D */
    }

    for (i=0; i<W2; i++)
    {
        st->syn_Overl_TDACFB[W2+i] = buf[W2+i]+buf[W1-1-W2-i];/* B-C */
    }

    for (i=0; i<W1; i++)
    {
        st->syn_Overl_TDACFB[i] *= w[W1-1-i];
    }

    st->tcx_cfg.tcx_curr_overlap_mode = FULL_OVERLAP;
    synth[-1] = pre_emph_buf;
    /* update memory for low band */
    lerp( st->syn_OverlFB, st->syn_Overl, st->L_frame/2, L_frame/2 );
    lerp( st->syn_Overl_TDACFB, st->syn_Overl_TDAC, st->L_frame/2, L_frame/2 );
    lerp( st->old_out, st->old_outLB, st->L_frame, L_frame );


    st->old_enr_LP = enr_1_Az( A_local, L_SUBFR );

    return;
}
