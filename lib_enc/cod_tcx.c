/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "cnst.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "prot.h"


/*-------------------------------------------------------------------*
* HBAutocorrelation()
*
*
*-------------------------------------------------------------------*/

void HBAutocorrelation(
    TCX_config *tcx_cfg,      /* input: configuration of TCX              */
    int left_overlap_mode,    /* input: overlap mode of left window half  */
    int right_overlap_mode,   /* input: overlap mode of right window half */
    float speech[],           /* input: synthesis                         */
    int L_frame,              /* input: frame length                      */
    float *r,                 /* output: autocorrelations vector          */
    int m                     /* input : order of LP filter               */
)
{
    int i, j, left_overlap, right_overlap;
    float s;
    float xn_buf[L_MDCT_OVLP_MAX+L_FRAME_PLUS+L_MDCT_OVLP_MAX];

    /*-----------------------------------------------------------*
     * Windowing                                                 *
     *-----------------------------------------------------------*/

    WindowSignal( tcx_cfg, tcx_cfg->tcx_offset, left_overlap_mode, right_overlap_mode,
                  &left_overlap, &right_overlap, speech, &L_frame, xn_buf, 0 );

    /*-----------------------------------------------------------*
     * Autocorrelation                                           *
     *-----------------------------------------------------------*/

    for (i = 0; i <= m; i++)
    {
        s = 0.0;

        for (j = 0; j < L_frame+(left_overlap+right_overlap)/2-i; j++)
        {
            s += xn_buf[j]*xn_buf[j+i];
        }
        r[i] = s;
    }

    if (r[0] < 100.0)
    {
        r[0] = 100.0;
    }

    return;

}

/*-------------------------------------------------------------------*
* TNSAnalysis()
*
*
*-------------------------------------------------------------------*/

void TNSAnalysis(
    TCX_config *tcx_cfg,  /* input: configuration of TCX */
    int L_frame,          /* input: frame length */
    int L_spec,
    const short tcxMode,  /* input: TCX mode for the frame/subframe - TCX20 | TCX10 | TCX 5 (meaning 2 x TCX 5) */
    int isAfterACELP,     /* input: Flag indicating if the last frame was ACELP. For the second TCX subframe it should be 0  */
    float spectrum[],     /* input: MDCT spectrum of the subframe */
    STnsData * pTnsData,  /* output: Tns data */
    int * pfUseTns,       /* output: Flag indicating if TNS is used */
    float* predictionGain
)
{
    float buff[8]; /* Buffer for the rearrangement of LF TCX5 */

    /* Init TNS */
    *pfUseTns = 0;

    if (tcx_cfg->fIsTNSAllowed)
    {

        tcx_cfg->pCurrentTnsConfig = &tcx_cfg->tnsConfig[tcxMode == TCX_20][isAfterACELP];
        L_spec = tcx_cfg->pCurrentTnsConfig->iFilterBorders[0];

        /*-----------------------------------------------------------*
         * Temporal Noise Shaping analysis                           *
         *-----------------------------------------------------------*/

        if (tcxMode == TCX_5)
        {
            /* rearrange LF sub-window lines prior to TNS analysis & filtering */
            if (L_spec < L_frame/2)
            {
                mvr2r(spectrum+8, spectrum+16, L_spec/2-8);
                mvr2r(spectrum+L_frame/4, spectrum+8, 8);
                mvr2r(spectrum+L_frame/4+8, spectrum+L_spec/2+8, L_spec/2-8);
            }
            else
            {
                mvr2r(spectrum+L_frame/4, buff, 8);
                mvr2r(spectrum+8, spectrum+16, L_frame/4-8);
                mvr2r(buff, spectrum+8, 8);
            }
        }

        *pfUseTns = DetectTnsFilt(tcx_cfg->pCurrentTnsConfig, spectrum, pTnsData, predictionGain);

        /* If TNS should be used then get the residual after applying it inplace in the spectrum */
        if (*pfUseTns)
        {
            ApplyTnsFilter(tcx_cfg->pCurrentTnsConfig, pTnsData, spectrum, 1);

        }

        if (tcxMode == TCX_5)
        {
            /* undo rearrangement of LF sub-window lines prior to TNS analysis */
            if (L_spec < L_frame/2)
            {
                mvr2r(spectrum+L_spec/2+8, spectrum+L_frame/4+8, L_spec/2-8);
                mvr2r(spectrum+8, spectrum+L_frame/4, 8);
                mvr2r(spectrum+16, spectrum+8, L_spec/2-8);
                set_zero(spectrum+L_spec/2, L_frame/4-L_spec/2);
                set_zero(spectrum+L_frame/4+L_spec/2, L_frame/4-L_spec/2);
            }
            else
            {
                mvr2r(spectrum+8, buff, 8);
                mvr2r(spectrum+16, spectrum+8, L_frame/4-8);
                mvr2r(buff, spectrum+L_frame/4, 8);
            }
        }
    }

    return;
}


/*-------------------------------------------------------------------*
* ShapeSpectrum()
*
*
*-------------------------------------------------------------------*/

void ShapeSpectrum(
    TCX_config *tcx_cfg,/*input: configuration of TCX                 */
    float A[],          /* input: quantized coefficients NxAz_q[M+1]  */
    float gainlpc[],    /* output: MDCT gains for the previous frame  */
    int L_frame_glob,   /* input: frame length                        */
    int L_spec,
    float spectrum[],   /* i/o: MDCT spectrum                         */
    int fUseTns,        /* output: Flag indicating if TNS is used     */
    Encoder_State *st
)
{
    int L_frame;
    int tcx_offset;
    float Ap[M+2];
    float gamma1;
    float max_low_pre = 0.f, max_high_pre = 0.f;

    /*-----------------------------------------------------------*
     * Init                                                      *
     *-----------------------------------------------------------*/

    /* Init lengths */
    L_frame = L_frame_glob;
    tcx_offset = tcx_cfg->tcx_offset;
    gamma1 = st->gamma;

    if (st->enableTcxLpc)
    {
        gamma1 = 1.0f;
    }

    /* if past frame is ACELP */

    if (st->last_core == 0)
    {
        L_frame += tcx_offset;
        L_spec += tcx_cfg->tcx_coded_lines >> 2;
        if(tcx_cfg->lfacNext<0)
        {
            L_frame -= tcx_cfg->lfacNext;
        }
    }

    tcxGetNoiseFillingTilt( A, L_frame, (st->total_brate >= ACELP_13k20 && !st->rf_mode), &st->noiseTiltFactor );

    /* Calculate Spectrum Flatness Measure for the TCX Concealment */
    if( st->enablePlcWaveadjust )
    {
        tcx_cfg->SFM2 = SFM_Cal(spectrum, min(200, L_frame));
    }

    if( (st->total_brate == ACELP_9k60 && st->bwidth == SWB) ||
            (st->total_brate == ACELP_13k20 && st->bwidth == SWB) )
    {
        int i;

        max_low_pre = 0.f;
        for (i = 0; i < L_frame; i++)
        {
            float tmp = fabs(spectrum[i]);
            if( tmp > max_low_pre )
            {
                max_low_pre = tmp;
            }
        }

        max_high_pre = 0.f;
        for (i = 0; i < L_spec - L_frame; i++)
        {
            float tmp = fabs(spectrum[L_frame + i]);
            if( tmp > max_high_pre )
            {
                max_high_pre = tmp;
            }
        }
    }

    /*-----------------------------------------------------------*
     * Pre-shaping in frequency domain using weighted LPC (Wz)   *
     *-----------------------------------------------------------*/
    weight_a(A, Ap, gamma1, M);

    lpc2mdct(Ap, M, gainlpc);

    mdct_preShaping(spectrum, L_frame, gainlpc);

    v_multc(spectrum+L_frame, 1.f/gainlpc[FDNS_NPTS-1], spectrum+L_frame, L_spec-L_frame);

    /* reduce the peaks in the IGF region, to make life of the core-coder easier... */
    if( (st->total_brate == ACELP_9k60 && st->bwidth == SWB) ||
            (st->total_brate == ACELP_13k20 && st->bwidth == SWB) )
    {
        int i, sf_width;
        int dist_low, dist_high;
        float max_fac;
        float max_low, max_low1, max_low2, max_high;

        if( st->tcx_lpc_shaped_ari )
        {
            max_fac = 1.5f;
        }
        else
        {
            max_fac = 3.f;
        }

        sf_width = L_frame / 2;

        max_low2 = 0.f;
        dist_low = 0;
        for (i = 0; i < sf_width; i++)
        {
            float tmp = fabs(spectrum[L_frame - 1 - i]);
            if( tmp > max_low2 )
            {
                max_low2 = tmp;
                dist_low = i;
            }
        }

        max_low1 = 0.f;
        for (i = 0; i < L_frame - sf_width; i++)
        {
            float tmp = fabs(spectrum[L_frame - sf_width - 1 - i]);
            if( tmp > max_low1 )
            {
                max_low1 = tmp;
            }
            if( tmp > max_low2 )
            {
                dist_low = sf_width + i;
            }
        }

        max_low = max(max_low1, max_low2);

        max_high = 0.f;
        dist_high = 0;
        for (i = 0; i < L_spec - L_frame; i++)
        {
            float tmp = fabs(spectrum[L_frame + i]);
            if( tmp > max_high)
            {
                max_high = tmp;
                dist_high = i;
            }
        }

        if( (4.f * dist_high * max_high > dist_low * max_low) &&  (16.f * max_low_pre > max_high_pre) && (max_high > max_fac * max_low2) )
        {
            float fac = max_fac * max_low2 / max_high;
            v_multc(spectrum + L_frame, fac, spectrum + L_frame, L_spec - L_frame);
        }

    }

    if (st->tcxonly && st->tcxltp && (st->tcxltp_gain > 0.0f) && !fUseTns )
    {
        PsychAdaptLowFreqEmph(spectrum, gainlpc);
    }

    return;
}


/*-------------------------------------------------------------------*
* QuantizeSpectrum()
*
*
*-------------------------------------------------------------------*/

void QuantizeSpectrum(
    TCX_config *tcx_cfg,/*input: configuration of TCX*/
    float A[],          /* input: quantized coefficients NxAz_q[M+1] */
    Word16 Aqind[],     /* input: frame-independent quantized coefficients (M+1) */
    float gainlpc[],    /* input: MDCT gains of the previous frame */
    float synth[],
    int L_frame_glob,   /* input: frame length             */
    int L_frameTCX_glob,
    int L_spec,
    int nb_bits,        /*input: bit budget*/
    int tcxonly,        /*input: only TCX flag*/
    float spectrum[],   /* i/o: MDCT spectrum, input is shaped MDCT spectrum */
    STnsData * pTnsData,/* input: Tns data */
    int fUseTns,        /* input: Flag indicating if TNS is used */
    int tnsSize,        /* input: number of tns parameters put into prm */
    LPD_state *LPDmem,  /*i/o: memories*/
    int prm[],          /* output: tcx parameters          */
    int frame_cnt,      /* input: frame counter in the super_frame */
    Encoder_State *st,
    CONTEXT_HM_CONFIG *hm_cfg
)
{
    int i, sqTargetBits, L_frame, tcx_offset, stop;
    int L_frameTCX;
    float fac_ns, ener, gain_tcx, gain_tcx_opt;
    float xn_buf[L_MDCT_OVLP_MAX+L_FRAME_PLUS+L_MDCT_OVLP_MAX];
    float x_orig[N_MAX];
    float sqGain = 1.0f;
    int sqBits;
    int sqBits_noStop;
    int overlap;
    int noiseFillingSize;
    int noiseTransWidth = MIN_NOISE_FILLING_HOLE;
    int nEncoded;
    int *sqQ;
    short LtpPitchLag;
    int ctxHmBits;
    int resQBits;
    int resQTargetBits = 0;
    short nf_seed = 0;
    int PeriodicityIndex, NumIndexBits;
    float *OriginalSpectrum;
    int nEncodedCtxHm, stopCtxHm, sqBitsCtxHm, Selector;
    int lastnz, lastnzCtxHm;
    float RelativeScore;
    int *signs;
    int signaling_bits;
    int *prm_ltp, *prm_tns, *prm_hm, *prm_lastnz, *prm_target;
    float SFM;
    float K, K2;
    float Aq_old[M+1];
    short aldo;       /* ALDO flag in current frame*/
    short nz;         /* non-zero length in ALDO window*/
    int maxNfCalcBw;

    /*-----------------------------------------------------------*
     * Init                                                      *
     *-----------------------------------------------------------*/

    /* Init lengths */
    L_frame = L_frame_glob;
    L_frameTCX = L_frameTCX_glob;
    overlap = tcx_cfg->tcx_mdct_window_length;
    aldo=0;
    nz=NS2SA(st->sr_core, N_ZERO_MDCT_NS);

    tcx_offset = tcx_cfg->tcx_offset;

    OriginalSpectrum = NULL;
    signs = NULL; /* silence warning */
    NumIndexBits = 0;
    sqBits = 0;
    ctxHmBits = 0;
    resQBits=0;
    prm_ltp = &prm[1+NOISE_FILL_RANGES];
    prm_tns = prm_ltp + LTPSIZE;
    prm_hm = prm_tns + tnsSize;
    prm_lastnz = prm_hm + 2;
    sqQ = prm_hm + NPRM_CTX_HM;

    /* if past frame is ACELP */

    if (st->last_core == 0)
    {
        tcx_cfg->last_aldo=0;

        L_frame += tcx_offset;
        L_frameTCX += tcx_cfg->tcx_offsetFB;
        L_spec += tcx_cfg->tcx_coded_lines >> 2;
        if(tcx_cfg->lfacNext<0)
        {
            L_frame -= tcx_cfg->lfacNext;
            L_frameTCX -= tcx_cfg->lfacNextFB;
            tcx_offset = tcx_cfg->lfacNext;
        }
        else
        {
            tcx_offset = 0;
        }
        st->noiseLevelMemory = 0;
    }


    lsp2a_stab( st->lsp_old, Aq_old, M );

    /* target bitrate for SQ */
    sqTargetBits = nb_bits-7-NBITS_NOISE_FILL_LEVEL;

    /*Unquantized spectrum here*/
    if(st->enablePlcWaveadjust)
    {

        SFM = SFM_Cal(spectrum, min(200, L_frame_glob));

        if (L_frame_glob <= 256)
        {
            K  = 0.4f;
            K2 = 0.1f;
        }
        else if (L_frame_glob == 320 || L_frame_glob == 512)
        {
            K  = 0.4f;
            K2 = 0.1f;
        }
        else /*FrameSize_Core == 640*/
        {
            K  = 0.35f;
            K2 = 0.04f;
        }


        if (SFM < K)
        {
            st->Tonal_SideInfo = 1;
        }
        else
        {
            st->Tonal_SideInfo = 0;
        }

        if (tcx_cfg->SFM2 < K2)
        {
            st->Tonal_SideInfo = 1;
        }
    }

    /* Save pre-shaped spectrum*/
    mvr2r(spectrum, x_orig, L_spec);

    /*-----------------------------------------------------------*
     * Bandwidth Limitation                                      *
     *-----------------------------------------------------------*/

    noiseFillingSize = L_spec;
    if (st->igf)
    {
        noiseFillingSize = (&st->hIGFEnc)->infoStartLine;
    }
    else
    {
        (&st->hIGFEnc)->infoStopLine = noiseFillingSize;
    }

    for (i=(&st->hIGFEnc)->infoStopLine; i < max(L_frame, L_frameTCX); i++)
    {
        spectrum[i] = 0.0f;
    }

    /*-----------------------------------------------------------*
     * Quantization                                              *
     *-----------------------------------------------------------*/

    if (!st->tcx_lpc_shaped_ari)
    {
        /* context based arithmetic coder */

        /* Fast estimation of the scalar quantizer step size */
        if (tcx_cfg->ctx_hm && (st->last_core != 0) )
        {

            LtpPitchLag = ((!tcxonly) && (st->tcxltp_pitch_int < st->L_frame)
                           ? ((2 * st->L_frame * st->pit_res_max) << kLtpHmFractionalResolution) / (st->tcxltp_pitch_int * st->pit_res_max + st->tcxltp_pitch_fr): -1);

            ++ctxHmBits;  /* ContextHM flag */
            --sqTargetBits; /* ContextHM flag */

            OriginalSpectrum = spectrum;

            PeriodicityIndex = SearchPeriodicityIndex( OriginalSpectrum, NULL, L_spec, sqTargetBits, LtpPitchLag,
                               st->tcxltp ? st->tcxltp_gain : -1.0f, &RelativeScore );

            ConfigureContextHm( L_spec, sqTargetBits, PeriodicityIndex, LtpPitchLag, hm_cfg );

            NumIndexBits = CountIndexBits( L_spec >= 256, PeriodicityIndex );

            /* Quantize original spectrum */
            sqGain = SQ_gain( OriginalSpectrum, (int)(LPDmem->tcx_target_bits_fac * (float)sqTargetBits), L_spec );

            tcx_scalar_quantization( OriginalSpectrum, sqQ, L_spec, sqGain, tcx_cfg->sq_rounding, st->memQuantZeros, tcxonly );

            /* Estimate original bitrate */
            stop = 0;
            sqBits = ACcontextMapping_encode2_estimate_no_mem_s17_LC( sqQ, L_spec, &lastnz, &nEncoded, sqTargetBits, &stop, NULL );

            /* Estimate context mapped bitrate */
            stopCtxHm = 0;

            /* Context Mapping */
            sqBitsCtxHm = ACcontextMapping_encode2_estimate_no_mem_s17_LC( sqQ, L_spec, &lastnzCtxHm, &nEncodedCtxHm,
                          sqTargetBits - NumIndexBits, &stopCtxHm, hm_cfg );

            /* Decide whether or not to use context mapping */

            Selector = max(stop, sqBits) - (max(stopCtxHm, sqBitsCtxHm) + NumIndexBits);

            if (Selector > 2 || (abs(Selector) <= 2 && kCtxHmOlRSThr < RelativeScore))
            {
                /* CtxHm is likely better */
                sqTargetBits -= NumIndexBits;
                ctxHmBits  += NumIndexBits;
                prm_hm[0] = 1;
                prm_hm[1] = PeriodicityIndex;
                *prm_lastnz = lastnzCtxHm;
                sqBits_noStop = sqBits = sqBitsCtxHm;
                nEncoded = nEncodedCtxHm;
                stop = stopCtxHm;
            }
            else
            {
                /* Original is better or not much difference */
                prm_hm[0] = 0;
                prm_hm[1] = PeriodicityIndex;
                *prm_lastnz = lastnz;
                PeriodicityIndex = -1;

                sqBits_noStop = sqBits;
            }

            if (stop != 0)
            {
                sqBits = stop;
            }
        }
        else
        {
            /* no context hm*/
            PeriodicityIndex = -1;

            sqGain = SQ_gain(spectrum, (int)(LPDmem->tcx_target_bits_fac * (float)sqTargetBits), L_spec);

            /* Quantize spectrum */
            tcx_scalar_quantization( spectrum, sqQ, L_spec, sqGain, tcx_cfg->sq_rounding, st->memQuantZeros, tcxonly );

            /* Estimate bitrate */
            stop = 0;
            sqBits_noStop = sqBits = ACcontextMapping_encode2_estimate_no_mem_s17_LC( sqQ, L_spec, prm_lastnz, &nEncoded, sqTargetBits, &stop, NULL);

            if(stop!=0)
            {
                sqBits=stop;
            }
        } /* end of if (ctx_hm) */

        /* Adjust correction factor */
        if ((L_spec & (L_spec - 1)) == 0)
        {
            /* power-of-2 */
            LPDmem->tcx_target_bits_fac *= (float)sqTargetBits / (float)(sqBits + 1);
        }
        else
        {
            LPDmem->tcx_target_bits_fac *= (float)sqTargetBits / (float)sqBits;
        }

        if (LPDmem->tcx_target_bits_fac>1.25)
        {
            LPDmem->tcx_target_bits_fac = 1.25;
        }
        if (LPDmem->tcx_target_bits_fac<0.75)
        {
            LPDmem->tcx_target_bits_fac = 0.75;
        }

        /* Refine quantizer step size with a rate-control-loop (optional) */
        sqBits = tcx_scalar_quantization_rateloop( spectrum, sqQ, L_spec, &sqGain, tcx_cfg->sq_rounding, st->memQuantZeros,
                 prm_lastnz, /* lastnz */ sqTargetBits, &nEncoded, &stop, sqBits_noStop, sqBits, tcx_cfg->tcxRateLoopOpt,
                 tcxonly, PeriodicityIndex >= 0 ? hm_cfg : NULL );

        if (ctxHmBits > 0)
        {
            /* Mapping tool is enabled */
            /* Truncate spectrum */
            for (i=nEncoded; i<L_spec; i++)
            {
                sqQ[i] = 0;
            }

            if (PeriodicityIndex >= 0)
            {
                /* Mapping is used */
                /* Estimate non-mapped bitrate */
                stopCtxHm = 1;

                sqBitsCtxHm = ACcontextMapping_encode2_estimate_no_mem_s17_LC(sqQ, L_spec, &lastnz, &nEncodedCtxHm, sqTargetBits, &stopCtxHm, NULL);

                /* Decide whether or not to revert mapping */
                Selector = sqBits - (sqBitsCtxHm + NumIndexBits);

                if (stopCtxHm == 0 && Selector > 0)
                {
                    /* Non-mapped is better */
                    sqTargetBits += NumIndexBits;
                    ctxHmBits  -= NumIndexBits;
                    prm_hm[0] = 0;
                    *prm_lastnz = lastnz;
                    PeriodicityIndex = -1;
                    sqBits_noStop = sqBits = sqBitsCtxHm;
                    nEncoded = nEncodedCtxHm;
                    stop = stopCtxHm;
                }
            }
            else
            {
                /* Mapping is not used */
                /* Estimate mapped bitrate */
                stopCtxHm = 1;

                sqBitsCtxHm = ACcontextMapping_encode2_estimate_no_mem_s17_LC( sqQ, L_spec, &lastnzCtxHm, &nEncodedCtxHm,
                              sqTargetBits - NumIndexBits, &stopCtxHm, hm_cfg );

                /* Decide whether or not to use mapping */
                Selector = sqBits - (sqBitsCtxHm + NumIndexBits);

                if (stopCtxHm == 0 && Selector > 0)
                {
                    /* Mapped is better */
                    sqTargetBits -= NumIndexBits;
                    ctxHmBits  += NumIndexBits;
                    prm_hm[0] = 1;
                    *prm_lastnz = lastnzCtxHm;
                    PeriodicityIndex = prm_hm[1];
                    sqBits_noStop = sqBits = sqBitsCtxHm;
                    nEncoded = nEncodedCtxHm;
                    stop = stopCtxHm;
                }
            }
        }

        /* Limit low sqGain for avoiding saturation of the gain quantizer*/
        if(sqGain < sqrt((float)NORM_MDCT_FACTOR / (float)L_spec))
        {
            sqGain = (float)sqrt((float)NORM_MDCT_FACTOR / (float)L_spec);

            tcx_scalar_quantization( spectrum, sqQ, L_spec, sqGain, tcx_cfg->sq_rounding, st->memQuantZeros, tcxonly );

            stop=1;

            sqBits = ACcontextMapping_encode2_estimate_no_mem_s17_LC( sqQ, L_spec, prm_lastnz, &nEncoded, sqTargetBits, &stop, PeriodicityIndex >= 0 ? hm_cfg : NULL);
        }

        /* Truncate spectrum (for CBR) */
        if ( stop )
        {
            for (i=nEncoded; i<L_spec; i++)
            {
                sqQ[i] = 0;
            }
        }
        /* Save quantized Values */

        for(i=0; i<L_spec; i++)
        {
            spectrum[i] = (float)sqQ[i];
            /* noise filling seed */
            nf_seed += (short)(abs(sqQ[i]) * i * 2);
        }

    }
    else
    {
        /* low rates: envelope based arithmetic coder */

        AdaptLowFreqEmph( spectrum, NULL, 0.f, 1, gainlpc, L_frame );

        prm_target = sqQ;
        sqQ = prm_target + 1;
        signs = hm_cfg->indexBuffer;

        LtpPitchLag = ((st->tcxltp_pitch_int < st->L_frame)
                       ? ((2 * st->L_frame * st->pit_res_max) << kLtpHmFractionalResolution) / (st->tcxltp_pitch_int * st->pit_res_max + st->tcxltp_pitch_fr): -1);

        tcx_arith_encode_envelope( spectrum, signs, L_frame, L_spec, st, Aqind, sqTargetBits, sqQ, st->last_core != ACELP_CORE,
                                   prm_hm, /* HM parameter area */ LtpPitchLag, &sqBits, &signaling_bits ,(st->bwidth > WB)?1:0 );

        sqTargetBits -= signaling_bits;
        *prm_target = sqTargetBits;

        /* Noise filling seed */
        for (i=0; i<noiseFillingSize; ++i)
        {
            nf_seed += (short)(abs((int)spectrum[i]) * i * 2);
        }
    }

    /*-----------------------------------------------------------*
     * Compute optimal TCX gain.                                 *
     *-----------------------------------------------------------*/

    /* initialize LF deemphasis factors in xn_buf */
    for (i = 0; i < L_spec; i++)
    {
        xn_buf[i] = 1.0f;
    }

    if (!tcxonly)
    {
        AdaptLowFreqDeemph( spectrum, st->tcx_lpc_shaped_ari, gainlpc, L_frame,  xn_buf /* LF deemphasis factors */ );
    }

    gain_tcx_opt = get_gain(x_orig, spectrum, L_spec, &ener);

    if (gain_tcx_opt <= 0.0f)
    {
        gain_tcx_opt = sqGain;
    }
    gain_tcx = gain_tcx_opt;

    /*-----------------------------------------------------------*
     * Quantize TCX gain                                         *
     *-----------------------------------------------------------*/

    /*  gain quantization here in case of VBR unvoiced coding; fixes problems of uninitialized global gain values */
    if (st->total_brate >= ACELP_13k20 && !st->rf_mode)
    {
        QuantizeGain(L_spec, &gain_tcx, &prm[0]);
    }

    /*-----------------------------------------------------------*
     * Residual Quantization                                     *
     *-----------------------------------------------------------*/

    if (tcx_cfg->resq)
    {
        resQTargetBits = sqTargetBits-sqBits;

        if (st->tcx_lpc_shaped_ari)
        {
            /* envelope based arithmetic coder */
            int *prm_resq;

            prm_resq = sqQ + sqTargetBits - resQTargetBits;

            resQBits = tcx_ari_res_Q_spec( x_orig, signs, spectrum, L_spec, gain_tcx, prm_resq, resQTargetBits,
                                           resQBits, tcx_cfg->sq_rounding, xn_buf /* LF deemphasis factors */ );

            /* Transmit zeros when there bits remain after RESQ */
            for (i=resQBits; i<resQTargetBits; ++i)
            {
                prm_resq[i] = 0;
            }
        }
        else
        {
            /* context based arithmetic coder */
            resQBits = tcx_res_Q_gain( gain_tcx_opt,&gain_tcx, sqQ+L_spec,resQTargetBits );

            resQBits = tcx_res_Q_spec( x_orig, spectrum, L_spec, gain_tcx, sqQ+L_spec, resQTargetBits,
                                       resQBits, tcx_cfg->sq_rounding, tcxonly ? NULL : xn_buf /* LF deemphasis factors */ );
        }

    }

    /*-----------------------------------------------------------*
     * ALFE tcx only bitrates                                    *
     *-----------------------------------------------------------*/

    if (st->tcxonly)
    {
        if (st->tcxltp && (st->tcxltp_gain > 0.0f) && !fUseTns)
        {
            PsychAdaptLowFreqDeemph( spectrum, gainlpc, NULL );
        }
    }

    /*-----------------------------------------------------------*
     * TCX SNR for Analysis purposes                             *
     *-----------------------------------------------------------*/




    maxNfCalcBw     = min(noiseFillingSize, (int)(st->measuredBwRatio * (float)L_frame + 0.5f));

    /*-----------------------------------------------------------*
     * Estimate and quantize noise factor                        *
     *-----------------------------------------------------------*/

    if (st->total_brate >= HQ_96k)
    {
        fac_ns = 0.0f;
        prm[1] = 0;
    }
    else
    {
        i = L_frame / ((st->total_brate >= ACELP_13k20 && !st->rf_mode) ? 6 : 8); /* noise filling start bin*/

        if (tcxonly)
        {
            noiseTransWidth = HOLE_SIZE_FROM_LTP(max(st->tcxltp_gain,(tcx_cfg->ctx_hm && st->last_core != 0) ? 0.3125f*prm_hm[0] : 0));

            if (L_frame == st->L_frame >> 1)
            {
                noiseTransWidth = 3;  /* minimum transition for noise filling in TCX-10 */
            }
        }

        tcx_noise_factor( x_orig, spectrum, i, maxNfCalcBw, noiseTransWidth, L_frame,
                          gain_tcx, st->noiseTiltFactor, &fac_ns, &prm[NOISE_FILL_RANGES] );

        /* hysteresis for very tonal passages (more stationary noise filling level) */

        if (prm[1] == 1)
        {
            st->noiseLevelMemory = 1 + abs(st->noiseLevelMemory);   /* update counter */
        }
        else
        {
            if ((prm[1] == 2) && (abs(st->noiseLevelMemory) > 5))
            {
                prm[1] = 1;   /* reduce noise filling level by one step */
                fac_ns = 0.75f / (1<<NBITS_NOISE_FILL_LEVEL);

                /* signal that noise level is changed by inverting sign of level memory */
                st->noiseLevelMemory = (st->noiseLevelMemory < 0) ? 5 : -1 - st->noiseLevelMemory;
            }
            else
            {
                st->noiseLevelMemory = 0;  /* reset memory since level is too different */
            }
        }
    } /* bitrate */


    /*-----------------------------------------------------------*
     * Internal TCX decoder                                      *
     *-----------------------------------------------------------*/

    /*-----------------------------------------------------------*
     * Noise Filling.                                            *
     *-----------------------------------------------------------*/

    /* Replication of ACELP formant enhancement for low rates */
    if ( st->total_brate <  ACELP_13k20 || st->rf_mode )
    {
        tcxFormantEnhancement(xn_buf, gainlpc, spectrum, L_frame);
    }

    if (fac_ns > 0.0f)
    {
        i = tcxGetNoiseFillingTilt( A, L_frame, (st->total_brate >= ACELP_13k20 && !st->rf_mode), &st->noiseTiltFactor);

        tcx_noise_filling( spectrum, nf_seed, i, noiseFillingSize, noiseTransWidth, L_frame, st->noiseTiltFactor, fac_ns, NULL );
    }

    if (st->total_brate < ACELP_13k20 || st->rf_mode)
    {
        /* partially recompute global gain (energy part), taking noise filling and formant enhancement into account */
        gain_tcx_opt = 1e-6f;
        for (i = 0; i < L_spec; i++)
        {
            gain_tcx_opt += spectrum[i] * spectrum[i];
        }
        gain_tcx *= (float)sqrt(ener / gain_tcx_opt);
        QuantizeGain(L_spec, &gain_tcx, &prm[0]);
    }

    /*end of noise filling*/

    /*-----------------------------------------------------------*
     * Noise shaping in frequency domain (1/Wz)                  *
     *-----------------------------------------------------------*/

    /* LPC gains already available */
    mdct_noiseShaping( spectrum, L_frame, gainlpc );

    /*-----------------------------------------------------------*
     * Apply gain                                                *
     *-----------------------------------------------------------*/

    if( st->tcx_cfg.coder_type == INACTIVE )
    {

        gain_tcx *= tcx_cfg->na_scale;
    }

    v_multc( spectrum, gain_tcx, spectrum, L_spec);

    stop = tcx_cfg->tcx_last_overlap_mode;   /* backup last TCX overlap mode */


    if ((L_frame == st->L_frame >> 1) && tcxonly)
    {
        int L = L_frame;

        if ((tcx_cfg->fIsTNSAllowed && fUseTns != 0) || (L_spec > L_frame)) L = L_spec;

        tcxInvertWindowGrouping( tcx_cfg, xn_buf, spectrum, L, fUseTns, st->last_core, stop, frame_cnt, 0);
    }

    /*-----------------------------------------------------------*
     * Temporal Noise Shaping Synthesis                          *
     *-----------------------------------------------------------*/

    if (tcx_cfg->fIsTNSAllowed)
    {
        SetTnsConfig(tcx_cfg, L_frame_glob == st->L_frame, (st->last_core == 0) && (frame_cnt == 0));

        /* Apply TNS to get the reconstructed signal */
        if (fUseTns != 0)
        {
            ApplyTnsFilter(tcx_cfg->pCurrentTnsConfig, pTnsData, spectrum, 0);

            if ((L_frame == st->L_frame >> 1) && (tcxonly))
            {
                if ((tcx_cfg->tcx_last_overlap_mode != FULL_OVERLAP) ||
                        ((tcx_cfg->tcx_curr_overlap_mode == FULL_OVERLAP) && (frame_cnt == 0) && (stop == 0))
                   )
                {
                    const int L_win = L_spec >> 1;

                    /* undo rearrangement of LF sub-window lines for TNS synthesis filter */
                    if (L_frame > L_spec)
                    {
                        assert(0);
                    }
                    else
                    {
                        mvr2r(spectrum+8, xn_buf, L_win);
                        mvr2r(xn_buf, spectrum+L_win, 8);
                        mvr2r(xn_buf+8, spectrum+8, L_win-8);
                    }
                }
            }
        }
    }

    /*-----------------------------------------------------------*
    * Compute inverse MDCT of spectrum[].                        *
    *-----------------------------------------------------------*/

    if ((L_frame == st->L_frame >> 1) && (tcxonly))
    {
        if (tcx_cfg->tcx_last_overlap_mode != FULL_OVERLAP)
        {
            /* minimum or half overlap, two transforms, grouping into one window */
            float win[(L_FRAME_PLUS+L_MDCT_OVLP_MAX)/2];
            const int L_win = L_frame >> 1;
            const int L_spec_TCX5 = max(L_frame, L_spec) >> 1;
            const int L_ola = (tcx_cfg->tcx_last_overlap_mode == MIN_OVERLAP) ? tcx_cfg->tcx_mdct_window_min_length : tcx_cfg->tcx_mdct_window_half_length;
            int w;

            set_f( win, 0, (L_FRAME_PLUS+L_MDCT_OVLP_MAX)/2 );
            set_zero(xn_buf, tcx_offset+(L_ola>>1));  /* zero left end of buffer */

            for (w = 0; w < 2; w++)
            {
                if (tcx_cfg->tcx_last_overlap_mode == MIN_OVERLAP)
                {
                    TCX_MDCT_Inverse(spectrum+w*L_spec_TCX5, win, L_ola, L_win-L_ola, L_ola);
                }
                else
                {
                    TCX_MDCT_Inverse(spectrum+w*L_spec_TCX5, win, L_ola, L_win-L_ola, L_ola);
                }

                tcx_windowing_synthesis_current_frame( win, tcx_cfg->tcx_aldo_window_2, tcx_cfg->tcx_mdct_window_half, tcx_cfg->tcx_mdct_window_minimum,
                                                       L_ola, tcx_cfg->tcx_mdct_window_half_length, tcx_cfg->tcx_mdct_window_min_length, w==0 && st->last_core==0,
                                                       (w > 0) || (w == 0 && stop == 2) ? MIN_OVERLAP : tcx_cfg->tcx_last_overlap_mode, LPDmem->acelp_zir,
                                                       st->LPDmem.Txnq, NULL, Aq_old, tcx_cfg->tcx_mdct_window_trans, L_win,
                                                       tcx_offset<0?-tcx_offset:0,
                                                       (w > 0) ? 1 : st->last_core, 0, 0 );

                if (w > 0)
                {
                    tcx_windowing_synthesis_past_frame( xn_buf+tcx_offset-(L_ola>>1)+w*L_win, tcx_cfg->tcx_aldo_window_1_trunc, tcx_cfg->tcx_mdct_window_half,
                                                        tcx_cfg->tcx_mdct_window_minimum, L_ola, tcx_cfg->tcx_mdct_window_half_length, tcx_cfg->tcx_mdct_window_min_length, 2 );
                }

                /* add part of current sub-window overlapping with previous window */
                v_add(win, xn_buf+tcx_offset-(L_ola>>1)+w*L_win, xn_buf+tcx_offset-(L_ola>>1)+w*L_win, L_ola);

                /* copy new sub-window region not overlapping with previous window */
                mvr2r(win+L_ola, xn_buf+tcx_offset+(L_ola>>1)+w*L_win, L_win);
            }

            /* To assure that no garbage values are copied to LPDmem->Txnq */
            set_zero(xn_buf+L_frame+tcx_offset+(L_ola>>1), overlap-tcx_offset-(L_ola>>1));
        }
        else if ((frame_cnt == 0) && (tcx_cfg->tcx_curr_overlap_mode == FULL_OVERLAP))
        {
            /* special overlap attempt, two transforms, grouping into one window */
            float win[(L_FRAME_PLUS+L_MDCT_OVLP_MAX)/2];
            const int L_win = L_frame >> 1;
            const int L_spec_TCX5 = max(L_frame, L_spec) >> 1;
            const int L_ola = tcx_cfg->tcx_mdct_window_min_length;
            int w;

            set_f( win, 0, (L_FRAME_PLUS+L_MDCT_OVLP_MAX)/2 );

            /* Resize overlap (affect only asymmetric window)*/
            overlap = st->tcx_cfg.tcx_mdct_window_delay;

            /* 1st TCX-5 window, special MDCT with minimum overlap on right side */
            TCX_MDCT_Inverse(spectrum, win+L_win, 0, L_win-(L_ola>>1), L_ola);

            /* copy new sub-window region not overlapping with previous window */
            mvr2r(win+L_win, xn_buf+(overlap>>1), L_win+(L_ola>>1));

            /* 2nd TCX-5 window, regular MDCT with minimum overlap on both sides */
            TCX_MDCT_Inverse(spectrum+L_spec_TCX5, win, L_ola, L_win-L_ola, L_ola);

            tcx_windowing_synthesis_current_frame( win, tcx_cfg->tcx_aldo_window_2, tcx_cfg->tcx_mdct_window_half, tcx_cfg->tcx_mdct_window_minimum,
                                                   L_ola, tcx_cfg->tcx_mdct_window_half_length, tcx_cfg->tcx_mdct_window_min_length, 0,  /* left_rect */
                                                   2,  /* left_mode */ LPDmem->acelp_zir, st->LPDmem.Txnq, NULL, Aq_old,
                                                   tcx_cfg->tcx_mdct_window_trans, L_win, tcx_offset<0?-tcx_offset:0,
                                                   1, /* not st->last_core */ 0 ,0 );

            tcx_windowing_synthesis_past_frame( xn_buf+(overlap>>1)+L_win-(L_ola>>1), tcx_cfg->tcx_aldo_window_1_trunc, tcx_cfg->tcx_mdct_window_half,
                                                tcx_cfg->tcx_mdct_window_minimum, L_ola, tcx_cfg->tcx_mdct_window_half_length, tcx_cfg->tcx_mdct_window_min_length, 2 );

            /* add part of current sub-window overlapping with previous window */
            v_add(win, xn_buf+(overlap>>1)+L_win-(L_ola>>1), xn_buf+(overlap>>1)+L_win-(L_ola>>1), L_ola);

            /* copy new sub-window region not overlapping with previous window */
            mvr2r(win+L_ola, xn_buf+(overlap>>1)+L_win+(L_ola>>1), L_win);

            /* extra folding-out on left side of win, for perfect reconstruction */
            for (w = (overlap>>1); w < overlap; w++)
            {
                xn_buf[overlap-1-w] = -1.0f * xn_buf[w];
            }

            tcx_windowing_synthesis_current_frame( xn_buf, tcx_cfg->tcx_aldo_window_2, tcx_cfg->tcx_mdct_window_half, tcx_cfg->tcx_mdct_window_minimum,
                                                   overlap, /*tcx_cfg->tcx_mdct_window_length*/ tcx_cfg->tcx_mdct_window_half_length, tcx_cfg->tcx_mdct_window_min_length,
                                                   st->last_core == ACELP_CORE, 0,  /* left_mode */ LPDmem->acelp_zir, st->LPDmem.Txnq, NULL,
                                                   Aq_old, tcx_cfg->tcx_mdct_window_trans, L_win, tcx_offset<0?-tcx_offset:0,
                                                   st->last_core, 0, 0 );
        }
        else
        {
            /* default, i.e. maximum overlap, single transform, no grouping */
            TCX_MDCT_Inverse(spectrum, xn_buf, overlap, L_frame-overlap, overlap);

            tcx_windowing_synthesis_current_frame( xn_buf, tcx_cfg->tcx_aldo_window_2, tcx_cfg->tcx_mdct_window_half, tcx_cfg->tcx_mdct_window_minimum,
                                                   overlap, /*tcx_cfg->tcx_mdct_window_length*/ tcx_cfg->tcx_mdct_window_half_length, tcx_cfg->tcx_mdct_window_min_length,
                                                   st->last_core == ACELP_CORE, (frame_cnt > 0) && (stop == 0) && (st->last_core!=0) ? 2 : stop, LPDmem->acelp_zir,
                                                   st->LPDmem.Txnq, NULL, Aq_old, tcx_cfg->tcx_mdct_window_trans, L_frame_glob >> 1,
                                                   tcx_offset<0?-tcx_offset:0,
                                                   st->last_core, 0, 0 );

        } /* tcx_last_overlap_mode != FULL_OVERLAP */
    }
    else
    {
        /* frame is TCX-20 or not TCX-only */

        if (st->tcx_cfg.tcx_last_overlap_mode != TRANSITION_OVERLAP)
        {
            float tmp[L_FRAME_PLUS];

            edct(spectrum, xn_buf+overlap/2+nz, L_frame);
            v_multc( xn_buf+overlap/2+nz, (float)sqrt((float)L_frame / NORM_MDCT_FACTOR), tmp, L_frame);

            window_ola( tmp, xn_buf, st->old_out, L_frame, tcx_cfg->tcx_last_overlap_mode, tcx_cfg->tcx_curr_overlap_mode, 0, 0, NULL );
            aldo = 1;
        }
        else
        {
            TCX_MDCT_Inverse(spectrum, xn_buf, overlap, L_frame-overlap, overlap);

            /*-----------------------------------------------------------*
             * Windowing, overlap and add                                *
             *-----------------------------------------------------------*/

            /* Window current frame */
            tcx_windowing_synthesis_current_frame( xn_buf, tcx_cfg->tcx_aldo_window_2, tcx_cfg->tcx_mdct_window_half, tcx_cfg->tcx_mdct_window_minimum,
                                                   overlap, /*tcx_cfg->tcx_mdct_window_length*/ tcx_cfg->tcx_mdct_window_half_length, tcx_cfg->tcx_mdct_window_min_length,
                                                   st->last_core == ACELP_CORE, tcx_cfg->tcx_last_overlap_mode, /*left mode*/ LPDmem->acelp_zir, st->LPDmem.Txnq,
                                                   NULL, Aq_old, tcx_cfg->tcx_mdct_window_trans, L_frame_glob >> 1, tcx_offset<0?-tcx_offset:0,
                                                   st->last_core, 0, 0 );
        }
    } /* TCX-20/TCX-10 and TCX-only */

    /* Window and overlap-add past frame if past frame is TCX */
    if ((st->last_core > 0) && (((L_frameTCX == st->L_frameTCX >> 1) && (st->tcxonly)) || (st->tcx_cfg.tcx_last_overlap_mode == TRANSITION_OVERLAP)) )
    {
        if (tcx_cfg->last_aldo)
        {
            for (i=0; i < overlap - tcx_cfg->tcx_mdct_window_min_length; i++)
            {
                xn_buf[i] += st->old_out[i+nz];
            }
            /* fade truncated ALDO window */
            for ( ; i < overlap; i++)
            {
                xn_buf[i] += st->old_out[i+nz] * tcx_cfg->tcx_mdct_window_minimum[overlap-1-i];
            }
        }
        else
        {
            if ((frame_cnt > 0) && (stop == 0) && (tcx_cfg->tcx_curr_overlap_mode == FULL_OVERLAP)&& (st->last_core!=0))
            {
                stop = 2;     /* use minimum overlap between the two TCX-10 windows */
            }

            tcx_windowing_synthesis_past_frame( LPDmem->Txnq, tcx_cfg->tcx_aldo_window_1_trunc, tcx_cfg->tcx_mdct_window_half,
                                                tcx_cfg->tcx_mdct_window_minimum, overlap, tcx_cfg->tcx_mdct_window_half_length, tcx_cfg->tcx_mdct_window_min_length,
                                                (stop == 0 || tcx_cfg->tcx_last_overlap_mode == MIN_OVERLAP) ? tcx_cfg->tcx_last_overlap_mode : stop );

            for (i=0; i<overlap; i++)
            {
                xn_buf[i] += LPDmem->Txnq[i];
            }
        }
    }

    if(!aldo && (((L_frameTCX == st->L_frameTCX >> 1)&&frame_cnt > 0) || L_frameTCX != (st->L_frameTCX >> 1)))
    {
        /*Compute windowed synthesis in case of switching to ALDO windows in next frame*/
        mvr2r(xn_buf+L_frame-nz, st->old_out, nz+overlap);
        set_zero(st->old_out+nz+overlap, nz);
        tcx_windowing_synthesis_past_frame( st->old_out+nz, tcx_cfg->tcx_aldo_window_1_trunc, tcx_cfg->tcx_mdct_window_half, tcx_cfg->tcx_mdct_window_minimum,
                                            overlap, tcx_cfg->tcx_mdct_window_half_length, tcx_cfg->tcx_mdct_window_min_length, tcx_cfg->tcx_curr_overlap_mode );

        if(tcx_cfg->tcx_curr_overlap_mode==FULL_OVERLAP)
        {
            for (i=0; i<nz; i++)
            {
                st->old_out[nz+overlap+i]=xn_buf[L_frame-1-i]*tcx_cfg->tcx_aldo_window_1_trunc[-1-i];
            }
            tcx_cfg->tcx_curr_overlap_mode=ALDO_WINDOW;
        }
    }
    tcx_cfg->last_aldo=aldo;

    /* Update Txnq */
    if (!tcx_cfg->last_aldo)
    {
        mvr2r(xn_buf+L_frame, LPDmem->Txnq, overlap);
    }

    /* Output */
    mvr2r( xn_buf+(overlap>>1)-tcx_offset, synth, L_frame_glob );

    /* Update L_frame_past */
    st->L_frame_past = L_frame;

    return;
}


/*-------------------------------------------------------------------*
* coder_tcx()
*
*
*-------------------------------------------------------------------*/

void coder_tcx(
    int n,
    TCX_config *tcx_cfg,  /*input: configuration of TCX*/
    float A[],         /* input: quantized coefficients NxAz_q[M+1] */
    Word16 Aqind[],    /* input: frame-independent quantized coefficients (M+1) */
    float synth[],
    int L_frame_glob,  /* input: frame length             */
    int L_frameTCX_glob,
    int L_spec,
    int nb_bits,       /*input: bit budget*/
    int tcxonly,       /*input: only TCX flag*/
    float spectrum[],  /* i/o: MDCT spectrum */
    LPD_state *LPDmem, /*i/o: memories*/
    int prm[],         /* output: tcx parameters          */
    Encoder_State *st,
    CONTEXT_HM_CONFIG *hm_cfg
)
{
    int L_frame;
    int left_overlap=-1, right_overlap=-1;
    int tnsSize = 0; /* number of tns parameters put into prm */
    int tnsBits = 0; /* number of tns bits in the frame */
    int ltpBits = 0;
    float gainlpc[FDNS_NPTS];
    float buf[N_MAX+L_MDCT_OVLP_MAX];
    float winMDST[N_MAX+L_MDCT_OVLP_MAX];
    float * win;
    float * powerSpec;

    powerSpec = win = buf; /* Share memory for windowed TD signal and for the power spectrum */

    L_frame = L_frameTCX_glob;

    /*-----------------------------------------------------------*
     * Windowing                                                 *
     *-----------------------------------------------------------*/
    if (st->tcx_cfg.tcx_last_overlap_mode == TRANSITION_OVERLAP)
    {
        WindowSignal( tcx_cfg, tcx_cfg->tcx_offsetFB, tcx_cfg->tcx_last_overlap_mode, tcx_cfg->tcx_curr_overlap_mode,
                      &left_overlap, &right_overlap, st->speech_TCX, &L_frame, win, 1 );

        /*-----------------------------------------------------------*
         * Compute MDCT for xn_buf[].                                *
         *-----------------------------------------------------------*/

        TCX_MDCT( win, spectrum, left_overlap, L_frame-(left_overlap+right_overlap)/2, right_overlap );

    }
    else
    {
        wtda( st->new_speech_TCX, win, NULL, tcx_cfg->tcx_last_overlap_mode, tcx_cfg->tcx_curr_overlap_mode, L_frame);

        WindowSignal( tcx_cfg, tcx_cfg->tcx_offsetFB, tcx_cfg->tcx_last_overlap_mode == ALDO_WINDOW ? FULL_OVERLAP : tcx_cfg->tcx_last_overlap_mode,
                      tcx_cfg->tcx_curr_overlap_mode == ALDO_WINDOW ? FULL_OVERLAP : tcx_cfg->tcx_curr_overlap_mode,
                      &left_overlap, &right_overlap, st->speech_TCX, &L_frame, winMDST, 1 );


        edct( win, spectrum, L_frame );

        v_multc(spectrum, (float)sqrt((float)NORM_MDCT_FACTOR / L_frame), spectrum, L_frame);
    }

    /*-----------------------------------------------------------*
     * Attenuate upper end of NB spectrum,                       *
     * to simulate ACELP behavior                                *
     *-----------------------------------------------------------*/

    if( st->narrowBand )
    {
        attenuateNbSpectrum(L_frame, spectrum);
    }

    /*-----------------------------------------------------------*
     * Compute noise-measure flags for spectrum filling          *
     * and quantization (0: tonal, 1: noise-like).               *
     * Detect low pass if present.                               *
     *-----------------------------------------------------------*/

    AnalyzePowerSpectrum( st, L_frame*st->L_frame/st->L_frameTCX, L_frame, left_overlap, right_overlap, spectrum,
                          (st->tcx_cfg.tcx_last_overlap_mode == TRANSITION_OVERLAP) ? win : winMDST, powerSpec );

    if (tcx_cfg->fIsTNSAllowed)
    {
        SetTnsConfig(tcx_cfg, L_frame_glob == st->L_frame, st->last_core == 0);
        TNSAnalysis(tcx_cfg, L_frame, L_spec, TCX_20, st->last_core == 0, spectrum, st->tnsData, st->fUseTns, &(&st->hIGFEnc)->tns_predictionGain );

    }
    else
    {
        st->fUseTns[0] = st->fUseTns[1] = 0;
    }

    if(st->igf)
    {
        ProcessIGF(&st->hIGFEnc, st, spectrum, powerSpec, 1, st->fUseTns[0], (st->last_core == ACELP_CORE), 0);
    }

    ShapeSpectrum( tcx_cfg, A, gainlpc, L_frame_glob, L_spec, spectrum, st->fUseTns[0], st );

    if(st->igf)
    {
        nb_bits -= st->hIGFEnc.infoTotalBitsPerFrameWritten;
    }

    if (tcx_cfg->fIsTNSAllowed)
    {
        EncodeTnsData(tcx_cfg->pCurrentTnsConfig, st->tnsData, prm+1+NOISE_FILL_RANGES+LTPSIZE, &tnsSize, &tnsBits);
    }

    QuantizeSpectrum( tcx_cfg, A, Aqind, gainlpc, synth, L_frame_glob, L_frameTCX_glob, L_spec,
                      nb_bits-tnsBits-ltpBits, tcxonly, spectrum,st->tnsData, st->fUseTns[0],
                      tnsSize, LPDmem, prm, n, st, hm_cfg );

    return;
}

/*-------------------------------------------------------------------*
* coder_tcx_post()
*
*
*-------------------------------------------------------------------*/

void coder_tcx_post(
    Encoder_State *st,
    LPD_state *LPDmem,
    TCX_config *tcx_cfg,
    float *synth,
    float *A,
    const float *Ai,
    float *wsig )
{
    float xn_buf[L_FRAME_MAX];

    /* TCX output */
    mvr2r( synth, xn_buf, st->L_frame );


    /*-----------------------------------------------------------*
     * Memory update                                             *
     *-----------------------------------------------------------*/

    /* Update LPDmem (Txnq,syn,syn_pe,old_exc,wsyn,Ai,Aq) */
    tcx_encoder_memory_update( wsig, xn_buf, st->L_frame, Ai, A, tcx_cfg->preemph_fac, LPDmem, st, M, synth );

    return;
}
