/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "prot.h"
#include "options.h"
#include "stat_com.h"
#include "cnst.h"


/*-----------------------------------------------------------------*
 * Local functions
 *-----------------------------------------------------------------*/

static void IMDCT( float *x, float *old_syn_overl, float *syn_Overl_TDAC, float *xn_buf, float *tcx_aldo_window_1_trunc, float *tcx_aldo_window_2,
                   float *tcx_mdct_window_half, float *tcx_mdct_window_minimum, float *tcx_mdct_window_trans, int tcx_mdct_window_half_length,
                   int tcx_mdct_window_min_length, int index, int left_rect, int tcx_offset, int overlap, int L_frame, int L_frameTCX, int L_spec_TCX5,
                   int L_frame_glob,
                   int frame_cnt, int bfi, float *old_out, short FB_flag, Decoder_State *st, int fullband, float *acelp_zir );

/*-----------------------------------------------------------------*
 * decoder_tcx()
 *
 *
 *-----------------------------------------------------------------*/

void decoder_tcx(
    TCX_config *tcx_cfg,  /* input: configuration of TCX         */
    int prm[],            /* input:  parameters                  */
    float A[],            /* input:  coefficients NxAz[M+1]      */
    Word16 Aind[],        /* input: frame-independent coefficients Az[M+1] */
    int L_frame_glob,     /* input:  frame length                */
    int L_frameTCX_glob,
    int L_spec,
    float synth[],        /* in/out: synth[-M..L_frame]          */
    float synthFB[],
    Decoder_State *st,    /* in/out: coder memory state          */
    const short coder_type,     /* input: coder type                   */
    int bfi,              /* input:  Bad frame indicator         */
    int frame_cnt,        /* input: frame counter in the super frame */
    float stab_fac        /* input: stability of isf             */
)
{
    short i, index, L_frame, tcx_offset;
    short L_frameTCX, tcx_offsetFB;
    short firstLine;
    float gain_tcx, fac_ns;
    float Ap[M+2];
    float x[N_MAX];
    int overlap;
    short overlapFB;
    short noiseFillingSize;
    short noiseTransWidth = MIN_NOISE_FILLING_HOLE;
    int tnsSize = 0; /* number of tns parameters put into prm   */
    int fUseTns = 0; /* flag that is set if TNS data is present */
    STnsData tnsData;
    short left_rect;
    float gainlpc2[FDNS_NPTS];
    float gamma1;
    float gamma;
    short nf_seed = 0;
    float gainCompensate = 1.f;
    float h1[L_FRAME_MAX/4+1];
    float mem[M];
    float tmp2;
    int arith_bits, signaling_bits;
    const int *prm_ltp, *prm_tns, *prm_hm, *prm_sqQ, *prm_target;
    float xn_buf[L_MDCT_OVLP_MAX+L_FRAME_PLUS+L_MDCT_OVLP_MAX];
    float xn_bufFB[L_MDCT_OVLP_MAX+L_FRAME_PLUS+L_MDCT_OVLP_MAX];
    float acelp_zir[L_FRAME_MAX/2];
    int sum_word32;
    short temp_concealment_method = 0;      /* to avoid compilation warnings */
    int infoIGFStartLine;

    prm_target = (const int *)NULL;       /* just to suppress MSVC warnigs */


    /*-----------------------------------------------------------------*
     * Initializations
     *-----------------------------------------------------------------*/

    /* Init lengths */
    overlap    = tcx_cfg->tcx_mdct_window_length;
    overlapFB  = tcx_cfg->tcx_mdct_window_lengthFB;
    tcx_offset = tcx_cfg->tcx_offset;
    tcx_offsetFB = tcx_cfg->tcx_offsetFB;
    gamma1     = st->gamma;

    if (st->enableTcxLpc)
    {
        gamma1 = 1.0f;
    }

    if( bfi )
    {
        /* PLC: [TCX: Memory update]
         * PLC: Init buffers */

        assert(st->L_frame_past > 0);
        L_frame    = st->L_frame_past;
        L_frameTCX = st->L_frameTCX_past;

        left_rect = st->prev_widow_left_rect;

        if( left_rect )
        {
            tcx_offset = tcx_cfg->lfacNext;
            tcx_offsetFB = tcx_cfg->lfacNextFB;
            L_spec += st->tcx_cfg.tcx_coded_lines >> 2;
        }
    }
    else
    {
        if ( frame_cnt == 0 && st->last_core == ACELP_CORE )
        {
            if (!st->prev_bfi)
            {
                tcx_cfg->last_aldo = 0;
            }

            /* if past frame is ACELP */
            L_frame = L_frame_glob + tcx_offset;
            L_frameTCX = L_frameTCX_glob + tcx_offsetFB;
            L_spec += st->tcx_cfg.tcx_coded_lines >> 2;

            assert(tcx_cfg->lfacNext<=0);
            L_frame     -= tcx_cfg->lfacNext;
            L_frameTCX  -= tcx_cfg->lfacNextFB;
            tcx_offset   = tcx_cfg->lfacNext;
            tcx_offsetFB = tcx_cfg->lfacNextFB;

            left_rect = 1;
            st->prev_widow_left_rect = 1;
        }
        else
        {
            L_frame = L_frame_glob;
            L_frameTCX = L_frameTCX_glob;
            left_rect = 0;
            st->prev_widow_left_rect = 0;
        }

        st->L_frame_past = L_frame;
        st->L_frameTCX_past = L_frameTCX;
    }

    if( (L_frame == st->L_frame >> 1) && (st->tcxonly) )
    {
        IGFDecUpdateInfo( &st->hIGFDec, IGF_GRID_LB_SHORT );
    }
    else
    {
        IGFDecUpdateInfo( &st->hIGFDec, (st->last_core == ACELP_CORE || (left_rect && st->bfi))?IGF_GRID_LB_TRAN:IGF_GRID_LB_NORM );
    }

    if( st->igf == 0 )
    {
        if( st->narrowBand == 0 )
        {
            /* minimum needed for output with sampling rates lower then the
               nominal sampling rate */
            infoIGFStartLine = min(L_frameTCX, L_frame);
        }
        else
        {
            infoIGFStartLine  = L_frameTCX;
        }
    }
    else
    {
        infoIGFStartLine = min(st->hIGFDec.infoIGFStartLine,L_frameTCX);
    }

    noiseFillingSize = L_spec;
    if( st->igf )
    {
        noiseFillingSize = st->hIGFDec.infoIGFStartLine;
    }

    prm_ltp = &prm[1+NOISE_FILL_RANGES];
    prm_tns = prm_ltp + LTPSIZE;

    /*-----------------------------------------------------------*
     * Read TCX parameters                                       *
     *-----------------------------------------------------------*/

    index = 0;

    if( !bfi )
    {
        index = prm[0];

        /* read noise level (fac_ns) */
        st->noise_filling_index = prm[1];
    }

    fac_ns = (float)st->noise_filling_index * 0.75f / (1<<NBITS_NOISE_FILL_LEVEL);

    /* read TNS data */
    if( !bfi && tcx_cfg->fIsTNSAllowed )
    {
        fUseTns = DecodeTnsData( tcx_cfg->pCurrentTnsConfig, prm_tns, &tnsSize, &tnsData );
    }
    else
    {
        fUseTns = 0;
    }

    prm_hm = prm_tns + tnsSize;
    prm_sqQ = prm_hm + NPRM_CTX_HM;

    /*-----------------------------------------------------------*
     * Spectrum data                                             *
     *-----------------------------------------------------------*/

    if( !bfi)
    {
        /*-----------------------------------------------------------*
         * Context HM                                                *
         *-----------------------------------------------------------*/

        if(tcx_cfg->ctx_hm && ( (st->last_core != ACELP_CORE) || (frame_cnt > 0) ) )
        {
            st->last_ctx_hm_enabled = prm_hm[0];
            {
                for (i = 0; i < L_spec; i++)
                {
                    /* no context harmonic model, copy MDCT coefficients to x */
                    x[i] = (float)prm_sqQ[i];
                }
            }
        }
        else  /* tcx_cfg->ctx_hm == 0 */
        {
            if( st->tcx_lpc_shaped_ari )  /* low rates: envelope based arithmetic coder */
            {
                prm_target = prm_sqQ;
                prm_sqQ = prm_target + 1;

                tcx_arith_decode_envelope( x, L_frame, L_spec, st, coder_type, Aind, st->tcxltp_gain, *prm_target, prm_sqQ, st->last_core != ACELP_CORE,
                                           prm_hm, /* HM parameter area */ st->tcx_hm_LtpPitchLag, &arith_bits, &signaling_bits ,(st->bwidth > WB)?1:0 );

                st->resQBits[frame_cnt] = *prm_target - arith_bits;

                /* Noise filling seed */
                for (i=0; i<noiseFillingSize; ++i)
                {
                    nf_seed += (short)(abs((int)x[i]) * i * 2);
                }
            }
            else  /* TCX-only: context based arithmetic coder */
            {
                for (i = 0; i < L_spec; i++)
                {
                    x[i] = (float)prm_sqQ[i];
                }

                for (i=L_spec ; i < L_frameTCX; i++)
                {
                    x[i] = 0.0f;
                }
            }

        }  /* else of if tcx_cfg->ctx_hm */

        for( i=L_spec ; i < max(L_frame, L_frameTCX); i++ )
        {
            x[i] = 0.0f;
        }

        /*-----------------------------------------------------------*
         * adaptive low frequency deemphasis.                        *
         *-----------------------------------------------------------*/

        weight_a( A, Ap, gamma1, M );
        lpc2mdct( Ap, M, gainlpc2 );

        /* initialize LF deemphasis factors in xn_buf */
        for( i = 0; i < max(L_spec, L_frameTCX); i++ )
        {
            xn_buf[i] = 1.0f;
        }

        if( !st->tcxonly )
        {
            AdaptLowFreqDeemph( x, st->tcx_lpc_shaped_ari, gainlpc2, L_frame, xn_buf /* LF deemphasis factors */ );
        }
    }

    st->damping = 0.f;

    if( bfi == 0 )
    {
        /*-----------------------------------------------------------*
         * Compute global gain                                       *
         *-----------------------------------------------------------*/

        gain_tcx = (float)pow(10.0f, index/28.0f) * (float)sqrt((float)NORM_MDCT_FACTOR / (float)L_spec);

        st->old_gaintcx_bfi = gain_tcx;

        st->cummulative_damping_tcx = 1.0f;
    }
    else /* bfi = 1 */
    {
        /* PLC: [TCX: Fade-out]
         * derivation of damping factor */
        if( st->use_partial_copy )
        {
            if( st->rf_frame_type == RF_TCXFD )
            {
                gain_tcx = (float)pow(10.0f, (int)st->old_gaintcx_bfi/28.0f) * (float)sqrt((float)NORM_MDCT_FACTOR / (float)L_spec);
                st->old_gaintcx_bfi = gain_tcx;
            }
            else
            {
                gain_tcx = st->old_gaintcx_bfi;
            }

            st->damping = 1;
        }
        else
        {
            gain_tcx = st->old_gaintcx_bfi;
            st->damping = Damping_fact( coder_type, st->nbLostCmpt, st->last_good, stab_fac, &(st->lp_gainp), st->last_core );
        }

        st->cummulative_damping_tcx *= st->damping;
    }

    if( bfi )
    {
        if( bfi && st->envWeighted )
        {
            gamma = st->gamma;
        }
        else
        {
            gamma = gamma1;
        }

        /* PLC: [TCX: Fade-out]
         * PLC: invert LPC weighting in case of PLC */
        if( bfi )
        {
            if (st->enableTcxLpc)
            {
                gamma = st->cummulative_damping_tcx * (st->gamma - 1) + 1;
            }
            else
            {
                gamma = st->cummulative_damping_tcx * (gamma1 - 1) + 1;
            }
        }

        weight_a( A, Ap, gamma, M );
        lpc2mdct( Ap, M, gainlpc2 );
    }

    tmp2 = 0;
    set_zero( h1, L_SUBFR+1 );
    set_zero( mem, M );
    h1[0] = 1.0f;
    syn_filt( Ap, M,h1, h1, L_SUBFR, mem, 0 ); /* impulse response of LPC     */
    deemph( h1, st->preemph_fac, L_SUBFR, &tmp2 );           /* impulse response of deemph  */

    /* impulse response level = gain introduced by synthesis+deemphasis */
    if( !bfi )
    {
        st->last_gain_syn_deemph = (float)sqrt(dotp( h1, h1, L_SUBFR) );
        /*for avoiding compiler warnings*/
        st->gainHelper = 1.f;
        st->stepCompensate=0.f;
    }
    else if (TCX_20_CORE == st->core || 1 == frame_cnt )
    {
        gainCompensate = st->last_gain_syn_deemph/(float)sqrt(dotp( h1, h1, L_SUBFR) );
        if (st->nbLostCmpt==1)
        {
            st->stepCompensate = (1.f - gainCompensate)/st->L_frame;
            st->gainHelper = 1.f;
        }
        else
        {
            st->stepCompensate = (st->last_concealed_gain_syn_deemph - gainCompensate)/st->L_frame;
            st->gainHelper = st->last_concealed_gain_syn_deemph;
        }
        st->last_concealed_gain_syn_deemph = gainCompensate;
    }

    /*-----------------------------------------------------------*
     * Residual inv. Q.                                          *
     *-----------------------------------------------------------*/

    if (!bfi && tcx_cfg->resq )
    {
        if( st->tcx_lpc_shaped_ari )
        {
            /* envelope based arithmetic coder */
            const int *prm_resq;
            prm_resq = prm_sqQ + *prm_target /* = targetBits */ - st->resQBits[frame_cnt];
            i = tcx_ari_res_invQ_spec( x, L_spec, prm_resq, st->resQBits[frame_cnt], 0, tcx_cfg->sq_rounding, xn_buf /* LF deemphasis factors */ );
        }
        else
        {
            /* context based arithmetic coder */
            i = tcx_res_invQ_gain(&gain_tcx, &prm_sqQ[L_spec], st->resQBits[frame_cnt] );
            tcx_res_invQ_spec( x, L_spec, &prm_sqQ[L_spec], st->resQBits[frame_cnt], i, tcx_cfg->sq_rounding, st->tcxonly ? NULL : xn_buf /* LF deemphasis factors */ );
        }
    }

    if( !bfi && st->tcxonly )
    {
        if( st->tcxltp && (st->tcxltp_gain > 0.0f) && !fUseTns )
        {
            PsychAdaptLowFreqDeemph(x, gainlpc2, NULL);
        }
    }

    if( !bfi && !st->tcxonly )
    {
        /* Replication of ACELP formant enhancement for low rates */
        if (st->bits_frame < 256)
        {
            tcxFormantEnhancement( xn_buf, gainlpc2, x, L_frame );
        }
    }

    /*-----------------------------------------------------------*
     * Add gain to the lpc gains                                 *
     *-----------------------------------------------------------*/

    if( st->VAD == 0 )
    {
        gain_tcx *= tcx_cfg->na_scale;
    }

    v_multc( gainlpc2, gain_tcx, gainlpc2, FDNS_NPTS);


    /*-----------------------------------------------------------*
     * Noise filling.                                            *
     *-----------------------------------------------------------*/

    if( !bfi && (fac_ns > 0.0f) )
    {
        float noiseTiltFactor;

        firstLine = tcxGetNoiseFillingTilt( A, L_frame, (st->bits_frame >= 256 && !st->rf_flag), &noiseTiltFactor );

        if( st->tcxonly )
        {
            noiseTransWidth = HOLE_SIZE_FROM_LTP(max(st->tcxltp_gain,(tcx_cfg->ctx_hm && st->last_core != 0) ? 0.3125f*st->last_ctx_hm_enabled : 0));

            if (L_frame == st->L_frame >> 1)
            {
                noiseTransWidth = 3;  /* minimum transition fading for noise filling in TCX-10 */
            }
        }

        if (!st->tcx_lpc_shaped_ari)
        {
            /* context based arithmetic coder */
            /* noise filling seed */
            for (i = 0; i < L_spec; i++)
            {
                nf_seed += (short)(abs(prm_sqQ[i]) * i * 2);
            }
        }


        tcx_noise_filling( x, nf_seed, firstLine, noiseFillingSize, noiseTransWidth,
                           L_frame, noiseTiltFactor, fac_ns, (st->igf)?st->hIGFDec.infoTCXNoise:NULL );
        st->seed_tcx_plc = nf_seed;
    }

    if( st->enablePlcWaveadjust )
    {
        if( bfi )
        {
            if( st->nbLostCmpt == 1 )
            {
                st->plcInfo.concealment_method = TCX_NONTONAL;
                /* tonal/non-tonal decision */
                if (st->plcInfo.Transient[0] == 1 && st->plcInfo.Transient[1] == 1 && st->plcInfo.Transient[2] == 1)
                {
                    sum_word32 = 0;

                    for (i = 9; i >= 0; i--)
                    {
                        sum_word32 += st->plcInfo.TCX_Tonality[i];
                    }

                    if (sum_word32 >= 6)
                    {
                        st->plcInfo.concealment_method = TCX_TONAL;
                    }
                }

                if (st->tonal_mdct_plc_active)
                {
                    st->plcInfo.concealment_method = TCX_TONAL;
                }
            }

            if (L_frameTCX > st->L_frameTCX )
            {
                st->plcInfo.concealment_method = TCX_TONAL;
            }

            temp_concealment_method = st->plcInfo.concealment_method;

            if (st->core == TCX_10_CORE)
            {
                temp_concealment_method = TCX_TONAL;
            }
        }

        /* get the starting location of the subframe in the frame */
        if (st->core == TCX_10_CORE)
        {
            st->plcInfo.subframe = frame_cnt*L_frameTCX_glob;
        }
    }

    /* PLC: [TCX: Tonal Concealment] */
    /* PLC: [TCX: Fade-out]
     * PLC: Fade out to white noise */
    if( !bfi )
    {
        TonalMDCTConceal_SaveFreqSignal( &st->tonalMDCTconceal, x, L_frameTCX, L_frame, gainlpc2 );
    }
    else
    {
        if( !st->enablePlcWaveadjust || (temp_concealment_method == TCX_TONAL))
        {
            /* set f to 1 to not fade out */
            /* set f to 0 to immediately switch to white noise */
            float f;
            float noiseTiltFactor;

            if (st->tcxonly)
            {
                f = 1.0f;
            }
            else
            {
                f = st->cummulative_damping_tcx;
            }

            if( (frame_cnt == 0) && (L_frameTCX == st->L_frameTCX >> 1)
                    && (st->tcxonly) && (!st->tonal_mdct_plc_active) && (st->nbLostCmpt == 1)
                    && (tcx_cfg->tcx_last_overlap_mode != FULL_OVERLAP)
                    && (tcx_cfg->tcx_curr_overlap_mode != FULL_OVERLAP) )
            {
                float E_2ndlast, E_last;

                E_2ndlast = E_last = EPSILON;
                for( i=0; i<infoIGFStartLine; i=i+2 )
                {
                    E_2ndlast += st->tonalMDCTconceal.lastBlockData.spectralData[i]  *st->tonalMDCTconceal.lastBlockData.spectralData[i];
                    E_last    += st->tonalMDCTconceal.lastBlockData.spectralData[i+1]*st->tonalMDCTconceal.lastBlockData.spectralData[i+1];
                }
                tmp2 = E_2ndlast/E_last;

                /* replace higher energy TCX5 frame by lower one to avoid energy fluctuation */
                if( tmp2 > 2 )
                {
                    for( i=0; i<infoIGFStartLine; i=i+2 )
                    {
                        st->tonalMDCTconceal.lastBlockData.spectralData[i] = st->tonalMDCTconceal.lastBlockData.spectralData[i+1];
                    }
                }
                else if( tmp2 < 0.5 )
                {
                    for( i=0; i<infoIGFStartLine; i=i+2 )
                    {
                        st->tonalMDCTconceal.lastBlockData.spectralData[i+1] = st->tonalMDCTconceal.lastBlockData.spectralData[i];
                    }
                }
            }

            noiseTiltFactor = 1.0f;

            tcxGetNoiseFillingTilt( A, L_frame, (st->bits_frame >= 256 && !st->rf_flag), &noiseTiltFactor );

            TonalMDCTConceal_InsertNoise( &st->tonalMDCTconceal, x, st->tonal_mdct_plc_active,
                                          &st->seed_tcx_plc, noiseTiltFactor, f, infoIGFStartLine );

        }
    }

    if( L_spec < L_frame )
    {
        set_zero( x+L_spec, L_frame-L_spec );
    }
    else if( L_spec > L_frameTCX )
    {
        set_zero( x+L_frameTCX, L_spec-L_frameTCX );
    }

    if( bfi && (!st->enablePlcWaveadjust || (temp_concealment_method == TCX_TONAL))
            && st->igf && (frame_cnt == 0) && (L_frameTCX == st->L_frameTCX >> 1)
            && (st->tcxonly) && (!st->tonal_mdct_plc_active) && (st->nbLostCmpt == 1)
            && (tcx_cfg->tcx_last_overlap_mode != FULL_OVERLAP)
            && (tcx_cfg->tcx_curr_overlap_mode != FULL_OVERLAP) )
    {
        IGFDecCopyLPCFlatSpectrum( &st->hIGFDec, x, IGF_GRID_LB_SHORT );
        mvi2i(st->hIGFDec.igfData.igf_curr_subframe[0][0], st->hIGFDec.igfData.igf_curr_subframe[1][0], IGF_MAX_SFB);
    }

    /*-----------------------------------------------------------*
     * Noise shaping in frequency domain (1/Wz)                  *
     *-----------------------------------------------------------*/

    if( st->igf && ! bfi )
    {
        if ((L_frame == st->L_frame >> 1) && (st->tcxonly))
        {
            IGFDecCopyLPCFlatSpectrum( &st->hIGFDec, x, IGF_GRID_LB_SHORT );
        }
        else
        {
            IGFDecCopyLPCFlatSpectrum( &st->hIGFDec, x, (st->last_core == ACELP_CORE)?IGF_GRID_LB_TRAN:IGF_GRID_LB_NORM );
        }
    }
    /* LPC gains already available */

    if( !st->enablePlcWaveadjust || !bfi || (temp_concealment_method == TCX_TONAL) )
    {
        mdct_noiseShaping( x, L_frame, gainlpc2 );

        if( !bfi )
        {
            v_multc( x+L_frame, gainlpc2[FDNS_NPTS-1], x+L_frame, L_spec-L_frame );
        }

        set_zero( x+L_spec, L_frameTCX-L_spec );
    }

    /* PLC: [TCX: Tonal Concealment] */
    if( bfi && st->tonal_mdct_plc_active )
    {
        TonalMDCTConceal_Apply( &st->tonalMDCTconceal, x );
    }

    TonalMDCTConceal_UpdateState( &st->tonalMDCTconceal, L_frameTCX, (st->tcxltp_last_gain_unmodified > 0) ? st->old_fpitch : 0,
                                  bfi, bfi && st->tonal_mdct_plc_active );


    if( st->enablePlcWaveadjust )
    {
        int core;
        core = st->core;

        /* spectrum concealment */
        if (bfi && temp_concealment_method == TCX_NONTONAL)
        {
            concealment_decode( core, x, &st->plcInfo );
        }

        /* update spectrum buffer, tonality flag, etc. */
        concealment_update( bfi, core, st->tonality_flag, x, &st->plcInfo );
    }

    /*-----------------------------------------------------------*
     * IGF                                                       *
     *-----------------------------------------------------------*/

    if (st->igf && !((L_frame == st->L_frame >> 1) && (st->tcxonly)))
    {
        /* copy low spectrum to IGF des buffer */
        st->hIGFDec.igfData.igfInfo.nfSeed = (short)(nf_seed * 31821L + 13849L);

        IGFDecApplyMono( &st->hIGFDec, x, (st->last_core == ACELP_CORE || (left_rect && bfi))?IGF_GRID_LB_TRAN:IGF_GRID_LB_NORM, bfi );

        /* IGF Decoder - All */
    }

    if (st->igf && ((L_frame == st->L_frame >> 1) && (st->tcxonly)))
    {
        /* copy low spectrum to IGF des buffer */
        st->hIGFDec.igfData.igfInfo.nfSeed = (short)(nf_seed * 31821L + 13849L);

        IGFDecApplyMono( &st->hIGFDec, x, IGF_GRID_LB_SHORT, bfi );

        /* IGF Decoder - All */
    }

    index = tcx_cfg->tcx_last_overlap_mode;  /* backup last TCX overlap mode */

    if( st->igf )
    {
        int proc = st->hIGFDec.flatteningTrigger;
        if( proc && fUseTns != 0 )
        {
            proc = 0;
        }

        if( proc )
        {
            short int startLine = st->hIGFDec.infoIGFStartLine;
            short int endLine = st->hIGFDec.infoIGFStopLine;
            float x_itf[N_MAX_TCX-IGF_START_MN];
            int j;

            const int* chk_sparse = st->hIGFDec.flag_sparse;
            const float* virtualSpec = st->hIGFDec.virtualSpec;

            const int maxOrder = 8;
            int curr_order = 0;
            float A[ITF_MAX_FILTER_ORDER+1];
            float predictionGain = 0;

            for (j = startLine; j < endLine; j++)
            {
                if( chk_sparse[j-IGF_START_MN] == 2 )
                {
                    x_itf[j-IGF_START_MN] = x[j];
                    x[j] = virtualSpec[j-IGF_START_MN];
                }
            }

            ITF_Detect( x+IGF_START_MN, startLine, endLine, maxOrder, A, &predictionGain, &curr_order );

            ITF_Apply(x, startLine, endLine, A, curr_order);

            for (j = startLine; j < endLine; j++)
            {
                if( chk_sparse[j-IGF_START_MN] == 2 )
                {
                    x[j] = x_itf[j-IGF_START_MN];
                }
            }
        }
        else
        {

        }
    }

    if( (L_frame == st->L_frame >> 1) && st->tcxonly )
    {
        int L = L_frameTCX;
        if( (tcx_cfg->fIsTNSAllowed && fUseTns != 0 && bfi!= 1) || (L_spec > L_frameTCX) )
        {
            L = L_spec;
        }
        tcxInvertWindowGrouping( tcx_cfg, xn_buf, x, L, fUseTns, st->last_core, index, frame_cnt, bfi );
    }

    /*-----------------------------------------------------------*
     * Temporal Noise Shaping Synthesis                          *
     *-----------------------------------------------------------*/
    if( tcx_cfg->fIsTNSAllowed && fUseTns != 0 && bfi!= 1 )
    {
        /* Apply TNS to get the reconstructed signal */

        SetTnsConfig(tcx_cfg, L_frame_glob == st->L_frame, (st->last_core == ACELP_CORE) && (frame_cnt == 0));

        ApplyTnsFilter(tcx_cfg->pCurrentTnsConfig, &tnsData, x, 0);

        if ((L_frame == st->L_frame >> 1) && (st->tcxonly))
        {
            if ((tcx_cfg->tcx_last_overlap_mode != FULL_OVERLAP) ||
                    ((tcx_cfg->tcx_curr_overlap_mode == FULL_OVERLAP) && (frame_cnt == 0) && (index == 0)) )
            {
                const int L_win = tcx_cfg->tnsConfig[0][0].iFilterBorders[0] >> 1;
                /* undo rearrangement of LF sub-window lines for TNS synthesis filtering */
                if (L_frameTCX > 2*L_win)
                {
                    const int L_win2 = L_frameTCX >> 1;
                    mvr2r( x + L_win + 8, x + L_win2 + 8, L_win - 8 );
                    mvr2r( x + 8, x + L_win2, 8 );
                    mvr2r( x + 16, x + 8, L_win - 8 );
                    set_zero( x + L_win, L_win2 - L_win );
                    set_zero( x + L_win2 + L_win, L_win2 - L_win );
                }
                else
                {
                    mvr2r(x+8, xn_buf, L_win );
                    mvr2r( xn_buf, x+L_win, 8 );
                    mvr2r( xn_buf+8, x+8, L_win-8 );
                }

            }
        }
    }


    /*-----------------------------------------------------------*
     * Prepare OLA buffer after waveadjustment.                  *
     * Compute inverse MDCT of x[].                              *
     *-----------------------------------------------------------*/

    mvr2r( x, xn_bufFB, max(L_spec, max(L_frame, L_frameTCX)) );

    if( st->igf )
    {
        set_zero(xn_bufFB+st->hIGFDec.infoIGFStartLine, L_frameTCX-st->hIGFDec.infoIGFStartLine);
    }

    IMDCT( xn_bufFB,
           st->syn_Overl,
           st->syn_Overl_TDAC,
           xn_buf,
           tcx_cfg->tcx_aldo_window_1_trunc,
           tcx_cfg->tcx_aldo_window_2,
           tcx_cfg->tcx_mdct_window_half,
           tcx_cfg->tcx_mdct_window_minimum,
           tcx_cfg->tcx_mdct_window_trans,
           tcx_cfg->tcx_mdct_window_half_length,
           tcx_cfg->tcx_mdct_window_min_length,
           index,
           left_rect,
           tcx_offset,
           overlap,
           L_frame,
           L_frameTCX,
           max(L_frameTCX, L_spec) >> 1,
           L_frame_glob,
           frame_cnt,
           bfi,
           st->old_outLB,
           0,
           st,
           0,
           acelp_zir );

    /* Generate additional comfort noise to mask potential coding artefacts */
    if( st->flag_cna)
    {
        generate_masking_noise_mdct(x, st->hFdCngDec->hFdCngCom);
    }

    IMDCT( x,
           st->syn_OverlFB,
           st->syn_Overl_TDACFB,
           xn_bufFB,
           tcx_cfg->tcx_aldo_window_1_FB_trunc,
           tcx_cfg->tcx_aldo_window_2_FB,
           tcx_cfg->tcx_mdct_window_halfFB,
           tcx_cfg->tcx_mdct_window_minimumFB,
           tcx_cfg->tcx_mdct_window_transFB,
           tcx_cfg->tcx_mdct_window_half_lengthFB,
           tcx_cfg->tcx_mdct_window_min_lengthFB,
           index,
           left_rect,
           tcx_offsetFB,
           overlapFB,
           L_frameTCX,
           L_frameTCX,
           max(L_frameTCX, L_spec) >> 1,
           L_frameTCX_glob,
           frame_cnt,
           bfi,
           st->old_out,
           1,
           st,
           FSCALE_DENOM * L_frameTCX_glob / L_frame_glob,
           acelp_zir );


    /* PLC: [TCX: Tonal Concealment] */
    if (!bfi)
    {

        st->second_last_tns_active = st->last_tns_active;
        st->last_tns_active = tcx_cfg->fIsTNSAllowed & fUseTns;
        st->tcxltp_third_last_pitch  = st->tcxltp_second_last_pitch;
        st->tcxltp_second_last_pitch = st->old_fpitch;
        st->old_fpitch = st->tcxltp_pitch_int + st->tcxltp_pitch_fr/(float)st->pit_res_max;
        st->old_fpitchFB = st->old_fpitch * (float)L_frameTCX / (float)L_frame;
    }

    /* Update old_syn_overl */
    if (!tcx_cfg->last_aldo)
    {
        mvr2r(xn_buf+L_frame, st->syn_Overl, overlap);
        mvr2r(xn_bufFB+L_frameTCX, st->syn_OverlFB, overlapFB);
    }

    /* Output */
    mvr2r( xn_buf+(overlap>>1)-tcx_offset, synth, L_frame_glob );
    mvr2r( xn_bufFB+(overlapFB>>1)-tcx_offsetFB, synthFB, L_frameTCX_glob );

    return;
}


/*-------------------------------------------------------------------*
 * decoder_tcx_post()
 *
 *
 *-------------------------------------------------------------------*/

void decoder_tcx_post(
    Decoder_State *st,
    float *synth,
    float *synthFB,
    float *A,
    int bfi
)
{
    int i;
    float level_syn, gainCNG = 0.0f, step;
    float xn_buf[L_FRAME_MAX];

    /* TCX output */
    mvr2r( synth, xn_buf, st->L_frame );

    /* first TCX frame after ACELP; overwrite ltp initialization done during acelp PLC */
    if( !bfi && st->prev_bfi && !st->last_core )
    {
        st->tcxltp_last_gain_unmodified = 0.0f;
    }

    if (bfi && !st->use_partial_copy)
    {
        /* run lpc gain compensation not for waveform adjustment */
        if (!st->enablePlcWaveadjust ||  st->plcInfo.concealment_method == TCX_TONAL )
        {
            float gainHelperFB     = st->gainHelper;
            float stepCompensateFB = st->stepCompensate * st->L_frame / st->L_frameTCX;

            for( i=0; i < st->L_frameTCX; i++ )
            {
                synthFB[i] *= gainHelperFB;
                gainHelperFB -= stepCompensateFB;
            }
        }

        for( i=0; i < st->L_frame; i++ )
        {
            xn_buf[i] *= st->gainHelper;
            st->gainHelper -= st->stepCompensate;
        }
    }

    /* PLC: [TCX: Fade-out]
     * PLC: estimate and update CNG energy */
    level_syn = (float)sqrt(( dotp(synthFB, synthFB, st->L_frameTCX)) / st->L_frameTCX);


    /* PLC: [TCX: Fade-out]
     * PLC: update or retrieve the background level */
    if( bfi == 0 && st->tcxonly && st->clas_dec == UNVOICED_CLAS )
    {
        minimumStatistics( st->NoiseLevelMemory_bfi, &st->NoiseLevelIndex_bfi, &st->CurrLevelIndex_bfi, &st->CngLevelBackgroundTrace_bfi,
                           &st->LastFrameLevel_bfi, level_syn, PLC_MIN_CNG_LEV, PLC_MIN_STAT_BUFF_SIZE);
    }

    /* PLC: [TCX: Fade-out]
     * PLC: fade-out in time domain */
    if( bfi )
    {
        float conceal_eof_gainFB;

        if( st->tcxonly )
        {
            gainCNG = st->CngLevelBackgroundTrace_bfi/(level_syn+0.01f);
        }
        else
        {
            gainCNG = st->cngTDLevel/(level_syn+0.01f);
        }
        if( st->nbLostCmpt == 1 )
        {
            st->conceal_eof_gain = 1.0f;
        }

        step = (st->conceal_eof_gain - ( st->conceal_eof_gain * st->damping + gainCNG * (1 - st->damping) )) / st->L_frame;
        {
            float stepFB = step * st->L_frame / st->L_frameTCX;
            conceal_eof_gainFB = st->conceal_eof_gain;

            for( i=0; i < st->L_frameTCX; i++ )
            {
                synthFB[i] *= conceal_eof_gainFB;
                conceal_eof_gainFB -= stepFB;
            }
        }

        for( i=0; i < st->L_frame; i++ )
        {
            xn_buf[i] *= st->conceal_eof_gain;
            st->conceal_eof_gain -= step;
        }

        /* run lpc gain compensation not for waveform adjustment */
        if( (!st->enablePlcWaveadjust || st->plcInfo.concealment_method == TCX_TONAL ) && !st->use_partial_copy )
        {
            st->plcInfo.recovery_gain = conceal_eof_gainFB * st->last_concealed_gain_syn_deemph;
        }
        else
        {
            st->plcInfo.recovery_gain = conceal_eof_gainFB;
        }
        st->plcInfo.step_concealgain = step * st->L_frame / st->L_frameTCX;
    }


    /*-----------------------------------------------------------*
     * Memory update                                             *
     *-----------------------------------------------------------*/

    /* Update synth, exc and old_Aq */
    tcx_decoder_memory_update( xn_buf, synth, st->L_frame, A, st, st->syn );

    /* PLC: [TCX: Memory update] */
    st->old_pitch_buf[0] = st->old_pitch_buf[st->nb_subfr];
    st->old_pitch_buf[1] = st->old_pitch_buf[st->nb_subfr+1];
    mvr2r( &st->old_pitch_buf[st->nb_subfr+2], &st->old_pitch_buf[2], st->nb_subfr );
    set_f( &st->old_pitch_buf[st->nb_subfr+2], st->old_fpitch, st->nb_subfr );
    st->bfi_pitch = st->old_fpitch;
    st->bfi_pitch_frame = st->L_frame;

    st->mem_pitch_gain[2*st->nb_subfr+1] = st->mem_pitch_gain[st->nb_subfr+1];
    st->mem_pitch_gain[2*st->nb_subfr]   = st->mem_pitch_gain[st->nb_subfr];

    for( i = 0; i < st->nb_subfr; i++ )
    {
        st->mem_pitch_gain[2*st->nb_subfr-1 - i] = st->mem_pitch_gain[st->nb_subfr-1 - i];
        st->mem_pitch_gain[st->nb_subfr-1 - i]   = st->tcxltp_last_gain_unmodified;
    }

    return;
}

/*-------------------------------------------------------------------*
 * IMDCT()
 *
 *
 *-------------------------------------------------------------------*/

static void IMDCT(
    float *x,
    float *old_syn_overl,
    float *syn_Overl_TDAC,
    float *xn_buf,
    float *tcx_aldo_window_1_trunc,
    float *tcx_aldo_window_2,
    float *tcx_mdct_window_half,
    float *tcx_mdct_window_minimum,
    float *tcx_mdct_window_trans,
    int tcx_mdct_window_half_length,
    int tcx_mdct_window_min_length,
    int index,
    int left_rect,
    int tcx_offset,
    int overlap,
    int L_frame,
    int L_frameTCX,
    int L_spec_TCX5,
    int L_frame_glob,
    int frame_cnt,
    int bfi,
    float *old_out,
    short FB_flag,
    Decoder_State *st,
    int fullbandScale,
    float *acelp_zir
)
{
    short i, nz, aldo;
    TCX_config *tcx_cfg = &st->tcx_cfg;

    aldo = 0;

    /* number of zero for ALDO windows*/
    nz = NS2SA(st->output_Fs, N_ZERO_MDCT_NS)*L_frame/L_frameTCX;

    if( (L_frameTCX == st->L_frameTCX >> 1) && (st->tcxonly) )
    {
        /* Mode decision in PLC

        last OL   curr OL   left TCX-10  right TCX-10
        -------------------------------------------------------------
            0      0         2x TCX-5* 1x  TCX-10
            0      2         1x TCX-10 1x  TCX-10
            0      3         1x TCX-10 1x  TCX-10
            2      0         2x TCX-5  1x  TCX-10
            2      2         2x TCX-5  2x  TCX-5
            2      3         2x TCX-5  2x  TCX-5
            3      0         2x TCX-5  1x  TCX-10
            3      2         2x TCX-5  2x  TCX-5
            3      3         2x TCX-5  2x  TCX-5
        */

        if( (!bfi && tcx_cfg->tcx_last_overlap_mode != FULL_OVERLAP) || (bfi && (tcx_cfg->tcx_last_overlap_mode != FULL_OVERLAP) && (tcx_cfg->tcx_curr_overlap_mode != FULL_OVERLAP)) )
        {
            /* minimum or half overlap, two transforms, grouping into one window */
            float win[(L_FRAME_PLUS+L_MDCT_OVLP_MAX)/2];
            const int L_win = L_frame >> 1;
            const int L_ola = (tcx_cfg->tcx_last_overlap_mode == MIN_OVERLAP) ? tcx_mdct_window_min_length : tcx_mdct_window_half_length;
            int w;

            set_f( win, 0, (L_FRAME_PLUS+L_MDCT_OVLP_MAX)/2 );
            set_zero(xn_buf, tcx_offset+(L_ola>>1));  /* zero left end of buffer */

            for (w = 0; w < 2; w++)
            {
                TCX_MDCT_Inverse(x+w*L_spec_TCX5, win, L_ola, L_win-L_ola, L_ola);

                tcx_windowing_synthesis_current_frame( win, tcx_aldo_window_2, tcx_mdct_window_half, tcx_mdct_window_minimum, L_ola, tcx_mdct_window_half_length, tcx_mdct_window_min_length,
                                                       (w > 0) ? 0 : left_rect, (w > 0) || (w == 0 && index == 2) ? MIN_OVERLAP : tcx_cfg->tcx_last_overlap_mode,
                                                       acelp_zir, st->old_syn_Overl, syn_Overl_TDAC, st->old_Aq_12_8, tcx_mdct_window_trans, L_win, tcx_offset<0?-tcx_offset:0,
                                                       (w > 0) || (frame_cnt > 0) ? 1 : st->last_core_bfi, (w > 0) || (frame_cnt > 0) ? 0 : st->last_is_cng, fullbandScale );

                if( w > 0 )
                {
                    tcx_windowing_synthesis_past_frame( xn_buf+tcx_offset-(L_ola>>1)+w*L_win, tcx_aldo_window_1_trunc, tcx_mdct_window_half,
                                                        tcx_mdct_window_minimum, L_ola, tcx_mdct_window_half_length, tcx_mdct_window_min_length, MIN_OVERLAP );
                }

                /* add part of current sub-window overlapping with previous window */
                v_add( win, xn_buf+tcx_offset-(L_ola>>1)+w*L_win, xn_buf+tcx_offset-(L_ola>>1)+w*L_win, L_ola );

                /* copy new sub-window region not overlapping with previous window */
                mvr2r( win+L_ola, xn_buf+tcx_offset+(L_ola>>1)+w*L_win, L_win );
            }

            /* To assure that no garbage values are passed to overlap */
            set_zero( xn_buf+L_frame+tcx_offset+(L_ola>>1), overlap-tcx_offset-(L_ola>>1) );
        }
        else if (!bfi && (frame_cnt == 0) && (tcx_cfg->tcx_curr_overlap_mode == FULL_OVERLAP))
        {
            /* special overlap attempt, two transforms, grouping into one window */
            float win[(L_FRAME_PLUS+L_MDCT_OVLP_MAX)/2];
            const int L_win = L_frame >> 1;
            const int L_ola = tcx_mdct_window_min_length;
            int w;

            set_f( win, 0, (L_FRAME_PLUS+L_MDCT_OVLP_MAX)/2 );

            /* 1st TCX-5 window, special MDCT with minimum overlap on right side */
            TCX_MDCT_Inverse(x, win+L_win, 0, L_win-(L_ola>>1), L_ola);

            set_zero(xn_buf, (overlap>>1));

            /* copy new sub-window region not overlapping with previous window */
            mvr2r(win+L_win, xn_buf+(overlap>>1), L_win+(L_ola>>1));

            /* 2nd TCX-5 window, regular MDCT with minimum overlap on both sides */
            TCX_MDCT_Inverse(x+L_spec_TCX5, win, L_ola, L_win-L_ola, L_ola);

            tcx_windowing_synthesis_current_frame( win, tcx_aldo_window_2, tcx_mdct_window_half, tcx_mdct_window_minimum, L_ola, tcx_mdct_window_half_length, tcx_mdct_window_min_length, 0,
                                                   /* left_rect */ MIN_OVERLAP,  /* left_mode */ acelp_zir, st->old_syn_Overl, syn_Overl_TDAC, st->old_Aq_12_8, tcx_mdct_window_trans, L_win,
                                                   tcx_offset<0?-tcx_offset:0,
                                                   1, /* st->last_mode_bfi */ 0, /* st->last_is_cng */ fullbandScale );

            tcx_windowing_synthesis_past_frame( xn_buf+(overlap>>1)+L_win-(L_ola>>1), tcx_aldo_window_1_trunc, tcx_mdct_window_half, tcx_mdct_window_minimum,
                                                L_ola, tcx_mdct_window_half_length, tcx_mdct_window_min_length, 2 );

            /* add part of current sub-window overlapping with previous window */
            v_add(win, xn_buf+(overlap>>1)+L_win-(L_ola>>1), xn_buf+(overlap>>1)+L_win-(L_ola>>1), L_ola);

            /* copy new sub-window region not overlapping with previous window */
            mvr2r(win+L_ola, xn_buf+(overlap>>1)+L_win+(L_ola>>1), L_win);

            /* extra folding-out on left side of win, for perfect reconstruction */
            for (w = (overlap>>1); w < overlap; w++)
            {
                xn_buf[overlap-1-w] = -1.0f * xn_buf[w];
            }

            tcx_windowing_synthesis_current_frame( xn_buf, tcx_aldo_window_2, tcx_mdct_window_half, tcx_mdct_window_minimum, overlap, tcx_mdct_window_half_length,
                                                   tcx_mdct_window_min_length, left_rect, 0,  /* left_mode */ acelp_zir, st->old_syn_Overl, syn_Overl_TDAC, st->old_Aq_12_8,
                                                   tcx_mdct_window_trans, 2*L_win, tcx_offset<0?-tcx_offset:0,
                                                   st->last_core_bfi, st->last_is_cng, fullbandScale );
        }
        else
        {
            /* default, i.e. maximum overlap, single transform, no grouping */

            TCX_MDCT_Inverse(x, xn_buf, overlap, L_frame-overlap, overlap);

            tcx_windowing_synthesis_current_frame( xn_buf, tcx_aldo_window_2, tcx_mdct_window_half, tcx_mdct_window_minimum, overlap, tcx_mdct_window_half_length,
                                                   tcx_mdct_window_min_length, left_rect, !bfi && (frame_cnt > 0) && (index == 0) && ( st->last_core != ACELP_CORE) ? MIN_OVERLAP : index,
                                                   acelp_zir, st->old_syn_Overl, syn_Overl_TDAC, st->old_Aq_12_8, tcx_mdct_window_trans, L_frame_glob >> 1,
                                                   tcx_offset<0?-tcx_offset:0,
                                                   (frame_cnt > 0 /*|| (st->last_con_tcx )*/) ? 1 : st->last_core_bfi, (frame_cnt > 0) ? 0 : st->last_is_cng, fullbandScale );

        } /* tcx_last_overlap_mode != FULL_OVERLAP */
    }
    else
    {
        /* frame is TCX-20 or not TCX-only */
        assert(frame_cnt == 0);
        if (st->tcx_cfg.tcx_last_overlap_mode != TRANSITION_OVERLAP)
        {
            float tmp[L_FRAME_PLUS];

            edct(x, xn_buf+overlap/2+nz, L_frame);
            v_multc( xn_buf+overlap/2+nz, (float)sqrt((float)L_frame / NORM_MDCT_FACTOR), tmp, L_frame);
            window_ola( tmp, xn_buf, old_out, L_frame, tcx_cfg->tcx_last_overlap_mode, tcx_cfg->tcx_curr_overlap_mode, 0, 0, NULL );
            aldo = 1;

        }
        else
        {
            TCX_MDCT_Inverse(x, xn_buf, overlap, L_frame-overlap, overlap);

            /*-----------------------------------------------------------*
             * Windowing, overlap and add                                *
             *-----------------------------------------------------------*/

            /* Window current frame */
            tcx_windowing_synthesis_current_frame( xn_buf, tcx_aldo_window_2, tcx_mdct_window_half, tcx_mdct_window_minimum, overlap, tcx_mdct_window_half_length, tcx_mdct_window_min_length,
                                                   left_rect, tcx_cfg->tcx_last_overlap_mode, acelp_zir, st->old_syn_Overl, syn_Overl_TDAC, st->old_Aq_12_8, tcx_mdct_window_trans, L_frame_glob >> 1,
                                                   tcx_offset<0?-tcx_offset:0,
                                                   st->last_core_bfi,
                                                   st->last_is_cng,
                                                   fullbandScale );
        }
    } /* TCX-20 and TCX-only */

    /* Window and overlap-add past frame if past frame is TCX */
    if( (frame_cnt != 0) || (st->last_core_bfi > ACELP_CORE) )
    {
        if (((L_frameTCX == st->L_frameTCX >> 1) && (st->tcxonly)) || (st->tcx_cfg.tcx_last_overlap_mode == TRANSITION_OVERLAP))
        {
            if (!bfi && (frame_cnt > 0) && (index == 0) && (tcx_cfg->tcx_curr_overlap_mode == FULL_OVERLAP) && (st->last_core != ACELP_CORE))
            {
                index = MIN_OVERLAP;     /* use minimum overlap between the two TCX-10 windows */
            }

            if( tcx_cfg->last_aldo )
            {
                for (i=0; i < overlap - tcx_mdct_window_min_length; i++)
                {
                    xn_buf[i+overlap/2-tcx_offset] += old_out[i+nz];
                }

                /* fade truncated ALDO window */
                for ( ; i < overlap; i++)
                {
                    xn_buf[i+overlap/2-tcx_offset] += old_out[i+nz] * tcx_mdct_window_minimum[overlap-1-i];
                }
            }
            else
            {
                tcx_windowing_synthesis_past_frame( old_syn_overl
                                                    ,
                                                    tcx_aldo_window_1_trunc, tcx_mdct_window_half, tcx_mdct_window_minimum,
                                                    overlap, tcx_mdct_window_half_length, tcx_mdct_window_min_length,
                                                    (index == 0 || tcx_cfg->tcx_last_overlap_mode == MIN_OVERLAP) ? tcx_cfg->tcx_last_overlap_mode : index );
                if (bfi)
                {
                    for (i=0; i<tcx_mdct_window_half_length; i++)
                    {
                        xn_buf[i+overlap/2-tcx_offset] += old_syn_overl[i]*tcx_mdct_window_half[tcx_mdct_window_half_length-1-i];
                    }
                }
                else if(!left_rect)
                {
                    for (i=0; i<overlap; i++)
                    {
                        xn_buf[i] += old_syn_overl[i];
                    }
                }
                else
                {
                    for (i=0; i<overlap; i++)
                    {
                        xn_buf[i+(overlap>>1)] += old_syn_overl[i];
                    }
                }
            }
        }
    }

    if( !aldo && (((L_frameTCX == st->L_frameTCX >> 1)&&frame_cnt > 0) || L_frameTCX != (st->L_frameTCX >> 1)) )
    {
        /*Compute windowed synthesis in case of switching to ALDO windows in next frame*/
        mvr2r( xn_buf+L_frame-nz, old_out, nz+overlap );
        set_zero( old_out+nz+overlap, nz );

        tcx_windowing_synthesis_past_frame( old_out+nz, tcx_aldo_window_1_trunc, tcx_mdct_window_half, tcx_mdct_window_minimum,
                                            overlap, tcx_mdct_window_half_length, tcx_mdct_window_min_length, tcx_cfg->tcx_curr_overlap_mode );

        /* If current overlap mode = FULL_OVERLAP -> ALDO_WINDOW */
        if( tcx_cfg->tcx_curr_overlap_mode == FULL_OVERLAP )
        {
            for (i=0; i<nz; i++)
            {
                old_out[nz+overlap+i]=xn_buf[L_frame-1-i]*tcx_aldo_window_1_trunc[-1-i];
            }
            aldo = 1;
        }
    }

    if( FB_flag )
    {
        tcx_cfg->last_aldo = aldo;
    }

    /* Smoothing between the ACELP PLC and TCX Transition frame. Using the shape of the half overlap window for the crossfading. */
    if( left_rect && (frame_cnt == 0) && (st->last_core_bfi == ACELP_CORE) && st->prev_bfi )
    {
        if( FB_flag )
        {
            for (i=0; i<tcx_mdct_window_half_length; i++)
            {
                xn_buf[i+overlap/2-tcx_offset] *= tcx_mdct_window_half[i];
                xn_buf[i+overlap/2-tcx_offset] += st->syn_OverlFB[i] * tcx_mdct_window_half[tcx_mdct_window_half_length-1-i] * tcx_mdct_window_half[tcx_mdct_window_half_length-1-i];
            }
        }
        else
        {
            for (i=0; i<tcx_mdct_window_half_length; i++)
            {
                xn_buf[i+overlap/2-tcx_offset] *= tcx_mdct_window_half[i];
                xn_buf[i+overlap/2-tcx_offset] += st->syn_Overl[i] * tcx_mdct_window_half[tcx_mdct_window_half_length-1-i] * tcx_mdct_window_half[tcx_mdct_window_half_length-1-i];
            }
        }
    }

    return;
}
