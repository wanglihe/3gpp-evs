/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"


/*-------------------------------------------------------------------*
* core_signal_analysis_high_bitrate()
*
*
*-------------------------------------------------------------------*/

void core_signal_analysis_high_bitrate(
    const float *new_samples,
    const short T_op[3],        /* i  : open-loop pitch values for quantiz. */
    const float voicing[3],     /* i  : open-loop pitch gains               */
    const short pitch[2],       /* i  : open-loop pitch @12.8kHz for adapt. lag windowing */
    float lsp_new[],
    float lsp_mid[],
    Encoder_State *st,
    int pTnsSize[],
    int pTnsBits[],
    int param_core[],
    int *ltpBits,
    int L_frame,
    int L_frameTCX
)
{
    const int last_overlap = st->tcx_cfg.tcx_last_overlap_mode;
    const int curr_overlap = st->tcx_cfg.tcx_curr_overlap_mode;
    int i, frameno;
    int L_subframe;
    int left_overlap=-1, right_overlap=-1, folding_offset;
    float buf[N_MAX+L_MDCT_OVLP_MAX]; /* Buffer for TCX20/TCX10 windowing */
    float mdstWin[N_MAX+L_MDCT_OVLP_MAX]; /* Buffer for MDST windowing */
    float * powerSpec;
    float * tcx20Win;
    float tcx5Win[N_TCX10_MAX/2+L_MDCT_OVLP_MAX]; /* Buffer for TCX5 windowing and interleaving. */
    float * interleaveBuf = tcx5Win;
    int nSubframes;
    short overlap_mode[3];
    short transform_type[2];
    float r[M+1];
    float A[M+1];
    float * lsp[2];
    const int tcx10SizeFB = 2*st->tcx_cfg.tcx5SizeFB;
    const int tcx5SizeFB  =   st->tcx_cfg.tcx5SizeFB;
    const int tcx10Size = 2*st->tcx_cfg.tcx5Size;
    short alw_pitch_lag_12k8[2], alw_pitch_lag_12k8_wc;
    float alw_voicing[2], alw_voicing_wc;

    powerSpec = tcx20Win = buf; /* Share memory for windowed TD signal and for the power spectrum */

    /*--------------------------------------------------------------*
    * Input Signal Processing: copy, HP filter, pre-emphasis
    *---------------------------------------------------------------*/

    /* Copy Samples */
    mvr2r( new_samples, st->new_speech_enc, L_frame );

    /*--------------------------------------------------------------*
    * TCX-LTP
    *---------------------------------------------------------------*/

    tcx_ltp_encode( st->tcxltp, st->tcxonly, st->tcxMode, L_frame, L_SUBFR, st->speech_enc+st->encoderLookahead_enc,
                    st->speech_ltp+st->encoderLookahead_enc, st->speech_enc+st->encoderLookahead_enc,
                    T_op[1], &param_core[1+NOISE_FILL_RANGES], ltpBits, &st->tcxltp_pitch_int, &st->tcxltp_pitch_fr,
                    &st->tcxltp_gain, &st->tcxltp_pitch_int_past, &st->tcxltp_pitch_fr_past, &st->tcxltp_gain_past,
                    &st->tcxltp_norm_corr_past, st->last_core, st->pit_min, st->pit_fr1, st->pit_fr2, st->pit_max,
                    st->pit_res_max, &st->transientDetection, (st->sr_core > 25600), NULL, M );

    mvr2r( st->speech_enc+st->encoderLookahead_enc, st->new_speech_enc_pe, L_frame );

    preemph( st->new_speech_enc_pe, st->preemph_fac, L_frame, &(st->mem_preemph_enc) );

    if( st->tcxMode == TCX_10 )
    {
        mvi2i( &param_core[1+NOISE_FILL_RANGES], &param_core[NPRM_DIV+1+NOISE_FILL_RANGES], LTPSIZE );
    }

    lsp[0] = lsp_new;
    lsp[1] = lsp_mid;

    /*-------------------------------------------------------------------------*
    * Decision matrix for the transform and overlap length
    *--------------------------------------------------------------------------*/

    alw_pitch_lag_12k8[0] = pitch[0];
    alw_pitch_lag_12k8[1] = pitch[1];
    alw_voicing[0] = voicing[0];
    alw_voicing[1] = voicing[1];
    alw_pitch_lag_12k8_wc = min(alw_pitch_lag_12k8[0], alw_pitch_lag_12k8[1]);
    alw_voicing_wc = max(alw_voicing[0], alw_voicing[1]);
    overlap_mode[0] = last_overlap; /* Overlap between the last and the current frame */

    if (st->tcxMode == TCX_20)
    {
        nSubframes = 1;
        transform_type[0] = TCX_20;
        overlap_mode[1] = curr_overlap; /* Overlap between the current and the next frame */
        alw_pitch_lag_12k8[0] = alw_pitch_lag_12k8_wc;
        alw_voicing[0] = alw_voicing_wc;
    }
    else
    {
        nSubframes = 2;
        if (curr_overlap == FULL_OVERLAP)
        {
            transform_type[0] = TCX_5;
            transform_type[1] = TCX_10;
            overlap_mode[1] = (last_overlap == HALF_OVERLAP) ? HALF_OVERLAP : MIN_OVERLAP; /* Overlap between 2nd and 3rd sub-frame */
        }
        else if (last_overlap == FULL_OVERLAP)
        {
            transform_type[0] = TCX_10;
            transform_type[1] = TCX_5;
            overlap_mode[1] = (curr_overlap == HALF_OVERLAP) ? HALF_OVERLAP : MIN_OVERLAP; /* Overlap between 1st and 2nd sub-frame */
        }
        else
        {
            transform_type[0] = transform_type[1] = TCX_5;
            overlap_mode[1] = (last_overlap == HALF_OVERLAP && curr_overlap == HALF_OVERLAP) ? HALF_OVERLAP : MIN_OVERLAP; /* Overlap between 2nd and 3rd sub-frame */
        }
        overlap_mode[2] = curr_overlap; /* Overlap between the current and the next frame */
    }

    if (transform_type[0] != TCX_20)
    {
        IGFEncResetTCX10BitCounter(&st->hIGFEnc);
    }

    for (frameno = 0; frameno < nSubframes; frameno++)
    {
        /*-------------------------------------------------------------------------*
        * Get MDCT output and TNS parameters. Apply TNS in the spectrum if needed
        *--------------------------------------------------------------------------*/

        L_subframe = L_frameTCX/nSubframes;

        if( (transform_type[frameno] == TCX_20) && (st->tcx_cfg.tcx_last_overlap_mode != TRANSITION_OVERLAP) )
        {
            wtda( st->new_speech_TCX, tcx20Win, NULL, overlap_mode[frameno], overlap_mode[frameno+1], L_frameTCX );

            /* Windowing of the 2xTCX5 subframes or 1xTCX10 or 1xTCX20 */
            WindowSignal( &st->tcx_cfg, st->tcx_cfg.tcx_offsetFB, overlap_mode[frameno] == ALDO_WINDOW ? FULL_OVERLAP : overlap_mode[frameno],
                          overlap_mode[frameno+1] == ALDO_WINDOW ? FULL_OVERLAP : overlap_mode[frameno+1], &left_overlap, &right_overlap,
                          &st->speech_TCX[frameno*tcx10SizeFB], &L_subframe, mdstWin, 1 );
        }
        else
        {
            /* Windowing of the 2xTCX5 subframes or 1xTCX10 or 1xTCX20 */
            WindowSignal( &st->tcx_cfg, st->tcx_cfg.tcx_offsetFB, overlap_mode[frameno], overlap_mode[frameno+1], &left_overlap, &right_overlap,
                          &st->speech_TCX[frameno*tcx10SizeFB], &L_subframe, tcx20Win, 1 );
        }

        if( transform_type[frameno] == TCX_5 )
        {
            /* Outter left folding */
            for (i = 0; i < left_overlap/2; i++)
            {
                tcx20Win[left_overlap/2+i] -= tcx20Win[left_overlap/2-1-i];
            }

            /* Outter right folding */
            for (i = 0; i < right_overlap/2; i++)
            {
                tcx20Win[L_subframe+left_overlap/2-1-i] += tcx20Win[L_subframe+left_overlap/2+i];
            }

            /* 2xTCX5 */
            L_subframe = tcx5SizeFB;
            folding_offset = left_overlap/2;

            for (i = 0; i < 2; i++)
            {
                WindowSignal( &st->tcx_cfg, folding_offset, i == 0 ? RECTANGULAR_OVERLAP : MIN_OVERLAP, i == 1 ? RECTANGULAR_OVERLAP : MIN_OVERLAP,
                              &left_overlap, &right_overlap, tcx20Win+i*tcx5SizeFB, &L_subframe, tcx5Win, 1 );

                TCX_MDCT( tcx5Win, st->spectrum[frameno]+i*tcx5SizeFB, left_overlap, L_subframe-(left_overlap+right_overlap)/2, right_overlap );
            }
        }
        else /* transform_type[frameno] != TCX_5 */
        {
            assert(transform_type[frameno] == TCX_10 || transform_type[frameno] == TCX_20);

            if( transform_type[frameno] == TCX_20 && st->tcx_cfg.tcx_last_overlap_mode != TRANSITION_OVERLAP )
            {
                edct(tcx20Win, st->spectrum[frameno], L_subframe);
                v_multc(st->spectrum[frameno], (float)sqrt((float)NORM_MDCT_FACTOR / L_subframe), st->spectrum[frameno], L_subframe);
            }
            else
            {
                /* TCX20/TCX10 */
                TCX_MDCT( tcx20Win, st->spectrum[frameno], left_overlap, L_subframe-(left_overlap+right_overlap)/2, right_overlap );
            }

            /* For TCX20 at bitrates up to 64 kbps we need the power spectrum */
            if ((st->tcxMode == TCX_20) && ((st->total_brate < HQ_96k) || st->igf))
            {
                /* Compute noise-measure flags for spectrum filling and quantization */
                AnalyzePowerSpectrum( st, L_subframe*st->L_frame/st->L_frameTCX, L_subframe, left_overlap, right_overlap, st->spectrum[frameno],
                                      ((st->tcxMode == TCX_20) && (st->tcx_cfg.tcx_last_overlap_mode != TRANSITION_OVERLAP)) ? mdstWin : tcx20Win, powerSpec );
            }
        }

        TNSAnalysis( &st->tcx_cfg, L_frameTCX, st->tcx_cfg.tcx_coded_lines, transform_type[frameno], (frameno == 0) && (st->last_core == 0),
                     st->spectrum[frameno], &st->tnsData[frameno], &st->fUseTns[frameno], NULL );

        EncodeTnsData(st->tcx_cfg.pCurrentTnsConfig, &st->tnsData[frameno],
                      param_core+frameno*NPRM_DIV+1+NOISE_FILL_RANGES+LTPSIZE, pTnsSize+frameno, pTnsBits+frameno);

        if( transform_type[frameno] == TCX_5 )
        {
            /* group sub-windows: interleave bins according to their frequencies */
            for( i = 0; i < tcx5SizeFB; i++ )
            {
                interleaveBuf[2*i] = st->spectrum[frameno][i];
                interleaveBuf[2*i+1] = st->spectrum[frameno][tcx5SizeFB+i];
            }

            mvr2r( interleaveBuf, st->spectrum[frameno], tcx10SizeFB );
        }

        /*--------------------------------------------------------------*
        * LPC analysis
        *---------------------------------------------------------------*/

        HBAutocorrelation( &st->tcx_cfg, overlap_mode[frameno]==ALDO_WINDOW?FULL_OVERLAP:overlap_mode[frameno],
                           overlap_mode[frameno+1]==ALDO_WINDOW?FULL_OVERLAP:overlap_mode[frameno+1],
                           &st->speech_enc_pe[frameno*tcx10Size], L_frame/nSubframes, r, M );

        adapt_lag_wind( r, M, alw_pitch_lag_12k8[frameno], alw_voicing[frameno], st->sr_core );

        lev_dur( A, r, M, NULL );

        a2lsp_stab( A, lsp[nSubframes-1-frameno], st->lspold_enc );

        if( st->igf )
        {
            ProcessIGF(&st->hIGFEnc, st, st->spectrum[frameno], powerSpec, transform_type[frameno] == TCX_20, st->fUseTns[frameno], (st->last_core == ACELP_CORE), frameno );
        }
    }

    /* Copy memory */
    mvr2r( lsp_new, st->lspold_enc, M );

    return;
}


