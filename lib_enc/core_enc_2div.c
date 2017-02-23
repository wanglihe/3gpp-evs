/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*-------------------------------------------------------------------*
 * core_encode_twodiv()
 *
 * Two-div core encoder
 *-------------------------------------------------------------------*/

void core_encode_twodiv(
    const float new_samples[],          /* i  : new samples                         */
    Encoder_State *st,                    /* i/o : coder memory state                 */
    const short coder_type,             /* i  : coding type                         */
    const short pitch[3],                /* i  : open-loop pitch values for quantiz. */
    const float voicing[3],             /* i  : open-loop pitch gains               */
    float Aw[NB_SUBFR16k*(M+1)]   /* i  : weighted A(z) unquant. for subframes*/
)
{
    short n;
    float lsp_new[M], lsp_mid[M];
    float lsf_q[M], lsp_q[M];
    float lspmid_q[M];
    float A_q[M+1];
    int param_lpc[NPRM_LPC_NEW];
    int nbits_lpc[2];
    int param_core[2*NPRM_DIV];
    int target_bits;
    float gainlpc[2][FDNS_NPTS];
    int tnsSize[2];   /* number of tns parameters put into prm */
    int tnsBits[2];   /* number of tns bits in the frame */
    int ltpBits;
    int bitsAvailable;
    int indexBuffer[2*((N_MAX/2)+1)];
    CONTEXT_HM_CONFIG hm_cfg[2];
    short bits_param_lpc[10], no_param_lpc;
    short i, T_op[3];


    hm_cfg[0].indexBuffer = &indexBuffer[0];
    hm_cfg[1].indexBuffer = &indexBuffer[N_MAX/2+1];

    set_i( tnsSize, 0, 2 );
    set_i( tnsBits, 0, 2 );
    ltpBits = 0;

    for( i = 0; i < 3; i++ )
    {
        T_op[i] = pitch[i];

        /* check minimum pitch for quantization */
        if( T_op[i] < PIT_MIN_SHORTER )
        {
            T_op[i] *= 2;
        }

        /* convert pitch values to core sampling-rate */
        if ( st->L_frame != L_FRAME )
        {
            T_op[i] = (short)(T_op[i] * (float)st->L_frame/(float)L_FRAME + 0.5f);
        }
    }


    /*--------------------------------------------------------------*
    * TCX20/TCX10 switching decision
    *---------------------------------------------------------------*/

    if ( st->tcxMode == TCX_10 )
    {
        st->core = TCX_10_CORE;
    }
    else if ( st->tcxMode == TCX_20 )
    {
        st->core = TCX_20_CORE;
    }


    /*--------------------------------------------------------------*
    * Core Signal Analysis: MDCT, TNS, LPC analysis
    *---------------------------------------------------------------*/

    core_signal_analysis_high_bitrate( new_samples, T_op, voicing, pitch, lsp_new, lsp_mid, st,
                                       tnsSize, tnsBits, param_core, &ltpBits, st->L_frame,st->L_frameTCX );

    /*--------------------------------------------------------------*
    * LPC Quantization
    *---------------------------------------------------------------*/
    lpc_quantization( st, st->core, st->lpcQuantization, st->lsf_old, lsp_new, lsp_mid, lsp_q, lsf_q,
                      lspmid_q, st->mem_MA, st->mem_AR, st->narrowBand, coder_type,
                      0, /*No acelp->no need to compute any mid-LPC*/
                      param_lpc, nbits_lpc, &(st->seed_acelp), st->sr_core, st->Bin_E, st->Bin_E_old,
                      bits_param_lpc, &no_param_lpc );

    /*--------------------------------------------------------------*
     * Rate switching
     *--------------------------------------------------------------*/

    if( st->rate_switching_reset )
    {
        mvr2r( lsp_q, st->lsp_old, M );
        mvr2r( lsf_q, st->lsf_old, M );
    }



    /*--------------------------------------------------------------*
    * Run Two TCX10
    *---------------------------------------------------------------*/

    if( st->core == TCX_10_CORE )
    {
        const int last_ace_mode = st->last_core;

        for (n = 0; n < 2; n++)
        {
            if ( n == 0 )
            {
                lsp2a_stab( lspmid_q, A_q, M );
            }
            else
            {
                lsp2a_stab(lsp_q, A_q, M);
            }

            /* Shape spectrum */
            ShapeSpectrum( &(st->tcx_cfg), A_q, gainlpc[n], st->L_frame/2, st->tcx_cfg.tcx_coded_lines/2, st->spectrum[n], st->fUseTns[n], st );

            st->last_core = st->core;
        }
        st->last_core = last_ace_mode;

        /* Calculate target bits */
        bitsAvailable = st->bits_frame_core - nbits_lpc[0] - nbits_lpc[1] - st->nb_bits_header_tcx;

        /* subtract bits for TCX overlap mode (1 bit: full, 2 bits: half or no overlap) */
        bitsAvailable -= (st->tcx_cfg.tcx_curr_overlap_mode == HALF_OVERLAP || st->tcx_cfg.tcx_curr_overlap_mode == MIN_OVERLAP) ? 2 : 1;
        bitsAvailable -= (&st->hIGFEnc)->infoTotalBitsWritten;

        st->measuredBwRatio = 1.f;

        for( n = 0; n < 2; n++ )
        {
            target_bits = (bitsAvailable + 1 - n)/2 - tnsBits[n];

            if( n == 0 )
            {
                target_bits -= ltpBits;
            }

            if(st->enablePlcWaveadjust && n)
            {
                target_bits  -= 1;
            }

            /* Run TCX10 encoder */
            QuantizeSpectrum( &(st->tcx_cfg), A_q, NULL, gainlpc[n], st->synth+n*st->L_frame/2, st->L_frame/2, st->L_frameTCX/2, st->tcx_cfg.tcx_coded_lines/2,
                              target_bits, st->tcxonly, st->spectrum[n], st->tnsData+n, st->fUseTns[n], tnsSize[n], &(st->LPDmem), param_core+n*NPRM_DIV, n, st, &hm_cfg[n] );

            /* Update tcx overlap mode */
            if( (n > 0) || !st->tcxonly )
            {
                st->tcx_cfg.tcx_last_overlap_mode = st->tcx_cfg.tcx_curr_overlap_mode;
            }


        }

        coder_tcx_post( st, &(st->LPDmem), &(st->tcx_cfg), st->synth, A_q, Aw, st->wspeech_enc );
    }

    /*--------------------------------------------------------------*
    * Run One TCX20
    *---------------------------------------------------------------*/

    if( st->core == TCX_20_CORE )
    {
        lsp2a_stab( lsp_q, A_q, M );

        ShapeSpectrum( &(st->tcx_cfg), A_q, gainlpc[0], st->L_frame, st->tcx_cfg.tcx_coded_lines, st->spectrum_long, st->fUseTns[0], st );

        st->measuredBwRatio = 1.f;

        /* Calculate target bits */
        target_bits = st->bits_frame_core - tnsBits[0] - nbits_lpc[0] - st->nb_bits_header_tcx - ltpBits;

        /* subtract bits for TCX overlap mode (1 bit: full, 2 bits: half or no overlap) */
        target_bits -= (st->tcx_cfg.tcx_curr_overlap_mode == HALF_OVERLAP || st->tcx_cfg.tcx_curr_overlap_mode == MIN_OVERLAP) ? 2 : 1;

        target_bits -= st->hIGFEnc.infoTotalBitsPerFrameWritten;

        if( st->enablePlcWaveadjust )
        {
            target_bits -= 1;
        }

        QuantizeSpectrum( &(st->tcx_cfg), A_q, NULL, gainlpc[0],st->synth, st->L_frame, st->L_frameTCX,st->tcx_cfg.tcx_coded_lines, target_bits,
                          st->tcxonly, st->spectrum_long, st->tnsData, st->fUseTns[0], tnsSize[0], &(st->LPDmem), param_core, 0, st, &hm_cfg[0]);

        coder_tcx_post( st, &(st->LPDmem), &(st->tcx_cfg), st->synth, A_q, Aw, st->wspeech_enc );

    }

    /* Update lsp/lsf memory */
    mvr2r( lsf_q, st->lsf_old, M );
    mvr2r( lsp_q, st->lsp_old, M );



    /*--------------------------------------------------------------*
    * Generate Bitstream
    *---------------------------------------------------------------*/

    enc_prm( coder_type, param_core, param_lpc, st,  st->L_frame, hm_cfg, bits_param_lpc, no_param_lpc );

    return;

}
