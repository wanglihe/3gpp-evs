/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"

/*-------------------------------------------------------------------*
 * core_coder_mode_switch()
 *
 *
 *-------------------------------------------------------------------*/

void core_coder_mode_switch(
    Encoder_State *st,
    int bandwidth_in,
    int bitrate
)
{
    int i, fscale, switchWB, sr_core;
    int bSwitchFromAmrwbIO;
    int bandwidth;

    bandwidth = bandwidth_in;

    switchWB = 0;
    bSwitchFromAmrwbIO = 0;
    if( st->last_core == AMR_WB_CORE )
    {
        bSwitchFromAmrwbIO = 1;
    }

    /* force active frame for the first frame when switching from high bitrates when dtx is enabled*/
    sr_core = getCoreSamplerateMode2(bitrate, bandwidth, st->rf_mode);
    fscale = sr2fscale(sr_core);

    if ( (bandwidth >= WB) && (fscale==(FSCALE_DENOM*16000)/12800) && (fscale == st->fscale) )
    {
        if ( ((bitrate>32000) && (st->tcxonly==0)) || ((bitrate<=32000) && (st->tcxonly==1)) )
        {
            switchWB = 1;
        }
    }
    if( st->last_codec_mode == MODE1 )
    {
        switchWB = 1;  /*force init when coming from MODE1*/
    }

    if( st->last_total_brate > ACELP_32k && st->total_brate <= ACELP_32k )
    {
        switchWB = 1;  /*force init when coming from MODE1*/
    }

    if( fscale == st->fscale && !bSwitchFromAmrwbIO && !switchWB )
    {

        st->total_brate = bitrate;
        st->sr_core = sr_core;
        st->L_frame = sr_core / 50;
        st->tcxonly = getTcxonly(st->total_brate);
        st->bits_frame_nominal = (int)( (float)st->L_frame/(float)st->fscale * (float)FSCALE_DENOM/128.0f * (float)st->total_brate/100.0f + 0.49f );
        st->igf = getIgfPresent(bitrate, bandwidth, st->rf_mode );

        /* switch IGF configuration */
        if (st->igf)
        {
            IGFEncSetMode( &st->hIGFEnc, st->total_brate, bandwidth, st->rf_mode );
        }

        st->tcx_cfg.tcx_coded_lines = getNumTcxCodedLines(bandwidth);
        st->tcx_cfg.bandwidth = getTcxBandwidth(bandwidth);
        st->tcx_cfg.tcxRateLoopOpt = (st->tcxonly) ? 2 : 0;
        st->tcx_cfg.ctx_hm = getCtxHm(st->total_brate, st->rf_mode );
        st->tcx_cfg.resq   = getResq(st->total_brate);
        st->tcx_lpc_shaped_ari = getTcxLpcShapedAri( st->total_brate, st->bwidth, st->rf_mode );

        st->tcx_cfg.tcxRateLoopOpt = (st->tcx_cfg.resq && !st->tcxonly) ? 1 : st->tcx_cfg.tcxRateLoopOpt;
        st->tcx_cfg.fIsTNSAllowed = getTnsAllowed(st->total_brate, st->igf);

        if (st->tcx_cfg.fIsTNSAllowed)
        {
            InitTnsConfigs( bwMode2fs[bandwidth], st->tcx_cfg.tcx_coded_lines, st->tcx_cfg.tnsConfig, (&st->hIGFEnc)->infoStopFrequency, st->total_brate);
        }

        if( bandwidth == NB )
        {
            st->narrowBand = 1;
            st->min_band = 1;
            st->max_band = 16;
        }
        else
        {
            st->narrowBand = 0;
            st->min_band = 0;
            st->max_band = 19;
        }

        for (i=0; i<FRAME_SIZE_NB; i++)
        {
            if (FrameSizeConfig[i].frame_bits==st->bits_frame_nominal)
            {
                st->frame_size_index = i;
                st->bits_frame = FrameSizeConfig[i].frame_bits;
                st->bits_frame_core = FrameSizeConfig[i].frame_net_bits;
                break;
            }
        }

        st->restrictedMode = getRestrictedMode( st->total_brate, 0 );

        core_coder_reconfig( st );
    }
    else
    {
        st->igf = getIgfPresent(bitrate, bandwidth, st->rf_mode );
        init_coder_ace_plus( st );
    }

    if( st->igf )
    {
        /* reset TBE */
        if( ( st->bwidth == WB && st->last_extl != WB_TBE ) ||
                ( st->bwidth == SWB && st->last_extl != SWB_TBE ) ||
                ( st->bwidth == FB && st->last_extl != FB_TBE ) )
        {
            TBEreset_enc( st, st->bwidth );
        }
        else
        {
            set_f( st->state_lpc_syn, 0.0f, LPC_SHB_ORDER );
            set_f( st->state_syn_shbexc, 0.0f, L_SHB_LAHEAD );
            set_f( st->mem_stp_swb, 0.0f, LPC_SHB_ORDER );
            set_f( st->mem_zero_swb, 0, LPC_SHB_ORDER );
            st->gain_prec_swb = 1.0f;
        }
    }

    if (st->envWeighted && !st->enableTcxLpc)
    {
        /* Unweight the envelope */
        E_LPC_lsp_unweight( st->lsp_old, st->lsp_old, st->lsf_old, 1.0f/st->gamma );
        st->envWeighted = 0;
    }

    if( bitrate >= HQ_48k )
    {
        st->enablePlcWaveadjust = 1;
    }
    else
    {
        st->enablePlcWaveadjust = 0;
    }

    if( (st->last_total_brate > HQ_32k) || (st->last_codec_mode == MODE1) )
    {
        st->glr_reset = 1;
    }

    return;
}
