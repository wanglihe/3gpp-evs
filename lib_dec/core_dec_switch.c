/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"


/*-------------------------------------------------------------*
 * mode_switch_decoder_LPD()
 *
 *
 *-------------------------------------------------------------*/

void mode_switch_decoder_LPD(
    Decoder_State *st,
    int bandwidth,
    int bitrate,
    int frame_size_index
)
{
    int fscale, switchWB, sr_core;
    int bSwitchFromAmrwbIO;
    int frame_size;

    switchWB = 0;
    bSwitchFromAmrwbIO = 0;
    if(st->last_core == AMR_WB_CORE)
    {
        bSwitchFromAmrwbIO = 1;
    }
    sr_core = getCoreSamplerateMode2(bitrate, bandwidth, st->rf_flag);
    fscale = sr2fscale(sr_core);

    /* set number of coded lines */
    st->tcx_cfg.tcx_coded_lines = getNumTcxCodedLines(bandwidth);

    if ( (bandwidth>=WB) && (fscale==(FSCALE_DENOM*16000)/12800) && (fscale == st->fscale) )
    {
        if ( ((bitrate > ACELP_32k) && (st->tcxonly==0)) || ((bitrate <= ACELP_32k) && (st->tcxonly==1)) )
        {
            switchWB = 1;
        }
    }

    if( st->last_L_frame > L_FRAME16k && st->total_brate <= ACELP_32k )
    {
        switchWB = 1;  /*force init when coming from MODE1*/
    }

    st->igf = getIgfPresent(bitrate, bandwidth, st->rf_flag );

    st->hIGFDec.infoIGFStopFreq = -1;
    if( st->igf )
    {
        /* switch IGF configuration */
        IGFDecSetMode( &st->hIGFDec, st->total_brate, bandwidth, -1, -1, st->rf_flag );
    }

    if( fscale != st->fscale || switchWB || bSwitchFromAmrwbIO || st->last_codec_mode == MODE1 || st->force_lpd_reset )
    {
        open_decoder_LPD( st, bitrate, bandwidth );
    }
    else
    {
        assert(fscale > (FSCALE_DENOM/2));
        st->fscale_old  = st->fscale;
        st->fscale      = fscale;
        st->L_frame = st->sr_core / 50;
        st->L_frameTCX = st->output_Fs / 50;

        st->tcx_cfg.ctx_hm = getCtxHm(bitrate, st->rf_flag);
        st->tcx_cfg.resq   = getResq(bitrate);

        st->tcx_lpc_shaped_ari = getTcxLpcShapedAri( bitrate, bandwidth, st->rf_flag );

        if ( bandwidth == NB )
        {
            st->narrowBand = 1;
        }
        else
        {
            st->narrowBand = 0;
        }
        st->TcxBandwidth = getTcxBandwidth(bandwidth);

        st->tcx_cfg.pCurrentTnsConfig = NULL;
        st->tcx_cfg.fIsTNSAllowed = getTnsAllowed(st->total_brate, st->igf);
        if( st->tcx_cfg.fIsTNSAllowed )
        {
            InitTnsConfigs( bwMode2fs[bandwidth], st->tcx_cfg.tcx_coded_lines, st->tcx_cfg.tnsConfig, st->hIGFDec.infoIGFStopFreq, st->total_brate );
        }
    }

    frame_size = FrameSizeConfig[frame_size_index].frame_net_bits;

    reconfig_decoder_LPD( st, frame_size, bandwidth, bitrate, st->last_L_frame );

    if (st->envWeighted && !st->enableTcxLpc)
    {
        mvr2r(st->lspold_uw, st->lsp_old, M);
        mvr2r(st->lsfold_uw, st->lsf_old, M);
        st->envWeighted = 0;
    }

    /* update PLC LSF memories */
    lsp2lsf( st->lsp_old, st->lsfoldbfi1, M, st->sr_core );
    mvr2r(st->lsfoldbfi1, st->lsfoldbfi0,M);
    mvr2r(st->lsfoldbfi1, st->lsf_adaptive_mean,M);

    if( st->igf )
    {
        /* reset TBE */
        if( ( st->bwidth == WB && st->last_extl != WB_TBE ) ||
                ( st->bwidth == SWB && st->last_extl != SWB_TBE ) ||
                ( st->bwidth == FB && st->last_extl != FB_TBE ) )
        {
            TBEreset_dec( st, st->bwidth );
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

    if( bandwidth == SWB && (st->total_brate == ACELP_16k40 || st->total_brate == ACELP_24k40) )
    {
        if(st->tec_tfa == 0)
        {
            set_zero(st->tecDec.loBuffer, MAX_TEC_SMOOTHING_DEG);
        }
        st->tec_tfa = 1;
    }
    else
    {
        st->tec_tfa = 0;
    }

    st->tec_flag = 0;
    st->tfa_flag = 0;

    /* needed in decoder to read the bitstream */
    if (
        ( bandwidth == FB  && bitrate == ACELP_24k40 ) ||
        ( bandwidth == WB  && bitrate == ACELP_24k40 ) ||
        ( bandwidth == SWB && bitrate == ACELP_24k40 )
    )
    {
        st->enableGplc = 1;
    }
    else
    {
        st->enableGplc = 0;
    }

    if( st->total_brate == ACELP_9k60 || st->total_brate == ACELP_16k40 || st->total_brate == ACELP_24k40 )
    {
        st->dec_glr = 1;
    }
    else
    {
        st->dec_glr = 0;
    }

    st->dec_glr_idx = 0;


    return;
}
