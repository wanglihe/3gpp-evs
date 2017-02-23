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
#include "rom_enc.h"


/*-----------------------------------------------------------------*
 * Funtion  core_coder_reconfig                                    *
 *          ~~~~~~~~~~~~~~~~~~~                                    *
 * - reconfig core coder when switching to another frame type      *
 *-----------------------------------------------------------------*/

void core_coder_reconfig(
    Encoder_State *st
)
{
    short bandwidth, i;

    /*Configuration of ACELP*/
    BITS_ALLOC_init_config_acelp( st->total_brate, st->narrowBand, st->nb_subfr, &(st->acelp_cfg) );

    /*Configuration of partial copy*/
    st->acelp_cfg_rf.mode_index = 1;
    st->acelp_cfg_rf.midLpc = 0;
    st->acelp_cfg_rf.midLpc_enable = 0;
    st->acelp_cfg_rf.pre_emphasis = 0;
    st->acelp_cfg_rf.formant_enh = 1;
    st->acelp_cfg_rf.formant_tilt = 1;
    st->acelp_cfg_rf.voice_tilt = 1;
    st->acelp_cfg_rf.formant_enh_num = FORMANT_SHARPENING_G1;
    st->acelp_cfg_rf.formant_enh_den = FORMANT_SHARPENING_G2;

    if( st->tcxonly )
    {
        st->nb_bits_header_tcx = 1+1;  /*TCX20/TCX10 + last_core*/
        st->nb_bits_header_tcx += 2;   /* Siganl class*/
    }
    else
    {
        st->nb_bits_header_ace = 1+2+1; /*TCX/ACELP+coder_type + last_core*/
        st->nb_bits_header_tcx = st->nb_bits_header_ace;

        if ( st->tcx_cfg.lfacNext<=0 )
        {
            st->nb_bits_header_ace--; /*No last_core*/
        }
    }

    /*Switch off  TCX or ACELP?*/
    if( st->sr_core==12800 )
    {
        st->acelpEnabled =  (st->restrictedMode & 1) == 1;
        st->tcx20Enabled = (st->restrictedMode & 2) == 2;
    }
    st->prevEnergyHF = st->currEnergyHF = 65535.0f; /* prevent block switch */


    /* TCX-LTP */
    st->tcxltp = getTcxLtp(st->sr_core);

    /*Use for 12.8 kHz sampling rate and low bitrates, the conventional pulse search->better SNR*/
    st->acelp_autocorr = 1;
    if( ((st->total_brate <= ACELP_9k60) && (st->sr_core == 12800)) )
    {
        st->acelp_autocorr = 0;
    }

    /*Get bandwidth mode*/
    if( st->narrowBand )
    {
        bandwidth = NB;
    }
    else if(st->sr_core<=16000)
    {
        bandwidth = WB;
    }
    else
    {
        bandwidth = SWB;
    }

    /*Scale TCX for non-active frames to adjust loudness with ACELP*/
    st->tcx_cfg.na_scale = 1.f;
    if( bandwidth < SWB && !(st->tcxonly) )
    {
        int scaleTableSize = sizeof (scaleTcxTable) / sizeof (scaleTcxTable[0]);
        for (i = 0 ; i < scaleTableSize ; i++)
        {

            if ( (bandwidth == scaleTcxTable[i].bwmode) &&
                    (st->total_brate >= scaleTcxTable[i].bitrateFrom) &&
                    (st->total_brate < scaleTcxTable[i].bitrateTo) )
            {
                if( st->rf_mode )
                {
                    i--;
                }
                st->tcx_cfg.na_scale=scaleTcxTable[i].scale;
                break;
            }
        }
    }

    st->enableTcxLpc = (st->lpcQuantization == 1) && (st->total_brate <= LOWRATE_TCXLPC_MAX_BR || st->rf_mode);

    if( st->ini_frame == 0 || st->last_codec_mode == MODE1 )
    {
        st->envWeighted = 0;
    }

    if( st->bwidth == SWB && (st->total_brate == ACELP_16k40 || st->total_brate == ACELP_24k40) )
    {
        if(st->tec_tfa == 0)
        {
            set_zero(st->tecEnc.loBuffer, MAX_TEC_SMOOTHING_DEG);
        }
        st->tec_tfa = 1;
    }
    else
    {
        st->tec_tfa = 0;
    }

    st->enablePlcWaveadjust = 0;
    if( st->total_brate >= HQ_48k )
    {
        st->enablePlcWaveadjust = 1;
    }

    st->glr = 0;
    if( st->total_brate == ACELP_9k60 || st->total_brate == ACELP_16k40 || st->total_brate == ACELP_24k40 )
    {
        st->glr = 1;
    }

    if( st->glr )
    {
        st->nb_bits_header_ace += G_LPC_RECOVERY_BITS;
    }

    if( st->bwidth == NB || st->bwidth == WB )
    {
        st->nmStartLine = startLineWB[min((sizeof(startLineWB)/sizeof(startLineWB[0])-1), (unsigned)max(0, (st->rf_mode==0) ? st->frame_size_index:st->frame_size_index-1))];
    }
    else /* (st->bwidth == SWB || st->bwidth == FB) */
    {
        st->nmStartLine = startLineSWB[min((sizeof(startLineSWB)/sizeof(startLineSWB[0])-1), (unsigned)max(3, (st->rf_mode==0) ? st->frame_size_index:st->frame_size_index-1) - 3)];
    }

    if( (st->total_brate < ACELP_24k40) && ( (st->total_brate > st->last_total_brate) || (st->last_codec_mode == MODE1) ) )
    {
        /* low-freq memQuantZeros must be reset partially if bitrate increased */
        set_i( st->memQuantZeros, 0, st->nmStartLine );
    }
    else
    {
        if( (st->total_brate >= ACELP_24k40) && (st->total_brate <= ACELP_32k) && (st->last_total_brate >= ACELP_13k20) && (st->last_total_brate < ACELP_24k40) )
        {
            set_i(st->memQuantZeros, 0, st->L_frame);
        }
    }


    return;
}

