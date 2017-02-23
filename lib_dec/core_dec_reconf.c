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
#include "rom_dec.h"

/*---------------------------------------------------------------------*
 * reconfig_decoder_LPD()
 *
 *
 *---------------------------------------------------------------------*/

void reconfig_decoder_LPD(
    Decoder_State *st,
    int bits_frame,
    int bandwidth,
    int bitrate,
    int L_frame_old
)
{
    short i;

    st->bits_frame = bits_frame;

    if( bandwidth == NB )
    {
        st->narrowBand = 1;
    }
    else if( bandwidth > NB )
    {
        st->narrowBand = 0;
    }

    BITS_ALLOC_init_config_acelp(bitrate, st->narrowBand, st->nb_subfr, &(st->acelp_cfg));

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

    st->flag_cna = getCnaPresent(bitrate, bandwidth);

    /* TCX-LTP */
    st->tcxltp = getTcxLtp(st->sr_core);

    /*Scale TCX for non-active frames to adjust loudness with ACELP*/
    st->tcx_cfg.na_scale = 1.f;

    if( bandwidth < SWB && !(st->tcxonly) )
    {
        int scaleTableSize = sizeof (scaleTcxTable) / sizeof (scaleTcxTable[0]);
        for (i = 0 ; i < scaleTableSize ; i++)
        {
            if ( (bandwidth == scaleTcxTable[i].bwmode) &&
                    (bitrate >= scaleTcxTable[i].bitrateFrom) &&
                    (bitrate < scaleTcxTable[i].bitrateTo) )
            {
                if( st->rf_flag )
                {
                    i--;
                }
                st->tcx_cfg.na_scale = scaleTcxTable[i].scale;
                break;
            }
        }
    }

    /*if its not the first frame resample overlap buffer to new sampling rate */
    if( st->ini_frame != 0 )
    {
        if( st->fscale!=st->fscale_old
                && ! (st->last_codec_mode == MODE1
                      && st->last_core == ACELP_CORE
                      && st->prev_bfi != 0))
            /* no resempling is needed here when recovering from mode 1
               acelp plc, since the buffers are already sampled with the
               correct sampling rate in open_decoder_LPD() */
        {
            unsigned short newLen;
            unsigned short oldLen;

            newLen = st->tcx_cfg.tcx_mdct_window_length;
            oldLen = st->tcx_cfg.tcx_mdct_window_length_old;

            if( (st->prev_bfi && st->last_core_bfi == ACELP_CORE) || st->last_core == ACELP_CORE )
            {
                newLen = st->L_frame/2;
                oldLen = L_frame_old/2;
            }

            lerp( st->old_syn_Overl, st->old_syn_Overl, newLen, oldLen );
            lerp( st->syn_Overl,     st->syn_Overl,     newLen, oldLen );

            if( st->prev_bfi && st->last_core_bfi == ACELP_CORE )
            {
                lerp( st->syn_Overl_TDAC, st->syn_Overl_TDAC, newLen, oldLen );
            }
        }

        if (st->L_frame <= L_FRAME16k)
        {
            if( st->last_L_frame <= L_FRAME16k )
            {
                if( st->L_frame!=st->last_L_frame )
                {
                    unsigned short oldLenClasBuff;
                    unsigned short newLenClasBuff;

                    if( st->L_frame > st->last_L_frame )
                    {
                        oldLenClasBuff = L_SYN_MEM_CLAS_ESTIM * st->last_L_frame/st->L_frame;
                        newLenClasBuff = L_SYN_MEM_CLAS_ESTIM;
                    }
                    else
                    {
                        oldLenClasBuff = L_SYN_MEM_CLAS_ESTIM;
                        newLenClasBuff = L_SYN_MEM_CLAS_ESTIM * st->L_frame/st->last_L_frame;
                    }
                    lerp( &st->mem_syn_clas_estim[L_SYN_MEM_CLAS_ESTIM-oldLenClasBuff], &st->mem_syn_clas_estim[L_SYN_MEM_CLAS_ESTIM-newLenClasBuff], newLenClasBuff, oldLenClasBuff );
                }
            }
            else
            {
                set_zero( st->mem_syn_clas_estim, L_SYN_MEM_CLAS_ESTIM );
            }
        }
    }

    st->enableTcxLpc = (st->numlpc == 1) && (st->lpcQuantization == 1) && (bitrate <= LOWRATE_TCXLPC_MAX_BR || st->rf_flag);

    if( st->ini_frame == 0 )
    {
        st->envWeighted = 0;
    }

    return;
}
