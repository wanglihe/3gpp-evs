/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "options.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * getTcxonly()
 *
 *
 *-------------------------------------------------------------------*/

short getTcxonly(
    const int bitrate
)
{
    short tcxonly = 0;

    if( bitrate > ACELP_32k )
    {
        tcxonly = 1;
    }

    return tcxonly;
}


/*-------------------------------------------------------------------*
 * getCtxHm()
 *
 *
 *-------------------------------------------------------------------*/

short getCtxHm(
    const int bitrate,
    const short rf_flag
)
{
    short ctx_hm = 0;
    if( (bitrate > LPC_SHAPED_ARI_MAX_RATE) && (bitrate <= 64000) && !rf_flag)
    {
        ctx_hm = 1;
    }

    return ctx_hm;
}

/*-------------------------------------------------------------------*
 * getResq()
 *
 *
 *-------------------------------------------------------------------*/

short getResq(
    const int bitrate
)
{
    short resq = 0;

    if(bitrate <= HQ_64k)
    {
        resq = 1;
    }

    return resq;
}

/*-------------------------------------------------------------------*
 * getTnsAllowed()
 *
 *
 *-------------------------------------------------------------------*/

short getTnsAllowed(
    const int bitrate,
    const short igf
)
{
    short tnsAllowed = 0;

    if( igf )
    {
        if( bitrate > HQ_16k40 )
        {
            tnsAllowed = 1;
        }
    }
    else if( bitrate > HQ_32k )
    {
        tnsAllowed = 1;
    }

    return tnsAllowed;
}

/*-------------------------------------------------------------------*
 * getRestrictedMode()
 *
 *
 *-------------------------------------------------------------------*/

short getRestrictedMode(
    const int bitrate,
    const short Opt_AMR_WB
)
{
    short restrictedMode = 3;

    if ( !Opt_AMR_WB && (bitrate > HQ_32k) )
    {
        restrictedMode = 6;
    }
    else if ( Opt_AMR_WB )
    {
        restrictedMode = 1;
    }

    return restrictedMode;
}

/*-------------------------------------------------------------------*
 * getMdctWindowLength()
 *
 *
 *-------------------------------------------------------------------*/

short getMdctWindowLength(
    const float fscale
)
{

    short mdctWindowLength;

    mdctWindowLength = (L_LOOK_12k8 * fscale)/FSCALE_DENOM;

    return mdctWindowLength;
}

/*-------------------------------------------------------------------*
 * sr2fscale()
 *
 *
 *-------------------------------------------------------------------*/

short sr2fscale(
    const int sr
)
{

    return (FSCALE_DENOM*sr)/12800;
}

/*-------------------------------------------------------------------*
 * getCoreSamplerateMode2()
 *
 *
 *-------------------------------------------------------------------*/

int getCoreSamplerateMode2(
    const int bitrate,
    const int bandwidth,
    const short rf_mode
)
{
    int sr_core = 0;

    if( bandwidth == NB )
    {
        sr_core = 12800;
    }
    else if( (bandwidth == WB && bitrate < ACELP_13k20) ||
             (bandwidth == SWB && bitrate <= ACELP_13k20) || (rf_mode == 1) )
    {
        sr_core = 12800;
    }
    else if( bandwidth == WB || (bandwidth == SWB && bitrate <= ACELP_32k)
             || (bandwidth == FB && bitrate <= ACELP_32k) )
    {
        sr_core = 16000;
    }
    else if( (bandwidth == SWB || bandwidth == FB) && bitrate <= HQ_64k)
    {
        sr_core = 25600;
    }
    else if( bandwidth == SWB || bandwidth == FB )
    {
        sr_core = 32000;
    }

    return sr_core;
}

/*-------------------------------------------------------------------*
 * getTcxBandwidth()
 *
 *
 *-------------------------------------------------------------------*/

float getTcxBandwidth(
    const int bandwidth
)
{
    float tcxBandwidth = 0.5f;

    if( bandwidth == NB )
    {
        tcxBandwidth = 0.3125f;
    }

    return tcxBandwidth;
}

/*-------------------------------------------------------------------*
 * getIgfPresent()
 *
 *
 *-------------------------------------------------------------------*/

short getIgfPresent(
    const int bitrate,
    const int bandwidth,
    const short rf_mode
)
{
    short igfPresent = 0;

    if( bandwidth == SWB && bitrate >= ACELP_9k60 && bitrate < HQ_96k )
    {
        igfPresent = 1;
    }
    else if( bandwidth == FB && bitrate >= ACELP_16k40 )
    {
        igfPresent = 1;
    }
    else if( bandwidth == WB && bitrate == ACELP_9k60 )
    {
        igfPresent = 1;
    }
    if( ((bandwidth == WB) || (bandwidth == SWB)) && (rf_mode == 1) && (bitrate == ACELP_13k20) )
    {
        igfPresent = 1;
    }

    return igfPresent;
}

/*-------------------------------------------------------------------*
 * getCnaPresent()
 *
 *
 *-------------------------------------------------------------------*/

short getCnaPresent(
    const int bitrate,
    const int bandwidth
)
{
    short flag_cna = 0;

    if( bandwidth == NB && bitrate <= ACELP_13k20 )
    {
        flag_cna = 1;
    }

    if( bandwidth == WB && bitrate <= ACELP_13k20 )
    {
        flag_cna = 1;
    }

    if( bandwidth == SWB && bitrate <= ACELP_13k20 )
    {
        flag_cna = 1;
    }

    return flag_cna;
}

/*-------------------------------------------------------------------*
 * getTcxLtp()
 *
 *
 *-------------------------------------------------------------------*/
short getTcxLtp(
    const int sr_core
)
{
    short tcxltp = 0;

    if ( (sr_core <= 25600) )
    {
        tcxltp = 1;
    }

    return tcxltp;
}

/*-------------------------------------------------------------------*
 * initPitchLagParameters()
 *
 *
 *-------------------------------------------------------------------*/

short initPitchLagParameters(
    const int sr_core,
    int *pit_min,
    int *pit_fr1,
    int *pit_fr1b,
    int *pit_fr2,
    int *pit_max
)
{

    short pit_res_max;


    if (sr_core==12800)
    {
        *pit_min = PIT_MIN_12k8;
        *pit_max = PIT_MAX_12k8;
        *pit_fr2 = PIT_FR2_12k8;
        *pit_fr1 = PIT_FR1_12k8;
        *pit_fr1b = PIT_FR1_8b_12k8;
        pit_res_max = 4;
    }
    else if (sr_core==16000)
    {
        *pit_min = PIT_MIN_16k;
        *pit_max = PIT16k_MAX;
        *pit_fr2 = PIT_FR2_16k;
        *pit_fr1 = PIT_FR1_16k;
        *pit_fr1b = PIT_FR1_8b_16k;
        pit_res_max = 6;
    }
    else if (sr_core==25600)
    {
        *pit_min = PIT_MIN_25k6;
        *pit_max = PIT_MAX_25k6;
        *pit_fr2 = PIT_FR2_25k6;
        *pit_fr1 = PIT_FR1_25k6;
        *pit_fr1b = PIT_FR1_8b_25k6;
        pit_res_max = 4;
    }
    else /* sr_core==32000 */
    {
        *pit_min = PIT_MIN_32k;
        *pit_max = PIT_MAX_32k;
        *pit_fr2 = PIT_FR2_32k;
        *pit_fr1 = PIT_FR1_32k;
        *pit_fr1b = PIT_FR1_8b_32k;
        pit_res_max = 6;
    }

    return pit_res_max;
}

/*-------------------------------------------------------------------*
 * getNumTcxCodedLines()
 *
 *
 *-------------------------------------------------------------------*/

short getNumTcxCodedLines(
    const short bwidth
)
{
    short tcx_coded_lines;

    switch (bwidth)
    {
    case NB:
        tcx_coded_lines = 160;
        break;
    case WB:
        tcx_coded_lines = 320;
        break;
    case SWB:
        tcx_coded_lines = 640;
        break;
    case FB:
        tcx_coded_lines = 960;
        break;
    default:
        tcx_coded_lines = 0;
        break;
    }

    return tcx_coded_lines;
}

/*-------------------------------------------------------------------*
 * getTcxLpcShapedAri()
 *
 *
 *-------------------------------------------------------------------*/

short getTcxLpcShapedAri(
    const int total_brate,
    const short bwidth,
    const short rf_mode
)
{
    short tcx_lpc_shaped_ari = 0;

    (void) bwidth;

    if( total_brate <= LPC_SHAPED_ARI_MAX_RATE || rf_mode )
    {
        tcx_lpc_shaped_ari = 1;
    }

    return tcx_lpc_shaped_ari;
}
