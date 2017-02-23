/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "cnst.h"
#include "stl.h"
#include "basop_util.h"
#include <assert.h>
#include <math.h>
#include "prot.h"
#include "rom_com.h"

/*-------------------------------------------------------------------*
 * DecodeIndex()
 *
 *
 *-------------------------------------------------------------------*/

int DecodeIndex(
    Decoder_State *st,
    int Bandwidth,
    int *PeriodicityIndex)
{
    if ((st->tcx_hm_LtpPitchLag > 0) && (st->tcxltp_gain > kLtpHmGainThr))
    {
        int LtpPitchIndex = ((st->tcx_hm_LtpPitchLag + (1 << (kLtpHmFractionalResolution - 1))) >> kLtpHmFractionalResolution) - 2;
        *PeriodicityIndex = kLtpHmFlag;
        *PeriodicityIndex |= get_next_indice(st, NumRatioBits[Bandwidth][LtpPitchIndex]);
        ++*PeriodicityIndex;
        *PeriodicityIndex |= LtpPitchIndex << 9;
        return NumRatioBits[Bandwidth][LtpPitchIndex];
    }
    else
    {
        *PeriodicityIndex = get_next_indice(st, 8);
        return 8;
    }
}


/*-------------------------------------------------------------------*
 * tcx_hm_dequantize_gain()
 *
 *
 *-------------------------------------------------------------------*/

static int tcx_hm_dequantize_gain(
    int coder_type,    /* i: GC/VC mode                            */
    int gain_idx,      /* i: quantization index                    */
    Word16 *gain       /* o: dequantized gain (Q11)                */
)
{
    assert(0 <= coder_type && coder_type <= 1);

    /* safety check in case of bit errors */
    if( !(0 <= gain_idx && gain_idx < (1 << kTcxHmNumGainBits)) )
    {
        *gain = 0;
        return 1;
    }

    *gain = qGains[coder_type][gain_idx];

    return 0;
}


/*-------------------------------------------------------------------*
 * tcx_hm_decode()
 *
 *
 *-------------------------------------------------------------------*/

void tcx_hm_decode(
    int L_frame,       /* i: number of spectral lines              */
    Word32 env[],      /* i/o: envelope shape (Q16)                */
    int targetBits,    /* i: target bit budget                     */
    int coder_type,    /* i: GC/VC mode                            */
    const int prm_hm[],/* i: HM parameters                         */
    short LtpPitchLag, /* i: LTP pitch lag or -1 if none           */
    float LtpGain,     /* i: LTP gain                              */
    int *hm_bits       /* o: bit consumption                       */
)
{
    int NumTargetBits, fract_res, lag;
    Word16 p[2*kTcxHmParabolaHalfWidth+1], gain;

    *hm_bits = 0;

    if( !(coder_type == VOICED || coder_type == GENERIC) )
    {
        /* A bit error was encountered */
        *hm_bits = -1;
        return;
    }

    NumTargetBits = CountIndexBits( L_frame >= 256, prm_hm[1]) + targetBits;

    if (coder_type == VOICED)
    {
        NumTargetBits += kTcxHmNumGainBits;
    }
    *hm_bits = NumTargetBits - targetBits + 1;

    /* Convert the index to lag */
    UnmapIndex( prm_hm[1], L_frame >= 256, LtpPitchLag, (NumTargetBits <= kSmallerLagsTargetBitsThreshold) || (L_frame < 256), &fract_res, &lag );

    /* Render the harmonic model */
    if( tcx_hm_render( lag, fract_res, LtpGain, p ) )
    {
        /* A bit error was encountered */
        *hm_bits = -1;
        return;
    }

    /* Dequantize gain */

    if( tcx_hm_dequantize_gain( coder_type==VOICED,prm_hm[2],&gain ) )
    {
        /* A bit error was encountered */
        *hm_bits = -1;
        return;
    }
    tcx_hm_modify_envelope( gain, lag, fract_res, p, env, L_frame );

    return;
}
