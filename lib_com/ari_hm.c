/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include <math.h>
#include "cnst.h"
#include "stl.h"
#include "basop_util.h"
#include "prot.h"
#include "rom_com.h"


/*-------------------------------------------------------------------*
 * UnmapIndex()
 *
 *
 *-------------------------------------------------------------------*/

void UnmapIndex(
    int PeriodicityIndex,
    int Bandwidth,
    short LtpPitchLag,
    int SmallerLags,
    int *FractionalResolution,
    int *Lag)
{
    if ((LtpPitchLag > 0) && (PeriodicityIndex & kLtpHmFlag))
    {
        int LtpPitchIndex, Multiplier;
        LtpPitchIndex = PeriodicityIndex >> 9;
        Multiplier = PeriodicityIndex & 0xff;
        assert(0 <= LtpPitchIndex && LtpPitchIndex <= 16);
        assert(1 <= Multiplier && Multiplier <= (1 << NumRatioBits[Bandwidth][LtpPitchIndex]));
        *FractionalResolution = kLtpHmFractionalResolution;
        *Lag = (LtpPitchLag * (int)(4 * Ratios[Bandwidth][LtpPitchIndex][Multiplier-1])) >> 2;
    }
    else
    {
        if (PeriodicityIndex < 16)
        {
            *FractionalResolution = 3;
            *Lag = PeriodicityIndex + GET_ADJ(0, 6);
        }
        else if (PeriodicityIndex < 80)
        {
            *FractionalResolution = 4;
            *Lag = PeriodicityIndex + GET_ADJ(16, 8);
        }
        else if (PeriodicityIndex < 208)
        {
            *FractionalResolution = 3;
            *Lag = PeriodicityIndex + GET_ADJ(80, 12);
        }
        else if (PeriodicityIndex < 224 || SmallerLags)
        {
            *FractionalResolution = 1;
            *Lag = PeriodicityIndex + GET_ADJ(208, 28);
        }
        else
        {
            *FractionalResolution = 0;
            *Lag = PeriodicityIndex + GET_ADJ(224, 188);
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * ConfigureContextHm()
 *
 *
 *-------------------------------------------------------------------*/

void ConfigureContextHm(
    int NumCoeffs,                      /* (I) Number of coefficients                         */
    int TargetBits,                     /* (I) Target bit budget (excl. Done flag)            */
    int PeriodicityIndex,               /* (I) Pitch related index                            */
    short LtpPitchLag,                  /* (I) TCX-LTP pitch in F.D.                          */
    CONTEXT_HM_CONFIG *hm_cfg           /* (O) Context-based harmonic model configuration     */
)
{
    int     Bandwidth, SmallerLags;
    int     i, Limit, Lag;
    int     j, Index, FractionalResolution;
    int     *tmp;

    Bandwidth = 0;
    if (NumCoeffs >= 256)
    {
        Bandwidth = 1;
    }

    SmallerLags = 0;

    if (TargetBits <= kSmallerLagsTargetBitsThreshold || Bandwidth == 0)
    {
        SmallerLags = 1;
    }

    UnmapIndex(PeriodicityIndex, Bandwidth, LtpPitchLag, SmallerLags, &FractionalResolution, &Lag );

    /* Set up and fill peakIndices */
    hm_cfg->peakIndices = hm_cfg->indexBuffer;
    tmp = hm_cfg->peakIndices;
    Limit = (NumCoeffs - 1) << FractionalResolution;

    for (i=Lag; i<Limit; i+=Lag)
    {
        Index = i >> FractionalResolution;
        *tmp++ = Index - 1;
        *tmp++ = Index;
        *tmp++ = Index + 1;
    }
    hm_cfg->numPeakIndices = tmp - hm_cfg->indexBuffer;

    /* Set up and fill holeIndices */
    hm_cfg->holeIndices = hm_cfg->indexBuffer + hm_cfg->numPeakIndices;
    tmp = hm_cfg->holeIndices;
    Index = 0;

    for (j=0; j<hm_cfg->numPeakIndices; j+=3)
    {
        for (; Index<hm_cfg->peakIndices[j]; ++Index)
        {
            *tmp++ = Index;
        }
        Index += 3; /* Skip the peak */
    }

    for (; Index<NumCoeffs; ++Index)
    {
        *tmp++ = Index;
    }
    hm_cfg->numHoleIndices = tmp - hm_cfg->holeIndices;
    /* Add extremal element signaling the end of the buffer */
    *tmp++ = NumCoeffs;

    return;
}


/*-------------------------------------------------------------------*
 * CountIndexBits()
 *
 *
 *-------------------------------------------------------------------*/

int CountIndexBits(
    int Bandwidth,
    int PeriodicityIndex)
{

    if (PeriodicityIndex & kLtpHmFlag)
    {
        int LtpPitchIndex = PeriodicityIndex >> 9;
        return NumRatioBits[Bandwidth][LtpPitchIndex];
    }
    return 8;
}


/*-------------------------------------------------------------------*
 * tcx_hm_render()
 *
 *
 *-------------------------------------------------------------------*/

int tcx_hm_render(
    int lag,           /* i: pitch lag                             */
    int fract_res,     /* i: fractional resolution of the lag      */
    float LtpGain,     /* i: LTP gain                              */
    Word16 p[]         /* o: harmonic model (Q13)                  */
)
{
    int k;
    Word32 f0, tmp32;
    Word16 height, PeakDeviation, tmp;

    /* Set up overall shape */
    (void)LtpGain;

    f0 = L_shl(lag, sub(15, fract_res)); /* Q15 */

    tmp32 = Mpy_32_16(f0, -26474);
    tmp32 = L_shr_r(BASOP_Util_InvLog2(L_shl(tmp32, 7)), 2);
    tmp32 = L_sub(603979776L, tmp32);
    tmp32 = L_add(L_add(tmp32, tmp32), Mpy_32_16(tmp32, 26214));
    height = round_fx(tmp32); /* Q13 */

    tmp32 = Mpy_32_16(f0, -18910);
    tmp32 = L_shr_r(BASOP_Util_InvLog2(L_shl(tmp32, 7)), 2);
    tmp32 = L_sub(1395864371L, tmp32);
    PeakDeviation = round_fx(tmp32); /* Q14 */

    IF( sub(13915,PeakDeviation) > 0 )
    {
        /* A bit error was encountered */
        return 1;
    }
    ELSE
    {
        tmp = div_s(13915, PeakDeviation);
        tmp = mult_r(tmp, tmp); /* Q15 */
    }

    /* Render the prototype peak */
    p[kTcxHmParabolaHalfWidth] = height;

    for (k=1; k<=kTcxHmParabolaHalfWidth; ++k)
    {
        p[kTcxHmParabolaHalfWidth+k] = round_fx(Mpy_32_16(BASOP_Util_InvLog2(L_shl(L_mult0(mult0(negate(k),k), tmp),10)), height));
    }
    /* Mirror */
    for (k=-kTcxHmParabolaHalfWidth; k<0; ++k)
    {
        p[kTcxHmParabolaHalfWidth+k] = p[kTcxHmParabolaHalfWidth-k];
    }

    return 0;
}


/*-------------------------------------------------------------------*
 * tcx_hm_modify_envelope()
 *
 *
 *-------------------------------------------------------------------*/

void tcx_hm_modify_envelope(
    Word16 gain,       /* i: HM gain (Q11)                         */
    int lag,
    int fract_res,
    Word16 p[],        /* i: harmonic model (Q13)                  */
    Word32 env[],      /* i/o: envelope (Q16)                      */
    int L_frame        /* i: number of spectral lines              */
)
{
    int k;
    int h, x;
    Word16 inv_shape[2*kTcxHmParabolaHalfWidth+1]; /* Q15 */

    if (gain == 0)
    {
        return;
    }

    for (k=0; k<2*kTcxHmParabolaHalfWidth+1; ++k)
    {
        inv_shape[k] = div_s(512, add(512, round_fx(L_mult(gain, p[k]))));
    }

    h = 1;
    k = lag >> fract_res;

    while (k <= L_frame + kTcxHmParabolaHalfWidth - 1)
    {

        for (x=max(0, k-kTcxHmParabolaHalfWidth); x<=min(k+kTcxHmParabolaHalfWidth, L_frame-1); ++x)
        {
            env[x] = Mpy_32_16(env[x], inv_shape[x-k+kTcxHmParabolaHalfWidth]);
        }
        ++h;
        k = (h * lag) >> fract_res;
    }

    return;
}
