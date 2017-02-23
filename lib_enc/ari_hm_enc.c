/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "cnst.h"
#include "rom_enc.h"
#include "stl.h"
#include "basop_util.h"
#include <assert.h>
#include <math.h>
#include "prot.h"
#include "rom_com.h"


/*-------------------------------------------------------------------*
 * EncodeIndex()
 *
 *
 *-------------------------------------------------------------------*/

int EncodeIndex(
    int Bandwidth,
    int PeriodicityIndex,
    Encoder_State *st
)
{

    if (PeriodicityIndex & kLtpHmFlag)
    {
        int LtpPitchIndex = PeriodicityIndex >> 9;
        assert(0 <= LtpPitchIndex && LtpPitchIndex <= 16);
        --PeriodicityIndex;
        assert((PeriodicityIndex & 0xff) < (1 << NumRatioBits[Bandwidth][LtpPitchIndex]));

        push_next_indice(st, PeriodicityIndex & 0xff, NumRatioBits[Bandwidth][LtpPitchIndex]);
        return NumRatioBits[Bandwidth][LtpPitchIndex];
    }
    else
    {
        push_next_indice(st, PeriodicityIndex, 8);
        return 8;
    }
}


/*-------------------------------------------------------------------*
 * GetWeight()
 *
 *
 *-------------------------------------------------------------------*/

static float GetWeight(int i)
{
    i = 3 * i - 2;

    return (float)(pow(i, 0.3) / pow(256 - 1, 0.3));
}


/*-------------------------------------------------------------------*
 * SearchPeriodicityIndex_Single()
 *
 *
 *-------------------------------------------------------------------*/

static float SearchPeriodicityIndex_Single(
    const float AbsMdct3[],
    int NumToConsider,
    int Lag,
    int FractionalResolution
)
{
    int   HighestMultiplier;
    float AbsMeanCurrent3;        /* Mean for BucketWidth == 3 */
    int   Limit;
    int   OldIndex, i;

    Limit = (NumToConsider - 1) << FractionalResolution;
    AbsMeanCurrent3 = 0;
    HighestMultiplier = 1;

    for (i=Lag; i<Limit; i+=Lag)
    {
        OldIndex = i >> FractionalResolution;
        AbsMeanCurrent3 += AbsMdct3[OldIndex] * GetWeight(HighestMultiplier);
        ++HighestMultiplier;
    }

    return AbsMeanCurrent3 / (HighestMultiplier - 1 + 0.00001f);
}


/*-------------------------------------------------------------------*
 * SearchPeriodicityIndex_Range()
 *
 *
 *-------------------------------------------------------------------*/

static void SearchPeriodicityIndex_Range(
    const float AbsMdct3[],
    int NumToConsider,
    int Lo,
    int Hi,
    int FractionalResolution,
    int Adj,
    int Spacing,
    int *PeriodicityIndex,
    float *Score
)
{
    int   Index, BestIndex;
    float CurrentScore, BestScore;
    int   B;

    BestScore = -1e30f;
    BestIndex = 0;

    for (Index = Lo; Index < Hi; Index += Spacing)
    {
        CurrentScore = SearchPeriodicityIndex_Single( AbsMdct3, NumToConsider, Index + Adj, FractionalResolution );

        if (CurrentScore > BestScore)
        {
            BestScore = CurrentScore;
            BestIndex = Index;
        }
    }

    if (BestScore > *Score)
    {
        *Score            = BestScore;
        *PeriodicityIndex = BestIndex;
    }


    B = BestIndex - (Spacing >> 1);
    B = max(Lo, B);

    for (Index = B; Index < BestIndex; ++Index)
    {
        CurrentScore = SearchPeriodicityIndex_Single( AbsMdct3, NumToConsider, Index + Adj, FractionalResolution );

        if (CurrentScore > *Score)
        {
            *Score            = CurrentScore;
            *PeriodicityIndex = Index;
        }
    }

    B = BestIndex + (Spacing >> 1);

    for (Index = BestIndex + 1; Index <= B; ++Index)
    {
        CurrentScore = SearchPeriodicityIndex_Single( AbsMdct3, NumToConsider, Index + Adj, FractionalResolution );

        if (CurrentScore > *Score)
        {
            *Score            = CurrentScore;
            *PeriodicityIndex = Index;
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * UnmapIndex()
 *
 *
 *-------------------------------------------------------------------*/

/* Returns: PeriodicityIndex */
int SearchPeriodicityIndex(
    const float Mdct[],                /* (I) Coefficients, Mdct[0..NumCoeffs-1]                      */
    const float UnfilteredMdct[],      /* (I) Unfiltered coefficients, UnfilteredMdct[0..NumCoeffs-1] */
    int NumCoeffs,                     /* (I) Number of coefficients                                  */
    int TargetBits,                    /* (I) Target bit budget (excl. Done flag)                     */
    short LtpPitchLag,
    float LtpGain,                     /* (I) LTP gain                                                */
    float *RelativeScore               /* (O) Energy concentration factor                             */
)
{
    float    AbsMdct3[MAX_LENGTH], A, B, C=0.f;
    int      i;
    int      MaxAt;
    float    Score;
    int      PeriodicityIndex = 0;
    int      NumToConsider = NumCoeffs;
    float    AbsTotal;


    Score = -1e30f;

    A = (float)fabs(Mdct[0]);
    B = (float)fabs(Mdct[1]);

    for (i = 1; i < NumToConsider - 3; i += 3)
    {
        C = (float)fabs(Mdct[i + 1]);
        AbsMdct3[i] = A + B + C;

        A = (float)fabs(Mdct[i + 2]);
        AbsMdct3[i + 1] = A + B + C;

        B = (float)fabs(Mdct[i + 3]);
        AbsMdct3[i + 2] = A + B + C;
    }

    if (i < NumToConsider - 1)
    {
        C = (float)fabs(Mdct[i + 1]);
        AbsMdct3[i] = A + B + C;
    }


    if (i + 1 < NumToConsider - 1)
    {
        A = (float)fabs(Mdct[i + 2]);
        AbsMdct3[i + 1] = A + B + C;
    }

    AbsTotal = 0.0f;

    if (UnfilteredMdct != NULL)
    {
        for (i = 0; i < NumToConsider; ++i)
        {
            AbsTotal += (float)fabs(UnfilteredMdct[i]);
        }
    }
    else
    {
        for (i = 1; i < NumToConsider - 1; i += 3)
        {
            AbsTotal += AbsMdct3[i];
        }
    }


    if ((LtpPitchLag > 0) && (LtpGain > kLtpHmGainThr))
    {
        int FractionalResolution = kLtpHmFractionalResolution;
        int Multiplier;
        int LtpPitchIndex;
        int Bandwidth;

        Bandwidth = NumCoeffs >= 256;
        LtpPitchIndex = ((LtpPitchLag + (1 << (kLtpHmFractionalResolution - 1))) >> kLtpHmFractionalResolution) - 2;
        assert(0 <= LtpPitchIndex && LtpPitchIndex <= 16);

        for (Multiplier = 1; Multiplier <= (1 << NumRatioBits[Bandwidth][LtpPitchIndex]); ++Multiplier)
        {
            float CurrentScore;
            int Lag;

            Lag = (LtpPitchLag * (int)(4 * Ratios[Bandwidth][LtpPitchIndex][Multiplier-1])) >> 2;

            if (Lag >= (4 << FractionalResolution) && (Lag <= ((NumToConsider-2) << FractionalResolution)))
            {
                CurrentScore = SearchPeriodicityIndex_Single( AbsMdct3, NumToConsider, Lag, FractionalResolution );

                if (CurrentScore > Score)
                {
                    Score = CurrentScore;
                    PeriodicityIndex = Multiplier | kLtpHmFlag;
                }
            }
        }
        PeriodicityIndex |= LtpPitchIndex << 9;
    }
    else
    {
        if (UnfilteredMdct != NULL)
        {
            MaxAt = 1;
            A = AbsMdct3[1];

            for (i = 4; i < NumToConsider - 1; i += 3)
            {

                if (AbsMdct3[i] > AbsMdct3[MaxAt])
                {
                    MaxAt = i;
                }
                A += AbsMdct3[i];
            }

            if (AbsMdct3[MaxAt] > A * 0.7f)
            {
                NumToConsider = min(NumToConsider, MaxAt + 4);
            }
        }

        SearchPeriodicityIndex_Range( AbsMdct3, NumToConsider, 0, 16, 3, GET_ADJ2(0, 6, 3), 4, &PeriodicityIndex, &Score );

        SearchPeriodicityIndex_Range( AbsMdct3, NumToConsider, 16, 80, 4, GET_ADJ2(16, 8, 4), 4, &PeriodicityIndex, &Score );

        SearchPeriodicityIndex_Range( AbsMdct3, NumToConsider, 80, 208, 3, GET_ADJ2(80, 12, 3), 4, &PeriodicityIndex, &Score );

        if (NumToConsider <= 128)
        {
            /* no long lags for band-limited MDCTs */

            SearchPeriodicityIndex_Range( AbsMdct3, NumToConsider, 208, 88 + NumToConsider, 0, GET_ADJ2(224, 188, 0), 1, &PeriodicityIndex, &Score );
        }
        else
        {

            if (TargetBits > kSmallerLagsTargetBitsThreshold && NumCoeffs >= 256)
            {
                SearchPeriodicityIndex_Range( AbsMdct3, NumToConsider, 208, 224, 1, GET_ADJ2(208, 28, 1), 1, &PeriodicityIndex, &Score );

                SearchPeriodicityIndex_Range( AbsMdct3, NumToConsider, 224, 256, 0, GET_ADJ2(224, 188, 0 ), 1, &PeriodicityIndex, &Score );
            }
            else
            {
                SearchPeriodicityIndex_Range( AbsMdct3, NumToConsider, 208, 256, 1, GET_ADJ2(208, 28, 1), 1, &PeriodicityIndex, &Score );
            }
        }
    }

    if (AbsTotal > 0)
    {
        *RelativeScore = Score/AbsTotal*(float)NumCoeffs;
    }
    else
    {
        *RelativeScore = 0;
    }


    return PeriodicityIndex;
}


/*-------------------------------------------------------------------*
 * PeakFilter()
 *
 *
 *-------------------------------------------------------------------*/

#define kPeakElevationThreshold 1.0f

static void PeakFilter(
    const float x[],                     /* (I) absolute spectrum                              */
    float y[],                           /* (O) filtered absolute spectrum, must not alias x[] */
    int L_frame                          /* (I) number of spectral lines                       */
)
{
    int flen, i;
    float a, m;

    flen = (L_frame >> 4);
    m = kPeakElevationThreshold / (float)(2*flen + 1);

    a = 0.0f;
    for (i=0; i<flen; ++i)
    {
        a += x[i];
    }

    for (i=0; i<flen; ++i)
    {
        y[i] = max(0.0f, x[i] - a*m);
        a += x[i+flen];
    }
    for (; i<L_frame-flen; ++i)
    {
        y[i] = max(0.0f, x[i] - a*m);
        a -= x[i-flen];
        a += x[i+flen];
    }

    for (; i<L_frame; ++i)
    {
        y[i] = max(0.0f, x[i] - a*m);
        a -= x[i-flen];
    }

    return;
}


/*-------------------------------------------------------------------*
 * tcx_hm_get_re()
 *
 *
 *-------------------------------------------------------------------*/

static float tcx_hm_get_re( /* Returns: RE error                 */
    const float x[],   /* i: absolute spectrum                     */
    Word16 gain,       /* i: HM gain (Q11)                         */
    int lag,
    int fract_res,
    Word16 p[],        /* i: harmonic model (Q13)                  */
    Word32 env[],      /* i: envelope (Q16)                        */
    int L_frame        /* i: number of spectral lines              */
)
{
    Word32 ne[N_MAX_ARI];
    float G, e;
    int i;

    /* Calculate new envelope with "gain" harmonic gain */
    for (i=0; i<L_frame; ++i)
    {
        ne[i] = env[i];
    }

    tcx_hm_modify_envelope( gain, lag, fract_res, p, ne, L_frame );

    /* Normalize */
    G = 0;

    for (i=0; i<L_frame; ++i)
    {
        G += x[i] * ne[i];
    }
    G = 1.0f / G;

    /* Calculate error */
    e = 0;

    for (i=0; i<L_frame; ++i)
    {
        e += (float)pow(x[i] * (ne[i] * G), 4);
    }

    return e;
}


/*-------------------------------------------------------------------*
 * tcx_hm_quantize_gain()
 *
 *
 *-------------------------------------------------------------------*/

static void tcx_hm_quantize_gain(
    const float x[],      /* i: absolute spectrum                     */
    Word32 env[],         /* i: envelope (Q16)                        */
    int lag,
    int fract_res,
    Word16 p[],           /* i: harmonic model (Q13)                  */
    int L_frame,          /* i: number of spectral lines              */
    int coder_type,       /* i: GC/VC mode                            */
    float relative_score, /* i: periodicity score                   */
    int *gain_idx,        /* o: quantization index                    */
    Word16 *gain          /* o: quantized harmonic model gain (Q11)   */
)
{
    int g;
    float be, e, pe;
    const float kLowPeriodicityThr[2] = { 0.5f, 0.2f };
    int s;

    assert(coder_type==VOICED || coder_type == GENERIC);

    s=0;
    if(coder_type==VOICED)
    {
        s=1;
    }

    *gain = 0;

    /* Disable the harmonic model if periodicity is very low */
    if (relative_score < kLowPeriodicityThr[s])
    {
        return;
    }

    be = tcx_hm_get_re(x, *gain, lag, fract_res, p, env, L_frame);

    if (coder_type == GENERIC)
    {
        e = tcx_hm_get_re(x, qGains[s][0], lag, fract_res, p, env, L_frame);
        pe = 1.05f;

        if (e * pe < be)
        {
            *gain_idx = 0;
            *gain     = qGains[s][0];
        }
    }
    else
    {
        /* Iterate over all possible gain values */
        for (g=0; g<(1 << kTcxHmNumGainBits); ++g)
        {

            e = tcx_hm_get_re(x, qGains[s][g], lag, fract_res, p, env, L_frame);

            /* Add bit penalty */
            pe = 1.0f;
            if (*gain == 0.0f)
            {
                pe = 1.05f;
            }

            /*  Minimum selection */
            if (e * pe < be)
            {
                be        = e;
                *gain_idx = g;
                *gain     = qGains[s][g];
            }
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * tcx_hm_analyse()
 *
 *
 *-------------------------------------------------------------------*/

void tcx_hm_analyse(
    const float abs_spectrum[], /* i: absolute spectrum            */
    int L_frame,       /* i: number of spectral lines              */
    Word32 env[],      /* i/o: envelope shape (Q16)                */
    int targetBits,    /* i: target bit budget                     */
    int coder_type,    /* i: GC/VC mode                            */
    int prm_hm[],      /* o: HM parameters                         */
    short LtpPitchLag, /* i: LTP pitch lag or -1 if none           */
    float LtpGain,     /* i: LTP gain                              */
    int *hm_bits       /* o: bit consumption                       */
)
{
    int lag, fract_res;
    float fspec[N_MAX_ARI], RelativeScore;
    Word16 p[2*kTcxHmParabolaHalfWidth+1], gain;

    /* Disable HM for non-GC,VC modes */
    if ((coder_type != VOICED) && (coder_type != GENERIC))
    {
        *hm_bits  = 0;
        prm_hm[0] = 0;

        return;
    }

    /* Bit consumption for the HM off case: 1 bit flag */
    *hm_bits = 1;

    /* Filter out noise and keep the peaks */
    PeakFilter(abs_spectrum, fspec, L_frame);

    /* Get the best lag index */
    prm_hm[1] = SearchPeriodicityIndex( fspec, abs_spectrum, L_frame, targetBits - *hm_bits, LtpPitchLag, LtpGain, &RelativeScore );

    /* Convert the index to lag */
    UnmapIndex( prm_hm[1], L_frame >= 256, LtpPitchLag, (targetBits - *hm_bits <= kSmallerLagsTargetBitsThreshold) || (L_frame < 256), &fract_res, &lag );

    /* Render harmonic model */
    tcx_hm_render( lag, fract_res, LtpGain, p );

    /* Calculate and quantize gain */
    gain = 0;

    tcx_hm_quantize_gain( abs_spectrum, env, lag, fract_res, p, L_frame, coder_type, RelativeScore, &prm_hm[2], &gain );

    /* Decision */
    if (gain > 0)
    {
        prm_hm[0] = 1; /* flag: on */

        *hm_bits += CountIndexBits( L_frame >= 256, prm_hm[1] );

        if (coder_type == VOICED)
        {
            *hm_bits += kTcxHmNumGainBits;
        }

        tcx_hm_modify_envelope( gain, lag, fract_res, p, env, L_frame );
    }
    else
    {
        prm_hm[0] = 0;  /* flag: off   */
        prm_hm[1] = -1; /* pitch index */
        prm_hm[2] = 0;  /* gain index  */
    }

    return;
}
