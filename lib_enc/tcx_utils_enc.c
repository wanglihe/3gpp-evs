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
#include "cnst.h"



/* compute noise-measure flags for spectrum filling and quantization (0: tonal, 1: noise-like) */
static void ComputeSpectrumNoiseMeasure(
    const float *powerSpec,
    int L_frame,
    int startLine,
    int resetMemory,
    int *noiseFlags,
    int lowpassLine
)
{
    int i, lastTone;
    float s;
    if (resetMemory)
    {
        for (i = 0; i < lowpassLine; i++)
        {
            noiseFlags[i] = 0;
        }
    }
    for (i = lowpassLine; i < L_frame; i++)
    {
        noiseFlags[i] = 1;
    }
    if (powerSpec && startLine+6<L_frame)
    {
        lastTone = 0;
        /* noise-measure flags for spectrum filling and quantization (0: tonal, 1: noise-like) */
        i = startLine - 1;
        s = powerSpec[i-7] + powerSpec[i-6] + powerSpec[i-5] +
            powerSpec[i-4] + powerSpec[i-3] + powerSpec[i-2] +
            powerSpec[i-1] + powerSpec[i  ] + powerSpec[i+1] +
            powerSpec[i+2] + powerSpec[i+3] + powerSpec[i+4] +
            powerSpec[i+5] + powerSpec[i+6] + powerSpec[i+7];
        for (i++; i < lowpassLine - 7; i++)
        {
            float c = powerSpec[i-1] + powerSpec[i] + powerSpec[i+1];                      /**/
            s += powerSpec[i+7] - powerSpec[i-8];
            if (s >= (1.75f - 0.5f * noiseFlags[i]) * c)
            {
                noiseFlags[i] = 1;
            }
            else
            {
                noiseFlags[i] = 0;
                lastTone = i;
            }
        }
        /* lower L_frame*startRatio lines are tonal (0), upper 7 lines are processed separately */
        for (; i < lowpassLine - 1; i++)
        {
            float c = powerSpec[i-1] + powerSpec[i] + powerSpec[i+1];                      /**/
            /* running sum can't be updated any more, just use the latest one */
            if (s >= (1.75f - 0.5f * noiseFlags[i]) * c)
            {
                noiseFlags[i] = 1;
            }
            else
            {
                noiseFlags[i] = 0;
                /* lastTone = i; */
            }
        }
        noiseFlags[i] = 1;   /* uppermost line is defined as noise-like (1) */

        if (lastTone > 0)    /* spread uppermost tonal line one line upward */
        {
            noiseFlags[lastTone+1] = 0;
        }
    }
}

static void detectLowpassFac(const float *powerSpec, int L_frame, int rectWin, float *pLpFac, int lowpassLine)
{
    int i;
    float threshold;

    threshold = 0.1f * 2*NORM_MDCT_FACTOR;
    if (rectWin)
    {
        /* compensate for bad side-lobe attenuation with asymmetric windows */
        threshold *= 2.f;
    }
    for (i = lowpassLine-1; i >= lowpassLine/2; i--)
    {
        if (powerSpec[i] > threshold)
        {
            break;
        }
    }
    *pLpFac =
        (0.3f * (*pLpFac)) +
        (0.7f * ((float) (i+1) / (float) L_frame));
}

/*-----------------------------------------------------------*
 * Compute noise-measure flags for spectrum filling          *
 * and quantization (0: tonal, 1: noise-like).               *
 * Detect low pass if present.                               *
 *-----------------------------------------------------------*/
void AnalyzePowerSpectrum(
    Encoder_State *st,              /* i/o: encoder states                                  */
    int L_frame,                    /* input: frame length                                  */
    int L_frameTCX,                 /* input: full band frame length                        */
    int left_overlap,               /* input: left overlap length                           */
    int right_overlap,              /* input: right overlap length                          */
    float const mdctSpectrum[],     /* input: MDCT spectrum                                 */
    float const signal[],           /* input: windowed signal corresponding to mdctSpectrum */
    float powerSpec[]               /* output: Power spectrum. Can point to signal          */
)
{
    int i, iStart, iEnd, lowpassLine;

    lowpassLine = L_frameTCX;
    {
        TCX_MDST(signal, powerSpec, left_overlap, L_frameTCX-(left_overlap+right_overlap)/2, right_overlap);
    }
    iStart = 0;
    iEnd = L_frameTCX;

    if(st->narrowBand)
    {
        attenuateNbSpectrum(L_frameTCX, powerSpec);
    }

    /* power spectrum: MDCT^2 + MDST^2 */
    for (i = iStart; i < iEnd; i++)
    {
        powerSpec[i] *= powerSpec[i];
        powerSpec[i] += mdctSpectrum[i] * mdctSpectrum[i];
    }
    ComputeSpectrumNoiseMeasure(powerSpec, L_frameTCX, st->nmStartLine*L_frame/st->L_frame, (st->L_frame*st->last_sr_core != st->L_frame_past*st->sr_core) || (st->last_core != TCX_20_CORE), st->memQuantZeros, lowpassLine);
    if( st->total_brate <= ACELP_24k40 )
    {
        lowpassLine = (int)(2.0f*st->tcx_cfg.bandwidth*L_frame);
        detectLowpassFac(powerSpec,  L_frame, (st->last_core == ACELP_CORE), &st->measuredBwRatio, lowpassLine);
    }
    else
    {
        st->measuredBwRatio = 1.0f;
    }
}

void mdct_preShaping(float x[], int lg, const float gains[])
{
    int i, j, k, l;
    float g;
    int m, n, k1, k2;

    j = 0;                                                                            /*  not counted, is included in ptr init */
    /* FDNS_NPTS = 64 !!! */
    k = lg/FDNS_NPTS;
    m = lg%FDNS_NPTS;
    if (m)
    {
        if ( m <= (FDNS_NPTS/2) )
        {
            n = FDNS_NPTS/m;
            k1 = k;
            k2 = k + 1;
        }
        else
        {
            n = FDNS_NPTS/(FDNS_NPTS-m);
            k1 = k + 1;
            k2 = k;
        }
        for (i=0; i<lg; )
        {
            if (j%n)
            {
                k = k1;
            }
            else
            {
                k = k2;
            }
            g = 1.f/gains[j++];

            /* Limit number of loops, if end is reached */
            k = min(k, lg-i);
            for (l=0; l < k; l++)
            {
                x[i++] *= g;
            }
        }
    }
    else
    {
        for (i=0; i<lg; )
        {
            g = 1.f/gains[j++];
            for (l=0; l < k; l++)
            {
                x[i++] *= g;
            }
        }
    }

    return;
}

void AdaptLowFreqEmph(float x[],
                      int xq[],
                      float invGain,
                      short tcx_lpc_shaped_ari,
                      const float lpcGains[],
                      const int lg
                     )
{
    int i, i_max, i_max_old;

    if(!tcx_lpc_shaped_ari)
    {
        /* 1. find first magnitude maximum in lower quarter of spectrum */
        invGain *= 2.0f;
        i_max = -1;
        for (i = 0; i < lg/4; i++)
        {
            if (((xq[i] <= -2) || (xq[i] >= 2)) &&
                    ((invGain * x[i] <= -3.625f) || (invGain * x[i] >= 3.625f)))
            {
                xq[i] += (xq[i] < 0) ? -2 : 2;
                i_max = i;
                break;
            }
        }
        /* 2. compress value range of all xq up to i_max: add two steps */
        for (i = 0; i < i_max; i++)
        {
            if (x[i] < 0.0f)
            {
                xq[i] = (int)(invGain * x[i] - 0.375f);
            }
            else
            {
                xq[i] = (int)(invGain * x[i] + 0.375f);
            }
        }
        /* 3. find first mag. maximum below i_max which is half as high */
        i_max_old = i_max;
        if (i_max_old > -1)
        {
            invGain *= 2.0f;
            i_max = -1;  /* reset first maximum, update inverse gain */
            for (i = 0; i < lg/4; i++)
            {
                if (((xq[i] <= -2) || (xq[i] >= 2)) &&
                        ((invGain * x[i] <= -3.625f) || (invGain * x[i] >= 3.625f)))
                {
                    xq[i] += (xq[i] < 0) ? -2 : 2;
                    i_max = i;
                    break;
                }
            }
        }
        /* 4. re-compress and quantize all xq up to half-height i_max+1 */
        for (i = 0; i < i_max; i++)
        {
            if (x[i] < 0.0f)
            {
                xq[i] = (int)(invGain * x[i] - 0.375f);
            }
            else
            {
                xq[i] = (int)(invGain * x[i] + 0.375f);
            }
        }
        /* 5. always compress 2 lines; lines could be at index 0 and 1! */
        if (i_max_old > -1)
        {
            invGain *= 0.5f;  /* reset inverse gain */
            if (i_max < i_max_old)
            {
                i_max = i_max_old;
            }
        }
        i = i_max + 1;
        if (x[i] < 0.0f)
        {
            xq[i] = (invGain * x[i] <= -3.625f) ? xq[i] - 2 : (int)(invGain * x[i] - 0.375f);
        }
        else
        {
            xq[i] = (invGain * x[i] >=  3.625f) ? xq[i] + 2 : (int)(invGain * x[i] + 0.375f);
        }
        i++;
        if (x[i] < 0.0f)
        {
            xq[i] = (invGain * x[i] <= -3.625f) ? xq[i] - 2 : (int)(invGain * x[i] - 0.375f);
        }
        else
        {
            xq[i] = (invGain * x[i] >=  3.625f) ? xq[i] + 2 : (int)(invGain * x[i] + 0.375f);
        }
    }
    else  /*if(!tcx_lpc_shaped_ari)*/
    {
        PsychAdaptLowFreqEmph(x, lpcGains);
    }/*if(!tcx_lpc_shaped_ari)*/
}

void PsychAdaptLowFreqEmph(float x[],
                           const float lpcGains[]
                          )
{
    int i;
    float max, fac, tmp;
    max = tmp = lpcGains[0];

    /* find minimum (tmp) and maximum (max) of LPC gains in low frequencies */
    for (i = 1; i < 9; i++)
    {
        if (tmp > lpcGains[i])
        {
            tmp = lpcGains[i];
        }
        if (max < lpcGains[i])
        {
            max = lpcGains[i];
        }
    }
    tmp *= 32.0f;
    if ((max < tmp) && (max > FLT_MIN))
    {
        fac = tmp = (float)pow(tmp / max, 0.0078125f);

        /* gradual boosting of lowest 32 bins; DC is boosted by (tmp/max)^1/4 */
        for (i = 31; i >= 0; i--)
        {
            x[i] *= fac;
            fac  *= tmp;
        }
    }

}

float SQ_gain(   /* output: SQ gain                   */
    float x[],     /* input:  vector to quantize        */
    int nbitsSQ,   /* input:  number of bits targeted   */
    int lg)        /* input:  vector size (2048 max)    */
{
    int    i, iter;
    float  ener, tmp, target, fac, offset;
    float  en[N_MAX/4];


    /* energy of quadruples with 9dB offset */
    for (i=0; i<lg; i+=4)
    {

        ener = 0.01f + x[i]*x[i] + x[i+1]*x[i+1] + x[i+2]*x[i+2] + x[i+3]*x[i+3];
        en[i>>2] = (float)log10(ener);  /* saves a MAC */
    }

    /* SQ scale: 4 bits / 6 dB per quadruple */
    target = 0.15f * (float)(nbitsSQ - (lg>>4));
    fac = 12.8f;
    offset = fac;

    /* find offset (0 to 128 dB with step of 0.125dB) */
    for (iter=0; iter<10; iter++)
    {
        fac *= 0.5f;
        offset -= fac;
        ener = 0.0f;

        for (i=0; i<lg/4; i++)
        {
            tmp = en[i] - offset;

            /* avoid SV with 1 bin of amp < 0.5f */
            if (tmp > 0.3f)
            {
                ener += tmp;

                /* if ener is above target -> break and increase offset */
                if (ener > target)
                {
                    offset += fac;
                    break;
                }
            }
        }
    }

    /* return gain */


    return (float)pow(10.0f, 0.45f + 0.5f*offset);
}

void tcx_scalar_quantization(
    float *x,                 /* i: input coefficients            */
    int *xq,                  /* o: quantized coefficients        */
    int L_frame,              /* i: frame length                  */
    float gain,               /* i: quantization gain             */
    float offset,             /* i: rounding offset (deadzone)    */
    int *memQuantZeros,       /* o: coefficients set to 0         */
    const int alfe_flag
)
{
    int i;
    float gainInv, rounding, rounding2;


    /* Init scalar quantizer */
    gainInv = 1.0f/gain;
    rounding = offset;
    rounding2 = -offset;

    for (i = L_frame - 1; (memQuantZeros[i]) && ((float)fabs(x[i]) * gainInv < 1.0f); i--)
    {
        xq[i] = 0;
    }
    for (; i >= 0; i--)
    {
        if (x[i]>0.f)
        {
            xq[i] = ((int) (rounding + x[i]*gainInv));
        }
        else
        {
            xq[i] = ((int) (rounding2 + x[i]*gainInv));
        }
    }

    /* don't instrument; BASOP code restricts to 16 bit by using 16 bit operators only*/
    for(i=0; i<L_frame; i++)
    {
        xq[i] = max(min(xq[i], 32767), -32768);
    }

    if (!alfe_flag)
    {
        AdaptLowFreqEmph( x,xq, gainInv, 0, NULL, L_frame );
    }

    return;
}

int tcx_scalar_quantization_rateloop(
    float *x,                 /* i  : input coefficients            */
    int *xq,                  /* o  : quantized coefficients        */
    int L_frame,              /* i  : frame length                  */
    float *gain,              /* i/o: quantization gain             */
    float offset,             /* i  : rounding offset (deadzone)    */
    int *memQuantZeros,       /* o  : coefficients set to 0         */
    int *lastnz_out,          /* i/o: last nonzero coeff index      */
    int target,               /* i  : target number of bits         */
    int *nEncoded,            /* o  : number of encoded coeff       */
    int *stop,                /* i/o: stop param                    */
    int sqBits_in_noStop,     /* i  : number of sqBits as determined in prev. quant. stage, w/o using stop mechanism (ie might exceed target bits) */
    int sqBits_in,            /* i  : number of sqBits as determined in prev. quant. stage, using stop mechanism (ie always <= target bits) */
    int tcxRateLoopOpt,       /* i  : turns on/off rateloop optimization */
    const int tcxonly,
    CONTEXT_HM_CONFIG *hm_cfg
)
{
    const int iter_max = 4;
    int sqBits, stopFlag;
    int ubfound,lbfound;
    float ub=0.f,lb=0.f;
    float shift;
    int iter;
    float sqGain;
    float w_lb, w_ub;
    const int kDampen = 10;
    int old_stopFlag;
    int old_nEncoded;
    int old_sqBits;
    float mod_adjust0, mod_adjust1;
    float inv_target;
    const float kMargin = 0.96f;
    int lastnz;



    /* Init */
    sqGain = *gain;
    stopFlag = *stop;
    ubfound = 0;
    lbfound = 0;
    shift = 0.25f;
    w_lb = 0.0f;
    w_ub = 0.0f;
    lastnz = *lastnz_out;

    old_stopFlag   = stopFlag;
    old_nEncoded   = *nEncoded;
    old_sqBits     = sqBits_in_noStop;

    sqBits      = sqBits_in;

    mod_adjust0 = max(1.0f, 2.3f - 0.0025f * target);
    mod_adjust1 = 1.0f / mod_adjust0;

    inv_target  = 1.0f/(float)target;

    /* Loop */
    for ( iter=0 ; iter<iter_max ; iter++ )
    {
        if(tcxRateLoopOpt == 2)
        {

            /* Ajust sqGain */
            if ( stopFlag )
            {
                lbfound = 1;
                lb = sqGain;
                w_lb = (float)(stopFlag-target+kDampen);
                if (ubfound)
                {
                    sqGain = (lb*w_ub + ub*w_lb)/(w_ub+w_lb);
                }
                else
                {
                    sqGain *= (1.0f+((float)(stopFlag/kMargin)*inv_target-1.0f)*mod_adjust0);
                }
            }
            else
            {
                ubfound = 1;
                ub = sqGain;
                w_ub = (float)(target-sqBits+kDampen);
                if (lbfound)
                {
                    sqGain = (lb*w_ub + ub*w_lb)/(w_ub+w_lb);
                }
                else
                {
                    sqGain *= (1.0f-(1.0f-(float)(sqBits*kMargin)*inv_target)*mod_adjust1);
                }
            }

        }
        else     /* tcxRateLoopOpt != 2 */
        {

            /* Ajust sqGain */
            if ( stopFlag )
            {
                lbfound = 1;
                lb = sqGain;
                if (ubfound)
                {
                    sqGain = (float)sqrt( lb * ub );
                }
                else
                {
                    sqGain = sqGain * (float)pow(10.0f, shift/10.0f);
                    shift *= 2.0f;
                }
            }
            else
            {
                ubfound = 1;
                ub = sqGain;
                if (lbfound)
                {
                    sqGain = (float)sqrt( lb * ub );
                }
                else
                {
                    sqGain = sqGain * (float)pow(10.0f, -shift/10.0f);
                    shift *= 2.0f;
                }
            }
        }

        /* Quantize spectrum */
        tcx_scalar_quantization( x, xq, L_frame, sqGain, offset, memQuantZeros, tcxonly );

        /* Estimate bitrate */
        if(tcxRateLoopOpt >= 1)
        {
            stopFlag = 0;
        }
        else
        {
            stopFlag = 1;
        }

        sqBits = ACcontextMapping_encode2_estimate_no_mem_s17_LC( xq, L_frame,
                 &lastnz,
                 nEncoded, target, &stopFlag, hm_cfg);

        if( tcxRateLoopOpt >= 1 )
        {
            if((*nEncoded>=old_nEncoded && (stopFlag>=old_stopFlag)) || (*nEncoded>old_nEncoded && (stopFlag==0 && old_stopFlag>0)) || (stopFlag==0 && old_stopFlag==0))
            {
                *gain = sqGain;
                old_nEncoded=*nEncoded;
                old_stopFlag=stopFlag;
                old_sqBits = sqBits;
                *lastnz_out = lastnz;
            }
        }
    } /* for ( iter=0 ; iter<iter_max ; iter++ ) */

    if( tcxRateLoopOpt >= 1 )
    {
        /* Quantize spectrum */
        tcx_scalar_quantization( x, xq, L_frame, *gain, offset, memQuantZeros, tcxonly );

        /* Output */
        *nEncoded = old_nEncoded;
        sqBits = old_sqBits;
        *stop  = old_stopFlag;
    }
    else
    {
        /* Output */
        *gain = sqGain;
        *stop = stopFlag;
        *lastnz_out = lastnz;
    }



    return sqBits;
}

void QuantizeGain(int n, float * pGain, int * pQuantizedGain)
{
    float ener, gain;
    int quantizedGain;

    ener = (float)sqrt((float)n / (float)NORM_MDCT_FACTOR);

    gain = *pGain * ener;

    assert(gain > 0);

    /* quantize gain with step of 0.714 dB */
    quantizedGain = (int)floor(0.5f + 28.0f * (float)log10(gain));

    if (quantizedGain < 0)
    {
        quantizedGain = 0;
    }
    if (quantizedGain > 127)
    {
        quantizedGain = 127;
    }

    *pQuantizedGain = quantizedGain;
    *pGain = (float)pow(10.0f, ((float)quantizedGain)/28.0f) / ener;
}

void tcx_noise_factor(
    float *x_orig,          /* i: unquantized mdct coefficients             */
    float *sqQ,             /* i: quantized mdct coefficients               */
    int iFirstLine,         /* i: first coefficient to be considered        */
    int lowpassLine,        /* i: last nonzero coefficients after low-pass  */
    int nTransWidth,        /* i: minimum size of hole to be checked        */
    int L_frame,            /* i: frame length                              */
    float gain_tcx,         /* i: tcx gain                                  */
    float tiltCompFactor,   /* i: LPC tilt compensation factor              */
    float *fac_ns,          /* o: noise factor                              */
    int *quantized_fac_ns   /* o: quantized noise factor                    */
)
{
    int i, k, win, segmentOffset;
    float inv_gain2, sqErrorNrg, n, tilt_factor, tmp;
    float att;  /* noise level attenuation factor for transient windows */


    /*Adjust noise filling level*/
    sqErrorNrg = 0.0f;
    n = 0.0f;
    /* max() */
    tilt_factor = 1.0f /(float)pow(max(0.375f, tiltCompFactor), 1.0f/(float)L_frame);    /* 1/(a^b) = a^-b */
    inv_gain2 = 1.0f / ((float)(nTransWidth * nTransWidth) * gain_tcx);

    /* find last nonzero line below iFirstLine, use it as start offset */
    for (i = iFirstLine; i > iFirstLine/2; i--)
    {
        if (sqQ[i] != 0.0f)
        {
            break;
        }
    }
    inv_gain2 *= (float)pow(tilt_factor, (float)i);

    segmentOffset = ++i;
    if (nTransWidth <= 3)
    {
        att = tmp = FLT_MIN;
        for (k = i & 0xFFFE; k < lowpassLine; k++)
        {
            att += x_orig[k] * x_orig[k]; /* even-index bins, left sub-win */
            k++;
            tmp += x_orig[k] * x_orig[k]; /* odd-index bins, right sub-win */
        }
        att = (float)sqrt((min(att, tmp)*2.0f) / (att + tmp));
    }
    else
    {
        att = 1.0f;
    }
    win = 0;
    for (; i < lowpassLine; i++)
    {
        inv_gain2 *= tilt_factor;
        if (sqQ[i] != 0)    /* current line is not zero, so reset pointers */
        {
            if (win > 0)   /* add segment sum to sum of segment magnitudes */
            {
                k = i - segmentOffset;
                if (nTransWidth <= 3)
                {
                    n += (k > 2 * nTransWidth - 4) ? (float)(k - nTransWidth + 1)
                         : (float)(k*k) * 0.28125f/nTransWidth;        /* table lookup instead of  */
                }
                else
                {
                    n += (k > 12) ? (float)k - 7.0f : (float)(k*k) * 0.03515625f;
                }
                for (k = segmentOffset; k < i-win; k++)
                {
                    sqErrorNrg += sqQ[k] * (float)nTransWidth;
                    sqQ[k] = 0;
                }
                for (; win > 0; win--)
                {
                    sqErrorNrg += sqQ[k] * (float)win;
                    sqQ[k++] = 0;
                }
            }
            segmentOffset = i + 1; /* new segment might start at next line */
        }
        else   /* current line is zero, so update pointers & segment sum */
        {
            if (win < nTransWidth)
            {
                win++;
            }
            /* update segment sum: magnitudes scaled by smoothing function */
            sqQ[i] = (float)fabs(x_orig[i]) * (float)win * inv_gain2;
        }
    }
    if (win > 0)    /* add last segment sum to sum of segment magnitudes */
    {
        k = i - segmentOffset;
        if (nTransWidth <= 3)
        {
            n += (k > 2 * nTransWidth - 4) ? (float)(k - nTransWidth + 1)
                 : (float)(k*k) * 0.28125f/nTransWidth;            /* table lookup instead of  */
        }
        else
        {
            n += (k > 12) ? (float)k - 7.0f : (float)(k*k) * 0.03515625f;
        }
        for (k = segmentOffset; k < i-win; k++)
        {
            sqErrorNrg += sqQ[k] * (float)nTransWidth;
            sqQ[k] = 0;
        }
        for (; win > 0; win--)
        {
            sqErrorNrg += sqQ[k] * (float)win;
            sqQ[k++] = 0;
        }
    }

    /* noise level factor: average of segment magnitudes of noise bins */
    if (n > 0.0f)
    {
        *fac_ns = (sqErrorNrg * att) / n;
    }
    else
    {
        *fac_ns = 0.0f;
    }

    /* quantize, dequantize noise level factor (range 0.09375 - 0.65625) */
    *quantized_fac_ns = (int)(0.5f + *fac_ns * 1.34375f*(1<<NBITS_NOISE_FILL_LEVEL));
    if (*quantized_fac_ns > (1<<NBITS_NOISE_FILL_LEVEL) - 1)
    {
        *quantized_fac_ns = (1<<NBITS_NOISE_FILL_LEVEL) - 1;
    }
    *fac_ns = (float)(*quantized_fac_ns) * 0.75f / (1<<NBITS_NOISE_FILL_LEVEL);

}

void tcx_encoder_memory_update(
    const float *wsig,      /* i : target weighted signal        */
    float *xn_buf,          /* i/o: mdct output buffer/time domain weigthed synthesis        */
    int L_frame_glob,       /* i: global frame length                         */
    const float *Ai,              /* i: Unquantized (interpolated) LPC coefficients */
    float *A,               /* i: Quantized LPC coefficients                  */
    float preemph_f,        /* i: preemphasis factor                          */
    LPD_state *LPDmem,      /* i/o: coder memory state                        */
    Encoder_State *st,
    int m,
    float *synthout
)
{
    float tmp;
    float buf[1+M+L_FRAME_PLUS];
    float *synth;



    /* Output synth */
    mvr2r(xn_buf, synthout, L_frame_glob);


    /* Update synth */
    synth = buf + 1 + m;
    mvr2r(LPDmem->syn, buf, 1+m);

    mvr2r(xn_buf, synth, L_frame_glob);
    mvr2r(synth+L_frame_glob-m-1, LPDmem->syn, 1+m);

    if(!st->tcxonly)
    {
        /* Update weighted synthesis */
        residu(Ai+(st->nb_subfr-1)*(M+1), M,synth+L_frame_glob-1, &tmp, 1);
        LPDmem->mem_w0=wsig[L_frame_glob-1]-tmp;
    }


    /* Emphasis of synth -> synth_pe */
    tmp = synth[-m-1];
    preemph(synth-m, preemph_f, m+L_frame_glob, &tmp);
    mvr2r(synth+L_frame_glob-m, LPDmem->mem_syn, m);
    mvr2r(synth+L_frame_glob-m, LPDmem->mem_syn2, m);
    mvr2r( synth+L_frame_glob-L_SYN_MEM, LPDmem->mem_syn_r, L_SYN_MEM);

    if ( !st->tcxonly || (L_frame_glob==L_FRAME16k))
    {
        /* Update excitation */
        if (L_frame_glob < L_EXC_MEM)
        {
            mvr2r( LPDmem->old_exc+(L_frame_glob), LPDmem->old_exc, L_EXC_MEM-(L_frame_glob) );
            residu(A, M,synth, LPDmem->old_exc+L_EXC_MEM-(L_frame_glob), (L_frame_glob));
        }
        else
        {
            residu(A, M,synth+(L_frame_glob)-L_EXC_MEM, LPDmem->old_exc, L_EXC_MEM);
        }
    }
}



/*---------------------------------------------------------------
 * Residual Quantization
 *--------------------------------------------------------------*/

/* Returns: number of bits used (including "bits") */
int tcx_ari_res_Q_spec(
    const float x_orig[], /* i: original spectrum                  */
    const int signs[],    /* i: signs (x_orig[.]<0)                */
    float x_Q[],          /* i/o: quantized spectrum               */
    int L_frame,          /* i: number of lines                    */
    float gain,           /* i: TCX gain                           */
    int prm[],            /* o: bit-stream                         */
    int target_bits,      /* i: number of bits available           */
    int bits,             /* i: number of bits used so far         */
    float deadzone,       /* i: quantizer deadzone                 */
    const float x_fac[]   /* i: spectrum post-quantization factors */
)
{
    int i, j, num_zeros;
    int zeros[L_FRAME_PLUS];
    float fac_m, fac_p, thres, sign, x_Q_m, x_Q_p;


    /* Limit the number of residual bits */
    target_bits = min(target_bits, NPRM_RESQ);


    /* Requantize the spectrum line-by-line */
    fac_m = deadzone * 0.5f;
    fac_p = 0.5f - fac_m;
    num_zeros = 0;
    for (i=0; i<L_frame; ++i)
    {
        if (bits >= target_bits)
        {
            /* no bits left */
            break;
        }
        if (x_Q[i] != 0)
        {
            sign = (1-2*signs[i])*x_fac[i];

            x_Q_m = x_Q[i] - sign*fac_m;
            x_Q_p = x_Q[i] + sign*fac_p;
            if (fabs(x_orig[i] - gain * x_Q_m) < fabs(x_orig[i] - gain * x_Q_p))   /* Decrease magnitude */
            {
                x_Q[i] = x_Q_m;
                prm[bits++] = 0;
            }
            else   /* Increase magnitude */
            {
                x_Q[i] = x_Q_p;
                prm[bits++] = 1;
            }
        }
        else
        {
            zeros[num_zeros++] = i;
        }
    }

    /* Requantize zeroed-lines of the spectrum */
    fac_p = (1.0f - deadzone)*0.33f;
    --target_bits; /* reserve 1 bit for the check below */
    for (j=0; j<num_zeros; ++j)
    {
        if (bits >= target_bits)
        {
            /* 1 or 0 bits left */
            break;
        }

        i = zeros[j];

        thres = fac_p * x_fac[i];
        if (fabs(x_orig[i]) > thres * gain)
        {
            prm[bits++] = 1;
            prm[bits++] = 1-signs[i];
            x_Q[i] = (2-4*signs[i]) * thres;
        }
        else
        {
            prm[bits++] = 0;
        }
    }

    return bits;
}

#define kMaxEstimatorOvershoot  5
#define kMaxEstimatorUndershoot 0


int tcx_res_Q_gain(
    float sqGain,
    float *gain_tcx,
    int *prm,
    int sqTargetBits
)
{
    int bits;
    float gain_reQ;

    /*Refine the gain quantization : Normal greedy gain coding */
    gain_reQ=*gain_tcx;
    for(bits=0; bits<TCX_RES_Q_BITS_GAIN; bits++)
    {
        if(sqGain<gain_reQ)
        {
            prm[bits]=0;
            gain_reQ*=gain_corr_inv_fac[bits];
        }
        else
        {
            prm[bits]=1;
            gain_reQ*=gain_corr_fac[bits];
        }
        if(bits<sqTargetBits)
        {
            *gain_tcx=gain_reQ;
        }
    }

    return(bits);
}


static void refine_0(float x_orig, float *x_Q, float sqGain, int *prm, int *bits, float sq_round, float lf_deemph_factor)
{
    float /*b,*/ thres;
    /* was  */
    /*b = x_orig/sqGain;*/
    thres = (1.0f-sq_round)*0.33f*lf_deemph_factor;
    if (x_orig > thres * sqGain)
    {
        /* was (b > thres) */
        prm[(*bits)++] = 1;
        prm[(*bits)++] = 1;
        *x_Q = 2.f*thres;
    }
    else if (x_orig < -thres * sqGain)
    {
        /* was (b < -thres) */
        prm[(*bits)++] = 1;
        prm[(*bits)++] = 0;
        *x_Q = -2.f*thres;
    }
    else
    {
        prm[(*bits)++] = 0;
    }

    return;
}


int tcx_res_Q_spec(
    float *x_orig,
    float *x_Q,
    int L_frame,
    float sqGain,
    int *prm,
    int sqTargetBits,
    int bits,
    float sq_round,
    const float lf_deemph_factors[]
)
{
    int i;
    float fac_m, fac_p;

    /* Limit the number of residual bits */
    sqTargetBits = min(sqTargetBits, NPRM_RESQ);

    /* Requantize the spectrum line-by-line */
    fac_p = 0.5f - sq_round * 0.5f;
    fac_m = sq_round * 0.5f;
    if (!lf_deemph_factors)
    {
        for (i = 0; i < L_frame; i++)
        {
            if (bits >= sqTargetBits-kMaxEstimatorUndershoot)
            {
                fac_m=fac_p=0;
                if (bits >= min(NPRM_RESQ,sqTargetBits+kMaxEstimatorOvershoot))
                {
                    break;
                }
            }

            if (x_Q[i] != 0.0f)
            {
                if(x_orig[i]<(sqGain)*x_Q[i])
                {
                    prm[bits++]=0;
                    x_Q[i]-=(x_Q[i]>0)?fac_m : fac_p;
                }
                else
                {
                    prm[bits++]=1;
                    x_Q[i]+=(x_Q[i]>0)?fac_p : fac_m;
                }
            }
        }
        sqTargetBits -= 2; /* Quantize zeroed lines of the spectrum */
        for (i = 0; (i < L_frame) && (bits < sqTargetBits); i++)
        {
            /* bits < sqTargetBits */
            if (x_Q[i] == 0.0f)
            {
                refine_0(x_orig[i], &x_Q[i], sqGain, prm, &bits, sq_round, 1.0f);  /*  inlined */
            }
        }
        /* Make sure that all possible bits are initialized */
        for (i = bits; i < NPRM_RESQ; i++)
        {
            prm[i] = 0;
        }
        return bits;
    }
    for (i = 0; i < L_frame; i++)
    {
        if (bits >= sqTargetBits-kMaxEstimatorUndershoot)
        {
            fac_m=fac_p=0;
            if (bits >= min(NPRM_RESQ,sqTargetBits+kMaxEstimatorOvershoot))
            {
                break;
            }
        }

        if (x_Q[i] != 0 && lf_deemph_factors[i] > 0.5f)
        {
            if(x_orig[i]<(sqGain)*x_Q[i])
            {
                prm[bits++]=0;
                x_Q[i]-=(x_Q[i]>0)?fac_m*lf_deemph_factors[i]:fac_p*lf_deemph_factors[i];
            }
            else
            {
                prm[bits++]=1;
                x_Q[i]+=(x_Q[i]>0)?fac_p*lf_deemph_factors[i]:fac_m*lf_deemph_factors[i];
            }
        }
    }

    /*Quantize zeroed-line of the spectrum*/
    for (i = 0; (i < L_frame) && (bits < (sqTargetBits-2)); i++)
    {
        /* For (bits >= (sqTargetBits-2)) */
        if (x_Q[i] == 0 && lf_deemph_factors[i] > 0.5f)
        {
            refine_0(x_orig[i], &x_Q[i], sqGain, prm, &bits, sq_round, lf_deemph_factors[i]);
        }
    }

    /*Be sure that every possible bits are initialized*/
    for (i = bits; i < NPRM_RESQ; i++)
    {
        prm[i]=0;
    }

    return bits;
}


void ProcessIGF(
    IGF_ENC_INSTANCE_HANDLE         const hInstance,          /**< in: instance handle of IGF Encoder */
    Encoder_State                        *st,                 /**< in: Encoder state */
    float                                *pMDCTSpectrum,      /**< in: MDCT spectrum */
    float                                *pPowerSpectrum,     /**< in: MDCT^2 + MDST^2 spectrum, or estimate */
    int                                   isTCX20,            /**< in: flag indicating if the input is TCX20 or TCX10/2xTCX5 */
    int                                   isTNSActive,        /**< in: flag indicating if the TNS is active */
    int                                   isTransition,       /**< in: flag indicating if the input is the transition from from ACELP to TCX20/TCX10 */
    int                                   frameno             /**< in: flag indicating index of current subframe */
)
{
    int igfGridIdx;
    int isIndepFlag;
    int bsBits;
    int pBsStart;

    isIndepFlag = 1;
    if (isTransition && isTCX20)
    {
        igfGridIdx = IGF_GRID_LB_TRAN;
    }
    else if (isTCX20)
    {
        igfGridIdx = IGF_GRID_LB_NORM;
    }
    else
    {
        /* It is short block */
        igfGridIdx = IGF_GRID_LB_SHORT;
        if (frameno == 1)
        {
            isIndepFlag = 0;
        }
    }

    IGFEncApplyMono(hInstance,                   /**< in: instance handle of IGF Encoder */
                    igfGridIdx,                  /**< in: IGF grid index */
                    st,                          /**< in: Encoder state */
                    pMDCTSpectrum,               /**< in: MDCT spectrum */
                    pPowerSpectrum,              /**< in: MDCT^2 + MDST^2 spectrum, or estimate */
                    isTCX20,                     /**< in: flag indicating if the input is TCX20 or TCX10/2xTCX5 */
                    isTNSActive,
                    (st->last_core == ACELP_CORE)
                   );                /**< in: flag indicating if the TNS is active */


    {
        const float tns_predGain = (&st->hIGFEnc)->tns_predictionGain;
        const short int startLine = (&st->hIGFEnc)->infoStartLine;
        const short int endLine = (&st->hIGFEnc)->infoStopLine;
        const int maxOrder = 8;
        const float* spec_before = (&st->hIGFEnc)->spec_be_igf;
        int curr_order = 0;
        float A[ITF_MAX_FILTER_ORDER+1];
        float predictionGain = 0;
        int* flatteningTrigger = &(&st->hIGFEnc)->flatteningTrigger;
        ITF_Detect( spec_before, startLine, endLine, maxOrder, A, &predictionGain, &curr_order );

        *flatteningTrigger = tns_predGain < 1.15 && predictionGain < 1.15;
    }

    pBsStart = st->next_ind;
    hInstance->infoTotalBitsPerFrameWritten = 0;
    IGFEncWriteBitstream( hInstance, isTCX20 ? NULL : st, &hInstance->infoTotalBitsPerFrameWritten, igfGridIdx, isIndepFlag );

    bsBits = st->next_ind - pBsStart;
    if (!isTCX20)
    {
        IGFEncConcatenateBitstream(hInstance, bsBits, &st->next_ind, &st->nb_bits_tot, st->ind_list);
    }

    return;
}


void attenuateNbSpectrum(int L_frame, float *spectrum)
{
    int i;
    int length = L_frame / 20;
    float att = (length == 8)?0.6f:0.66f;

    for(i=0; i<length; i++)
    {
        spectrum[L_frame - length + i] *= att;
        att *= att;
    }

    return;
}
