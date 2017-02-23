/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/


#include "cnst.h"
#include "prot.h"
#include "cnst.h"
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <assert.h>


#define NB_PULSES_MAX                         15

/*-------------------------------------------------------------------------
*
* Perform resynchronisation of the last glottal pulse in voiced lost frame
*
*------------------------------------------------------------------------*/


/** Get the location of the minimum energy in the given signal.
  * @param x Input signal.
  * @param length The length of the input signal.
  * @param filterLength the length of the filter length used for the energy calculation.
  * @returns Index of the position of the minimum energy, that is the position i where filter(x[i-filterLength/2],...,x[i+(filterLength-filterLength/2)-1]) is at maximum.
  */
static int GetMinimumPosition(float const * x, int length, int filterLength)
{
    int iMinEnergyPos, center, i;
    float energy, minEnergy;



    filterLength = min(filterLength, length);
    center = filterLength/2;
    iMinEnergyPos = center;
    if (filterLength > 0)
    {
        minEnergy = sum2_f(x, filterLength);
        energy = 0;
        center += 1;                                                                           /* To avoid adding 1 in the loop */
        for (i = 0; i < length-filterLength; i++)
        {
            energy  -= x[i]*x[i];
            energy  += x[i+filterLength]*x[i+filterLength];
            if (energy < 0)
            {
                minEnergy += energy;
                energy = 0;
                iMinEnergyPos = i+center;
            }
        }
    }

    return iMinEnergyPos;
}

/** Get the location of the maximum peak in the given signal.
  * @param x Input signal.
  * @param length The length of the input signal.
  * @returns Index of the position of the maximum peak, that is the position i where abs(x[i]) has it's maximum.
  */
static int FindMaxPeak(float const * x, int length)
{
    int iMax, i;
    float maxVal;



    iMax = 0;
    maxVal = 0.0f;
    for (i = 0; i < length; i++)
    {
        float const absVal = (float)fabs(x[i]);
        if (absVal > maxVal)
        {
            maxVal = absVal;
            iMax = i;
        }
    }

    return iMax;
}

static void AddSamples(float const * old_exc, float * new_exc, int L_frame, int n_samples_to_add, int const min_pos[], int const points_by_pos[], int nb_min)
{
    float * pt_dest;
    float const * pt_src;
    int last_min_pos, i, j;



    pt_dest = new_exc;
    pt_src = old_exc;
    last_min_pos = 0;
    for (i = 0; i < nb_min; i++)
    {
        float ftmp;
        /* Copy section */
        for (j = min_pos[i] - last_min_pos; j > 0; j--)
        {
            *pt_dest++ = *pt_src++;
        }
        /* Add some samples */
        ftmp = -(*pt_src/20);
        for (j = 0; j < points_by_pos[i];  j++)
        {
            *pt_dest++ = ftmp;
            ftmp = -ftmp;
        }
        /* Prepare for the next loop iteration */
        last_min_pos = min_pos[i];
    }
    /* Copy remaining length */
    for (j = L_frame-n_samples_to_add-last_min_pos; j > 0; j--)
    {
        *pt_dest++ = *pt_src++;
    }

}

static void RemoveSamples(float const * old_exc, float * new_exc, int L_frame, int n_samples_to_add, int const min_pos[], int const points_by_pos[], int nb_min)
{
    float * pt_dest;
    float const * pt_src;
    int last_min_pos, i, j;


    pt_dest = new_exc+L_frame;
    last_min_pos = L_frame-n_samples_to_add;
    for(i = nb_min-1; i >= 0; i--)
    {
        /* Compute len to copy */
        /* Copy section, removing some samples */
        pt_src = old_exc+last_min_pos;
        for (j = last_min_pos - (min_pos[i]+points_by_pos[i]); j > 0; j--)
        {
            *--pt_dest = *--pt_src;
        }
        /* Prepare for the next loop iteration */
        last_min_pos = min_pos[i];
    }
    /* Copy remaining length */
    pt_src = old_exc+last_min_pos;
    for (j = last_min_pos; j > 0; j--)
    {
        *--pt_dest = *--pt_src;
    }

}

/** Resynchronize glottal pulse positions of the signal in src_exc and store it in dst_exc.
  * src_exc holds on call the harmonic part of the signal with the constant pitch, constructed by repeating the last pitch cycle of length pitchStart.
  * dst_exc holds on return the harmonic part of the signal with the pitch changing from pitchStart to pitchEnd.
  * src_exc and dst_exc can overlap, but src_exc < dst_exc must be fullfiled.
  * @param src_exc Input excitation buffer.
  * @param dst_exc Output excitation buffer.
  * @param nFrameLength Length of the frame, that is the length of the valid data in the excitation buffer on return.
  * @param nSubframes Number of subframes in the excitation buffer. nFrameLength must be divisible by nSubframes.
  * @param pitchStart Pitch at the end of the last frame.
  * @param pitchEnd Pitch at the end of the current frame.
  */
void PulseResynchronization(float const * src_exc, float * dst_exc, int nFrameLength, int nSubframes, float pitchStart, float pitchEnd)
{
    int T0, i, k;
    float pitchDelta, samplesDelta, perCycleDeltaDelta, cycleDelta, freqStart, fractionalLeft, absPitchDiff;
    int roundedPitchStart, nSamplesDelta, nSamplesDeltaRemain, iMinPos1, iMinPos[NB_PULSES_MAX+1], iDeltaSamples[NB_PULSES_MAX+1], maxDeltaSamples, roundedCycleDelta;



    assert((nFrameLength > 0) && (nFrameLength > pitchStart) && (nSubframes > 1) && (nSubframes <= 5) && (nFrameLength%nSubframes == 0) && (pitchStart > 0) && (pitchEnd > 0) && (pitchEnd/pitchStart > 1-2.0f/(nSubframes+1)) && (src_exc != NULL) && (dst_exc != NULL) && (src_exc < dst_exc));

    pitchDelta = (pitchEnd - pitchStart)/nSubframes;
    roundedPitchStart = (int)(pitchStart+0.5f);
    freqStart = 1.0f/roundedPitchStart;

    /* Calculate number of samples to be removed (if negative) or added (if positive) */
    samplesDelta = 0.5f*pitchDelta*nFrameLength*(nSubframes+1)*freqStart;
    samplesDelta -= nFrameLength*(1.0f-pitchStart*freqStart);
    /* To have enough samples in the buffer of length nFrameLength*(nSubframes+1)/nSubframes, pitchEnd/pitchEnd must be bigger than (nSubframes-1)/(nSubframes+1)=1-2/(nSubframes+1) */
    /* Thus nSubframes must be bigger than 1 */
    nSamplesDelta = (int)floor(samplesDelta+0.5f);
    nSamplesDeltaRemain = abs(nSamplesDelta);
    /* Find the location of the glottal pulse */
    T0 = FindMaxPeak(src_exc, roundedPitchStart);
    /* Get the index of the last pulse in the resynchronized frame */
    k = (int)ceil((nFrameLength-nSamplesDelta-T0)*freqStart - 1);
    if ((k >= 0) && (k+1 <= NB_PULSES_MAX))
    {
        absPitchDiff = (float)fabs(roundedPitchStart-pitchEnd);
        /* Calculate the delta of the samples to be added/removed between consecutive cycles */
        perCycleDeltaDelta = (absPitchDiff*(nFrameLength-samplesDelta) - (float)fabs(samplesDelta)*roundedPitchStart)
                             / ((k+1)*(T0+0.5f*k*roundedPitchStart));
        /* Calculate the integer number of samples to be added/removed in each pitch cycle */
        cycleDelta = max(0, (absPitchDiff-(k+1)*perCycleDeltaDelta)*T0*freqStart);
        roundedCycleDelta = (int)(cycleDelta);
        iDeltaSamples[0] = roundedCycleDelta;
        fractionalLeft = cycleDelta-roundedCycleDelta;
        nSamplesDeltaRemain -= roundedCycleDelta;
        for (i = 1; i <= k; i++)
        {
            cycleDelta = (absPitchDiff-(k+1-i)*perCycleDeltaDelta) + fractionalLeft;
            cycleDelta = max(0, cycleDelta);
            /* Make sure that the number of samples increases */
            if (roundedCycleDelta > cycleDelta)
            {
                iDeltaSamples[i] = roundedCycleDelta;
                roundedCycleDelta = (int)(cycleDelta);
                iDeltaSamples[i-1] = roundedCycleDelta;

            }
            else
            {
                roundedCycleDelta = (int)(cycleDelta);
                iDeltaSamples[i] = roundedCycleDelta;
            }
            fractionalLeft = cycleDelta-roundedCycleDelta;
            nSamplesDeltaRemain -= roundedCycleDelta;
        }
        iDeltaSamples[k+1] = max(0, nSamplesDeltaRemain);
        maxDeltaSamples = max(iDeltaSamples[k], iDeltaSamples[k+1]);

        /* Find the location of the minimum energy between the first two pulses */
        iMinPos1 = T0+GetMinimumPosition(src_exc+T0, min(roundedPitchStart, (nSubframes+1)*nFrameLength/nSubframes-T0), maxDeltaSamples);
        if (nSamplesDelta < 0)
        {
            /* Find the location of the minimum energy before the first pulse */
            if (iMinPos1 > roundedPitchStart + iDeltaSamples[0]/2)
            {
                iMinPos[0] = iMinPos1 - roundedPitchStart - iDeltaSamples[0]/2;
            }
            else
            {
                iMinPos[0] = GetMinimumPosition(src_exc, T0, iDeltaSamples[0]) - iDeltaSamples[0]/2;
            }
            /* Find the location of the minimum energy between the pulses */
            for (i = 1; i <= k; i++)
            {
                iMinPos[i] = iMinPos1 + (i-1)*roundedPitchStart - iDeltaSamples[i]/2;
            }
            /* Find the location of the minimum energy after the last pulse */
            if (iMinPos1 + k*roundedPitchStart + iDeltaSamples[k+1] - iDeltaSamples[k+1]/2 < nFrameLength-nSamplesDelta)
            {
                iMinPos[k+1] = iMinPos1 + k*roundedPitchStart - iDeltaSamples[k+1]/2;
            }
            else
            {
                iMinPos[k+1] = T0+k*roundedPitchStart
                               + GetMinimumPosition(src_exc+T0+k*roundedPitchStart, nFrameLength-nSamplesDelta-(T0+k*roundedPitchStart), iDeltaSamples[k+1])
                               - iDeltaSamples[k+1]/2;
            }
            if (iMinPos[k+1]+iDeltaSamples[k+1] > nFrameLength-nSamplesDelta)
            {
                iDeltaSamples[k] += iMinPos[k+1]+iDeltaSamples[k+1] - (nFrameLength-nSamplesDelta);
                iDeltaSamples[k+1] = nFrameLength-nSamplesDelta-iMinPos[k+1];
            }
            /* Remove samples at the given positions */
            RemoveSamples(src_exc, dst_exc, nFrameLength, nSamplesDelta, iMinPos, iDeltaSamples, k+2);
        }
        else
        {
            /* Find the location of the minimum energy before the first pulse */
            if (iMinPos1 > roundedPitchStart)
            {
                iMinPos[0] = iMinPos1 - roundedPitchStart;
            }
            else
            {
                iMinPos[0] = GetMinimumPosition(src_exc, T0, iDeltaSamples[0]);
            }
            /* Find the location of the minimum energy between the pulses */
            for (i = 1; i <= k; i++)
            {
                iMinPos[i] = iMinPos1;
                iMinPos1 += roundedPitchStart;
            }
            /* Find the location of the minimum energy after the last pulse */
            if (iMinPos1 < nFrameLength-nSamplesDelta)
            {
                iMinPos[k+1] = iMinPos1;
            }
            else
            {
                iMinPos[k+1] = T0+k*roundedPitchStart
                               + GetMinimumPosition(src_exc+T0+k*roundedPitchStart, nFrameLength-nSamplesDelta-(T0+k*roundedPitchStart), iDeltaSamples[k+1]);
            }
            if (iMinPos[k+1]+iDeltaSamples[k+1] > nFrameLength-nSamplesDelta)
            {
                iDeltaSamples[k] += iMinPos[k+1]+iDeltaSamples[k+1] - (nFrameLength-nSamplesDelta);
                iDeltaSamples[k+1] = nFrameLength-nSamplesDelta-iMinPos[k+1];
            }
            /* Add samples at the given positions */
            AddSamples(src_exc, dst_exc, nFrameLength, nSamplesDelta, iMinPos, iDeltaSamples, k+2);
        }
    }

}


