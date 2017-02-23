/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#define _USE_MATH_DEFINES

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "options.h"
#include "prot.h"
#include "cnst.h"
#include "stat_com.h"



/***********************************************************************************/
/* forward declaration for local functions, see implementation at end of this file */
/***********************************************************************************/

static void calcPseudoSpec(float const * mdctSpec, unsigned int nSamples, float floorPowerSpectrum, float * powerSpec);
static void getEnvelope(unsigned int nSamples, float const * powerSpec, float F0, float * envelope, float * smoothedSpectrum);
static void GetF0(unsigned int const nSamples, unsigned int const nSamplesCore, float const * const powerSpectrum, float const pitchLag, float * const pOrigF0, float * const pF0);
static void findStrongestHarmonics(unsigned int nSamples, float const * powerSpectrum, float F0, unsigned int nTotalHarmonics, unsigned int * pHarmonicIndexes, unsigned int * pnHarmonics);
static void CorrectF0(unsigned int const * pHarmonicIndexes, unsigned int nHarmonics, float * pF0);
static void findCandidates(unsigned int nSamples, const float * MDCTSpectrum, float * thresholdModificationNew, float floorPowerSpectrum /* i: lower limit for power spectrum bins  */ );
static void modifyThreshold(int i, float F0, float threshold, float * thresholdModification);
static void modifyThresholds(float F0, float origF0, float * thresholdModification);
static void RefineThresholdsUsingPitch(unsigned int nSamples, unsigned int nSamplesCore, float const * powerSpectrum, float lastPitchLag, float currentPitchLag, float * pF0, float * thresholdModification);
static void findTonalComponents(unsigned short int * indexOfTonalPeak, unsigned short int * lowerIndex, unsigned short int * upperIndex, unsigned int *numIndexes, unsigned int nSamples, const float * powerSpectrum, float F0, float const * thresholdModification);

/*******************************************************/
/*-------------- public functions -------------------- */
/*******************************************************/

/* Detect tonal components in the lastMDCTSpectrum, use
 * secondLastPowerSpectrum for the precise location of the peaks and
 * store them in indexOfTonalPeak.  Updates lowerIndex, upperIndex,
 * pNumIndexes accordingly. */
void DetectTonalComponents(unsigned short int indexOfTonalPeak[],
                           unsigned short int lowerIndex[],
                           unsigned short int upperIndex[],
                           unsigned int * pNumIndexes,
                           float lastPitchLag, float currentPitchLag,
                           float const lastMDCTSpectrum[],
                           float const scaleFactors[],
                           float const secondLastPowerSpectrum[],
                           unsigned int nSamples
                           ,unsigned int nSamplesCore
                           ,float floorPowerSpectrum        /* i: lower limit for power spectrum bins  */
                          )
{
    float F0;
    float thresholdModification[L_FRAME_MAX];
    float * const pScaledMdctSpectrum = thresholdModification;


    /* Convert from 16 bit to 32 bit */
    mvr2r(lastMDCTSpectrum, pScaledMdctSpectrum, nSamples);
    mdct_noiseShaping(pScaledMdctSpectrum, nSamplesCore, scaleFactors);
    v_multc(pScaledMdctSpectrum + nSamplesCore, scaleFactors[FDNS_NPTS-1], pScaledMdctSpectrum + nSamplesCore, nSamples - nSamplesCore );

    /* Find peak candidates in the last frame. */
    findCandidates(nSamples, pScaledMdctSpectrum, thresholdModification, floorPowerSpectrum );

    /* Refine peak candidates using the pitch information */
    RefineThresholdsUsingPitch(nSamples, nSamplesCore, secondLastPowerSpectrum, lastPitchLag, currentPitchLag, &F0, thresholdModification);
    /* Find peaks in the second last frame */
    findTonalComponents(indexOfTonalPeak, lowerIndex, upperIndex,  pNumIndexes, nSamples, secondLastPowerSpectrum, F0, thresholdModification);

}

/* When called, the tonal components are already stored in
 * indexOfTonalPeak.  Detect tonal components in the lastMDCTSpectrum,
 * use secondLastPowerSpectrum for the precise location of the peaks and
 * then keep in indexOfTonalPeak only the tonal components that are
 * again detected Updates indexOfTonalPeak, lowerIndex, upperIndex,
 * phaseDiff, phases, pNumIndexes accordingly. */
void RefineTonalComponents(unsigned short int indexOfTonalPeak[],
                           unsigned short int lowerIndex[],
                           unsigned short int upperIndex[],
                           float phaseDiff[],
                           float phases[], unsigned int * pNumIndexes,
                           float lastPitchLag, float currentPitchLag,
                           float const lastMDCTSpectrum[],
                           float const scaleFactors[],
                           float const secondLastPowerSpectrum[],
                           unsigned int nSamples
                           ,unsigned int nSamplesCore
                           ,float floorPowerSpectrum        /* i: lower limit for power spectrum bins  */
                          )
{
    unsigned short int newIndexOfTonalPeak[MAX_NUMBER_OF_IDX];
    unsigned short int newLowerIndex[MAX_NUMBER_OF_IDX];
    unsigned short int newUpperIndex[MAX_NUMBER_OF_IDX];
    unsigned int newNumIndexes, nPreservedPeaks;
    unsigned int iNew, iOld, j;
    float * pOldPhase, * pNewPhase;

    DetectTonalComponents(newIndexOfTonalPeak,
                          newLowerIndex,
                          newUpperIndex,
                          &newNumIndexes,
                          lastPitchLag,
                          currentPitchLag,
                          lastMDCTSpectrum,
                          scaleFactors,
                          secondLastPowerSpectrum,
                          nSamples,
                          nSamplesCore,
                          floorPowerSpectrum );

    nPreservedPeaks = 0;
    iNew = 0;
    pOldPhase = phases;
    pNewPhase = phases;
    for (iOld = 0; iOld < *pNumIndexes; iOld++)
    {
        /* We don't want that the old peak index is at the border of the new peak region, that is why >= newUpperIndex and > newLowerIndex */
        while ((iNew < newNumIndexes) && (indexOfTonalPeak[iOld] >= newUpperIndex[iNew]))
        {
            ++iNew;
        }
        if ((iNew < newNumIndexes) && (indexOfTonalPeak[iOld] > newLowerIndex[iNew]))
        {
            newIndexOfTonalPeak[nPreservedPeaks] = indexOfTonalPeak[iOld];
            newLowerIndex[nPreservedPeaks] = lowerIndex[iOld];
            newUpperIndex[nPreservedPeaks] = upperIndex[iOld];
            phaseDiff[nPreservedPeaks] = phaseDiff[iOld];
            for (j = lowerIndex[iOld]; j <= upperIndex[iOld]; j++)
            {
                *pNewPhase++ = *pOldPhase++;
            }
            ++nPreservedPeaks;
        }
        else
        {
            pOldPhase += upperIndex[iOld]-lowerIndex[iOld]+1;
        }
    }
    for (iNew = 0; iNew < nPreservedPeaks; iNew++)
    {
        indexOfTonalPeak[iNew] = newIndexOfTonalPeak[iNew];
        lowerIndex[iNew] = newLowerIndex[iNew];
        upperIndex[iNew] = newUpperIndex[iNew];
    }
    *pNumIndexes = nPreservedPeaks;

}

/*****************************************************
---------------- private functions -------------------
******************************************************/

static void calcPseudoSpec(float const * mdctSpec,
                           unsigned int nSamples,
                           float floorPowerSpectrum,
                           float * powerSpec)
{
    unsigned int k;
    float x;

    for (k = 1; k <= nSamples-2; k++)
    {
        x = (mdctSpec[k+1] - mdctSpec[k-1]); /* An MDST estimate */
        x = mdctSpec[k] * mdctSpec[k] + x * x;
        powerSpec[k] = max(floorPowerSpectrum, x);
    }
    powerSpec[0] = 0.5f*powerSpec[1];
    powerSpec[nSamples-1] = 0.5f*powerSpec[nSamples-2];

}

static void getEnvelope(unsigned int nSamples,
                        float const * powerSpec,
                        float F0,
                        float * envelope,
                        float * smoothedSpectrum)
{
    unsigned int nFilterLength, nHalfFilterLength, nSecondHalfFilterLength, n1, n2;
    unsigned int i;
    double sum; /* double is required to avoid precision problems in the optimized version. */


    if (F0 == 0)
    {
        nFilterLength = 15;
    }
    else if (F0 <= 10.0f)
    {
        nFilterLength = 11;
    }
    else if (F0 >= 22.0f)
    {
        /* For F0 >= 22 peak is isolated well enough with the filter length of 23.
           This case is however not triggered due to the limit of pit_min,
           but the line is left for security reasons. */
        nFilterLength = 23;
    }
    else
    {
        nFilterLength = 1+2*(int)(F0/2);
    }
    nHalfFilterLength = nFilterLength/2;
    n1 = nHalfFilterLength+1;
    nSecondHalfFilterLength = nFilterLength-nHalfFilterLength;
    n2 = nSecondHalfFilterLength-1;

    assert((nFilterLength >= 7) && (nFilterLength <= 23) && (nFilterLength %2 == 1));

    sum = 0;
    for (i = 0; i < n2; i++)
    {
        sum += LEVEL_ABOVE_ENVELOPE*powerSpec[i];
    }
    /* No need for PTR_INIT for powerSpec[i+n2] as we continue from the previous loop */
    for (i = 0; i < n1; i++)
    {
        sum += LEVEL_ABOVE_ENVELOPE*powerSpec[i+n2];
        /* 1/(i+nSecondHalfFilterLength) needs to be stored in a table for each filter length */
        envelope[i] = (float)sum / (i+nSecondHalfFilterLength);
    }
    sum /= nFilterLength;                                                          /* 1/nFilterLength is from a table */
    /* No need for PTR_INIT for powerSpec[i+n2] as we continue from the previous loop */
    for (i = n1; i < nSamples-n2; i++)
    {
        sum += LEVEL_ABOVE_ENVELOPE*(powerSpec[i+n2]-powerSpec[i-n1])/nFilterLength; /* LEVEL_ABOVE_ENVELOPE/nFilterLength is from a table */
        envelope[i] = (float)sum;
    }
    sum *= nFilterLength;
    /* No need for PTR_INIT for powerSpec[i-n1] as we continue from the previous loop */
    for (i = nSamples-n2; i < nSamples; i++)
    {
        sum -= LEVEL_ABOVE_ENVELOPE*powerSpec[i-n1];;
        /* 1/(nSamples-(i-nHalfFilterLength)) needs to be stored in a table for each filter length */
        envelope[i] = (float)sum / (nSamples-(i-nHalfFilterLength));
    }
    for (i = 1; i < nSamples-1; i++)
    {
        smoothedSpectrum[i] = 0.75f*powerSpec[i-1]+powerSpec[i]+0.75f*powerSpec[i+1];
    }
    smoothedSpectrum[0] = powerSpec[0]+0.75f*powerSpec[1];
    smoothedSpectrum[nSamples-1] = 0.75f*powerSpec[nSamples-2]+powerSpec[nSamples-1];


}

static void GetF0(unsigned int const nSamples,            /*i*/
                  unsigned int const nSamplesCore,        /*i*/
                  float const * const powerSpectrum,      /*i*/
                  float const pitchLag,                   /*i*/
                  float * const pOrigF0,                  /*i/o*/
                  float * const pF0)                      /*i/o*/
{
    float halfPitchLag;
    unsigned int rgiStrongHarmonics[MAX_PEAKS_FROM_PITCH];
    unsigned int nTotalHarmonics, nStrongHarmonics;


    assert(LAST_HARMONIC_POS_TO_CHECK <= nSamplesCore);

    /* Use only F0 >= 100 Hz */
    if ((pitchLag > 0) && (pitchLag <= 0.5f*nSamplesCore))
    {
        halfPitchLag = 0.5f*pitchLag;
        *pF0 = nSamplesCore/halfPitchLag;
        *pOrigF0 = *pF0;
        if (nSamples < 2*LAST_HARMONIC_POS_TO_CHECK)
        {
            nTotalHarmonics = (unsigned int)(nSamples/(*pF0));
        }
        else
        {
            nTotalHarmonics = (int)(2*LAST_HARMONIC_POS_TO_CHECK/(*pF0)); /* For correcting F0 we go 2 times the last harmonic position that will be used */
        }

        /* Get in rgiStrongHarmonics all i for which i*F0 are the strongest harmonics */
        findStrongestHarmonics(nSamples, powerSpectrum, *pF0, nTotalHarmonics, rgiStrongHarmonics, &nStrongHarmonics);
        CorrectF0(rgiStrongHarmonics, nStrongHarmonics, pF0);
    }
    else
    {
        *pF0 = 0;
        *pOrigF0 = 0;
    }

}

/* This is very fast with almost ordered vectors. */
static void sort(unsigned int *in, unsigned int n)
{
    int i;
    unsigned int j, tempr;

    for (i = n-2; i >= 0; i--)
    {
        tempr = in[i];
        for (j = i+1; (j < n) && (tempr > in[j]); j++ )
        {
            in[j-1] = in[j];
            /* (tempr > r[j]) */
        }
        in[j-1] = tempr;
    }

}

static void findStrongestHarmonics(unsigned int nSamples, float const * powerSpectrum, float F0, unsigned int nTotalHarmonics, unsigned int * pHarmonicIndexes, unsigned int * pnHarmonics)
{
    float peaks[MAX_PEAKS_FROM_PITCH], smallestPeak;
    unsigned int nPeaksToCheck, nPeaks, iSmallestPeak;
    unsigned int i, l, k;
    (void)nSamples;

    nPeaks = 0;
    iSmallestPeak = 0;
    smallestPeak = FLT_MAX;
    nPeaksToCheck = min(nTotalHarmonics, MAX_PEAKS_FROM_PITCH+1);
    for (i = 1; i < nPeaksToCheck; i++)
    {
        float newPeak;
        k = (int)(i*F0);
        assert(k > 0 && k < 2*LAST_HARMONIC_POS_TO_CHECK && k < nSamples);
        newPeak = powerSpectrum[k];
        peaks[nPeaks] = newPeak;
        pHarmonicIndexes[nPeaks] = i;
        if (newPeak <= smallestPeak)
        {
            iSmallestPeak = nPeaks;
            smallestPeak = newPeak;
        }
        ++nPeaks;
    }
    for (; i < nTotalHarmonics; i++)
    {
        float newPeak;
        k = (int)(i*F0);
        assert(k > 0 && k < 2*LAST_HARMONIC_POS_TO_CHECK && k < nSamples);
        newPeak = powerSpectrum[k];
        if (newPeak > smallestPeak)
        {
            peaks[iSmallestPeak] = newPeak;
            pHarmonicIndexes[iSmallestPeak] = i;
            smallestPeak = newPeak;
            for (l = 0; l < MAX_PEAKS_FROM_PITCH; l++)
            {
                if (peaks[l] <= smallestPeak)
                {
                    iSmallestPeak = l;
                    smallestPeak = peaks[l];
                }
            }
        }
    }
    sort(pHarmonicIndexes, nPeaks);
    *pnHarmonics = nPeaks;

}

/* Use new F0, for which harmonics are most common in pHarmonicIndexes */
static void CorrectF0(unsigned int const * pHarmonicIndexes,
                      unsigned int nHarmonics,
                      float * pF0)
{
    unsigned int i;
    float F0;
    unsigned int diff[MAX_PEAKS_FROM_PITCH-1], sortedDiff[MAX_PEAKS_FROM_PITCH-1];
    unsigned int iMostCommonDiff, nMostCommonDiff, nSameDiff, iMult;

    F0 = *pF0;
    if (F0 > 0 && nHarmonics > 0)
    {
        for (i = 0; i < nHarmonics-1; i++)
        {
            diff[i] = pHarmonicIndexes[i+1]-pHarmonicIndexes[i];
            sortedDiff[i] = diff[i];
        }
        sort(sortedDiff, nHarmonics-1);
        iMostCommonDiff = sortedDiff[0];
        nSameDiff = 1;
        i = 1;
        if (sortedDiff[0]*pHarmonicIndexes[0] == 1)
        {
            /* Find how many distances between peaks have length 1 */
            for (; i < nHarmonics-1; i++)
            {
                if (sortedDiff[i] == 1)
                {
                    ++nSameDiff;
                }
            }
        }
        nMostCommonDiff = nSameDiff;

        /* If there are at least 3 distances between peaks with length 1 and if the 1st harmonic is in pHarmonicIndexes then keep the original F0 */
        /* Otherwise find the most common distance between peaks */
        if (nSameDiff < 3)
        {

            /* Find the most common difference */
            for (i = nSameDiff; i < nHarmonics-1; i++)
            {
                if (sortedDiff[i] == sortedDiff[i-1])
                {
                    ++nSameDiff;
                }
                else
                {
                    if (nSameDiff > nMostCommonDiff)
                    {
                        nMostCommonDiff = nSameDiff;
                        iMostCommonDiff = sortedDiff[i-1];
                    }
                    else
                    {
                        if ((nSameDiff == nMostCommonDiff) && (abs((int)iMostCommonDiff-(int)pHarmonicIndexes[0]) > abs((int)sortedDiff[i-1]-(int)pHarmonicIndexes[0])))
                        {
                            nMostCommonDiff = nSameDiff;
                            iMostCommonDiff = sortedDiff[i-1];
                        }
                    }
                    nSameDiff = 1;
                }
            }
            if (nSameDiff > nMostCommonDiff)
            {
                nMostCommonDiff = nSameDiff;
                iMostCommonDiff = sortedDiff[nHarmonics-2];
            }
        }

        /* If there are enough peaks at the same distance */
        if (nMostCommonDiff >= MAX_PEAKS_FROM_PITCH/2)
        {
            iMult = 1;
            for (i = 0; i < nHarmonics-1; i++)
            {
                if (diff[i] == iMostCommonDiff)
                {
                    iMult = pHarmonicIndexes[i];
                    break;
                }
                /* for rare cases of octave mismatch or missing harmonics */
                if ((i < nHarmonics-2) && (diff[i] == diff[i+1]) && (diff[i]+diff[i+1] == iMostCommonDiff))
                {
                    iMult = pHarmonicIndexes[i];
                    break;
                }
            }

            /* If the real F0 is much higher than the original F0 from the pitch */
            if (iMult <= 3)
            {
                /* Use iMostCommonDiff, because the lowest pHarmonicIndexes[i] (which is equal to iMult) may not correspond to the new F0, but to it's multiple */
                F0 = iMostCommonDiff*F0;
            }
            else
            {
                F0 = 0;
            }
        }
        /* Otherwise if there are at least 3 distances between peaks with length 1 and if the 1st harmonic is in pHarmonicIndexes then keep the original F0 */
        /* Otherwise don't use F0 */
        else if ((iMostCommonDiff > 1) || (nMostCommonDiff < 3))
        {
            /* Not enough peaks at the same distance => don't use the pitch. */
            F0 = 0;
        }
        *pF0 = F0;
    }

}

static void modifyThreshold(int i,
                            float F0,
                            float threshold,
                            float * thresholdModification)
{
    float harmonic, fractional, twoTimesFract;
    int k;

    harmonic = i*F0;
    k = (int)harmonic;
    fractional = harmonic - k; /* Fractional part of the i*F0 */
    twoTimesFract = 2.0f*fractional; /* threshold if the centar of the peek is between k-1 and k, threshold+2 if the centar of the peek is between k and k+1 */
    thresholdModification[k]   = threshold;
    thresholdModification[k-1] = threshold +        twoTimesFract;
    thresholdModification[k+1] = threshold + 2.0f - twoTimesFract;

}

static void modifyThresholds(float F0,
                             float origF0,
                             float * thresholdModification)
{
    int i, nHarmonics;



    if ((F0 == 0) && (origF0 > 0))
    {
        nHarmonics = min(MAX_PEAKS_FROM_PITCH, (int)(LAST_HARMONIC_POS_TO_CHECK/origF0));
        for (i = 1; i <= nHarmonics; i++)
        {
            modifyThreshold(i, origF0, 0.7f, thresholdModification);
        }
    }
    else if ((F0 > 0) && (origF0 > 0))
    {
        nHarmonics = min(MAX_PEAKS_FROM_PITCH, (int)(LAST_HARMONIC_POS_TO_CHECK/F0));

        for (i = (int)(F0/origF0+0.5f); i > 0; i--)
        {
            modifyThreshold(i, origF0, 0.35f, thresholdModification);
        }
        for (i = 1; i <= nHarmonics; i++)
        {
            modifyThreshold(i, F0, 0.35f, thresholdModification);
        }
    }

}

static void findCandidates(unsigned int  nSamples,
                           const float * MDCTSpectrum,
                           float       * thresholdModificationNew
                           ,float floorPowerSpectrum             /* i: lower limit for power spectrum bins  */
                          )
{
    float powerSpectrum[L_FRAME_MAX];
    float envelope[L_FRAME_MAX];
    float smoothedSpectrum[L_FRAME_MAX];
    unsigned int upperIdx, lowerIdx;
    unsigned int k, j;

    calcPseudoSpec(MDCTSpectrum, nSamples, floorPowerSpectrum, powerSpectrum);
    getEnvelope(nSamples, powerSpectrum, 0, envelope, smoothedSpectrum);

    set_f( thresholdModificationNew, UNREACHABLE_THRESHOLD, nSamples);
    for (k = GROUP_LENGTH/2 ; k <= nSamples - (GROUP_LENGTH-GROUP_LENGTH/2); k++)
    {
        if (smoothedSpectrum[k] > envelope[k])
        {
            /* The check that bin at k is bigger than bins at k-1 and k+1 is needed to avoid deadlocks when the thresholds are low. */
            /* It removes some true peaks, especially if non weighted sum is used for the smoothed spectrum. */
            float biggerNeighbor;
            biggerNeighbor = max(powerSpectrum[k-1], powerSpectrum[k+1]);
            if (powerSpectrum[k] >= biggerNeighbor)
            {
                /* Find the right foot */
                for (upperIdx = k+1; upperIdx < nSamples-1; upperIdx++)
                {
                    if (powerSpectrum[upperIdx] < powerSpectrum[upperIdx+1])
                    {
                        /* Side lobes may increase for certain ammount */
                        if (ALLOWED_SIDE_LOBE_FLUCTUATION*powerSpectrum[upperIdx] < powerSpectrum[upperIdx+1])
                        {
                            break;
                        }
                        /* Check for further decrease after a side lobe increase */
                        for (j = upperIdx+1; j < nSamples-1; j++)
                        {
                            if (powerSpectrum[j] < ALLOWED_SIDE_LOBE_FLUCTUATION*powerSpectrum[j+1])
                            {
                                break;
                            }
                        }
                        /* Side lobe increase must be 2 times smaller than the decrease to the foot */
                        /* Eq. to 2.0f*powerSpectrum[lowerIdx-1]/powerSpectrum[lowerIdx] > powerSpectrum[lowerIdx]/powerSpectrum[j] */
                        if (2.0f*powerSpectrum[upperIdx+1]*powerSpectrum[j] > powerSpectrum[upperIdx]*powerSpectrum[upperIdx])
                        {
                            break;
                        }
                        upperIdx = j-1;                                                     /* Reinitialize pointers due to the change of upperIdx */
                    }
                }
                /* left foot */
                for (lowerIdx = k-1; lowerIdx > 0; lowerIdx--)
                {
                    if (powerSpectrum[lowerIdx] < powerSpectrum[lowerIdx-1])
                    {
                        /* Side lobes may increase for certain ammount */
                        if (ALLOWED_SIDE_LOBE_FLUCTUATION*powerSpectrum[lowerIdx] < powerSpectrum[lowerIdx-1])
                        {
                            break;
                        }
                        /* Check for further decrease after a side lobe increase */
                        for (j = lowerIdx-1; j > 0; j--)
                        {
                            if (powerSpectrum[j] < ALLOWED_SIDE_LOBE_FLUCTUATION*powerSpectrum[j-1])
                            {
                                break;
                            }
                        }
                        /* Side lobe increase must be 2 times smaller than the decrease to the foot */
                        /* Eq. to 2.0f*powerSpectrum[lowerIdx-1]/powerSpectrum[lowerIdx] > powerSpectrum[lowerIdx]/powerSpectrum[j] */
                        if (2.0f*powerSpectrum[lowerIdx-1]*powerSpectrum[j] > powerSpectrum[lowerIdx]*powerSpectrum[lowerIdx])
                        {
                            break;
                        }
                        lowerIdx = j+1;                                                     /* Reinitialize pointers due to the change of lowerIdx */
                    }
                }

                /* Check if there is a bigger peak up to the next peak foot */
                for (j = max(GROUP_LENGTH/2, lowerIdx); j <= min(upperIdx, nSamples - (GROUP_LENGTH-GROUP_LENGTH/2)); j++)
                {
                    if (powerSpectrum[j] > powerSpectrum[k])
                    {
                        k = j;                                                              /* PTR_INIT for powerSpectrum[k] */
                    }
                }

                /* Modify thresholds for the following frame */
                for (j = k-1; j < k+2; j++)
                {
                    if (smoothedSpectrum[j] > envelope[j])
                    {
                        thresholdModificationNew[j] = SMALL_THRESHOLD;
                    }
                    else
                    {
                        thresholdModificationNew[j] = BIG_THRESHOLD;
                    }
                }
                /* Jump to the next foot of the peak. */
                k = upperIdx;                                                           /* Reinitialize pointers due to the change of k */
            }
        }
    }

}

static void RefineThresholdsUsingPitch(unsigned int nSamples,
                                       unsigned int nSamplesCore,
                                       float const * powerSpectrum,
                                       float lastPitchLag,
                                       float currentPitchLag,
                                       float * pF0,
                                       float * thresholdModification)
{
    int pitchIsStable;
    float origF0;

    pitchIsStable = (fabs(lastPitchLag-currentPitchLag) < 0.25f);
    if (pitchIsStable)
    {
        GetF0(nSamples,
              nSamplesCore,
              powerSpectrum, lastPitchLag, &origF0, pF0);
        modifyThresholds(*pF0, origF0, thresholdModification);
    }
    else
    {
        *pF0 = 0;
    }

}

static void findTonalComponents(unsigned short int * indexOfTonalPeak, /* OUT */
                                unsigned short int * lowerIndex,       /* OUT */
                                unsigned short int * upperIndex,       /* OUT */
                                unsigned int * numIndexes,             /* OUT */
                                unsigned int nSamples,                 /* IN */
                                const float * powerSpectrum,           /* IN */
                                float F0,                              /* IN */
                                float const * thresholdModification)   /* IN */
{
    float envelope[L_FRAME_MAX];
    float smoothedSpectrum[L_FRAME_MAX];

    unsigned int nrOfFIS;
    unsigned int upperIdx, lowerIdx, lowerBound;
    unsigned int k, j;

    getEnvelope(nSamples, powerSpectrum, F0, envelope, smoothedSpectrum);
    nrOfFIS = 0;
    lowerBound = 0;
    for (k = GROUP_LENGTH/2 ; k <= nSamples - (GROUP_LENGTH-GROUP_LENGTH/2); k++)
    {
        if (smoothedSpectrum[k] > envelope[k]*thresholdModification[k])
        {
            /* The check that bin at k is bigger than bins at k-1 and k+1 is needed to avoid deadlocks when the thresholds are low. */
            /* It removes some true peaks, especially if non weighted sum is used for the smoothed spectrum. */
            float biggerNeighbor;
            biggerNeighbor = max(powerSpectrum[k-1], powerSpectrum[k+1]);
            if (powerSpectrum[k] >= biggerNeighbor)
            {
                /* Find the right foot */
                for (upperIdx = k+1; upperIdx < nSamples-1; upperIdx++)
                {
                    if (powerSpectrum[upperIdx] < powerSpectrum[upperIdx+1])
                    {
                        /* Side lobes may increase for certain ammount */
                        if (ALLOWED_SIDE_LOBE_FLUCTUATION*powerSpectrum[upperIdx] < powerSpectrum[upperIdx+1])
                        {
                            break;
                        }
                        /* Check for further decrease after a side lobe increase */
                        for (j = upperIdx+1; j < nSamples-1; j++)
                        {
                            if (powerSpectrum[j] < ALLOWED_SIDE_LOBE_FLUCTUATION*powerSpectrum[j+1])
                            {
                                break;
                            }
                        }
                        /* Side lobe increase must be 2 times smaller than the decrease to the foot */
                        /* Eq. to 2.0f*powerSpectrum[lowerIdx-1]/powerSpectrum[lowerIdx] > powerSpectrum[lowerIdx]/powerSpectrum[j] */
                        if (2.0f*powerSpectrum[upperIdx+1]*powerSpectrum[j] > powerSpectrum[upperIdx]*powerSpectrum[upperIdx])
                        {
                            break;
                        }
                        upperIdx = j-1;                                                     /* Reinitialize pointers due to the change of upperIdx */
                    }
                }
                /* left foot */
                for (lowerIdx = k-1; lowerIdx > lowerBound; lowerIdx--)
                {
                    if (powerSpectrum[lowerIdx] < powerSpectrum[lowerIdx-1])
                    {
                        /* Side lobes may increase for certain ammount */
                        if (ALLOWED_SIDE_LOBE_FLUCTUATION*powerSpectrum[lowerIdx] < powerSpectrum[lowerIdx-1])
                        {
                            break;
                        }
                        /* Check for further decrease after a side lobe increase */
                        for (j = lowerIdx-1; j > 0; j--)
                        {
                            if (powerSpectrum[j] < ALLOWED_SIDE_LOBE_FLUCTUATION*powerSpectrum[j-1])
                            {
                                break;
                            }
                        }
                        /* Side lobe increase must be 2 times smaller than the decrease to the foot */
                        /* Eq. to 2.0f*powerSpectrum[lowerIdx-1]/powerSpectrum[lowerIdx] > powerSpectrum[lowerIdx]/powerSpectrum[j] */
                        if (2.0f*powerSpectrum[lowerIdx-1]*powerSpectrum[j] > powerSpectrum[lowerIdx]*powerSpectrum[lowerIdx])
                        {
                            break;
                        }
                        lowerIdx = j+1;                                                     /* Reinitialize pointers due to the change of lowerIdx */
                    }
                }
                lowerBound = upperIdx;

                /* Check if there is a bigger peak up to the next peak foot */
                for (j = max(GROUP_LENGTH/2, lowerIdx); j <= min(upperIdx, nSamples - (GROUP_LENGTH-GROUP_LENGTH/2)); j++)
                {
                    if (powerSpectrum[j] > powerSpectrum[k])
                    {
                        k = j;                                                              /* PTR_INIT for powerSpectrum[k] */
                    }
                }

                assert((nrOfFIS == 0) || (indexOfTonalPeak[nrOfFIS-1] < k));
                lowerIndex[nrOfFIS] = k-GROUP_LENGTH/2;
                upperIndex[nrOfFIS] = k+(GROUP_LENGTH-GROUP_LENGTH/2-1);
                if ((nrOfFIS > 0) && (lowerIndex[nrOfFIS] <= upperIndex[nrOfFIS-1]))
                {
                    int m = (k+indexOfTonalPeak[nrOfFIS-1])/2;
                    upperIndex[nrOfFIS-1] = m;
                    lowerIndex[nrOfFIS] = m+1;
                }
                indexOfTonalPeak[nrOfFIS++] = k;
                if (nrOfFIS == MAX_NUMBER_OF_IDX)
                {
                    break;
                }
                /* Jump to the next foot of the peak. */
                k = upperIdx;                                                           /* Reinitialize pointers due to the change of k */
            }
        }
    }
    *numIndexes = nrOfFIS;

}

