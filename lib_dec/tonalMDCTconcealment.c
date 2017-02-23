/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#define _USE_MATH_DEFINES

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "options.h"
#include "prot.h"




/************************************************************************************/
/* forward declarations for local functions, see implementation at end of this file */
/************************************************************************************/

static void CalcMDXT(TonalMDCTConcealPtr         const self,
                     char                        const type,
                     float               const * const timeSignal,
                     float                     * const mdxtOutput);

static void CalcPowerSpec(float const * mdctSpec,
                          float const * mdstSpec,
                          unsigned int nSamples,
                          float floorPowerSpectrum,
                          float * powerSpec);

static void CalcPowerSpecAndDetectTonalComponents(TonalMDCTConcealPtr const self,
        float secondLastMDST[],
        float secondLastMDCT[],
        float const pitchLag);

static void FindPhases(TonalMDCTConcealPtr const self,
                       float const secondLastMDCT[],
                       float const secondLastMDST[]);

static void FindPhaseDifferences(TonalMDCTConcealPtr const self,
                                 float powerSpectrum[]);

/*******************************************************/
/*-------------- public functions -------------------- */
/*******************************************************/

TONALMDCTCONCEAL_ERROR TonalMDCTConceal_Init( TonalMDCTConcealPtr self,
        unsigned int nSamples,
        unsigned int nSamplesCore,
        unsigned int nScaleFactors,
        TCX_config * tcx_cfg
                                            )
{
    if (nSamples > L_FRAME_MAX || nScaleFactors > FDNS_NPTS)
    {
        assert(nSamples <= L_FRAME_MAX);
        assert(nScaleFactors <= FDNS_NPTS);
        return TONALMDCTCONCEAL_NSAMPLES_LARGER_THAN_MAXBLOCKSIZE;
    }
    assert((self->nScaleFactors == nScaleFactors) || (self->nSamples != nSamples)); /* If nSamples doesn't change then also nScaleFactors must stay the same */

    self->tcx_cfg = tcx_cfg;
    self->lastBlockData.spectralData       = self->spectralDataBuffers[0];
    self->secondLastBlockData.spectralData = self->spectralDataBuffers[1];
    self->secondLastPowerSpectrum = self->secondLastBlockData.spectralData;
    self->lastBlockData.scaleFactors       = self->scaleFactorsBuffers[0];
    self->secondLastBlockData.scaleFactors = self->scaleFactorsBuffers[1];
    self->lastBlockData.blockIsValid       = 0;
    self->secondLastBlockData.blockIsValid = 0;
    self->nSamples = 0;
    self->nScaleFactors = 0;

    self->lastBlockData.blockIsConcealed       = 0;
    self->secondLastBlockData.blockIsConcealed = 0;
    self->pTCI = (TonalComponentsInfo *)self->timeDataBuffer;

    self->lastPitchLag = 0;

    if (self->nSamples != nSamples)
    {
        self->secondLastBlockData.blockIsValid = 0;
        self->lastBlockData.blockIsValid = 0;
    }
    self->nSamples = nSamples;
    self->nSamplesCore = nSamplesCore;
    self->nScaleFactors = nScaleFactors;

    /* Offset the pointer to the end of buffer, so that pTCI is not destroyed when
       new time samples are stored in lastPcmOut */
    /* just the second half of the second last pcm output is needed */
    self->secondLastPcmOut        = &self->timeDataBuffer[(3*L_FRAME_MAX)/2-(3*min(L_FRAME_MAX, nSamples))/2];
    self->lastPcmOut              = &self->timeDataBuffer[(3*L_FRAME_MAX)/2-min(L_FRAME_MAX, nSamples)];

    /* If the second last frame was lost and concealed with tonal PLC, we
       reuse saved TonalComponentsInfo and don't update pcm buffers */
    assert(sizeof(*self->pTCI) <= (self->lastPcmOut-self->timeDataBuffer)*sizeof(self->timeDataBuffer[0]));

    return TONALMDCTCONCEAL_OK;
}


TONALMDCTCONCEAL_ERROR TonalMDCTConceal_SaveFreqSignal( TonalMDCTConcealPtr self,
        float const *mdctSpectrum,
        unsigned int nNewSamples,
        unsigned int nNewSamplesCore,
        float const *scaleFactors
                                                      )
{
    float * temp;
    int nOldSamples;

    assert(nNewSamples > 0 && nNewSamples <= 2*L_FRAME_MAX);

    /* Avoid overwriting self->secondLastPowerSpectrum stored in spectralData,
       because it is needed if the second last and the current frame are lost
       and concealed using the Tonal MDCT PLC */
    if (!self->lastBlockData.tonalConcealmentActive || (self->lastBlockData.nSamples != nNewSamples))
    {
        if (nNewSamples <= L_FRAME_MAX)
        {
            /* Shift the buffers */
            temp = self->secondLastBlockData.spectralData; /* Save the pointer */
            self->secondLastBlockData.spectralData = self->lastBlockData.spectralData;
            self->lastBlockData.spectralData = temp;
            temp = self->secondLastBlockData.scaleFactors;
            self->secondLastBlockData.scaleFactors = self->lastBlockData.scaleFactors;
            self->lastBlockData.scaleFactors = temp;
        }
        else
        {
            /* Order the buffers so that even transition frame can fit in if written into the first buffer */
            self->lastBlockData.spectralData       = self->spectralDataBuffers[0];
            self->secondLastBlockData.spectralData = self->spectralDataBuffers[1];
            self->lastBlockData.scaleFactors       = self->scaleFactorsBuffers[0];
            self->secondLastBlockData.scaleFactors = self->scaleFactorsBuffers[1];
        }
        nOldSamples = self->lastBlockData.nSamples;
        self->lastBlockData.nSamples = nNewSamples;
        self->secondLastBlockData.nSamples = nOldSamples;
        nOldSamples = self->lastBlockData.nSamplesCore;
        self->lastBlockData.nSamplesCore = nNewSamplesCore;
        self->secondLastBlockData.nSamplesCore = nOldSamples;
    }
    if ((nNewSamples > 0) && (nNewSamples <= 2*L_FRAME_MAX))
    {
        /* Store new data */
        mvr2r(mdctSpectrum, self->lastBlockData.spectralData, nNewSamples);
        mvr2r(scaleFactors, self->lastBlockData.scaleFactors, self->nScaleFactors);
    }
    return TONALMDCTCONCEAL_OK;
}


TONALMDCTCONCEAL_ERROR TonalMDCTConceal_UpdateState(TonalMDCTConcealPtr self,
        int nNewSamples,
        float pitchLag,
        int badBlock,
        int tonalConcealmentActive
                                                   )
{
    int newBlockIsValid;


    assert(!(!badBlock && tonalConcealmentActive));

    if (badBlock)
    {
        newBlockIsValid = self->lastBlockData.blockIsValid;
    }
    else
    {
        newBlockIsValid = (nNewSamples <= 2*L_FRAME_MAX) && (nNewSamples > 0);
    }

    /* Shift old state */
    self->secondLastBlockData.blockIsConcealed = self->lastBlockData.blockIsConcealed;
    self->secondLastBlockData.blockIsValid = self->lastBlockData.blockIsValid;
    self->secondLastBlockData.tonalConcealmentActive = self->lastBlockData.tonalConcealmentActive;

    /* Store new state */
    self->lastBlockData.blockIsConcealed = badBlock;
    self->lastBlockData.blockIsValid = newBlockIsValid;
    self->lastBlockData.tonalConcealmentActive = tonalConcealmentActive;
    self->lastPitchLag = pitchLag;

    return TONALMDCTCONCEAL_OK;
}


static void FindPhases(TonalMDCTConcealPtr const self, float const secondLastMDCT[], float const secondLastMDST[])
{
    unsigned int i;
    int l;
    float * pCurrentPhase;


    pCurrentPhase = self->pTCI->phase_currentFramePredicted;
    /* for each index/index group */
    for( i = 0; i < self->pTCI->numIndexes; i++)
    {
        for (l = self->pTCI->lowerIndex[i]; l <= self->pTCI->upperIndex[i]; l++)
        {
            *pCurrentPhase++ = (float)atan2(secondLastMDST[l], secondLastMDCT[l]);
        }
    }

    return;
}


static void FindPhaseDifferences(TonalMDCTConcealPtr const self, float powerSpectrum[])
{
    static float const bandwidth = 7.0f;
    float const m = (float)cos(EVS_PI/bandwidth);
    float const s = (float)cos((3*EVS_PI)/bandwidth);
    float const n = (float)sin(EVS_PI/bandwidth);
    float const j = (float)sin((3*EVS_PI)/bandwidth);
    static float const G = (float)(1.0/(2*1.36));
    static float const maxRatio = 44.8f; /* Maximum ratio |ODFT[k-1]|/|ODFT[k+1]| is 16.5 dB, that is maximum ratio (for fractional = 0) is (cos(EVS_PI/bandwidth)/cos(3PI/bandwidth))^1.36 */

    unsigned int i, k;
    float odft_left, odft_right;
    float * phaseDiff;
    float fractional;
    float Q, a;

    phaseDiff = self->pTCI->phaseDiff;
    for (i = 0; i < self->pTCI->numIndexes; i++)
    {
        k = self->pTCI->indexOfTonalPeak[i];
        odft_left  = powerSpectrum[k-1];
        odft_right = powerSpectrum[k+1];
        if (odft_left >= maxRatio*odft_right)
        {
            a = (float)tan(0.0f*EVS_PI/bandwidth);
        }
        else
        {
            if (odft_right >= maxRatio*odft_left)
            {
                a = (float)tan(2.0f*EVS_PI/bandwidth);
            }
            else
            {
                Q = (float)pow(odft_left/odft_right, G);
                a = (m - Q * s) / (n + Q * j);
            }
        }
        fractional = (float)atan(a) * (bandwidth/2.0f);
        assert((fractional >= 0) && (fractional <= EVS_PI + 1.192092896e-07F));
        phaseDiff[i] = fractional + EVS_PI*(k%4);
    }

    return;
}


static void CalcPowerSpecAndDetectTonalComponents(TonalMDCTConcealPtr const self,
        float secondLastMDST[],
        float secondLastMDCT[],
        float const pitchLag)
{
    unsigned int nSamples;
    unsigned int i;
    float floorPowerSpectrum; /* Minimum significant value of a spectral line in the power spectrum */
    float powerSpectrum[L_FRAME_MAX]; /* 32 bits are required */
    float invScaleFactors[FDNS_NPTS];

    nSamples = self->nNonZeroSamples;

    /* It is taken into account that the MDCT is not normalized. */
    floorPowerSpectrum = self->nSamples*self->nSamples/400.0f;
    CalcPowerSpec(secondLastMDCT, secondLastMDST, nSamples, floorPowerSpectrum, powerSpectrum );

    /* This setting to minimal level is required because the power spectrum is used in the threshold adaptation using the pitch up to self->nSamples. */
    set_f(powerSpectrum+nSamples, floorPowerSpectrum, self->nSamples-nSamples);
    /* this setting to zero is needed since the FDNS needs to be called
       with self->nSamplesCore; it relevant only for nb; it has no effect
       to the output, but memory checker may complain otherwise due to the
       usage of uninitialized values */
    if (self->nSamplesCore > self->nSamples)
    {
        set_zero(powerSpectrum+self->nSamples, self->nSamplesCore-self->nSamples);
    }
    DetectTonalComponents(self->pTCI->indexOfTonalPeak,
                          self->pTCI->lowerIndex,
                          self->pTCI->upperIndex,
                          &self->pTCI->numIndexes,
                          self->lastPitchLag,
                          pitchLag,
                          self->lastBlockData.spectralData,
                          self->lastBlockData.scaleFactors,
                          powerSpectrum,
                          nSamples,
                          self->nSamplesCore,
                          floorPowerSpectrum );
    FindPhases(self, secondLastMDCT, secondLastMDST);
    FindPhaseDifferences(self, powerSpectrum);
    if (self->pTCI->numIndexes > 0)
    {
        self->secondLastPowerSpectrum = self->secondLastBlockData.spectralData;
        for ( i=0; i<nSamples; i++)
        {
            powerSpectrum[i] = (float) sqrt(powerSpectrum[i]);
        }
        for (i = 0; i < self->nScaleFactors; i++)
        {
            invScaleFactors[i] = 1.0f/self->secondLastBlockData.scaleFactors[i];
        }
        mdct_noiseShaping(powerSpectrum, self->nSamplesCore, invScaleFactors);
        v_multc( powerSpectrum + self->nSamplesCore, invScaleFactors[FDNS_NPTS-1], powerSpectrum + self->nSamplesCore, self->nSamples - self->nSamplesCore);
        mvr2r( powerSpectrum, self->secondLastPowerSpectrum, self->nSamples); /* 16 bits are now enough for storing the power spectrum */
    }

    return;
}


static void CalcMDXT(TonalMDCTConcealPtr         const self,
                     char                        const type,
                     float               const * const timeSignal,
                     float                     * const mdxtOutput)
{
    float windowedTimeSignal[L_FRAME_PLUS+2*L_MDCT_OVLP_MAX];
    int left_overlap, right_overlap, L_frame;

    L_frame = self->nSamples;
    WindowSignal( self->tcx_cfg, self->tcx_cfg->tcx_offsetFB, FULL_OVERLAP, FULL_OVERLAP,
                  &left_overlap, &right_overlap, timeSignal, &L_frame, windowedTimeSignal, 1 );
    if (type == 'S')
    {
        TCX_MDST( windowedTimeSignal, mdxtOutput, left_overlap, L_frame - (left_overlap+right_overlap)/2, right_overlap );
    }
    else
    {
        TCX_MDCT( windowedTimeSignal, mdxtOutput, left_overlap, L_frame - (left_overlap+right_overlap)/2, right_overlap );
    }

    return;
}


TONALMDCTCONCEAL_ERROR TonalMDCTConceal_Detect( TonalMDCTConcealPtr const self,
        float const pitchLag,
        int * const numIndices)
{
    float secondLastMDST[L_FRAME_MAX]; /* 32 bits are required */
    float secondLastMDCT[L_FRAME_MAX]; /* 32 bits are required */
    float * powerSpectrum = secondLastMDST;
    unsigned int nSamples;
    unsigned int i;

    nSamples = self->nSamples;
    if (self->lastBlockData.blockIsValid
            && self->secondLastBlockData.blockIsValid
            && (self->lastBlockData.nSamples == nSamples)
            && (self->secondLastBlockData.nSamples == nSamples)
            && (!self->secondLastBlockData.blockIsConcealed
                || self->secondLastBlockData.tonalConcealmentActive
                || (pitchLag != 0))
            /* Safety if the second last frame was concealed and tonal concealment was inactive */
       )
    {
        if (!self->lastBlockData.blockIsConcealed)
        {
            if (!self->secondLastBlockData.tonalConcealmentActive)
            {
                CalcMDXT(self, 'S', self->secondLastPcmOut, secondLastMDST);
                CalcMDXT(self, 'C', self->secondLastPcmOut, secondLastMDCT);
                self->nNonZeroSamples = 0;
                for (i = 0; i < self->nSamples; i++)
                {
                    if (self->secondLastBlockData.spectralData[i] != 0)
                    {
                        self->nNonZeroSamples = i;
                    }
                }
                /* 23 is the maximum length of the MA filter in getEnvelope */
                self->nNonZeroSamples = min(self->nSamples, self->nNonZeroSamples+23);   ;
                CalcPowerSpecAndDetectTonalComponents(self, secondLastMDST, secondLastMDCT, pitchLag);
            }
            else
            {
                /* If the second last frame was also lost, it is expected that pastTimeSignal could hold a bit different signal (e.g. including fade-out) from the one stored in TonalMDCTConceal_SaveTimeSignal. */
                /* That is why we reuse the already stored information about the concealed spectrum in the second last frame */
                nSamples = self->nNonZeroSamples;
                mvr2r(self->secondLastPowerSpectrum, powerSpectrum, nSamples); /* Convert from 16 bits to 32 bits */
                mdct_noiseShaping(powerSpectrum, self->nSamplesCore, self->secondLastBlockData.scaleFactors);
                v_multc(powerSpectrum + self->nSamplesCore,
                        self->secondLastBlockData.scaleFactors[FDNS_NPTS-1],
                        powerSpectrum + self->nSamplesCore,
                        nSamples - self->nSamplesCore);
                v_mult(powerSpectrum, powerSpectrum, powerSpectrum, nSamples);
                RefineTonalComponents(self->pTCI->indexOfTonalPeak,
                                      self->pTCI->lowerIndex,
                                      self->pTCI->upperIndex,
                                      self->pTCI->phaseDiff,
                                      self->pTCI->phase_currentFramePredicted,
                                      &self->pTCI->numIndexes,
                                      self->lastPitchLag,
                                      pitchLag,
                                      self->lastBlockData.spectralData,
                                      self->lastBlockData.scaleFactors,
                                      powerSpectrum,
                                      nSamples,
                                      self->nSamplesCore,
                                      self->nSamples*self->nSamples/400.0f /* floorPowerSpectrum */ );
            }
        }
    }
    else
    {
        self->pTCI->numIndexes = 0;
    }
    *numIndices = self->pTCI->numIndexes;


    return TONALMDCTCONCEAL_OK;
}


TONALMDCTCONCEAL_ERROR TonalMDCTConceal_InsertNoise( TonalMDCTConcealPtr self,      /*IN */
        float* mdctSpectrum,           /*OUT*/
        int    tonalConcealmentActive,
        short* pSeed,
        float  tiltCompFactor,
        float  crossfadeGain,
        int    crossOverFreq)
{
    unsigned int i;
    Word16 rnd;
    float g, nrgNoiseInLastFrame, nrgWhiteNoise, tiltFactor, tilt;


    g = 1.0f-crossfadeGain;
    if (!self->lastBlockData.blockIsConcealed)
    {
        rnd = 1977;
    }
    else
    {
        rnd = *pSeed;
    }
    if (!self->lastBlockData.blockIsValid)
    {
        /* may just become active if the very first frame is lost */
        set_f( mdctSpectrum, 0.0f, self->nSamples);
    }
    else
    {
        /* based on what is done in tcx_noise_filling() */
        tiltFactor = (float)pow(max(0.375f, tiltCompFactor), 1.0f/self->lastBlockData.nSamples);
        tilt = 1.0f;
        nrgNoiseInLastFrame = nrgWhiteNoise = 0.0f;
        if (!tonalConcealmentActive)
        {
            for (i = 0; i < (unsigned int)crossOverFreq; i++)
            {
                float const x = self->lastBlockData.spectralData[i];
                float y;
                rnd = (Word16) (rnd * 31821L + 13849L);
                y = tilt*rnd;
                nrgNoiseInLastFrame += x*x;
                nrgWhiteNoise += y*y;
                mdctSpectrum[i] = y;
                tilt *= tiltFactor;
            }
            if (nrgWhiteNoise > 0)
            {
                g *= (float)sqrt(nrgNoiseInLastFrame/nrgWhiteNoise);
            }
            for (i = 0; i < (unsigned int)crossOverFreq; i++)
            {
                float const x = self->lastBlockData.spectralData[i];
                float const y = mdctSpectrum[i];

                if (y > 0)
                {
                    mdctSpectrum[i] = g*y + crossfadeGain*x;
                }
                else
                {
                    mdctSpectrum[i] = g*y - crossfadeGain*x;
                }
            }

            for (i = (unsigned int)crossOverFreq; i < self->lastBlockData.nSamples; i++)
            {
                mdctSpectrum[i] = self->lastBlockData.spectralData[i];
            }
        }
        else
        {
            unsigned int l;
            assert(self->pTCI->numIndexes > 0);
            for (l = 0; l < self->pTCI->lowerIndex[0]; l++)
            {
                float const x = self->lastBlockData.spectralData[l];
                float y;
                rnd = (Word16) (rnd * 31821L + 13849L);
                y = tilt*rnd;
                nrgNoiseInLastFrame += x*x;
                nrgWhiteNoise += y*y;
                mdctSpectrum[l] = y;
                tilt *= tiltFactor;
            }
            for (i = 1; i < self->pTCI->numIndexes; i++)
            {
                tilt *= (float)pow(tiltFactor, self->pTCI->upperIndex[i-1]-self->pTCI->lowerIndex[i-1]+1);
                for (l = self->pTCI->upperIndex[i-1]+1; l < self->pTCI->lowerIndex[i]; l++)
                {
                    float const x = self->lastBlockData.spectralData[l];
                    float y;
                    rnd = (Word16) (rnd * 31821L + 13849L);
                    y = tilt*rnd;
                    nrgNoiseInLastFrame += x*x;
                    nrgWhiteNoise += y*y;
                    mdctSpectrum[l] = y;
                    tilt *= tiltFactor;
                }
            }
            tilt *= (float)pow(tiltFactor, self->pTCI->upperIndex[self->pTCI->numIndexes-1]-self->pTCI->lowerIndex[self->pTCI->numIndexes-1]+1);

            for (l = self->pTCI->upperIndex[self->pTCI->numIndexes-1]+1; l < (unsigned int)crossOverFreq; l++)
            {
                float const x = self->lastBlockData.spectralData[l];
                float y;
                rnd = (Word16) (rnd * 31821L + 13849L);
                y = tilt*rnd;
                nrgNoiseInLastFrame += x*x;
                nrgWhiteNoise += y*y;
                mdctSpectrum[l] = y;
                tilt *= tiltFactor;
            }
            if (nrgWhiteNoise > 0)
            {
                g *= (float)sqrt(nrgNoiseInLastFrame/nrgWhiteNoise);
            }
            for (l = 0; l < self->pTCI->lowerIndex[0]; l++)
            {
                float const x = self->lastBlockData.spectralData[l];
                float const y = mdctSpectrum[l];

                if (y > 0)
                {
                    mdctSpectrum[l] = g*y + crossfadeGain*x;
                }
                else
                {
                    mdctSpectrum[l] = g*y - crossfadeGain*x;
                }
            }
            for (i = 1; i < self->pTCI->numIndexes; i++)
            {
                for (l = self->pTCI->upperIndex[i-1]+1; l < self->pTCI->lowerIndex[i]; l++)
                {
                    float const x = self->lastBlockData.spectralData[l];
                    float const y = mdctSpectrum[l];

                    if (y > 0)
                    {
                        mdctSpectrum[l] = g*y + crossfadeGain*x;
                    }
                    else
                    {
                        mdctSpectrum[l] = g*y - crossfadeGain*x;
                    }
                }
            }
            /* initialize bins of tonal components with zero: basically not
               necessary, but currently the whole spectrum is rescaled in
               mdct_noiseShaping() and then there would be a processing of
               uninitialized values */
            for (i = 0; i < self->pTCI->numIndexes; i++)
            {
                for (l = self->pTCI->lowerIndex[i]; l <= self->pTCI->upperIndex[i]; l++)
                {
                    mdctSpectrum[l] = 0;

                }
            }
            for (l = self->pTCI->upperIndex[self->pTCI->numIndexes-1]+1; l < (unsigned int)crossOverFreq; l++)
            {
                float const x = self->lastBlockData.spectralData[l];
                float const y = mdctSpectrum[l];

                if (y > 0)
                {
                    mdctSpectrum[l] = g*y + crossfadeGain*x;
                }
                else
                {
                    mdctSpectrum[l] = g*y - crossfadeGain*x;
                }
            }
            for (l = (unsigned int)crossOverFreq; l < self->lastBlockData.nSamples; l++)
            {
                mdctSpectrum[l] = self->lastBlockData.spectralData[l];
            }
        }
    }

    *pSeed = rnd;

    return TONALMDCTCONCEAL_OK;
}


TONALMDCTCONCEAL_ERROR TonalMDCTConceal_Apply(TonalMDCTConcealPtr self,     /*IN */
        float *mdctSpectrum           /*OUT */
                                             )
{
    unsigned int i, l;
    float * phaseDiff, * pCurrentPhase;
    float phaseToAdd;
    float powerSpectrum[L_FRAME_MAX];
    unsigned int nSamples;


    if (self->lastBlockData.blockIsValid & self->secondLastBlockData.blockIsValid)
    {
        assert(self->pTCI->numIndexes > 0);

        nSamples = self->nNonZeroSamples;
        assert(self->pTCI->upperIndex[self->pTCI->numIndexes-1] < nSamples);
        mvr2r(self->secondLastPowerSpectrum, powerSpectrum, nSamples); /* Convert from 16 bits to 32 bits */
        mdct_noiseShaping(powerSpectrum, self->nSamplesCore, self->secondLastBlockData.scaleFactors);
        v_multc( powerSpectrum + self->nSamplesCore, self->secondLastBlockData.scaleFactors[FDNS_NPTS-1], powerSpectrum + self->nSamplesCore, nSamples - self->nSamplesCore );
        phaseDiff = self->pTCI->phaseDiff; /* if multiple frame loss occurs use the phase from the last frame and continue rotating */
        pCurrentPhase = self->pTCI->phase_currentFramePredicted;
        if (!self->lastBlockData.blockIsConcealed)
        {
            if (self->secondLastBlockData.tonalConcealmentActive)
            {
                self->nFramesLost += 1;
            }
            else
            {
                self->nFramesLost = 1.5;
            }
        }
        /* for each index group */
        for (i = 0; i < self->pTCI->numIndexes; i++)
        {
            phaseToAdd = self->nFramesLost*phaseDiff[i];
            /* Move phaseToAdd to range -EVS_PI..EVS_PI */
            while (phaseToAdd > EVS_PI)
            {
                phaseToAdd -= 2*EVS_PI;
            }
            while (phaseToAdd < -EVS_PI)
            {
                /* should never occur in flt - kept for safety reasons */
                phaseToAdd += 2*EVS_PI;
            }
            for (l = self->pTCI->lowerIndex[i]; l <= self->pTCI->upperIndex[i]; l++)
            {
                float const currentPhase = (*pCurrentPhase++) + phaseToAdd;              /* *pCurrentPhase and phaseToAdd are in range -EVS_PI..EVS_PI */
                mdctSpectrum[l] = (float)cos(currentPhase) * powerSpectrum[l];
            }
        }
    }
    self->nFramesLost++;

    return TONALMDCTCONCEAL_OK;
}


TONALMDCTCONCEAL_ERROR TonalMDCTConceal_SaveTimeSignal( TonalMDCTConcealPtr self,
        float*              timeSignal,
        unsigned int        nNewSamples
                                                      )
{
    if (nNewSamples == self->nSamples)
    {
        assert(nNewSamples <= L_FRAME_MAX);
        if (!self->secondLastBlockData.tonalConcealmentActive)
        {
            mvr2r(self->lastPcmOut + self->nSamples/2, self->secondLastPcmOut, self->nSamples/2);
        }
        mvr2r(timeSignal, self->lastPcmOut, self->nSamples);
    }
    return TONALMDCTCONCEAL_OK;
}


static void CalcPowerSpec(float const * mdctSpec, float const * mdstSpec,
                          unsigned int nSamples, float floorPowerSpectrum,
                          float * powerSpec)
{
    unsigned int k;
    float x;

    for (k = 1; k <= nSamples-2; k++)
    {
        x = mdctSpec[k] * mdctSpec[k] + mdstSpec[k] * mdstSpec[k];
        powerSpec[k] = max(floorPowerSpectrum, x);
    }
    powerSpec[0] = 0.5f*powerSpec[1];
    powerSpec[nSamples-1] = 0.5f*powerSpec[nSamples-2];


    return;
}

