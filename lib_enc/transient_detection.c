/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "stat_enc.h"
#include "cnst.h"
#include "prot.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>



#define MIN_BLOCK_ENERGY 107.37f


/************************************************/
/*                                              */
/*        Internal functions prototypes         */
/*                                              */
/************************************************/

static void InitDelayBuffer(int nFrameLength, int nDelay, DelayBuffer * pDelayBuffer);
static void InitSubblockEnergies(int nFrameLength, int nDelay, DelayBuffer * pDelayBuffer, SubblockEnergies * pSubblockEnergies);
static void InitTransientDetector(SubblockEnergies * pSubblockEnergies, int nDelay, int nSubblocksToCheck,
                                  TCheckSubblocksForAttack pCheckSubblocksForAttack,
                                  float attackRatioThreshold, TransientDetector * pTransientDetector);
static void UpdateDelayBuffer(float const * input, int nSamplesAvailable, DelayBuffer * pDelayBuffer);
static void HighPassFilter(float const * input, int length, float * pFirState1, float * pFirState2, float * output);
static void UpdateSubblockEnergies(float const * input, int nSamplesAvailable, SubblockEnergies * pSubblockEnergies);
static void CalculateSubblockEnergies(float const * input, int nSamplesAvailable, SubblockEnergies * pSubblockEnergies);
static void RunTransientDetector(TransientDetector * pTransientDetector);

/************************************************/
/*                                              */
/*  Functions that define transient detectors   */
/*                                              */
/************************************************/

/** TCX decision.
  * Check if there is an attack in a subblock. Version for TCX Long/Short decision.
  * See TCheckSubblocksForAttack for definition of parameters.
  * It is assumed that the delay of MDCT overlap was not taken into account, so that the last subblock corresponds to the newest input subblock.
  */
static void GetAttackForTCXDecision(float const * pSubblockNrg, float const * pAccSubblockNrg, int nSubblocks, int nPastSubblocks, float attackRatioThreshold, int * pbIsAttackPresent, int * pAttackIndex)
{
    int i;
    int bIsAttackPresent, attackIndex;

    (void)nPastSubblocks;
    (void)nSubblocks;
    assert(nSubblocks >= NSUBBLOCKS);
    assert(nPastSubblocks >= 2);

    bIsAttackPresent = FALSE;
    attackIndex = -1;
    /* Search for the last attack in the subblocks */
    if ((pSubblockNrg[-1] > pAccSubblockNrg[-1] * attackRatioThreshold)
            || (pSubblockNrg[-2] > pAccSubblockNrg[-2] * attackRatioThreshold)
       )
    {
        bIsAttackPresent = TRUE;
        attackIndex = 0;
    }
    /* PTR_INIT for pSubblockNrg[i-1] */
    for (i = 0; i < NSUBBLOCKS; i++)
    {
        if (pSubblockNrg[i] > pAccSubblockNrg[i] * attackRatioThreshold)
        {
            if (i < 6)
            {
                bIsAttackPresent = TRUE;
            }
            if ((attackIndex != 2) && (attackIndex != 6))
            {
                attackIndex = i;
                if ((pSubblockNrg[i] < pAccSubblockNrg[i] * 1.125f * attackRatioThreshold) && (i == 2 || i == 6))
                {
                    attackIndex++;  /* avoid minimum overlap to prevent clicks */
                }
            }
        }
        else     /* no attack, but set index anyway in case of strong energy increase */
        {
            if ((pSubblockNrg[i] > pSubblockNrg[i-1] * 1.5f * attackRatioThreshold) &&
                    (pSubblockNrg[i] > pSubblockNrg[i-2] * 1.5f * attackRatioThreshold))
            {
                if ((attackIndex != 2) && (attackIndex != 6))
                {
                    attackIndex = i;

                    if (((pSubblockNrg[i] < pSubblockNrg[i-1] * 2.0f * attackRatioThreshold) ||
                            (pSubblockNrg[i] < pSubblockNrg[i-2] * 2.0f * attackRatioThreshold)) && (i == 2 || i == 6))
                    {
                        attackIndex++;  /* avoid minimum overlap to prevent clicks */
                    }
                }
            }
        }
    }
    /* avoid post-echos on click sounds (very short transients) due to TNS aliasing */
    if (attackIndex == 4)
    {
        attackIndex = 7;
    }
    else if (attackIndex == 5)
    {
        attackIndex = 6;
    }
    *pAttackIndex = attackIndex;
    *pbIsAttackPresent = bIsAttackPresent;
}

/** Initialize TCX transient detector.
  * See InitTransientDetector for definition of parameters.
  */
static void InitTCXTransientDetector(int nDelay, SubblockEnergies * pSubblockEnergies, TransientDetector * pTransientDetector)
{
    InitTransientDetector(pSubblockEnergies, nDelay, NSUBBLOCKS, GetAttackForTCXDecision, 8.5f, pTransientDetector);
}

/************************************************/
/*                                              */
/*              Interface functions             */
/*                                              */
/************************************************/

void InitTransientDetection(int nFrameLength,
                            int nTCXDelay,
                            TransientDetection * pTransientDetection)
{
    /* Init the delay buffer. */
    InitDelayBuffer(nFrameLength, nTCXDelay, &pTransientDetection->delayBuffer);
    /* Init a subblock energies buffer used for the TCX Short/Long decision. */
    InitSubblockEnergies(nFrameLength, nTCXDelay, &pTransientDetection->delayBuffer, &pTransientDetection->subblockEnergies);
    /* Init the TCX Short/Long transient detector. */
    InitTCXTransientDetector(nTCXDelay, &pTransientDetection->subblockEnergies, &pTransientDetection->transientDetector);
    /* We need two past subblocks for the TCX TD and NSUBBLOCKS+1 for the temporal flatness measure for the TCX LTP. */
    pTransientDetection->transientDetector.pSubblockEnergies->nDelay += NSUBBLOCKS+1;
}

float GetTCXAvgTemporalFlatnessMeasure(struct TransientDetection const * pTransientDetection, int nCurrentSubblocks, int nPrevSubblocks)
{
    int i;
    TransientDetector const * pTransientDetector = &pTransientDetection->transientDetector;
    SubblockEnergies const * pSubblockEnergies = pTransientDetector->pSubblockEnergies;
    int const nDelay = pTransientDetector->nDelay;                                                              /*  */
    int const nRelativeDelay = pSubblockEnergies->nDelay - nDelay;                                              /*   */
    float const * pSubblockNrgChange = NULL;
    float sumTempFlatness;
    int const nTotBlocks = nCurrentSubblocks+nPrevSubblocks;                                                    /*  */
    /* Initialization */
    assert(nTotBlocks > 0);
    sumTempFlatness = 0.0f;

    assert((nPrevSubblocks <= nRelativeDelay) && (nCurrentSubblocks <= NSUBBLOCKS+nDelay));
    pSubblockNrgChange = &pSubblockEnergies->subblockNrgChange[nRelativeDelay-nPrevSubblocks];
    for (i = 0; i < nTotBlocks; i++)
    {
        sumTempFlatness += pSubblockNrgChange[i];
    }
    return sumTempFlatness / (float)nTotBlocks;
}

float GetTCXMaxenergyChange(struct TransientDetection const * pTransientDetection,
                            const int isTCX10,
                            const int nCurrentSubblocks, const int nPrevSubblocks)
{
    int i;
    TransientDetector const * pTransientDetector = &pTransientDetection->transientDetector;
    SubblockEnergies const * pSubblockEnergies = pTransientDetector->pSubblockEnergies;
    int const nDelay = pTransientDetector->nDelay;                                                              /*  */
    int const nRelativeDelay = pSubblockEnergies->nDelay - nDelay;                                              /*   */
    float const * pSubblockNrgChange = NULL;
    float maxEnergyChange;
    int nTotBlocks = nCurrentSubblocks+nPrevSubblocks;                                                          /*  */
    /* Initialization */
    assert(nTotBlocks > 0);
    maxEnergyChange = 0.0f;

    assert((nPrevSubblocks <= nRelativeDelay) && (nCurrentSubblocks <= NSUBBLOCKS+nDelay));
    pSubblockNrgChange = &pSubblockEnergies->subblockNrgChange[nRelativeDelay-nPrevSubblocks];
    if (pTransientDetector->bIsAttackPresent || isTCX10)    /* frame is TCX-10 */
    {
        float const * pSubblockNrg = &pSubblockEnergies->subblockNrg[nRelativeDelay-nPrevSubblocks];
        float nrgMin, nrgMax = pSubblockNrg[0];
        int idxMax = 0;
        /* find subblock with maximum energy */
        for (i = 1; i < nTotBlocks; i++)
        {
            if (nrgMax < pSubblockNrg[i])
            {
                nrgMax = pSubblockNrg[i];
                idxMax = i;
            }
        }
        nrgMin = nrgMax;
        /* find minimum energy after maximum */
        for (i = idxMax + 1; i < nTotBlocks; i++)
        {
            if (nrgMin > pSubblockNrg[i])
            {
                nrgMin = pSubblockNrg[i];
            }
        }
        /* lower maxEnergyChange if energy doesn't decrease much after energy peak */
        if (nrgMin > 0.375f * nrgMax)
        {
            nTotBlocks = idxMax - 3;
        }
    }
    for (i = 0; i < nTotBlocks; i++)
    {
        maxEnergyChange = max(maxEnergyChange, pSubblockNrgChange[i]);
    }

    return maxEnergyChange;
}

/*---------------------------------------------------------------*
 * RunTransientDetection()
 *
 * Time Domain Transient Detector for TCX
 *---------------------------------------------------------------*/

void RunTransientDetection(
    float const * input,
    int nSamplesAvailable,
    TransientDetection * pTransientDetection
)
{
    float filteredInput[L_FRAME_MAX];
    SubblockEnergies * pSubblockEnergies = &pTransientDetection->subblockEnergies;                              /*  */
    TransientDetector * pTransientDetector = &pTransientDetection->transientDetector;                           /*  */

    assert((input != NULL) && (pTransientDetection != NULL) && (pSubblockEnergies != NULL) && (pTransientDetector != NULL));
    /* Variable initializations */
    HighPassFilter(input, nSamplesAvailable, &pSubblockEnergies->firState1, &pSubblockEnergies->firState2, filteredInput);

    /* Update subblock energies. */
    UpdateSubblockEnergies(filteredInput, nSamplesAvailable, pSubblockEnergies);

    /* Run transient detectors. */
    RunTransientDetector(pTransientDetector);

    /* Update the delay buffer. */
    UpdateDelayBuffer(filteredInput, nSamplesAvailable, &pTransientDetection->delayBuffer);

    return;
}


void SetTCXModeInfo(Encoder_State *st,                  /* i/o: encoder state structure                  */
                    TransientDetection const * pTransientDetection,
                    short * tcxModeOverlap
                   )
{
    if( st->codec_mode == MODE2 )
    {
        assert(pTransientDetection != NULL);

        /* determine window sequence (1 long or 2 short windows) */
        if (st->tcx10Enabled && st->tcx20Enabled)
        {
            /* window switching based on transient detector output */
            if (
                ((pTransientDetection->transientDetector.bIsAttackPresent)
                 || (st->currEnergyHF > st->prevEnergyHF * 39.0f)) &&
                ((st->last_core != ACELP_CORE) && (st->last_core != AMR_WB_CORE)))
            {
                st->tcxMode = TCX_10;
            }
            else
            {
                st->tcxMode = TCX_20;
            }
        }
        else
        {
            /* window selection (non-adaptive) based on flags only */
            if (st->tcx10Enabled)
            {
                st->tcxMode = TCX_10;
            }
            else if (st->tcx20Enabled)
            {
                st->tcxMode = TCX_20;
            }
            else
            {
                st->tcxMode = NO_TCX;
            }
        }

        /* set the left window overlap */

        if (st->last_core == ACELP_CORE || st->last_core == AMR_WB_CORE)
        {
            st->tcx_cfg.tcx_last_overlap_mode = TRANSITION_OVERLAP;
        }
        else if ((st->tcxMode == TCX_10) && (st->tcx_cfg.tcx_curr_overlap_mode == ALDO_WINDOW))
        {
            st->tcx_cfg.tcx_last_overlap_mode = FULL_OVERLAP;
        }
        else
        {
            st->tcx_cfg.tcx_last_overlap_mode = st->tcx_cfg.tcx_curr_overlap_mode;
        }

        /* determine the right window overlap */
        if( st->tcxMode == TCX_10 )
        {
            if (pTransientDetection->transientDetector.attackIndex < 0)
            {
                *tcxModeOverlap = HALF_OVERLAP;
            }
            else
            {
                *tcxModeOverlap = pTransientDetection->transientDetector.attackIndex % 4;
                if (*tcxModeOverlap == 1)
                {
                    *tcxModeOverlap = FULL_OVERLAP;
                }
            }
        }
        else if (st->tcxMode == TCX_20)
        {
            if (pTransientDetection->transientDetector.attackIndex == 7)
            {
                *tcxModeOverlap = HALF_OVERLAP;
            }
            else if (pTransientDetection->transientDetector.attackIndex == 6)
            {
                *tcxModeOverlap = MIN_OVERLAP;
            }
            else
            {
                *tcxModeOverlap = ALDO_WINDOW;
            }
        }
        else
        {
            /* NO_TCX */
            *tcxModeOverlap = TRANSITION_OVERLAP;
        }

        /* for the ACELP -> TCX transition frames use full right window overlap */

        if ((st->tcx_cfg.tcx_last_overlap_mode == TRANSITION_OVERLAP) && (*tcxModeOverlap == ALDO_WINDOW))
        {
            *tcxModeOverlap = FULL_OVERLAP;
        }
    }

    return;
}

/************************************************/
/*                                              */
/*              Internal functions              */
/*                                              */
/************************************************/

static void InitDelayBuffer(int nFrameLength, int nDelay, DelayBuffer * pDelayBuffer)
{
    int const nMaxBuffSize = sizeof(pDelayBuffer->buffer)/sizeof(pDelayBuffer->buffer[0]);

    assert((nFrameLength > NSUBBLOCKS) && (nFrameLength % NSUBBLOCKS == 0) && (nDelay >= 0) && (pDelayBuffer != NULL));
    pDelayBuffer->nSubblockSize = nFrameLength/NSUBBLOCKS;
    assert(pDelayBuffer->nSubblockSize <= nMaxBuffSize);
    set_f(pDelayBuffer->buffer, 0.0f, nMaxBuffSize);
    pDelayBuffer->nDelay = nDelay % pDelayBuffer->nSubblockSize;
    assert(pDelayBuffer->nDelay <= nMaxBuffSize);

    return;
}


static void InitSubblockEnergies(int nFrameLength, int nDelay, DelayBuffer * pDelayBuffer, SubblockEnergies * pSubblockEnergies)
{
    int const nMaxBuffSize = sizeof(pSubblockEnergies->subblockNrg)/sizeof(pSubblockEnergies->subblockNrg[0]);
    (void)nFrameLength;

    assert((pDelayBuffer != NULL) && (pSubblockEnergies != NULL) && (pDelayBuffer->nSubblockSize * NSUBBLOCKS == nFrameLength) && (pDelayBuffer->nSubblockSize > 0));

    set_f(pSubblockEnergies->subblockNrg,    MIN_BLOCK_ENERGY, nMaxBuffSize);
    set_f(pSubblockEnergies->accSubblockNrg, MIN_BLOCK_ENERGY, nMaxBuffSize+1);
    set_f(pSubblockEnergies->subblockNrgChange, 1.0f, nMaxBuffSize);
    pSubblockEnergies->nDelay = nDelay / pDelayBuffer->nSubblockSize;
    assert(pSubblockEnergies->nDelay < nMaxBuffSize);
    pSubblockEnergies->nPartialDelay = nDelay % pDelayBuffer->nSubblockSize;
    pSubblockEnergies->facAccSubblockNrg = 0.8125f; /* Energy accumulation factor */
    pSubblockEnergies->firState1 = 0.0f;
    pSubblockEnergies->firState2 = 0.0f;

    pSubblockEnergies->pDelayBuffer = pDelayBuffer;
    pDelayBuffer->nDelay = max(pDelayBuffer->nDelay, pSubblockEnergies->nPartialDelay);

    return;
}

/** Init transient detector.
  * Fills TransientDetector structure with sensible content and enable it.
  * @param pSubblockEnergies Subblock energies used in this transient detector.
  * @param nDelay Delay for this transient detector.
  * @param nSubblocksToCheck Number of subblocks to check in this transient detector.
  * @param pCheckSubblockForAttack Attack detection function for this transient detector.
  * @param pSetAttackPosition Function for finalizing this transient detector.
  * @param attackRatioThreshold Attack ratio threshold.
  * @param pTransientDetector Structure to be initialized.
  */
static void InitTransientDetector(SubblockEnergies * pSubblockEnergies, int nDelay, int nSubblocksToCheck,
                                  TCheckSubblocksForAttack pCheckSubblocksForAttack,
                                  float attackRatioThreshold, TransientDetector * pTransientDetector)
{
    int const nMaxBuffSize = sizeof(pSubblockEnergies->subblockNrg)/sizeof(pSubblockEnergies->subblockNrg[0]);
    (void)nMaxBuffSize;
    assert((pSubblockEnergies != NULL) && (pSubblockEnergies->pDelayBuffer != NULL) && (pTransientDetector != NULL) && (pSubblockEnergies->pDelayBuffer->nSubblockSize != 0));
    pTransientDetector->pSubblockEnergies = pSubblockEnergies;
    pTransientDetector->nDelay = (nDelay - pSubblockEnergies->nPartialDelay) / pSubblockEnergies->pDelayBuffer->nSubblockSize;
    assert(nDelay == pTransientDetector->nDelay * pSubblockEnergies->pDelayBuffer->nSubblockSize + pSubblockEnergies->nPartialDelay);
    assert(pTransientDetector->nDelay < nMaxBuffSize);
    pSubblockEnergies->nDelay = max(pSubblockEnergies->nDelay, pTransientDetector->nDelay);
    assert(nSubblocksToCheck <= NSUBBLOCKS + pTransientDetector->nDelay);
    pTransientDetector->nSubblocksToCheck = nSubblocksToCheck;
    pTransientDetector->CheckSubblocksForAttack = pCheckSubblocksForAttack;
    pTransientDetector->attackRatioThreshold = attackRatioThreshold;
    pTransientDetector->bIsAttackPresent = FALSE;
    pTransientDetector->attackIndex = -1;

    return;
}

/* This function should be inlined and WMOPS instrumentation takes that into account, meaning that all references are considered as local variables */
static float InlineFilter(float inValue, float firState1, float firState2)
{
    return 0.375f * inValue - 0.5f * firState1 + 0.125f * firState2;
}


static void HighPassFilter(float const * input, int length, float * pFirState1, float * pFirState2, float * output)
{
    int i;
    output[0] = InlineFilter(input[0], *pFirState1, *pFirState2);
    output[1] = InlineFilter(input[1], input[0], *pFirState1);
    for (i = 2; i < length; i++)
    {
        output[i] = InlineFilter(input[i], input[i-1], input[i-2]);
    }
    /* update filter states: shift time samples through delay line */
    *pFirState2 = input[length-2];
    *pFirState1 = input[length-1];

    return;
}


static void RunTransientDetector(TransientDetector * pTransientDetector)
{
    float const attackRatioThreshold = pTransientDetector->attackRatioThreshold;                                /*  */
    SubblockEnergies const * pSubblockEnergies = pTransientDetector->pSubblockEnergies;                         /*  */
    int const nDelay = pTransientDetector->nDelay;                                                              /*  */
    int const nRelativeDelay = pSubblockEnergies->nDelay - nDelay;                                              /*   */
    float const * pSubblockNrg = &pSubblockEnergies->subblockNrg[nRelativeDelay];                               /*  */
    float const * pAccSubblockNrg = &pSubblockEnergies->accSubblockNrg[nRelativeDelay];                         /*  */

    assert((pTransientDetector->CheckSubblocksForAttack != NULL));
    /* Variable initialization */
    pTransientDetector->CheckSubblocksForAttack(pSubblockNrg, pAccSubblockNrg,
            NSUBBLOCKS+nDelay, nRelativeDelay,
            attackRatioThreshold,
            &pTransientDetector->bIsAttackPresent, &pTransientDetector->attackIndex);

    return;
}

static void UpdateDelayBuffer(float const * input, int nSamplesAvailable, DelayBuffer * pDelayBuffer)
{
    int i;
    int const nDelay = pDelayBuffer->nDelay;
    assert((nDelay >= 0) && (nDelay <= (int)sizeof(pDelayBuffer->buffer)/(int)sizeof(pDelayBuffer->buffer[0])));
    assert(nSamplesAvailable <= NSUBBLOCKS*pDelayBuffer->nSubblockSize);
    /* If this is not the last frame */
    if (nSamplesAvailable == NSUBBLOCKS*pDelayBuffer->nSubblockSize)
    {
        /* Store the newest samples into the delay buffer */
        for (i = 0; i < nDelay; i++)
        {
            pDelayBuffer->buffer[i] = input[i+nSamplesAvailable-nDelay];
        }
    }

    return;
}


static void UpdateSubblockEnergies(float const * input, int nSamplesAvailable, SubblockEnergies * pSubblockEnergies)
{
    int i;
    assert((pSubblockEnergies->nDelay >= 0) && (pSubblockEnergies->nDelay+NSUBBLOCKS <= (int)sizeof(pSubblockEnergies->subblockNrg)/(int)sizeof(pSubblockEnergies->subblockNrg[0])));
    assert(pSubblockEnergies->nPartialDelay <= pSubblockEnergies->pDelayBuffer->nDelay);
    /* At least one block delay is required when subblock energy change is required */
    assert(pSubblockEnergies->nDelay >= 1);

    /* Shift old subblock energies */
    for (i = 0; i < pSubblockEnergies->nDelay; i++)
    {
        pSubblockEnergies->subblockNrg[i] = pSubblockEnergies->subblockNrg[i+NSUBBLOCKS];
        pSubblockEnergies->accSubblockNrg[i] = pSubblockEnergies->accSubblockNrg[i+NSUBBLOCKS];
        pSubblockEnergies->subblockNrgChange[i] = pSubblockEnergies->subblockNrgChange[i+NSUBBLOCKS];
    }

    /* Compute filtered subblock energies for the new samples */
    CalculateSubblockEnergies(input, nSamplesAvailable, pSubblockEnergies);

    return;
}


/* This function should be inlined and WMOPS instrumentation takes that into account, meaning that all references are considered as local variables */
static void UpdatedAndStoreAccWindowNrg(float newWindowNrgF, float * pAccSubblockNrg, float facAccSubblockNrg, float * pOutAccWindowNrgF)
{
    /* Store the accumulated energy */
    *pOutAccWindowNrgF = *pAccSubblockNrg;
    /* Update the accumulated energy: maximum of the current and the accumulated energy */
    *pAccSubblockNrg *= facAccSubblockNrg;
    if (newWindowNrgF > *pAccSubblockNrg)
    {
        *pAccSubblockNrg = newWindowNrgF;
    }

    return;
}


static void CalculateSubblockEnergies(float const * input, int nSamplesAvailable, SubblockEnergies * pSubblockEnergies)
{
    DelayBuffer * pDelayBuffer = pSubblockEnergies->pDelayBuffer;                                               /*  */
    int const nSubblockSize = pDelayBuffer->nSubblockSize;                                                      /*  */
    int const nDelay = pSubblockEnergies->nDelay;                                                               /*  */
    int const nPartialDelay = pSubblockEnergies->nPartialDelay;                                                 /*  */
    float const * delayBuffer = &pDelayBuffer->buffer[pDelayBuffer->nDelay - nPartialDelay];                    /*   */
    float const facAccSubblockNrg = pSubblockEnergies->facAccSubblockNrg;                                       /*  */
    float * pSubblockNrg = &pSubblockEnergies->subblockNrg[nDelay];                                             /*  */
    float * pAccSubblockNrg = &pSubblockEnergies->accSubblockNrg[nDelay];                                       /*  */
    float * pSubblockNrgChange = &pSubblockEnergies->subblockNrgChange[nDelay];                                 /*  */
    float * pAccSubblockTmp;
    int nWindows;
    int i, w, k;
    /* Variable initializations */
    nWindows = (nSamplesAvailable+nPartialDelay)/nSubblockSize;
    pAccSubblockTmp = &pAccSubblockNrg[nWindows];

    set_f(pSubblockNrg, MIN_BLOCK_ENERGY, NSUBBLOCKS);
    if (nWindows > 0)
    {
        /* Process left over samples from the previous frame. */
        for (k = 0; k < nPartialDelay; k++)
        {
            pSubblockNrg[0] += delayBuffer[k] * delayBuffer[k];
        }
        k = 0;
        /* Process new samples in the 0. subblock. */
        for (i = nPartialDelay; i < nSubblockSize; i++, k++)
        {
            pSubblockNrg[0] += input[k] * input[k];
        }
        /* Set accumulated subblock energy at this point. */
        UpdatedAndStoreAccWindowNrg(pSubblockNrg[0], pAccSubblockTmp, facAccSubblockNrg, &pAccSubblockNrg[0]);
        for (w = 1; w < nWindows; w++)
        {
            /* Process new samples in the w. subblock. */
            for (i = 0; i < nSubblockSize; i++, k++)
            {
                pSubblockNrg[w] += input[k] * input[k];
            }
            /* Set accumulated subblock energy at this point. */
            UpdatedAndStoreAccWindowNrg(pSubblockNrg[w], pAccSubblockTmp, facAccSubblockNrg, &pAccSubblockNrg[w]);
        }
        /* Calculate energy change for each block. */
        for (w = 0; w < nWindows; w++)
        {
            if (pSubblockNrg[w] > pSubblockNrg[w-1])
            {
                pSubblockNrgChange[w] = pSubblockNrg[w]/pSubblockNrg[w-1];
            }
            else
            {
                pSubblockNrgChange[w] = pSubblockNrg[w-1]/pSubblockNrg[w];
            }
        }
    }

    return;
}
