/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

#include <memory.h>
#include <math.h>
#include <assert.h>


#define HLM_MIN_NRG (32768.0f * 2*NORM_MDCT_FACTOR / (640*640))

static float Autocorrelation(float const x[], int n, int lag)
{
    return dotp(x, x+lag, n-lag);
}

/** Get TNS filter parameters from autocorrelation.
  *
  * @param rxx Autocorrelation function/coefficients.
  * @param maxOrder Maximum filter order/number of coefficients.
  * @param pTnsFilter Pointer to the output filter.
  */
static void GetFilterParameters(float const rxx[], int maxOrder, STnsFilter * pTnsFilter);

/** Quantization for reflection coefficients.
  *
  * @param parCoeff input reflection coefficients.
  * @param index output quantized values.
  * @param order number of coefficients/values.
  */
static void Parcor2Index(float const parCoeff[], int index[], int order);

/** Linear prediction analysis/synthesis filter definition.
  * @param order filter order.
  * @param parCoeff filter (PARCOR) coefficients.
  * @param state state of the filter. Must be at least of 'order' size.
  * @param x the current input value.
  * @return the output of the filter.
  */
typedef float (* TLinearPredictionFilter)(int order, float const parCoeff[], float * state, float x);

/********************************/
/*      Interface functions     */
/********************************/

#define MAX_SUBDIVISIONS 3

int DetectTnsFilt(STnsConfig const * pTnsConfig,
                  float const pSpectrum[],
                  STnsData * pTnsData
                  , float* predictionGain
                 )
{
    float norms[TNS_MAX_NUM_OF_FILTERS][MAX_SUBDIVISIONS];
    int iFilter = 0;
    short i;

    for( i=0; i<TNS_MAX_NUM_OF_FILTERS; i++ )
    {
        set_f( norms[i], 0, MAX_SUBDIVISIONS );
    }
    /* WMOPS: All initializations are either for safety or static (tables) and thus not to be counted */
    ResetTnsData(pTnsData);
    if (pTnsConfig->maxOrder <= 0)
    {
        return 0;
    }

    /* Calculate norms for each spectrum part */
    for (iFilter = 0; iFilter < pTnsConfig->nMaxFilters; iFilter++)
    {
        int const idx0 = pTnsConfig->iFilterBorders[iFilter+1];
        int const idx1 = pTnsConfig->iFilterBorders[iFilter];
        int const nSubdivisions = pTnsConfig->pTnsParameters[iFilter].nSubdivisions;
        int iSubdivisions;
        /* Variable initialization */
        assert(pTnsConfig->pTnsParameters[iFilter].nSubdivisions <= MAX_SUBDIVISIONS);
        for (iSubdivisions = 0; iSubdivisions < nSubdivisions; iSubdivisions++)
        {
            int const iStartLine = idx0 + (idx1-idx0)*iSubdivisions/nSubdivisions;
            int const iEndLine = idx0 + (idx1-idx0)*(iSubdivisions+1)/nSubdivisions;
            /* Variable initialization */
            norms[iFilter][iSubdivisions] = sum2_f(pSpectrum+iStartLine, iEndLine-iStartLine);                    /* No need for indirect, PTR_INIT used instead */

        }
    }
    /* Calculate normalized autocorrelation for spectrum subdivision and get TNS filter parameters based on it */
    for (iFilter = 0; iFilter < pTnsConfig->nMaxFilters; iFilter++)
    {
        float rxx[TNS_MAX_FILTER_ORDER+1];
        int const idx0 = pTnsConfig->iFilterBorders[iFilter+1];
        int const idx1 = pTnsConfig->iFilterBorders[iFilter];
        int const spectrumLength = idx1 - idx0;
        STnsFilter * const pFilter = pTnsData->filter + iFilter;
        int const nSubdivisions = pTnsConfig->pTnsParameters[iFilter].nSubdivisions;
        int iSubdivisions;

        set_f( rxx, 0, TNS_MAX_FILTER_ORDER+1 );            /* WMOPS: This initialization is required */
        /* Variable initialization */
        for (iSubdivisions = 0; (iSubdivisions < nSubdivisions) && (norms[iFilter][iSubdivisions] > HLM_MIN_NRG); iSubdivisions++)
        {
            float const fac = 1.0f/norms[iFilter][iSubdivisions];
            int const iStartLine = idx0 + spectrumLength*iSubdivisions/nSubdivisions;
            int const iEndLine = idx0 + spectrumLength*(iSubdivisions+1)/nSubdivisions;
            float const * pWindow = tnsAcfWindow;
            float const * pEndWindowCoeff = &tnsAcfWindow[sizeof(tnsAcfWindow)/sizeof(tnsAcfWindow[0])]; /* WMOPS: Not counted as it is used only in the assertion */
            int lag;
            (void)pEndWindowCoeff;
            /* For additional loop condition */
            /* Variable initialization */
            for (lag = 1; lag <= pTnsConfig->maxOrder; lag++)
            {
                rxx[lag] += fac * (*pWindow) * Autocorrelation(pSpectrum+iStartLine, iEndLine-iStartLine, lag);
                pWindow++;
            }
        }
        if (iSubdivisions == nSubdivisions) /* meaning there is no subdivision with low energy */
        {
            rxx[0] = (float)pTnsConfig->pTnsParameters[iFilter].nSubdivisions;
            pFilter->spectrumLength = spectrumLength;
            /* Limit the maximum order to spectrum length/4 */
            GetFilterParameters(rxx, min(pTnsConfig->maxOrder, pFilter->spectrumLength/4), pFilter);
        }
    }

    if (predictionGain)
    {
        assert(pTnsConfig->nMaxFilters == 1);
        *predictionGain = pTnsData->filter->predictionGain;
    }

    /* We check the filter's decisions in the opposite direction */
    for (iFilter = pTnsConfig->nMaxFilters-1; iFilter >= 0; iFilter--)
    {
        STnsFilter * const pFilter = pTnsData->filter + iFilter;
        struct TnsParameters const * const pTnsParameters = pTnsConfig->pTnsParameters + iFilter;

        /* TNS decision function */
        if ((pFilter->predictionGain > pTnsParameters->minPredictionGain)
                || (pFilter->avgSqrCoef > pTnsParameters->minAvgSqrCoef))
        {
            ++pTnsData->nFilters;
        }
        else if (pTnsData->nFilters > 0) /* If a previous filter is turned on */
        {
            /* Since TNS filter of order 0 is not allowed we haved to signal in the stream filter of order 1 with the 0th coefficient equal to 0 */
            ClearTnsFilterCoefficients(pFilter);
            pFilter->order = 1;
            ++pTnsData->nFilters;
        }
        else
        {
            ClearTnsFilterCoefficients(pFilter);
        }
    }
    return (pTnsData->nFilters > 0) ? 1 : 0;
}

TNS_ERROR EncodeTnsData(STnsConfig const * pTnsConfig, STnsData const * pTnsData, int * stream, int * pnSize, int * pnBits)
{
    *pnSize = 0;
    *pnBits = 0;
    if (pTnsConfig->nMaxFilters > 1)
    {
        if (pTnsConfig->iFilterBorders[0] < 512)
        {
            GetParameters(&tnsEnabledSWBTCX10BitMap, 1, pTnsData, &stream, pnSize, pnBits);
        }
        else
        {
            GetParameters(&tnsEnabledSWBTCX20BitMap, 1, pTnsData, &stream, pnSize, pnBits);
        }
    }
    else
    {
        GetParameters(&tnsEnabledWBTCX20BitMap, 1, pTnsData, &stream, pnSize, pnBits);
    }
    return TNS_NO_ERROR;
}

TNS_ERROR WriteTnsData(STnsConfig const * pTnsConfig, int const * stream, int * pnSize, Encoder_State *st, int * pnBits)
{
    if (pTnsConfig->nMaxFilters > 1)
    {
        if (pTnsConfig->iFilterBorders[0] < 512)
        {
            WriteToBitstream(&tnsEnabledSWBTCX10BitMap, 1, &stream, pnSize, st, pnBits);
        }
        else
        {
            WriteToBitstream(&tnsEnabledSWBTCX20BitMap, 1, &stream, pnSize, st, pnBits);
        }
    }
    else
    {
        WriteToBitstream(&tnsEnabledWBTCX20BitMap, 1, &stream, pnSize, st, pnBits);
    }
    return TNS_NO_ERROR;
}

/*********************************************************************************************/
/*  Definitions of functions used in the mapping between TNS parameters and a bitstream.     */
/*********************************************************************************************/

/* Helper functions for hufmann table coding */

/********************************/
/*      Private functions       */
/********************************/

/** Autocorrelation to parcor coefficients.
  * Conversion of autocorrelation to parcor/reflection coefficients.
  * @param input Autocorrelation function/coefficients.
  * @param parCoeff output filter (PARCOR) coefficients.
  * @param order filter order.
  * @return prediction gain.
  */
static float AutoToParcor(const float input[], float parCoeff[], int order)
{
    int i, j;
    float tmp, tmp2;
    float workBuffer[2*TNS_MAX_FILTER_ORDER];
    float * const pWorkBuffer = &workBuffer[order]; /* temp pointer */                                          /* WMOPS: No need for counting as it can be realized as a separate buffer */

    for(i=0; i<order; i++)
    {
        workBuffer[i] = input[i];
        pWorkBuffer[i] = input[i+1];
    }
    for(i=0; i<order; i++)
    {
        if (workBuffer[0] < 1.0f/65536.0f)
        {
            tmp = 0;
        }
        else
        {
            tmp = -pWorkBuffer[i]/workBuffer[0];
        }
        /*
          compensate for calculation inaccuracies
          limit reflection coefs to ]-1,1[
        */
        tmp = min(0.999f, max(-0.999f, tmp ));

        parCoeff[i] = tmp;
        for(j=i; j<order; j++)
        {
            tmp2 = pWorkBuffer[j] + tmp * workBuffer[j-i];
            workBuffer[j-i] += tmp * pWorkBuffer[j];
            pWorkBuffer[j] = tmp2;
        }
    }
    return ( (input[0] + 1e-30f) / (workBuffer[0] + 1e-30f) );
}

static void GetFilterParameters(const float rxx[], int maxOrder, STnsFilter * pTnsFilter)
{
    int i;
    float parCoeff[TNS_MAX_FILTER_ORDER];
    float const * values = tnsCoeff4;
    int * indexes = pTnsFilter->coefIndex;
    /* Variable initialization */
    /* compute TNS filter in lattice (ParCor) form with LeRoux-Gueguen algorithm */
    pTnsFilter->predictionGain = AutoToParcor(rxx, parCoeff, maxOrder);
    /* non-linear quantization of TNS lattice coefficients with given resolution */
    Parcor2Index(parCoeff, indexes, maxOrder);

    /* reduce filter order by truncating trailing zeros */
    i = maxOrder - 1;
    while ((i >= 0) && (indexes[i] == 0))
    {
        --i;
        /* Loop condition */
    }
    pTnsFilter->order = i + 1;

    /* compute avg(coef*coef) */
    pTnsFilter->avgSqrCoef = 0;
    for (i = pTnsFilter->order-1; i >= 0; i--)
    {
        float const value = values[indexes[i]+INDEX_SHIFT];
        /* Variable initialization */
        pTnsFilter->avgSqrCoef += value * value;

    }
    pTnsFilter->avgSqrCoef /= maxOrder;
}

static void Parcor2Index(const float parCoeff[], int index[], int order)
{
    int const nValues = 1 << TNS_COEF_RES;                                                                      /* WMOPS: This is just a constant */
    float const * values = tnsCoeff4;
    int i;
    int iIndex;
    float x;
    for (i = 0; i < order; i++)
    {
        iIndex = 1;
        x = parCoeff[i];
        /* Variable initialization */
        assert((x >= -1.0f) && (x <= 1.0f));
        while ((iIndex < nValues) && (x > 0.5f*(values[iIndex-1]+values[iIndex])))
        {
            ++iIndex;
            /* Loop condition */
        }
        index[i] = (iIndex-1) - INDEX_SHIFT;
    }
}

