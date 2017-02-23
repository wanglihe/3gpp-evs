/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <memory.h>
#include <math.h>
#include <assert.h>
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"




/** Inverse quantization for reflection coefficients.
  *
  * @param index input quantized values.
  * @param parCoeff output reflection coefficients.
  * @param order number of coefficients/values.
  */
static void Index2Parcor(int const index[], float parCoeff[], int order);

/** Linear prediction analysis/synthesis filter definition.
  * @param order filter order.
  * @param parCoeff filter (PARCOR) coefficients.
  * @param state state of the filter. Must be at least of 'order' size.
  * @param x the current input value.
  * @return the output of the filter.
  */
typedef float (* TLinearPredictionFilter)(int order, float const parCoeff[], float * state, float x);

/** Linear prediction analysis filter.
  * See TLinearPredictionFilter for details.
  */
static float FIRLattice(int order, const float *parCoeff, float *state, float x);

/** Linear prediction synthesis filter.
  * See TLinearPredictionFilter for details.
  */
static float IIRLattice(int order, const float *parCoeff, float *state, float x);

/** TNS analysis/synthesis filter.
  * @param spectrum input spectrum values.
  * @param numOfLines number of lines in the spectrum.
  * @param parCoeff filter (PARCOR) coefficients.
  * @param order filter order.
  * @param filter function that implements filtering.
    By this function it is defined whether analysis or synthesis is performed.
  * @param output filtered output spectrum values.
    Inplace operation is supported, so it can be equal to spectrum.
  */
static void TnsFilter(float const spectrum[], int numOfLines,
                      float const parCoeff[], int order,
                      TLinearPredictionFilter filter, float * state,
                      float output[]);


/********************************/
/*      Interface functions     */
/********************************/

TNS_ERROR InitTnsConfiguration(
    int nSampleRate,
    int frameLength,
    STnsConfig * pTnsConfig,
    int igfStopFreq,
    int bitrate
)
{
    int iFilter = 0;
    short int * startLineFilter = &pTnsConfig->iFilterBorders[1];

    /* Sanity checks */
    assert((nSampleRate > 0) && (frameLength > 0) && (pTnsConfig != NULL));
    if ((nSampleRate <= 0) || (frameLength <= 0) || (pTnsConfig == NULL))
        return TNS_FATAL_ERROR;
    /* WMOPS: All initializations are either for safety or static (tables) and thus not to be counted */
    /* Initialize TNS filter flag and maximum order */
    pTnsConfig->maxOrder    = TNS_MAX_FILTER_ORDER;
    if (bitrate <= ACELP_32k)
    {
        pTnsConfig->nMaxFilters = sizeof(tnsParametersIGF32kHz_LowBR)/sizeof(tnsParametersIGF32kHz_LowBR[0]);

        pTnsConfig->pTnsParameters = tnsParametersIGF32kHz_LowBR;
    }
    else
    {
        if (nSampleRate > 32000 && nSampleRate == 100 * frameLength)
        {
            pTnsConfig->nMaxFilters = sizeof(tnsParameters48kHz_grouped)/sizeof(tnsParameters48kHz_grouped[0]);
            pTnsConfig->pTnsParameters = tnsParameters48kHz_grouped;
        }
        else if (nSampleRate > 16000)
        {
            pTnsConfig->nMaxFilters = sizeof(tnsParameters32kHz)/sizeof(tnsParameters32kHz[0]);
            if (nSampleRate == 100 * frameLength)    /* sub-frame length is <= 10 ms */
            {
                pTnsConfig->pTnsParameters = tnsParameters32kHz_grouped;
            }
            else
            {
                pTnsConfig->pTnsParameters = tnsParameters32kHz;
            }
        }
        else
        {
            if (nSampleRate == 100 * frameLength)    /* sub-frame length is <= 10 ms */
            {
                pTnsConfig->nMaxFilters = sizeof(tnsParameters16kHz_grouped)/sizeof(tnsParameters16kHz_grouped[0]);
                pTnsConfig->pTnsParameters = tnsParameters16kHz_grouped;
            }
            else
            {
                pTnsConfig->nMaxFilters = sizeof(tnsParameters16kHz)/sizeof(tnsParameters16kHz[0]);
                pTnsConfig->pTnsParameters = tnsParameters16kHz;
            }
        }
    }

    assert(pTnsConfig->nMaxFilters <= TNS_MAX_NUM_OF_FILTERS);

    /* Set starting MDCT line for each filter based on the starting frequencies from the TNS table */
    for (iFilter = 0; iFilter < pTnsConfig->nMaxFilters; iFilter++)
    {
        assert(pTnsConfig->pTnsParameters[iFilter].startLineFrequency < 0.5f*nSampleRate);
        startLineFilter[iFilter] = (frameLength * 2 * pTnsConfig->pTnsParameters[iFilter].startLineFrequency) / nSampleRate;
    }
    if (igfStopFreq > 0)
    {
        pTnsConfig->iFilterBorders[0] = (frameLength * 2 * igfStopFreq) / nSampleRate;
    }
    else
    {
        pTnsConfig->iFilterBorders[0] = frameLength;
    }

    return TNS_NO_ERROR;
}

TNS_ERROR ApplyTnsFilter(
    STnsConfig const * pTnsConfig,
    STnsData const * pTnsData,
    float spectrum[],
    int fIsAnalysis)
{
    TLinearPredictionFilter filter;
    float state[TNS_MAX_FILTER_ORDER];
    int iFilter;
    int stopLine, startLine;
    short int const * pBorders;
    filter = fIsAnalysis ? FIRLattice : IIRLattice;
    set_f(state, 0, TNS_MAX_FILTER_ORDER);
    pBorders = pTnsConfig->iFilterBorders;
    for (iFilter = pTnsConfig->nMaxFilters-1; iFilter >= 0; iFilter--)
    {
        float parCoeff[TNS_MAX_FILTER_ORDER];
        STnsFilter const * const pFilter = &pTnsData->filter[iFilter];

        set_f( parCoeff, 0, TNS_MAX_FILTER_ORDER );                                                              /* WMOPS: Initialization for safety and thus not to be counted */
        stopLine = pBorders[iFilter];
        startLine = pBorders[iFilter+1];

        Index2Parcor(pFilter->coefIndex, parCoeff, pFilter->order);

        TnsFilter(&spectrum[startLine], stopLine-startLine,
                  parCoeff, pFilter->order,
                  filter, state,
                  &spectrum[startLine]);

    }
    return (pTnsData->nFilters < 0) ? TNS_FATAL_ERROR : TNS_NO_ERROR;
}

/*********************************************************************************************/
/*  Definitions of functions used in the mapping between TNS parameters and a bitstream.     */
/*********************************************************************************************/

/* Helper functions for hufmann table coding */

/** Get number of bits from a Huffman table.
  * The table must be sorted by values.
  */
static int GetBitsFromTable(int value, const Coding codes[], int nSize)
{
    (void)nSize;
    assert((value >= 0) && (value < nSize) && (nSize >= 0) && (nSize <= 256));
    return codes[value].nBits;
}

/** Get the code for a value from a Huffman table.
  * The table must be sorted by values.
  */
static int EncodeUsingTable(int value, const Coding codes[], int nSize)
{
    (void)nSize;
    assert((value >= 0) && (value < nSize) && (nSize >= 0) && (nSize <= 256));
    return codes[value].code;
}

/** Decode a value from a bitstream using a Huffman table. */
static int DecodeUsingTable(Decoder_State *st, int * pValue, const Coding codes[], int nSize)
{
    unsigned short int code = 0;
    unsigned char nBits = 0;
    unsigned char valueIndex = nSize;

    assert((nSize >= 0) && (nSize <= 256));
    /* Variable initialization. All are required! */
    while (valueIndex == nSize)
    {
        code = (code << 1) + (unsigned short int)get_next_indice(st, 1);
        ++nBits;
        if (nBits > nSize || nBits > 16)
        {
            st->BER_detect = 1;
            *pValue = 0;
            return -1;
        }
        for (valueIndex = 0; valueIndex < nSize; valueIndex++)
        {
            if (codes[valueIndex].nBits == nBits)
            {
                if (codes[valueIndex].code == code)
                    break;
            }
        }
        /* Loop condition */
    }

    if (valueIndex < nSize)
    {
        *pValue = codes[valueIndex].value;
    }
    else
    {
        st->BER_detect = 1;
        *pValue = 0;
        return -1;
    }

    return nBits;
}

/* TNS filter coefficients */

void const * GetTnsFilterCoeff(void const * p, int index, int * pValue)
{
    *pValue = ((int const *)p)[index] + INDEX_SHIFT;
    return NULL;
}

void * SetTnsFilterCoeff(void * p, int index, int value)
{
    ((int *)p)[index] = value - INDEX_SHIFT;
    return NULL;
}

int GetSWBTCX20TnsFilterCoeffBits(int value, int index)
{
    assert((index >= 0) && (index < nTnsCoeffTables));
    return GetBitsFromTable(value, codesTnsCoeffSWBTCX20[index], nTnsCoeffCodes);
}

int EncodeSWBTCX20TnsFilterCoeff(int value, int index)
{
    assert((index >= 0) && (index < nTnsCoeffTables));
    return EncodeUsingTable(value, codesTnsCoeffSWBTCX20[index], nTnsCoeffCodes);
}

int DecodeSWBTCX20TnsFilterCoeff(Decoder_State *st, int index, int * pValue)
{
    assert((index >= 0) && (index < nTnsCoeffTables));
    return DecodeUsingTable(st, pValue, codesTnsCoeffSWBTCX20[index], nTnsCoeffCodes);
}

int GetSWBTCX10TnsFilterCoeffBits(int value, int index)
{
    assert((index >= 0) && (index < nTnsCoeffTables));
    return GetBitsFromTable(value, codesTnsCoeffSWBTCX10[index], nTnsCoeffCodes);
}

int EncodeSWBTCX10TnsFilterCoeff(int value, int index)
{
    assert((index >= 0) && (index < nTnsCoeffTables));
    return EncodeUsingTable(value, codesTnsCoeffSWBTCX10[index], nTnsCoeffCodes);
}

int DecodeSWBTCX10TnsFilterCoeff(Decoder_State *st, int index, int * pValue)
{
    assert((index >= 0) && (index < nTnsCoeffTables));
    return DecodeUsingTable(st, pValue, codesTnsCoeffSWBTCX10[index], nTnsCoeffCodes);
}

int GetWBTCX20TnsFilterCoeffBits(int value, int index)
{
    assert((index >= 0) && (index < nTnsCoeffTables));
    return GetBitsFromTable(value, codesTnsCoeffWBTCX20[index], nTnsCoeffCodes);
}

int EncodeWBTCX20TnsFilterCoeff(int value, int index)
{
    assert((index >= 0) && (index < nTnsCoeffTables));
    return EncodeUsingTable(value, codesTnsCoeffWBTCX20[index], nTnsCoeffCodes);
}

int DecodeWBTCX20TnsFilterCoeff(Decoder_State *st, int index, int * pValue)
{
    assert((index >= 0) && (index < nTnsCoeffTables));
    return DecodeUsingTable(st, pValue, codesTnsCoeffWBTCX20[index], nTnsCoeffCodes);
}


/* TNS filter order */

void const * GetTnsFilterOrder(void const * p, int index, int * pValue)
{
    *pValue = ((STnsFilter const *)p)[index].order;
    return ((STnsFilter const *)p)[index].coefIndex;
}

void * SetTnsFilterOrder(void * p, int index, int value)
{
    ((STnsFilter *)p)[index].order = value;
    return ((STnsFilter *)p)[index].coefIndex;
}

int GetTnsFilterOrderBitsSWBTCX20(int value, int index)
{
    (void)index;
    return GetBitsFromTable(value-1, codesTnsOrderTCX20, nTnsOrderCodes);
}

int EncodeTnsFilterOrderSWBTCX20(int value, int index)
{
    (void)index;
    return EncodeUsingTable(value-1, codesTnsOrderTCX20, nTnsOrderCodes);
}

int DecodeTnsFilterOrderSWBTCX20(Decoder_State *st, int index, int * pValue)
{
    (void)index;
    return DecodeUsingTable(st, pValue, codesTnsOrderTCX20, nTnsOrderCodes);
}

int GetTnsFilterOrderBitsSWBTCX10(int value, int index)
{
    (void)index;
    return GetBitsFromTable(value-1, codesTnsOrderTCX10, nTnsOrderCodes);
}

int EncodeTnsFilterOrderSWBTCX10(int value, int index)
{
    (void)index;
    return EncodeUsingTable(value-1, codesTnsOrderTCX10, nTnsOrderCodes);
}

int DecodeTnsFilterOrderSWBTCX10(Decoder_State *st, int index, int * pValue)
{
    (void)index;
    return DecodeUsingTable(st, pValue, codesTnsOrderTCX10, nTnsOrderCodes);
}

int GetTnsFilterOrderBits(int value, int index)
{
    (void)index;
    return GetBitsFromTable(value-1, codesTnsOrder, nTnsOrderCodes);
}

int EncodeTnsFilterOrder(int value, int index)
{
    (void)index;
    return EncodeUsingTable(value-1, codesTnsOrder, nTnsOrderCodes);
}

int DecodeTnsFilterOrder(Decoder_State *st, int index, int * pValue)
{
    (void)index;
    return DecodeUsingTable(st, pValue, codesTnsOrder, nTnsOrderCodes);
}

/* Number of TNS filters */

void const * GetNumOfTnsFilters(void const * p, int index, int * pValue)
{
    *pValue = ((STnsData const *)p)[index].nFilters;
    return ((STnsData const *)p)[index].filter;
}

void * SetNumOfTnsFilters(void * p, int index, int value)
{
    ((STnsData *)p)[index].nFilters = value;
    return ((STnsData *)p)[index].filter;
}

/* TNS enabled/disabled flag */

void const * GetTnsEnabled(void const * p, int index, int * pValue)
{
    *pValue = ((STnsData const *)p)[index].nFilters > 0 ? 1 : 0;
    return NULL;
}

void * SetTnsEnabled(void * p, int index, int value)
{
    (void)p,(void)index,(void)value;
    return NULL;
}

void const * GetTnsEnabledSingleFilter(void const * p, int index, int * pValue)
{
    *pValue = ((STnsData const *)p)[index].nFilters > 0 ? 1 : 0;
    return ((STnsData const *)p)[index].filter;
}

void * SetTnsEnabledSingleFilter(void * p, int index, int value)
{
    ((STnsData *)p)[index].nFilters = value;
    return ((STnsData *)p)[index].filter;
}

/********************************/
/*      Private functions       */
/********************************/

void ResetTnsData(STnsData * pTnsData)
{
    unsigned int iFilter;

    pTnsData->nFilters = 0;
    for (iFilter = 0; iFilter < sizeof(pTnsData->filter)/sizeof(pTnsData->filter[0]); iFilter++)
    {
        STnsFilter * const pTnsFilter = &pTnsData->filter[iFilter];
        pTnsFilter->spectrumLength = 0;
        pTnsFilter->predictionGain  = 1.0f;
        pTnsFilter->avgSqrCoef = 0;
        ClearTnsFilterCoefficients(pTnsFilter);
    }

}

void ClearTnsFilterCoefficients(STnsFilter * pTnsFilter)
{
    pTnsFilter->order = 0;
    set_i(pTnsFilter->coefIndex, 0, (int)sizeof(pTnsFilter->coefIndex)/(int)sizeof(pTnsFilter->coefIndex[0]));
}

static void Index2Parcor(const int index[], float parCoeff[], int order)
{
    float const * values = tnsCoeff4;
    int i;
    for (i = 0; i < order; i++)
    {
        parCoeff[i] = values[index[i] + INDEX_SHIFT];
    }

}

static float FIRLattice(int order, const float *parCoeff, float *state, float x)
{
    int i;
    float tmpSave;
    tmpSave = x;
    for (i = 0; i < order-1; i++)
    {
        float const tmp = parCoeff[i] * x + state[i];
        /* Variable initialization */
        x += parCoeff[i] * state[i];
        state[i] = tmpSave;
        tmpSave = tmp;
    }

    /* last stage: only need half operations */
    x += parCoeff[order-1] * state[order-1];
    state[order-1] = tmpSave;
    return x;
}

static float IIRLattice(int order, const float *parCoeff, float *state, float x)
{
    int i;
    /* first stage: no need to calculate state[order-1] */
    x -= parCoeff[order-1] * state[order-1];
    for (i = order-2; i >= 0; i--)
    {
        x -= parCoeff[i] * state[i];
        state[i+1] = parCoeff[i] * x + state[i];
    }

    state[0] = x;
    return x;
}

static void TnsFilter(float const spectrum[], int numOfLines,
                      float const parCoeff[], int order,
                      TLinearPredictionFilter filter, float * state,
                      float output[])
{
    int j;

    assert((order >= 0) && (order <= TNS_MAX_FILTER_ORDER));
    assert((numOfLines > 0) || ((numOfLines == 0) && (order == 0)));
    if (order == 0)
    {
        if ((spectrum != output) && (numOfLines > 0))
        {
            mvr2r(spectrum, output, numOfLines);
        }
    }
    else
    {
        {
            for (j = 0; j < numOfLines; j++)
            {
                output[j] = filter(order, parCoeff, state, spectrum[j]);
            }
        }
    }
}


static void ITF_TnsFilter(float const spectrum[], int numOfLines,
                          const float A[], int order,
                          float output[])
{
    int j;
    assert((order >= 0) && (order <= ITF_MAX_FILTER_ORDER));
    assert((numOfLines > 0) || ((numOfLines == 0) && (order == 0)));
    if (order == 0)
    {
        if ((spectrum != output) && (numOfLines > 0))
        {
            mvr2r(spectrum, output, numOfLines);
        }
    }
    else
    {
        int i;
        float buf[ITF_MAX_FILTER_ORDER + N_MAX];
        float* p = buf + ITF_MAX_FILTER_ORDER;
        set_f(buf, 0, ITF_MAX_FILTER_ORDER);
        mvr2r(spectrum, p, numOfLines);
        for (j = 0; j < numOfLines; j++, p++)
        {
            output[j] = p[0];
            for (i = 1; i < order; i++)
            {
                output[j] += A[i] * p[-i];
            }
        }
    }
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
static float ITF_AutoToLPcoef(const float input[],  float a[], int order)
{
    int i, j;
    float tmp, tmp2;
    float workBuffer[2*ITF_MAX_FILTER_ORDER];
    float parCoeff[ITF_MAX_FILTER_ORDER];
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

    {
        /* Convert ParCor / reflection coefficients to LPC */
        a[0] = 1.0f;
        a[1] = parCoeff[0];

        for(i=1; i<order; i++)
        {
            for(j=0; j<i/2; j++)
            {
                tmp = a[j+1];
                a[j+1] += parCoeff[i] * a[i-1-j+1];
                a[i-1-j+1] += parCoeff[i] * tmp;
            }
            if (i & 1)
            {
                a[j+1] += parCoeff[i] * a[j+1];
            }

            a[i+1] = parCoeff[i];
        }
    }
    return ( (input[0] + 1e-30f) / (workBuffer[0] + 1e-30f) );
}


TNS_ERROR ITF_Apply(float spectrum[],
                    short int startLine,
                    short int stopLine,
                    const float* A,
                    int order)
{


    ITF_TnsFilter(&spectrum[startLine], stopLine-startLine,
                  A, order,
                  &spectrum[startLine]);

    return TNS_NO_ERROR;
}

#define HLM_MIN_NRG 32768.0f

static float Autocorrelation(float const x[], int n, int lag)
{
    return dotp(x, x+lag, n-lag);
}

/** Linear prediction analysis/synthesis filter definition.
  * @param order filter order.
  * @param parCoeff filter (PARCOR) coefficients.
  * @param state state of the filter. Must be at least of 'order' size.
  * @param x the current input value.
  * @return the output of the filter.
  */

/********************************/
/*      Interface functions     */
/********************************/

#define MAX_SUBDIVISIONS 3

int ITF_Detect(float const pSpectrum[],
               short int startLine, short int stopLine, int maxOrder,
               float* A, float* predictionGain, int* curr_order)
{
    float norms[MAX_SUBDIVISIONS] = {0};
    int const nSubdivisions = MAX_SUBDIVISIONS;
    int const idx0 = startLine;
    int const idx1 = stopLine;
    int iSubdivisions;
    if (maxOrder <= 0)
    {
        return 0;
    }

    /* Calculate norms for each spectrum part */
    for (iSubdivisions = 0; iSubdivisions < nSubdivisions; iSubdivisions++)
    {
        int const iStartLine = idx0 + (idx1-idx0)*iSubdivisions/nSubdivisions;
        int const iEndLine = idx0 + (idx1-idx0)*(iSubdivisions+1)/nSubdivisions;
        /* Variable initialization */
        norms[iSubdivisions] = sum2_f(pSpectrum+iStartLine-IGF_START_MN, iEndLine-iStartLine);                      /* No need for indirect, PTR_INIT used instead */
    }

    /* Calculate normalized autocorrelation for spectrum subdivision and get TNS filter parameters based on it */
    {
        float rxx[ITF_MAX_FILTER_ORDER+1] = {0};                                                                /* WMOPS: This initialization is required */
        int const spectrumLength = idx1 - idx0;
        /* Variable initialization */
        for (iSubdivisions = 0; (iSubdivisions < nSubdivisions) && (norms[iSubdivisions] > HLM_MIN_NRG); iSubdivisions++)
        {
            float const fac = 1.0f/norms[iSubdivisions];
            int const iStartLine = idx0 + spectrumLength*iSubdivisions/nSubdivisions;
            int const iEndLine = idx0 + spectrumLength*(iSubdivisions+1)/nSubdivisions;
            float const * pWindow = tnsAcfWindow;
            float const * pEndWindowCoeff = &tnsAcfWindow[sizeof(tnsAcfWindow)/sizeof(tnsAcfWindow[0])]; /* WMOPS: Not counted as it is used only in the assertion */
            int lag;
            (void)pEndWindowCoeff;
            /* For additional loop condition */
            /* Variable initialization */
            for (lag = 1; lag <= maxOrder; lag++)
            {
                rxx[lag] += fac * (*pWindow) * Autocorrelation(pSpectrum+iStartLine-IGF_START_MN, iEndLine-iStartLine, lag);
                pWindow++;
            }
        }
        *predictionGain = 0;
        if (iSubdivisions == nSubdivisions) /* meaning there is no subdivision with low energy */
        {
            rxx[0] = (float)nSubdivisions;

            /* Limit the maximum order to spectrum length/4 */
            *predictionGain = ITF_AutoToLPcoef(rxx, A, min(maxOrder, spectrumLength/4) );

            *curr_order = maxOrder;
        }
    }

    return 1;
}

