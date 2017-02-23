/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "cnst.h"
#include "rom_com.h"
#include "prot.h"
#include <memory.h>
#include <math.h>
#include <assert.h>



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

TNS_ERROR ReadTnsData(STnsConfig const * pTnsConfig, Decoder_State * st, int * pnBits, int * stream, int * pnSize)
{
    int start_bit_pos;
    start_bit_pos = st->next_bit_pos;
    if (pTnsConfig->nMaxFilters > 1)
    {
        if (pTnsConfig->iFilterBorders[0] < 512)
        {
            ReadFromBitstream(&tnsEnabledSWBTCX10BitMap, 1, st, &stream, pnSize);
        }
        else
        {
            ReadFromBitstream(&tnsEnabledSWBTCX20BitMap, 1, st, &stream, pnSize);
        }
    }
    else
    {
        ReadFromBitstream(&tnsEnabledWBTCX20BitMap, 1, st, &stream, pnSize);
    }

    *pnBits = st->next_bit_pos - start_bit_pos;
    return TNS_NO_ERROR;
}

int DecodeTnsData(STnsConfig const * pTnsConfig, int const * stream, int * pnSize, STnsData * pTnsData)
{
    ResetTnsData(pTnsData);
    if (pTnsConfig->nMaxFilters > 1)
    {
        if (pTnsConfig->iFilterBorders[0] < 512)
        {
            SetParameters(&tnsEnabledSWBTCX10BitMap, 1, pTnsData, &stream, pnSize);
        }
        else
        {
            SetParameters(&tnsEnabledSWBTCX20BitMap, 1, pTnsData, &stream, pnSize);
        }
    }
    else
    {
        SetParameters(&tnsEnabledWBTCX20BitMap, 1, pTnsData, &stream, pnSize);
    }
    return (pTnsData->nFilters > 0) ? TRUE : FALSE;
}

