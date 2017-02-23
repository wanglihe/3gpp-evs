/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "stat_com.h"
#include "prot.h"
#include <assert.h>




/********************************/
/*      Helper functions        */
/********************************/

/** Put nBits long encoded value from *pStream into bitstream. Using the function EncodeValue for encoding. */
static int PutIntoBitstream(int const ** pStream, TEncodeValue EncodeValue, int index, Encoder_State *st, int nBits)
{
    int const value = *(*pStream)++;
    int const codedValue = EncodeValue(value, index);
    /* Variable initialization */
    push_next_indice(st, codedValue, nBits);

    return value;
}

/** Get nBits long value from bitstream into *pStream. */
static int GetFromBitstream(Decoder_State *st, TDecodeValue DecodeValue, int index, int nFixedBits, int ** pStream)
{
    int value = 0;
    if (DecodeValue != NULL)
    {
        DecodeValue(st, index, &value);
    }
    else
    {
        value = get_next_indice(st, nFixedBits);
    }
    *(*pStream)++ = value;

    return value;
}

static int FixedWidthEncoding(int value, int index)
{
    (void)index;
    return value;
}

/********************************/
/*      Interface functions     */
/********************************/

void GetParameters(ParamsBitMap const * paramsBitMap, int nArrayLength, void const * pParameter, int ** pStream, int * pnSize, int * pnBits)
{
    int index;
    int iParam, nParams;
    int value;
    void const * pSubStruct;


    assert((paramsBitMap != NULL) && (nArrayLength > 0) && (pParameter != NULL) && (pStream != NULL) && (pnSize != NULL) && (pnBits != NULL));
    nParams = paramsBitMap->nParams;
    for (index = 0; index < nArrayLength; index++)
    {
        for (iParam = 0; iParam < nParams; iParam++)
        {
            ParamBitMap const * const param = & paramsBitMap->params[iParam];                                       /* WMOPS: Just a shortcut */

            pSubStruct = param->GetParamValue(pParameter, index, &value);
            /* If a function for encoding/decoding value is defined than it should take care of 0 */
            if (param->fZeroAllowed || (param->EncodeValue != NULL))
            {
                *(*pStream)++ = value;
            }
            else
            {
                *(*pStream)++ = value - 1;
            }
            ++*pnSize;
            *pnBits += (param->nBits != 0) ? param->nBits : param->GetNumberOfBits(value, index);
            if ((param->pSubParamBitMap != NULL) && (value > 0))
            {
                GetParameters(param->pSubParamBitMap, value, (pSubStruct != NULL) ? pSubStruct : pParameter, pStream, pnSize, pnBits);
            }
        }
    }
}

void SetParameters(ParamsBitMap const * paramsBitMap, int nArrayLength, void * pParameter, int const ** pStream, int * pnSize)
{
    int index;
    int iParam, nParams;
    int value;
    void * pSubStruct;
    assert((paramsBitMap != NULL) && (nArrayLength > 0) && (pParameter != NULL) && (pStream != NULL) && (pnSize != NULL));
    nParams = paramsBitMap->nParams;
    for (index = 0; index < nArrayLength; index++)
    {
        for (iParam = 0; iParam < nParams; iParam++)
        {
            ParamBitMap const * const param = & paramsBitMap->params[iParam];                                       /* WMOPS: Just a shortcut */
            /* If a function for encoding/decoding value is defined than it should take care of 0 */

            value = *(*pStream)++ + (param->fZeroAllowed || (param->EncodeValue != NULL) ? 0 : 1);
            pSubStruct = param->SetParamValue(pParameter, index, value);
            ++*pnSize;
            if ((param->pSubParamBitMap != NULL) && (value > 0))
            {
                SetParameters(param->pSubParamBitMap, value, (pSubStruct != NULL) ? pSubStruct : pParameter, pStream, pnSize);
            }
        }
    }
}

void WriteToBitstream(ParamsBitMap const * paramsBitMap, int nArrayLength, int const ** pStream, int * pnSize, Encoder_State *st, int * pnBits)
{
    int index;
    int iParam, nParams;
    assert((paramsBitMap != NULL) && (nArrayLength > 0) && (pStream != NULL) && (pnSize != NULL) && (st != NULL) && (pnBits != NULL));
    nParams = paramsBitMap->nParams;
    for (index = 0; index < nArrayLength; index++)
    {
        for (iParam = 0; iParam < nParams; iParam++)
        {
            ParamBitMap const * const param = & paramsBitMap->params[iParam];                                       /* WMOPS: Just a shortcut */
            int nBits;
            /* If a function for encoding/decoding value is defined than it should take care of 0 */
            int fShiftValue;
            TEncodeValue EncodeValue;
            int value;

            nBits = (param->nBits != 0) ? param->nBits : param->GetNumberOfBits(**pStream, index);
            fShiftValue = !param->fZeroAllowed && (param->EncodeValue == NULL);
            EncodeValue = (param->EncodeValue == NULL) ? &FixedWidthEncoding : param->EncodeValue;
            value = PutIntoBitstream(pStream, EncodeValue, index, st, nBits) + (fShiftValue ? 1 : 0);
            ++*pnSize;
            *pnBits += nBits;
            if ((param->pSubParamBitMap != NULL) && (value > 0))
            {
                WriteToBitstream(param->pSubParamBitMap, value, pStream, pnSize, st, pnBits);
            }
        }
    }
}

void ReadFromBitstream(ParamsBitMap const * paramsBitMap, int nArrayLength, Decoder_State *st, int ** pStream, int * pnSize)
{
    int index;
    int iParam, nParams;
    int fShiftValue;
    int value;
    assert((paramsBitMap != NULL) && (nArrayLength > 0) && (pStream != NULL) && (pnSize != NULL) && (st != NULL));
    nParams = paramsBitMap->nParams;
    for (index = 0; index < nArrayLength; index++)
    {
        for (iParam = 0; iParam < nParams; iParam++)
        {
            ParamBitMap const * param = & paramsBitMap->params[iParam];                                             /* WMOPS: Just a shortcut */
            /* If a function for encoding/decoding value is defined than it should take care of 0 */

            fShiftValue = !param->fZeroAllowed && (param->EncodeValue == NULL);
            value = GetFromBitstream(st, param->DecodeValue, index, param->nBits, pStream) + (fShiftValue ? 1 : 0);
            if ((param->pSubParamBitMap != NULL) && (value > 0))
            {
                ReadFromBitstream(param->pSubParamBitMap, value, st, pStream, pnSize);
            }
        }
    }
    *pnSize += nParams*nArrayLength;
}
