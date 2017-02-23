/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "options.h"
#include "prot.h"
#include "cnst.h"
#include "stat_enc.h"

/**********************************************************************/ /*
write single bit to stream
**************************************************************************/
static void IGF_write_bit(Encoder_State                                *st,                 /**< in:     | encoder state structure  */
                          int                                          *bitCount,           /**< in/out: | bit counter              */
                          int                                           bit                 /**< in/out: | bit counter              */
                         )
{
    IGFCommonFuncsWriteSerialBit(st, bitCount, bit);

    return;
}

/**********************************************************************/ /*
write bits to stream
**************************************************************************/
static void IGF_write_bits(Encoder_State                               *st,                 /**< in:     | encoder state structure  */
                           int                                         *bitCount,           /**< in/out: | bit counter              */
                           int                                          value,              /**< in/out: | bit counter              */
                           int                                          bits                /**< in:     | number of bits           */
                          )
{
    while (bits--)
    {
        IGF_write_bit(st, bitCount, ((value & (1 << bits)) == 0) ? 0 : 1);
    }

    return;
}


/**********************************************************************/ /*
envelope estimation
**************************************************************************/
static void IGF_CalculateEnvelope(const IGF_ENC_INSTANCE_HANDLE         hInstance,          /**< in:     | instance handle of IGF Encoder                     */
                                  float                                *pMDCTSpectrum,      /**< in:     | MDCT spectrum                                      */
                                  float                                *pPowerSpectrum,     /**< in:     | MDCT^2 + MDST^2 spectrum, or estimate              */
                                  const int                             igfGridIdx          /**< in:     | IGF grid index                                     */
                                 )
{
    IGF_ENC_PRIVATE_DATA_HANDLE hPrivateData;
    H_IGF_GRID hGrid;
    int    *swb_offset;
    int     sfb;                  /* this is the actual scalefactor band */
    int     width;                /* this is width in subbands of the actual scalefactor band */
    int     tile_idx;
    int     strt_cpy;
    float   gain;                 /* the gain which has to be applied to the source tile to get the destination energy */
    int     sb;
    float   sfbEnergyR;
    float   sfbEnergyC;           /* the energy of the destination region of the tile */
    float   sfbEnergyTileR;
    float   sfbEnergyTileC;       /* the energy of the destination region of the tile */

    mvr2r(pMDCTSpectrum + IGF_START_MN, hInstance->spec_be_igf, hInstance->infoStopLine-IGF_START_MN);
    hPrivateData = &hInstance->igfData;
    hGrid        = &hPrivateData->igfInfo.grid[(int)igfGridIdx];
    swb_offset   = hGrid->swb_offset;


    for (tile_idx = 0; tile_idx < hGrid->nTiles; tile_idx++)
    {
        strt_cpy = hGrid->sbWrap[tile_idx];

        for (sfb = hGrid->sfbWrap[tile_idx]; sfb < hGrid->sfbWrap[tile_idx + 1]; sfb++)
        {


            width           = swb_offset[sfb + 1] - swb_offset[sfb];
            sfbEnergyTileR  = FLT_EPSILON;
            sfbEnergyTileC  = FLT_EPSILON;
            sfbEnergyC      = FLT_EPSILON;

            if (pPowerSpectrum)
            {

                for (sb = swb_offset[sfb]; sb < swb_offset[sfb + 1]; sb++)
                {
                    sfbEnergyC     += pPowerSpectrum[sb];
                    sfbEnergyTileR += pMDCTSpectrum[strt_cpy] * pMDCTSpectrum[strt_cpy];
                    sfbEnergyTileC += pPowerSpectrum[strt_cpy];


                    strt_cpy++;
                }

                sfbEnergyTileR /= width;
                gain = (float)(sfbEnergyTileR * (sfbEnergyC / sfbEnergyTileC));

            }
            else
            {
                sfbEnergyR  = FLT_EPSILON + sum2_f(pMDCTSpectrum + swb_offset[ sfb ], width) / width;
                gain        = (float)(sfbEnergyR);

            }

            gain = 0.5f + (float) ((2.885390081777927f * log(gain) + 16.f));
            gain = min(gain, 91.f);  /* 13+15+63, see arithocde encode residual */
            gain = max(gain,  0.f);

            hPrivateData->igfScfQuantized[sfb] = (int)(gain);
        }
    }

    return;
}

/**********************************************************************/ /*
writes IGF SCF values
**************************************************************************/
static int IGF_WriteEnvelope(                                                               /**< out:    | number of bits writen                                                        */
    const IGF_ENC_INSTANCE_HANDLE              hInstance,          /**< in:     | instance handle of IGF Encoder                                               */
    Encoder_State                             *st,                 /**< in:     | encoder state                                                                */
    int                                       *pBitOffset,         /**< in:     | ptr to bitOffset counter                                                     */
    const int                                  igfGridIdx,         /**< in:     | igf grid index see declaration of IGF_GRID_IDX for details                   */
    const int                                  isIndepFlag,        /**< in:     | if 1 frame is independent, 0 = frame is coded with data from previous frame  */
    int                                       *igfAllZero          /**< in:     | returns 1 if all IGF scfs are zero, else 0                                   */
)
{
    IGF_ENC_PRIVATE_DATA_HANDLE hPrivateData;
    H_IGF_GRID hGrid;
    int totBitCount;
    int startBitCount;
    int sfb;

    startBitCount = *pBitOffset;
    totBitCount   = 0;
    *igfAllZero   = 1;
    hPrivateData  = &hInstance->igfData;
    hGrid         = &hPrivateData->igfInfo.grid[igfGridIdx];

    for (sfb = hGrid->startSfb; sfb < hGrid->stopSfb; sfb++)
    {
        if (hPrivateData->igfScfQuantized[sfb]!=0)
        {
            *igfAllZero = 0;
            break;
        }
    }

    if (*igfAllZero)
    {
        IGF_write_bit(st, pBitOffset, 1);
        if (NULL == st)
        {
            IGFSCFEncoderSaveContextState(&hPrivateData->hIGFSCFArithEnc);
        }
        IGFSCFEncoderReset(&hPrivateData->hIGFSCFArithEnc);
        if (NULL == st)
        {
            IGFSCFEncoderRestoreContextState(&hPrivateData->hIGFSCFArithEnc);
        }
    }
    else
    {
        IGF_write_bit(st, pBitOffset, 0);

        if (NULL == st)
        {
            IGFSCFEncoderSaveContextState(&hPrivateData->hIGFSCFArithEnc);
        }

        *pBitOffset = IGFSCFEncoderEncode( &hPrivateData->hIGFSCFArithEnc,st, *pBitOffset, &hPrivateData->igfScfQuantized[hGrid->startSfb], isIndepFlag, (NULL != st) );
        if (NULL == st)
        {
            IGFSCFEncoderRestoreContextState(&hPrivateData->hIGFSCFArithEnc);
        }
    }
    totBitCount = *pBitOffset - startBitCount;

    return totBitCount;
}

/**********************************************************************/ /*
identifies significant spectral content
**************************************************************************/
static float IGF_ErodeSpectrum(                                                             /**< out:    | highPassEnergy                 */
    const IGF_ENC_INSTANCE_HANDLE      hInstance,                /**< in:     | instance handle of IGF Encoder */
    float                             *pSpectrum,                /**< in/out: | MDCT spectrum                  */
    float                             *pPowerSpectrum,           /**< in/out: | power spectrum                 */
    const int                          igfGridIdx                /**< in: Q0  | IGF grid index                 */
)
{
    IGF_ENC_PRIVATE_DATA_HANDLE hPrivateData;
    H_IGF_GRID hGrid;
    int i;
    int igfBgn;
    int igfEnd;
    float highPassEner;
    float lastLine;
    float nextLine;
    float factor;
    int *igfScaleF;
    int startSfb;
    int stopSfb;
    int  *swb_offset;
    int tmp;
    int sfb;
    int line;

    highPassEner = 0.f;
    hPrivateData = &hInstance->igfData;
    hGrid        = &hPrivateData->igfInfo.grid[(int)igfGridIdx];
    igfBgn       = hGrid->startLine;
    igfEnd       = hGrid->stopLine;
    startSfb     = hGrid->startSfb;
    stopSfb      = hGrid->stopSfb;
    swb_offset   = hGrid->swb_offset;
    igfScaleF    = hPrivateData->igfScfQuantized;

    if (NULL == pPowerSpectrum)
    {
        for (i = igfBgn; i< hGrid->infoGranuleLen; i++)
        {
            pSpectrum[i] = 0.f;
        }
        return 0;
    }

    if (igfBgn > 0)
    {
        for (i = 0; i < igfBgn; i++)
        {
            highPassEner += (float)i * pPowerSpectrum[i];
        }

        if ( hPrivateData->igfInfo.bitRateIndex == IGF_BITRATE_SWB_9600 ||
                hPrivateData->igfInfo.bitRateIndex == IGF_BITRATE_RF_SWB_13200 ||
                hPrivateData->igfInfo.bitRateIndex == IGF_BITRATE_SWB_13200 )
        {
            factor    = 1.f;
        }
        else
        {
            factor    = 2.f;
        }
        highPassEner /= igfBgn * factor; /* for 9.6kbs use 1.f */
        lastLine      = pSpectrum[i-1];
        nextLine      = (pPowerSpectrum[i-1] < highPassEner) ? 0.0f : pSpectrum[i];

        for (/*i*/; i < igfEnd - 1; i++)
        {
            if (pPowerSpectrum[i] < highPassEner)
            {
                lastLine       = pSpectrum[i];
                pSpectrum[i]   = nextLine;
                nextLine       = 0.0f;
            }
            else
            {
                pSpectrum[i-1] = lastLine;
                lastLine       = pSpectrum[i];
                nextLine       = pSpectrum[i+1];
            }
        }

        /* i == igfEnd - 1 */
        if (pPowerSpectrum[i] < highPassEner)
        {
            pSpectrum[i] = 0.f;
        }
    }

    /* delete spectrum above igfEnd: */
    for (i = igfEnd; i < hGrid->infoGranuleLen; i++)
    {
        pSpectrum[i] = 0.f;
        pPowerSpectrum[i] = 0.f;
    }

    if (NULL != pPowerSpectrum)
    {
        for (sfb = startSfb; sfb < stopSfb; sfb++)
        {
            tmp = 0;
            for (line = swb_offset[sfb]; line < swb_offset[sfb+1]; line++)
            {
                if (pSpectrum[line] != 0.f) tmp++;
            }
            if(tmp && igfScaleF[sfb]) igfScaleF[sfb]--;
        }
    }

    return highPassEner;
}

/**********************************************************************/ /*
crest factor calculation
**************************************************************************/
static float IGF_getCrest(                                                                  /**< out:    | crest factor                 */
    const float                                  *powerSpectrum,      /**< in:     | power spectrum               */
    const int                                     start,              /**< in:     | start subband index          */
    const int                                     stop                /**< in:     | stop subband index           */
)
{
    int i;
    int x;
    int x_eff   = 0;
    int x_max   = 0;
    float crest = 1.f;

    for (i = start; i < stop; i++)
    {
        x      = max(0, (int)(log(max(FLT_MIN, powerSpectrum[i])) * INV_LOG_2));
        x_eff += x * x;

        if (x > x_max)
        {
            x_max = x;
        }
    }

    x_eff /= (stop - start);

    if (x_eff > 0 && x_max > 0)
    {
        crest = max(1.f, (float)x_max/sqrt(x_eff));
    }

    return crest;
}

/*************************************************************************
calculates spectral flatness measurment
**************************************************************************/
static float IGF_getSFM(                                                                    /**< out:    | SFM value              */
    const float                                    *powerSpectrum,      /**< in:     | energies               */
    const int                                       start,              /**< in:     | start subband index    */
    const int                                       stop                /**< in:     | stop subband index     */
)
{
    int n;
    int i;
    int num;
    float denom;
    float numf;
    float tmp;
    float sfm;

    num   = 0.f;
    denom = 1.f;
    sfm   = 1.f;

    for (i = start; i < stop; i++)
    {
        tmp    = powerSpectrum[i];
        n      = max(0, (int)(log(max(FLT_MIN, tmp)) * INV_LOG_2));
        num   += n;
        denom += tmp;
    }

    numf   = (float)num / (float)(stop - start);
    denom /= (float)(stop - start);

    if (denom != 0)
    {
        sfm = min(((float)pow(2.0, numf + 0.5f) / denom), 1.0f);
    }

    return sfm;
}

/**********************************************************************/ /*
calculates the IGF whitening levels by SFM and crest
**************************************************************************/
static void IGF_Whitening(const IGF_ENC_INSTANCE_HANDLE                 hInstance,          /**< in:     | instance handle of IGF Encoder               */
                          float                                        *powerSpectrum,      /**< in: Q31 | MDCT/MDST power spectrum                     */
                          const int                                     igfGridIdx,         /**< in: Q0  | IGF grid index                               */
                          int                                           isTransient,        /**< in: Q0  | boolean, indicating if transient is detected */
                          int                                           last_core_acelp     /**< in: Q0  | indictaor if last frame was acelp coded      */
                         )
{
    IGF_ENC_PRIVATE_DATA_HANDLE hPrivateData;
    H_IGF_GRID hGrid;
    int p;
    float tmp;
    float SFM;

    hPrivateData = &hInstance->igfData;
    hGrid        = &hPrivateData->igfInfo.grid[(int)igfGridIdx];

    if (igfGridIdx != IGF_GRID_LB_NORM)
    {
        for (p = 0; p < hGrid->nTiles; p++)
        {
            /* reset filter */
            hPrivateData->prevSFM_FIR[p] = 0.f;
            hPrivateData->prevSFM_IIR[p] = 0.f;

            /* preset values: */
            hPrivateData->igfCurrWhiteningLevel[p] = IGF_WHITENING_OFF;
        }
    }

    for (p = 0; p < IGF_MAX_TILES; p++)
    {
        /* update prev data: */
        hPrivateData->igfPrevWhiteningLevel[p] = hPrivateData->igfCurrWhiteningLevel[p];
        /* preset values: */
        hPrivateData->igfCurrWhiteningLevel[p] = IGF_WHITENING_OFF;
    }

    if (!(isTransient || hPrivateData->wasTransient))
    {
        if (powerSpectrum)
        {
            for (p = 0; p < hGrid->nTiles; p++)
            {
                tmp = IGF_getSFM(powerSpectrum, hGrid->tile[p], hGrid->tile[p + 1]) / IGF_getCrest(powerSpectrum, hGrid->tile[p], hGrid->tile[p + 1]);

                if(last_core_acelp || hPrivateData->wasTransient)
                {
                    hPrivateData->prevSFM_FIR[p] = hPrivateData->prevSFM_IIR[p] = tmp;
                }

                SFM = tmp + hPrivateData->prevSFM_FIR[p] + 0.5f * hPrivateData->prevSFM_IIR[p];
                SFM = min(2.7f, SFM);

                hPrivateData->prevSFM_FIR[p] = tmp;
                hPrivateData->prevSFM_IIR[p] = SFM;

                if (SFM > hGrid->whiteningThreshold[1][p])
                {
                    hPrivateData->igfCurrWhiteningLevel[p] = IGF_WHITENING_STRONG;
                }
                else if (SFM > hGrid->whiteningThreshold[0][p])
                {
                    hPrivateData->igfCurrWhiteningLevel[p] = IGF_WHITENING_MID;
                }
                else
                {
                    hPrivateData->igfCurrWhiteningLevel[p] = IGF_WHITENING_OFF;
                }
            }

            switch (hPrivateData->igfInfo.bitRateIndex)
            {
            case IGF_BITRATE_WB_9600:
            case IGF_BITRATE_RF_WB_13200:
            case IGF_BITRATE_RF_SWB_13200:
            case IGF_BITRATE_SWB_9600:
            case IGF_BITRATE_SWB_16400:
            case IGF_BITRATE_SWB_24400:
            case IGF_BITRATE_SWB_32000:
            case IGF_BITRATE_FB_16400:
            case IGF_BITRATE_FB_24400:
            case IGF_BITRATE_FB_32000:
                hPrivateData->igfCurrWhiteningLevel[hGrid->nTiles - 1] = hPrivateData->igfCurrWhiteningLevel[hGrid->nTiles - 2];
                break;
            default:
                break;
            }
        }
        else
        {
            for (p = 0; p < hGrid->nTiles; p++)
            {
                hPrivateData->igfCurrWhiteningLevel[p] = IGF_WHITENING_MID;
            }
        }
    }
    else
    {
        /* reset filter */
        for (p = 0; p < IGF_MAX_TILES; p++)
        {
            hPrivateData->prevSFM_FIR[p] = 0.f;
            hPrivateData->prevSFM_IIR[p] = 0.f;
        }
    }

    hPrivateData->wasTransient = isTransient;

    return;
}


/**********************************************************************/ /*
write whitening levels into bitstream
**************************************************************************/
static int IGF_WriteWhiteningTile(                                                          /**< out:    | number of bits written     */
    Encoder_State                        *st,                 /**< in:     | encoder state handle       */
    int                                  *pBitOffset,         /**< in:     | ptr to bitOffset counter   */
    int                                   whiteningLevel      /**< in: Q0  | whitening levels to write  */
)
{
    int totBitCount;
    int startBitCount;

    totBitCount   = 0;
    startBitCount = *pBitOffset;

    if (whiteningLevel == IGF_WHITENING_MID)
    {
        IGF_write_bits(st, pBitOffset, 0, 1);
    }
    else
    {
        IGF_write_bits(st, pBitOffset, 1, 1);
        if(whiteningLevel == IGF_WHITENING_OFF)
        {
            IGF_write_bits(st, pBitOffset, 0, 1);
        }
        else
        {
            IGF_write_bits(st, pBitOffset, 1, 1);
        }
    }
    totBitCount = *pBitOffset - startBitCount;

    return totBitCount;
}

/**********************************************************************/ /*
writes the whitening levels
**************************************************************************/
static int IGF_WriteWhiteningLevels(                                                        /**< out: Q0 | total number of bits written                                                 */
    const IGF_ENC_INSTANCE_HANDLE       hInstance,          /**< in:     | instance handle of IGF encoder                                               */
    Encoder_State                      *st,                 /**< in:     | encoder state                                                                */
    int                                *pBitOffset,         /**< in:     | ptr to bitOffset counter                                                     */
    const int                           igfGridIdx,         /**< in: Q0  | igf grid index see declaration of IGF_GRID_IDX for details                   */
    const int                           isIndepFlag         /**< in: Q0  | if 1 frame is independent, 0 = frame is coded with data from previous frame  */
)
{
    IGF_ENC_PRIVATE_DATA_HANDLE hPrivateData;
    H_IGF_GRID hGrid;
    int p;
    int nTiles;
    int totBitCount;
    int isSame;
    int startBitCount;

    totBitCount   = 0;
    isSame        = 1;
    startBitCount = *pBitOffset;
    hPrivateData  = &hInstance->igfData;
    hGrid         = &hPrivateData->igfInfo.grid[igfGridIdx];
    nTiles        = hGrid->nTiles;

    if (isIndepFlag)
    {
        isSame = 0;
    }
    else
    {
        for (p = 0; p < nTiles ; p++)
        {
            if (hPrivateData->igfCurrWhiteningLevel[p] != hPrivateData->igfPrevWhiteningLevel[p])
            {
                isSame = 0;
                break;
            }
        }
    }
    if (isSame)
    {
        IGF_write_bits(st, pBitOffset, 1, 1);
    }
    else
    {
        if (!isIndepFlag)
        {
            IGF_write_bits(st, pBitOffset, 0, 1);
        }
        IGF_WriteWhiteningTile(st, pBitOffset, hPrivateData->igfCurrWhiteningLevel[0]);
        for (p = 1; p < nTiles ; p++)
        {
            isSame = 1;
            if (hPrivateData->igfCurrWhiteningLevel[p] != hPrivateData->igfCurrWhiteningLevel[p-1])
            {
                isSame = 0;
                break;
            }
        }
        if (!isSame)
        {
            IGF_write_bits(st, pBitOffset, 1, 1);
            for (p = 1; p < nTiles ; p++)
            {
                IGF_WriteWhiteningTile(st, pBitOffset, hPrivateData->igfCurrWhiteningLevel[p]);
            }
        }
        else
        {
            IGF_write_bits(st, pBitOffset, 0, 1);
        }
    }
    totBitCount = *pBitOffset - startBitCount;

    return totBitCount;
}

/**********************************************************************/ /*
write flattening trigger
**************************************************************************/
static int IGF_WriteFlatteningTrigger(                                                      /**< out:    | number of bits written         */
    const IGF_ENC_INSTANCE_HANDLE     hInstance,          /**< in:     | instance handle of IGF Encoder */
    Encoder_State                    *st,                 /**< in:     | encoder state                  */
    int                              *pBitOffset          /**< in:     | ptr to bitOffset counter       */
)
{
    int totBitCount;
    int startBitCount;
    int flatteningTrigger;

    totBitCount       = 0;
    startBitCount     = *pBitOffset;
    flatteningTrigger = hInstance->flatteningTrigger;

    IGF_write_bits(st, pBitOffset, flatteningTrigger, 1);

    totBitCount = *pBitOffset - startBitCount;

    return totBitCount;
}

/**********************************************************************/ /*
updates the start/stop frequency of IGF according to igfGridIdx
**************************************************************************/
static void IGF_UpdateInfo(const IGF_ENC_INSTANCE_HANDLE                hInstance,          /**< in:     | instance handle of IGF Encoder */
                           const int                                    igfGridIdx          /**< in:     | IGF grid index                 */
                          )
{
    IGF_ENC_PRIVATE_DATA_HANDLE hPrivateData;
    H_IGF_GRID hGrid;

    hPrivateData                  = &hInstance->igfData;
    hGrid                         = &hPrivateData->igfInfo.grid[igfGridIdx];
    hInstance->infoStartFrequency = hGrid->startFrequency;
    hInstance->infoStopFrequency  = hGrid->stopFrequency;
    hInstance->infoStartLine      = hGrid->startLine;
    hInstance->infoStopLine       = hGrid->stopLine;

    return;
}

/**********************************************************************/ /*
IGF bitsream writer
**************************************************************************/
int IGFEncWriteBitstream(                                                                /**< out:    | number of bits written per frame                                             */
    const IGF_ENC_INSTANCE_HANDLE               hInstance,          /**< in:     | instance handle of IGF Encoder                                               */
    void                                       *st,                 /**< in:     | encoder state                                                                */
    int                                        *pBitOffset,         /**< in:     | ptr to bitOffset counter                                                     */
    const int                                   igfGridIdx,         /**< in:     | igf grid index see declaration of IGF_GRID_IDX for details                   */
    const int                                   isIndepFlag         /**< in:     | if 1 frame is independent, 0 = frame is coded with data from previous frame  */
)
{
    int igfAllZero;
    int startBitCount = *pBitOffset;

    hInstance->infoTotalBitsPerFrameWritten = 0;
    hInstance->infoFrameCount++;
    if (isIndepFlag)
    {
        hInstance->infoTotalBitsWritten = 0;
    }

    IGF_WriteEnvelope(hInstance,                        /* i: instance handle of IGF Encoder                                              */
                      st,                               /* i: encoder state                                                               */
                      pBitOffset,                       /* i: ptr to bitOffset counter                                                    */
                      igfGridIdx,                       /* i: igf grid index see definition of IGF_GRID_IDX for details                   */
                      isIndepFlag,                      /* i: if 1 frame is independent, 0 = frame is coded with data from previous frame */
                      &igfAllZero);                      /* o: *igfAllZero                                                                 */

    IGF_WriteWhiteningLevels(hInstance,                 /* i: instance handle of IGF Encoder                                              */
                             st,                        /* i: encoder state                                                               */
                             pBitOffset,                /* i: ptr to bitOffset counter                                                    */
                             igfGridIdx,                /* i: igf grid index see definition of IGF_GRID_IDX for details                   */
                             isIndepFlag);              /* i: if 1 frame is independent, 0 = frame is coded with data from previous frame */

    IGF_WriteFlatteningTrigger(hInstance,               /* i: instance handle of IGF Encoder                                              */
                               st,                      /* i: encoder state                                                               */
                               pBitOffset);             /* i: ptr to bitOffset counter                                                    */

    hInstance->infoTotalBitsPerFrameWritten  = (*pBitOffset - startBitCount);
    hInstance->infoTotalBitsWritten         += hInstance->infoTotalBitsPerFrameWritten;

    return hInstance->infoTotalBitsPerFrameWritten;
}

/**********************************************************************/ /*
sets the IGF mode according to given bitrate
**************************************************************************/
void IGFEncSetMode(
    const IGF_ENC_INSTANCE_HANDLE    hInstance,          /**< in:     | instance handle of IGF Encoder */
    const int                        bitRate,            /**< in:     | encoder bitrate                */
    const int                        mode,               /**< in:     | encoder bandwidth mode         */
    int const                        rf_mode             /**< in:       flag to signal the RF mode     */
)
{
    short i;
    IGF_ENC_PRIVATE_DATA_HANDLE hPrivateData;

    hPrivateData = &hInstance->igfData;
    hPrivateData->igfBitstreamBits = 0;
    set_i(hPrivateData->igfScfQuantized, 0, IGF_MAX_SFB);
    set_i(hPrivateData->igfCurrWhiteningLevel, 0, IGF_MAX_TILES);
    set_i(hPrivateData->igfPrevWhiteningLevel, 0, IGF_MAX_TILES);
    for( i=0; i<BITBUFSIZE/8; i++ )
    {
        hPrivateData->igfBitstream[i] = 0;
    }
    hPrivateData->wasTransient     = 0;
    set_f(hPrivateData->prevSFM_FIR, 0, IGF_MAX_TILES);
    set_f(hPrivateData->prevSFM_IIR, 0, IGF_MAX_TILES);


    if (IGFCommonFuncsIGFConfiguration(bitRate, mode, &hPrivateData->igfInfo, rf_mode))
    {
        IGFSCFEncoderOpen( &hPrivateData->hIGFSCFArithEnc, hPrivateData->igfInfo.grid[0].stopSfb - hPrivateData->igfInfo.grid[0].startSfb, bitRate, mode, rf_mode );

        hInstance->infoSamplingRate   = hPrivateData->igfInfo.sampleRate;
        hInstance->infoStartFrequency = hPrivateData->igfInfo.grid[0].startFrequency;
        hInstance->infoStopFrequency  = hPrivateData->igfInfo.grid[0].stopFrequency;
        hInstance->infoStartLine      = hPrivateData->igfInfo.grid[0].startLine;
        hInstance->infoStopLine       = hPrivateData->igfInfo.grid[0].stopLine;
    }
    else
    {
        /* IGF configuration failed -> error! */
        hInstance->infoSamplingRate   = 0;
        hInstance->infoStartFrequency = -1;
        hInstance->infoStopFrequency  = -1;
        hInstance->infoStartLine      = -1;
        hInstance->infoStopLine       = -1;
        fprintf(stderr,"IGFEncSetMode: initialization error!\n");
    }

    /* reset remaining variables */
    hInstance->infoTotalBitsWritten         = 0;
    hInstance->infoTotalBitsPerFrameWritten = 0;
    hInstance->flatteningTrigger            = 0;
    hInstance->tns_predictionGain           = 0;
    set_f(hInstance->spec_be_igf, 0, N_MAX_TCX-IGF_START_MN);
    return;
}

/**********************************************************************/ /*
IGF bitsream concatenation for TCX10 modes
**************************************************************************/
void IGFEncConcatenateBitstream(const IGF_ENC_INSTANCE_HANDLE        hInstance,          /**< in:     | instance handle of IGF Encoder                 */
                                short                                bsBits,             /**< in:     | number of IGF bits written to list of indices  */
                                short                               *next_ind,           /**< in/out: | pointer to actual bit indice                   */
                                short                               *nb_bits,            /**< in/out: | total number of bits already written           */
                                Indice                              *ind_list            /**< in:     | pointer to list of indices                     */
                               )
{
    IGF_ENC_PRIVATE_DATA_HANDLE hPrivateData;

    hPrivateData = &hInstance->igfData;
    *next_ind   -= bsBits;

    indices_to_serial_generic(
        &ind_list[*next_ind],
        bsBits,
        hPrivateData->igfBitstream,
        &hPrivateData->igfBitstreamBits
    );

    *nb_bits                       -= hInstance->infoTotalBitsPerFrameWritten;

    return;
}

/**********************************************************************/ /*
IGF reset bitsream bit counter for TCX10 modes
**************************************************************************/
void IGFEncResetTCX10BitCounter(const IGF_ENC_INSTANCE_HANDLE        hInstance           /**< in:     | instance handle of IGF Encoder */
                               )
{
    IGF_ENC_PRIVATE_DATA_HANDLE hPrivateData;

    hPrivateData                    = &hInstance->igfData;
    hPrivateData->igfBitstreamBits  = 0;
    hInstance->infoTotalBitsWritten = 0;
}

/**********************************************************************/ /*
IGF write concatenated bitsream for TCX10 modes
**************************************************************************/
int IGFEncWriteConcatenatedBitstream(                                                    /**< out:    | total number of bits written   */
    const IGF_ENC_INSTANCE_HANDLE    hInstance,         /**< in:     | instance handle of IGF Encoder */
    void                             *st                /**< in:     | encoder state                  */
)
{
    IGF_ENC_PRIVATE_DATA_HANDLE hPrivateData;
    int i;
    int bitsLeft;
    UWord8 *pBitstream;

    hPrivateData = &hInstance->igfData;
    pBitstream   = hPrivateData->igfBitstream;

    for (i = 0; i < hPrivateData->igfBitstreamBits >> 3; i++)
    {
        push_next_indice(st, pBitstream[i], 8);
    }

    bitsLeft = hPrivateData->igfBitstreamBits & 0x7;
    if(bitsLeft > 0)
    {
        push_next_indice(st, pBitstream[i] >> (8 - bitsLeft), bitsLeft);
    }

    return hInstance->infoTotalBitsWritten;
}

/**********************************************************************/ /*
apply the IGF encoder, main encoder interface
**************************************************************************/
void IGFEncApplyMono(const IGF_ENC_INSTANCE_HANDLE                   hInstance,          /**< in:     | instance handle of IGF Encoder                         */
                     const int                                       igfGridIdx,         /**< in:     | IGF grid index                                         */
                     Encoder_State                                  *st,                 /**< in:     | Encoder state                                          */
                     float                                          *pMDCTSpectrum,      /**< in:     | MDCT spectrum                                          */
                     float                                          *pPowerSpectrum,     /**< in:     | MDCT^2 + MDST^2 spectrum, or estimate                  */
                     int                                             isTCX20,            /**< in:     | flag indicating if the input is TCX20 or TCX10/2xTCX5  */
                     int                                             isTNSActive,        /**< in:     | flag indicating if the TNS is active                   */
                     int                                             last_core_acelp     /**< in:     | indictaor if last frame was acelp coded                */
                    )
{
    float *pPowerSpectrumParameter;                           /* If it is NULL it informs a function that specific handling is needed */

    pPowerSpectrumParameter = !isTNSActive && isTCX20 ? pPowerSpectrum : NULL;

    IGF_UpdateInfo(hInstance,                                 /* i: instance handle of IGF Encoder            */
                   igfGridIdx);                               /* i: IGF grid index                            */

    IGF_CalculateEnvelope(hInstance,                          /* i: instance handle of IGF Encoder            */
                          pMDCTSpectrum,                      /* i: MDCT spectrum                             */
                          pPowerSpectrumParameter,            /* i: MDCT^2 + MDST^2 spectrum, or estimate     */
                          igfGridIdx                          /* i: IGF grid index                            */
                         );


    IGF_Whitening(hInstance,                                  /* i: instance handle of IGF Encoder            */
                  isTCX20 ? pPowerSpectrum : NULL,            /* i: MDCT^2 + MDST^2 spectrum, or estimate     */
                  igfGridIdx,                                 /* i: IGF grid index                            */
                  (st->transientDetection.transientDetector.bIsAttackPresent == 1),
                  last_core_acelp                             /* i: last frame was acelp indicator            */
                 );

    pPowerSpectrumParameter = isTCX20 ? pPowerSpectrum : NULL;

    IGF_ErodeSpectrum(                                        /* o: highpass energy                           */
        hInstance,                              /* i: instance handle of IGF Encoder            */
        pMDCTSpectrum,                          /* i: MDCT spectrum                             */
        pPowerSpectrumParameter,                /* i: MDCT^2 + MDST^2 spectrum, or estimate     */
        igfGridIdx);                            /* i: IGF grid index                            */
}

