/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"

/*---------------------------------------------------------------------*
 * IGF_ApplyTransFac()
 *
 *
 *---------------------------------------------------------------------*/

static int IGF_ApplyTransFac(   /**< out:    | multiplication factor                                                        */
    const int   val,              /**< in: Q15 | input value for multiplication, Q15                                          */
    const float transFac          /**< in: Q14 | multiplicator for variable val, Q14: 1.25f=0x5000, 1.0f=0x4000, 0.5f=0x2000  */
)
{
    int ret = val;

    if (transFac!=1.f)
    {
        ret = round_f(val*transFac);
        ret +=(ret & 1);
    }

    return ret;
}


/*---------------------------------------------------------------------*
 * IGF_MapBitRateToIndex()
 *
 * maps a given bitrate to the IGF_BITRATE index
 *---------------------------------------------------------------------*/

static short IGF_MapBitRateToIndex(       /**< out: Q0 | return bit rate index  */
    int bitRate,                            /**< in:     | bitrate                */
    int mode,                               /**< in:     | bandwidth mode         */
    int rf_mode                             /**< in:     | flag to signal the RF mode */
)
{
    short bitRateIndex = IGF_BITRATE_UNKNOWN;

    switch (mode)
    {
    case IGF_MODE_WB:
        switch (bitRate)
        {
        case 13200:
            if (rf_mode == 1)
            {
                bitRateIndex = IGF_BITRATE_RF_WB_13200;
            }
            break;
        case 9600:
            bitRateIndex = IGF_BITRATE_WB_9600;
            break;
        default:
            assert(0);
        }
        break;
    case IGF_MODE_SWB:
        switch (bitRate)
        {
        case  9600:
            bitRateIndex = IGF_BITRATE_SWB_9600;
            break;
        case 13200:
            bitRateIndex = IGF_BITRATE_SWB_13200;
            if (rf_mode == 1)
            {
                bitRateIndex = IGF_BITRATE_RF_SWB_13200;
            }
            break;
        case 16400:
            bitRateIndex = IGF_BITRATE_SWB_16400;
            break;
        case 24400:
            bitRateIndex = IGF_BITRATE_SWB_24400;
            break;
        case 32000:
            bitRateIndex = IGF_BITRATE_SWB_32000;
            break;
        case 48000:
            bitRateIndex = IGF_BITRATE_SWB_48000;
            break;
        default:
            assert(0);
        }
        break;
    case IGF_MODE_FB:
        switch (bitRate)
        {
        case 16400:
            bitRateIndex = IGF_BITRATE_FB_16400;
            break;
        case 24400:
            bitRateIndex = IGF_BITRATE_FB_24400;
            break;
        case 32000:
            bitRateIndex = IGF_BITRATE_FB_32000;
            break;
        case 48000:
            bitRateIndex = IGF_BITRATE_FB_48000;
            break;
        case 96000:
            bitRateIndex = IGF_BITRATE_FB_96000;
            break;
        case 128000:
            bitRateIndex = IGF_BITRATE_FB_128000;
            break;
        default:
            assert(0);
        }
        break;
    default:
        assert(0);
    }

    return bitRateIndex;
}


/*---------------------------------------------------------------------*
 * IGF_gridSetUp()
 *
 * IGF grid setup
 *---------------------------------------------------------------------*/

static void IGF_gridSetUp(
    H_IGF_GRID  hGrid,              /**< out:    | IGF grid handle                                                    */
    short       bitRateIndex,       /**< in:     | IGF bitrate index                                                  */
    int         sampleRate,         /**< in:     | sample rate                                                        */
    int         frameLength,        /**< in:     | frame length                                                       */
    float       transFac,           /**< in:     | transFac                                                           */
    int         igfMinFq            /**< in:     | IGF minimum frequency indicating lower start frequency for copy up */
)
{
    int t;
    int sfb;
    int swb_offset_len;
    int wrp_sfb;
    const int *swb_offset;
    float bandwidth;

    /* inits */
    swb_offset     = NULL;
    swb_offset_len = 0;

    switch (bitRateIndex)
    {
    case IGF_BITRATE_WB_9600:
    case IGF_BITRATE_RF_WB_13200:
    case IGF_BITRATE_SWB_9600:
    case IGF_BITRATE_RF_SWB_13200:
    case IGF_BITRATE_SWB_13200:
    case IGF_BITRATE_SWB_16400:
    case IGF_BITRATE_SWB_24400:
    case IGF_BITRATE_SWB_32000:
    case IGF_BITRATE_SWB_48000:
        swb_offset     = &swb_offset_LB_new[bitRateIndex][1];
        swb_offset_len = swb_offset_LB_new[bitRateIndex][0];
        mvr2r(&igf_whitening_TH[(int)bitRateIndex][0][0], &hGrid->whiteningThreshold[0][0], IGF_MAX_TILES * 2);
        break;
    case IGF_BITRATE_FB_16400:
    case IGF_BITRATE_FB_24400:
    case IGF_BITRATE_FB_32000:
        swb_offset     = &swb_offset_LB_new[bitRateIndex][1];
        swb_offset_len = swb_offset_LB_new[bitRateIndex][0];
        mvr2r(&igf_whitening_TH[(int)bitRateIndex][0][0], &hGrid->whiteningThreshold[0][0], IGF_MAX_TILES * 2);
        break;
    case IGF_BITRATE_FB_48000:
    case IGF_BITRATE_FB_96000:
    case IGF_BITRATE_FB_128000:
        swb_offset     = &swb_offset_LB_new[bitRateIndex][1];
        swb_offset_len = swb_offset_LB_new[bitRateIndex][0];
        mvr2r(&igf_whitening_TH[(int)bitRateIndex][0][0], &hGrid->whiteningThreshold[0][0], IGF_MAX_TILES * 2);
        break;
    case IGF_BITRATE_UNKNOWN:
    default:
        assert(0);
    }

    for(sfb = 0; sfb < swb_offset_len; sfb++)
    {
        hGrid->swb_offset[sfb] = IGF_ApplyTransFac(swb_offset[sfb], transFac);
    }

    hGrid->infoIsRefined     = 0;
    frameLength              = IGF_ApplyTransFac(frameLength, transFac);
    bandwidth                = (float)sampleRate / 2.0f / (float)frameLength;
    hGrid->swb_offset_len    = swb_offset_len;
    hGrid->startSfb          = 0;
    hGrid->stopSfb           = hGrid->swb_offset_len-1;
    hGrid->startLine         = hGrid->swb_offset[ hGrid->startSfb ];
    hGrid->stopLine          = hGrid->swb_offset[ hGrid->stopSfb ];
    hGrid->startFrequency    = round_f(bandwidth * hGrid->startLine);
    hGrid->stopFrequency     = round_f(bandwidth * hGrid->stopLine);
    hGrid->minSrcSubband     = round_f((igfMinFq * (frameLength)) / (sampleRate >> 1));
    hGrid->minSrcSubband    += hGrid->minSrcSubband % 2;
    hGrid->minSrcFrequency   = round_f(bandwidth * hGrid->minSrcSubband);
    hGrid->infoGranuleLen    = frameLength;
    hGrid->infoTransFac      = transFac;
    hGrid->sfbWrap[0]        = 0;
    hGrid->tile[0]           = hGrid->startLine;

    switch (bitRateIndex)
    {
    /* SWB 13200 */
    case IGF_BITRATE_WB_9600:
        hGrid->nTiles           = 2;
        wrp_sfb                 = 2;

        /*1st*/
        hGrid->sfbWrap[0+1]     = wrp_sfb;
        hGrid->sbWrap[0]        = hGrid->minSrcSubband;
        hGrid->tile[0+1]        = hGrid->swb_offset[wrp_sfb];

        /*2nd*/
        hGrid->sfbWrap[1+1]     = hGrid->stopSfb;
        hGrid->sbWrap[1]        = hGrid->minSrcSubband;
        hGrid->tile[1+1]        = hGrid->swb_offset[hGrid->stopSfb];

        break;
    case IGF_BITRATE_RF_WB_13200:
        hGrid->nTiles           = 2;
        wrp_sfb                 = 2;

        /*1st*/
        hGrid->sfbWrap[0+1]     = wrp_sfb;
        hGrid->sbWrap[0]        = hGrid->minSrcSubband;
        hGrid->tile[0+1]        = hGrid->swb_offset[wrp_sfb];

        /*2nd*/
        hGrid->sfbWrap[1+1]     = hGrid->stopSfb;
        hGrid->sbWrap[1]        = hGrid->minSrcSubband;
        hGrid->tile[1+1]        = hGrid->swb_offset[hGrid->stopSfb];

        break;
    case IGF_BITRATE_SWB_9600:
        hGrid->nTiles           = 3;
        wrp_sfb                 = 1;

        /*1st*/
        hGrid->sfbWrap[0+1]     = wrp_sfb;
        hGrid->sbWrap[0]        = hGrid->minSrcSubband;
        hGrid->tile[0+1]        = hGrid->swb_offset[wrp_sfb];

        /*2nd*/
        wrp_sfb                 = 2;
        hGrid->sfbWrap[1+1]     = wrp_sfb;
        hGrid->sbWrap[1]        = hGrid->minSrcSubband + IGF_ApplyTransFac(32, transFac);
        hGrid->tile[1+1]        = hGrid->swb_offset[wrp_sfb];

        /*3rd*/
        hGrid->sfbWrap[2+1]     = hGrid->stopSfb;
        hGrid->sbWrap[2]        = hGrid->minSrcSubband + IGF_ApplyTransFac(46, transFac);
        hGrid->tile[2+1]        = hGrid->swb_offset[hGrid->stopSfb];

        break;
    case IGF_BITRATE_RF_SWB_13200:
        hGrid->nTiles           = 3;
        wrp_sfb                 = 1;

        /*1st*/
        hGrid->sfbWrap[0+1]     = wrp_sfb;
        hGrid->sbWrap[0]        = hGrid->minSrcSubband;
        hGrid->tile[0+1]        = hGrid->swb_offset[wrp_sfb];

        /*2nd*/
        wrp_sfb                 = 2;
        hGrid->sfbWrap[1+1]     = wrp_sfb;
        hGrid->sbWrap[1]        = hGrid->minSrcSubband + IGF_ApplyTransFac(32, transFac);
        hGrid->tile[1+1]        = hGrid->swb_offset[wrp_sfb];

        /*3rd*/
        hGrid->sfbWrap[2+1]     = hGrid->stopSfb;
        hGrid->sbWrap[2]        = hGrid->minSrcSubband + IGF_ApplyTransFac(46, transFac);
        hGrid->tile[2+1]        = hGrid->swb_offset[hGrid->stopSfb];
        break;
    case IGF_BITRATE_SWB_13200:
        hGrid->nTiles           = 2;
        wrp_sfb                 = 4;

        /*1st*/
        hGrid->sfbWrap[0+1]     = wrp_sfb;
        hGrid->sbWrap[0]        = hGrid->minSrcSubband;
        hGrid->tile[0+1]        = hGrid->swb_offset[wrp_sfb];

        /*2nd*/
        hGrid->sfbWrap[1+1]     = hGrid->stopSfb;
        hGrid->sbWrap[1]        = hGrid->minSrcSubband + IGF_ApplyTransFac(32, transFac);
        hGrid->tile[1+1]        = hGrid->swb_offset[hGrid->stopSfb];

        break;
    case IGF_BITRATE_SWB_16400:
        hGrid->nTiles           = 3;
        wrp_sfb                 = 4;

        /*1st*/
        hGrid->sfbWrap[0+1]     = wrp_sfb;
        hGrid->sbWrap[0]        = hGrid->minSrcSubband;
        hGrid->tile[0+1]        = hGrid->swb_offset[wrp_sfb];

        /*2nd*/
        hGrid->sfbWrap[1+1]     = 6;
        hGrid->sbWrap[1]        = hGrid->minSrcSubband + IGF_ApplyTransFac(48, transFac);
        hGrid->tile[1+1]        = hGrid->swb_offset[6];

        /*3nd*/
        hGrid->sfbWrap[2+1]     = hGrid->stopSfb;
        hGrid->sbWrap[2]        = hGrid->minSrcSubband + IGF_ApplyTransFac(64, transFac);
        hGrid->tile[2+1]        = hGrid->swb_offset[hGrid->stopSfb];

        break;
    case IGF_BITRATE_SWB_24400:
    case IGF_BITRATE_SWB_32000:
        hGrid->nTiles           = 3;
        wrp_sfb                 = 4;

        /*1st*/
        hGrid->sfbWrap[0+1]     = wrp_sfb;
        hGrid->sbWrap[0]        = hGrid->minSrcSubband;
        hGrid->tile[0+1]        = hGrid->swb_offset[wrp_sfb];

        /*2nd*/
        hGrid->sfbWrap[1+1]     = 7;
        hGrid->sbWrap[1]        = hGrid->minSrcSubband + IGF_ApplyTransFac(32, transFac);
        hGrid->tile[1+1]        = hGrid->swb_offset[7];

        /*3nd*/
        hGrid->sfbWrap[2+1]     = hGrid->stopSfb;
        hGrid->sbWrap[2]        = hGrid->minSrcSubband + IGF_ApplyTransFac(64, transFac);
        hGrid->tile[2+1]        = hGrid->swb_offset[hGrid->stopSfb];
        break;
    case IGF_BITRATE_SWB_48000:
        hGrid->nTiles           = 1;
        wrp_sfb                 = hGrid->stopSfb;

        /*1st*/
        hGrid->sfbWrap[0+1]     = hGrid->stopSfb;
        hGrid->sbWrap[0]        = 2*hGrid->startLine - hGrid->stopLine;
        hGrid->tile[0+1]        = hGrid->swb_offset[hGrid->stopSfb];

        break;
    case IGF_BITRATE_FB_16400:
        hGrid->nTiles           = 3;
        wrp_sfb                 = 4;

        /*1st*/
        hGrid->sfbWrap[0+1]     = wrp_sfb;
        hGrid->sbWrap[0]        = hGrid->minSrcSubband;
        hGrid->tile[0+1]        = hGrid->swb_offset[wrp_sfb];
        wrp_sfb                 = 7;

        /*2nd*/
        hGrid->sfbWrap[1+1]     = wrp_sfb;
        hGrid->sbWrap[1]        = hGrid->minSrcSubband;
        hGrid->tile[1+1]        = hGrid->swb_offset[wrp_sfb];

        /*3nd*/
        hGrid->sfbWrap[2+1]     = hGrid->stopSfb;
        hGrid->sbWrap[2]        = hGrid->minSrcSubband;
        hGrid->tile[2+1]        = hGrid->swb_offset[hGrid->stopSfb];

        break;
    case IGF_BITRATE_FB_24400:
    case IGF_BITRATE_FB_32000:
        hGrid->nTiles           = 4;
        wrp_sfb                 = 4;

        /*1st*/
        hGrid->sfbWrap[0+1]     = wrp_sfb;
        hGrid->sbWrap[0]        = hGrid->minSrcSubband;
        hGrid->tile[0+1]        = hGrid->swb_offset[wrp_sfb];
        wrp_sfb                 = 6;

        /*2nd*/
        hGrid->sfbWrap[1+1]     = wrp_sfb;
        hGrid->sbWrap[1]        = hGrid->minSrcSubband + IGF_ApplyTransFac(32, transFac);
        hGrid->tile[1+1]        = hGrid->swb_offset[wrp_sfb];
        wrp_sfb                 = 9;

        /*3nd*/
        hGrid->sfbWrap[2+1]     = wrp_sfb;
        hGrid->sbWrap[2]        = hGrid->minSrcSubband;
        hGrid->tile[2+1]        = hGrid->swb_offset[wrp_sfb];

        /*4nd*/
        hGrid->sfbWrap[3+1]     = hGrid->stopSfb;
        hGrid->sbWrap[3]        = hGrid->minSrcSubband + (hGrid->swb_offset[9] - hGrid->swb_offset[8]);
        hGrid->tile[3+1]        = hGrid->swb_offset[hGrid->stopSfb];
        break;
    case IGF_BITRATE_FB_48000:
    case IGF_BITRATE_FB_96000:
    case IGF_BITRATE_FB_128000:
        hGrid->nTiles           = 1;

        /*1st*/
        hGrid->sfbWrap[0+1]     = hGrid->stopSfb;
        hGrid->sbWrap[0]        = 2*hGrid->startLine - hGrid->stopLine;
        hGrid->tile[0+1]        = hGrid->swb_offset[hGrid->stopSfb];

        break;
    default:
        assert(0);
    }/*switch*/


    /* adapt level envelope: */
    switch(bitRateIndex)
    {
    case IGF_BITRATE_RF_WB_13200:
    case IGF_BITRATE_WB_9600:
        hGrid->gFactor = 0.800f;
        hGrid->fFactor = 0.70f;
        hGrid->lFactor = 0.60f;
        break;
    case IGF_BITRATE_SWB_13200:
    case IGF_BITRATE_FB_16400:
    case IGF_BITRATE_SWB_16400:
        hGrid->gFactor = 0.930f;
        hGrid->fFactor = 0.20f;
        hGrid->lFactor = 0.85f;
        break;
    case IGF_BITRATE_FB_24400:
    case IGF_BITRATE_SWB_24400:
    case IGF_BITRATE_FB_32000:
    case IGF_BITRATE_SWB_32000:
        hGrid->gFactor = 0.965f;
        hGrid->fFactor = 0.20f;
        hGrid->lFactor = 0.85f;
        break;
    case IGF_BITRATE_FB_48000:
    case IGF_BITRATE_SWB_48000:
        hGrid->gFactor = 1.000f;
        hGrid->fFactor = 0.20f;
        hGrid->lFactor = 1.000f;
        break;
    case IGF_BITRATE_SWB_9600:
    case IGF_BITRATE_RF_SWB_13200:
    default:
        hGrid->gFactor = 1.000f;
        hGrid->fFactor = 0.00f;
        hGrid->lFactor = 1.000f;
    }

    for ( t = hGrid->nTiles+1; t < IGF_MAX_TILES; t++)
    {
        hGrid->tile[t] = 0;
        hGrid->sbWrap[t - 1] = 0;
        hGrid->sfbWrap[t]    = 0;
    }

    return;
}


/*---------------------------------------------------------------------*
 * IGFCommonFuncsWriteSerialBit()
 *
 * write bits to bitstream
 *---------------------------------------------------------------------*/

void IGFCommonFuncsWriteSerialBit(
    void *st,                 /**< in:     | encoder/decoder state structure  */
    int  *pBitOffset,         /**< out:    | bit offset                       */
    int   bit                 /**< in:     | value of bit                     */
)
{
    if (st)
    {
        push_next_indice(st, bit, 1);
    }
    (*pBitOffset)++;

    return;
}


/*---------------------------------------------------------------------*
 * IGFCommonFuncsIGFConfiguration()
 *
 * changes the IGF configuration
 *---------------------------------------------------------------------*/

int IGFCommonFuncsIGFConfiguration( /**< out:    | error value: 0 -> error, 1 -> ok   */
    int        bitRate,               /**< in:     | bitrate in bs e.g. 9600 for 9.6kbs */
    int        mode,                  /**< in:     | bandwidth mode                     */
    H_IGF_INFO hIGFInfo,              /**< out:    | IGF info handle                    */
    int        rf_mode                /**< in:       flag to signal the RF mode         */
)
{
    H_IGF_GRID hGrid;
    int retValue;
    int sampleRate;
    int frameLength;
    int igfMinFq;
    int maxHopsize;

    retValue = 0;     /* bitrate index is unknown -> error! */

    /* interface call for reading in settings */
    hIGFInfo->bitRateIndex = IGF_MapBitRateToIndex( bitRate, mode, rf_mode );

    if (hIGFInfo->bitRateIndex != IGF_BITRATE_UNKNOWN)
    {
        retValue = 1; /* no error */

        /* mapping to local values */
        sampleRate  = igfMode[(int)hIGFInfo->bitRateIndex].sampleRate;
        frameLength = igfMode[(int)hIGFInfo->bitRateIndex].frameLength;
        igfMinFq    = igfMode[(int)hIGFInfo->bitRateIndex].igfMinFq;
        maxHopsize  = igfMode[(int)hIGFInfo->bitRateIndex].maxHopsize;

        /* basic information */
        hIGFInfo->sampleRate  = sampleRate;
        hIGFInfo->frameLength = frameLength;
        hIGFInfo->maxHopsize  = maxHopsize;
        hIGFInfo->nfSeed      = 0;

        /* set up regular IGF grid for TCX 20  (transfac = 1.f) */
        hGrid                 = &hIGFInfo->grid[IGF_GRID_LB_NORM];
        IGF_gridSetUp(hGrid, hIGFInfo->bitRateIndex, sampleRate, frameLength, 1.00f, igfMinFq );

        /* set up IGF grid for CELP->TCX 20 transitions (transfac = 1.25) */
        hGrid                 = &hIGFInfo->grid[IGF_GRID_LB_TRAN];
        IGF_gridSetUp( hGrid, hIGFInfo->bitRateIndex, sampleRate, frameLength, 1.25f, igfMinFq );

        /* set up IGF grid for TCX 10 (transfac = 0.5) */
        hGrid                 = &hIGFInfo->grid[IGF_GRID_LB_SHORT];
        IGF_gridSetUp( hGrid, hIGFInfo->bitRateIndex, sampleRate, frameLength, 0.50f, igfMinFq );
    }

    return retValue;
}


/*---------------------------------------------------------------------*
 * IGFCommonFuncsIGFGetCFTables()
 *
 * selects cumulative frequency tables and offsets for the IGF SCF arithmetic coder
 *---------------------------------------------------------------------*/

int IGFCommonFuncsIGFGetCFTables(            /**< out:    | error value: 0 -> error, 1 -> ok     */
    int                    bitRate,            /**< in:     | bitrate in bs e.g. 9600 for 9.6kbs   */
    int                    mode,               /**< in:     | bandwidth mode                       */
    int                    rf_mode,            /**< in:     | flag to signal the RF mode           */
    const unsigned short **cf_se00,            /**< out:    | CF table for t == 0 and f == 0       */
    const unsigned short **cf_se01,            /**< out:    | CF table for t == 0 and f == 1       */
    short                 *cf_off_se01,        /**< out:    | offset for CF table above            */
    const unsigned short **cf_se02,            /**< out:    | CF tables for t == 0 and f >= 2      */
    const short          **cf_off_se02,        /**< out:    | offsets for CF tables above          */
    const unsigned short **cf_se10,            /**< out:    | CF table for t == 1 and f == 0       */
    short                 *cf_off_se10,        /**< out:    | offset for CF table above            */
    const unsigned short **cf_se11,            /**< out:    | CF tables for t == 1 and f >= 1      */
    const short          **cf_off_se11         /**< out:    | offsets for CF tables above          */
)
{
    int retValue;
    short bitRateIndex;

    retValue     = 0;     /* bitrate index is unknown -> error! */
    bitRateIndex = IGF_MapBitRateToIndex( bitRate, mode, rf_mode );

    if(bitRateIndex != IGF_BITRATE_UNKNOWN)
    {
        retValue = 1; /* no error */

        switch(bitRateIndex)
        {
        case IGF_BITRATE_WB_9600:
        case IGF_BITRATE_RF_WB_13200:
        case IGF_BITRATE_SWB_9600:
        case IGF_BITRATE_SWB_13200:
        case IGF_BITRATE_RF_SWB_13200:
        case IGF_BITRATE_SWB_16400:
        case IGF_BITRATE_SWB_24400:
        case IGF_BITRATE_SWB_32000:
        case IGF_BITRATE_SWB_48000:
            *cf_se00      = cf_se00_tab;
            *cf_se01      = cf_se01_tab[bitRateIndex];
            *cf_off_se01  = cf_off_se01_tab[bitRateIndex];
            *cf_se02      = &cf_se02_tab[bitRateIndex][0][0];
            *cf_off_se02  = &cf_off_se02_tab[bitRateIndex][0];
            *cf_se10      = &cf_se10_tab[0];
            *cf_off_se10  = cf_off_se10_tab;
            *cf_se11      = &cf_se11_tab[0][0][0];
            *cf_off_se11  = &cf_off_se11_tab[0][0];
            break;
        case IGF_BITRATE_FB_16400:
        case IGF_BITRATE_FB_24400:
        case IGF_BITRATE_FB_32000:
            bitRateIndex  = bitRateIndex-IGF_BITRATE_FB_16400+IGF_BITRATE_SWB_16400;
            *cf_se00      = cf_se00_tab;
            *cf_se01      = cf_se01_tab[bitRateIndex];
            *cf_off_se01  = cf_off_se01_tab[bitRateIndex];
            *cf_se02      = &cf_se02_tab[bitRateIndex][0][0];
            *cf_off_se02  = &cf_off_se02_tab[bitRateIndex][0];
            *cf_se10      = &cf_se10_tab[0];
            *cf_off_se10  = cf_off_se10_tab;
            *cf_se11      = &cf_se11_tab[0][0][0];
            *cf_off_se11  = &cf_off_se11_tab[0][0];
            break;
        case IGF_BITRATE_FB_48000:
            bitRateIndex  = bitRateIndex-IGF_BITRATE_FB_48000+IGF_BITRATE_SWB_48000;
            *cf_se00      = cf_se00_tab;
            *cf_se01      = cf_se01_tab[bitRateIndex];
            *cf_off_se01  = cf_off_se01_tab[bitRateIndex];
            *cf_se02      = &cf_se02_tab[bitRateIndex][0][0];
            *cf_off_se02  = &cf_off_se02_tab[bitRateIndex][0];
            *cf_se10      = &cf_se10_tab[0];
            *cf_off_se10  = cf_off_se10_tab;
            *cf_se11      = &cf_se11_tab[0][0][0];
            *cf_off_se11  = &cf_off_se11_tab[0][0];
            break;
        case IGF_BITRATE_FB_96000:
        case IGF_BITRATE_FB_128000:
            bitRateIndex  = IGF_BITRATE_SWB_48000;
            *cf_se00      = cf_se00_tab;
            *cf_se01      = cf_se01_tab[bitRateIndex];
            *cf_off_se01  = cf_off_se01_tab[bitRateIndex];
            *cf_se02      = &cf_se02_tab[bitRateIndex][0][0];
            *cf_off_se02  = &cf_off_se02_tab[bitRateIndex][0];
            *cf_se10      = &cf_se10_tab[0];
            *cf_off_se10  = cf_off_se10_tab;
            *cf_se11      = &cf_se11_tab[0][0][0];
            *cf_off_se11  = &cf_off_se11_tab[0][0];
            break;
        case IGF_BITRATE_UNKNOWN:
        default:
            assert(0);
        }
    }

    return retValue;
}
