/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <float.h>
#include "options.h"
#include "prot.h"
#include "cnst.h"
#include "stat_dec.h"

/**********************************************************************/ /*
measures TCX noise
**************************************************************************/
static int IGF_replaceTCXNoise_1(                                                           /**< out:    | number of noise bands      */
    const float                           *in,                 /**< in:     | MDCT spectrum              */
    const unsigned char                   *TCXNoise,           /**< in:     | tcx noise indicator vector */
    const int                              start,              /**< in:     | start MDCT subband index   */
    const int                              stop,               /**< in:     | stop MDCT subband index    */
    float                                 *totalNoiseNrg       /**< out:    | measured noise energy      */
)
{
    int sb;
    int noise;
    float nE;
    float val;

    noise = 0;
    nE    = FLT_MIN;

    for (sb = start; sb <  stop; sb++)
    {
        if (TCXNoise[sb])
        {
            val = in[sb];
            nE += val * val;
            noise++;
        }
    }

    *totalNoiseNrg = nE;

    return noise;
}

/**********************************************************************/ /*
replaces TCX noise
**************************************************************************/
static void IGF_replaceTCXNoise_2(float                                *in,                 /**< in/out: | MDCT spectrum                */
                                  const unsigned char                  *TCXNoise,           /**< in:     | tcx noise indicator vector   */
                                  const int                             start,              /**< in:     | start MDCT subband index     */
                                  const int                             stop,               /**< in:     | stop MDCT subband index      */
                                  float                                 totalNoiseNrg,      /**< in:     | measured noise energy        */
                                  short                                *nfSeed              /**< in:     | random generator noise seed  */
                                 )
{
    int sb;
    float rE;
    float g;
    float val;

    rE = FLT_MIN;
    for (sb = start; sb <  stop; sb ++)
    {
        if (TCXNoise[sb])
        {
            val    = (float)own_random(nfSeed);
            in[sb] = val;
            rE    += val * val;
        }
    }

    g = (float)sqrt(totalNoiseNrg/rE);

    for (sb = start; sb <  stop; sb ++)
    {
        if (TCXNoise[sb])
        {
            in[sb] *= g;
        }
    }
}

/**********************************************************************/ /*
reads whitening levels
**************************************************************************/
static void IGF_decode_whitening_level(Decoder_State                   *st,                 /**< in:     | decoder state                    */
                                       IGF_DEC_PRIVATE_DATA_HANDLE      hPrivateData,       /**< in:     | instance handle of IGF Deccoder  */
                                       const int                        p                   /**< in:     | tile index, p = [0, 3]           */
                                      )
{
    int tmp;

    tmp = get_next_indice(st, 1);

    if (tmp == 1)
    {
        tmp = get_next_indice(st, 1);
        if (tmp == 1)
        {
            hPrivateData->currWhiteningLevel[p] = IGF_WHITENING_STRONG;
        }
        else
        {
            hPrivateData->currWhiteningLevel[p] = IGF_WHITENING_OFF;
        }
    }
    else
    {
        hPrivateData->currWhiteningLevel[p]   = IGF_WHITENING_MID;
    }
}

/**********************************************************************/ /*
reads flattening trigger
**************************************************************************/
static void IGF_decode_temp_flattening_trigger(Decoder_State           *st,                 /**< in:     | decoder state                   */
        IGF_DEC_INSTANCE_HANDLE  hInstance           /**< in:     | instance handle of IGF Deccoder */
                                              )
{
    hInstance->flatteningTrigger = get_next_indice(st, 1);
}

/**********************************************************************/ /*
square the MDCT spectrum
**************************************************************************/
static void IGF_getMDCTSquare(const int                                 startLine,          /**< in:     | start MDCT subband index       */
                              const int                                 stopLine,           /**< in:     | stop  MDCT subband index       */
                              const float                              *pSpectralData,      /**< in:     | MDCT spectrum                  */
                              float                                    *pSpecDataSqaure     /**< out:    | Squared MDCT spectrum          */
                             )
{
    int i;

    for (i = startLine; i < stopLine; i++)
    {
        pSpecDataSqaure[i] = pSpectralData[i] * pSpectralData[i];
    }
}

/**********************************************************************/ /*
calculate energy per SFB
**************************************************************************/
static void IGF_calcSfbEnergy(const int                                 startSfb,           /**< in:     | start sfb index                                        */
                              const int                                 stopSfb,            /**< in:     | stop  sfb index                                        */
                              const int                                *swb_offset,         /**< in:     | IGF swb offset table                                   */
                              const float                              *pPowerSpectrum,     /**< in:     | power spectrum                                         */
                              float                                    *sfbEnergy           /**< out:    | SFB energies, will be initialized inside this function */
                             )
{
    int sfb;
    int line;

    for (sfb = startSfb; sfb < stopSfb; sfb++)
    {
        sfbEnergy[sfb] = 0.f;

        for (line = swb_offset[sfb]; line < swb_offset[sfb + 1]; line++)
        {
            sfbEnergy[sfb] += pPowerSpectrum[line];
        }
    }
}

/**********************************************************************/ /*
set power spectrum values to zero, needed for energy calculation
**************************************************************************/
static void IGF_setLinesToZero(const int                                startLine,          /**< in:     | start MDCT subband index       */
                               const int                                stopLine,           /**< in:     | stop  MDCT subband index       */
                               const float                             *pSpectralData,      /**< in:     | original MDCT spectrum         */
                               float                                   *squareSpecIGF       /**< in/out: | prepared IGF energy spectrum   */
                              )
{
    int i;

    for (i = startLine; i < stopLine; i++)
    {
        if (pSpectralData[i] != 0.f)
        {
            squareSpecIGF[i] = 0.f;
        }
    }
}

/**********************************************************************/ /*
prepare IGF spectrum
**************************************************************************/
static void IGF_prep(IGF_DEC_PRIVATE_DATA_HANDLE                        hPrivateData,       /**< in:     | IGF private data handle                              */
                     const int                                          igfGridIdx,         /**< in:     | in case of CELP->TCX switching, use 1.25 framelength */
                     const unsigned char                               *TCXNoise,           /**< in:     | TCX noise vector                                     */
                     float                                             *igf_spec,           /**< in:     | prepared IGF spectrum                                */
                     float                                             *src_spec            /**< in:     | source spectrum                                      */
                    )
{
    H_IGF_GRID hGrid;
    H_IGF_INFO hInfo;
    int i;
    int tb;
    int sfb;
    int strt_cpy;
    int tile_idx;
    int *swb_offset;
    float *sel_spec;

    hInfo       = &hPrivateData->igfInfo;
    hGrid       = &hPrivateData->igfInfo.grid[igfGridIdx];
    swb_offset  = hGrid->swb_offset;

    for (tile_idx = 0; tile_idx < hGrid->nTiles; tile_idx ++)
    {
        strt_cpy = hGrid->sbWrap[tile_idx];

        if (IGF_WHITENING_STRONG == hPrivateData->currWhiteningLevel[tile_idx])
        {
            float abs_sum;
            abs_sum = 0.f;
            for(i = strt_cpy; i < hGrid->startLine; i++)
            {
                abs_sum += (float) fabs(src_spec[i]);
            }

            tb = swb_offset[hGrid->sfbWrap[tile_idx]];

            if( abs_sum > 0.f )
            {
                for (i = strt_cpy; i < hGrid->startLine; i++)
                {
                    igf_spec[tb++] = own_random(&hInfo->nfSeed);
                }
            }
            else
            {
                for (i = strt_cpy; i < hGrid->startLine; i++)
                {
                    igf_spec[tb++] = 0.f;
                }
            }
        }
        else
        {
            if (IGF_WHITENING_MID == hPrivateData->currWhiteningLevel[tile_idx])
            {
                if (hPrivateData->n_noise_bands)
                {
                    IGF_replaceTCXNoise_2(igf_spec,
                                          TCXNoise,
                                          hGrid->minSrcSubband,
                                          hGrid->startLine,
                                          hPrivateData->totalNoiseNrg,
                                          &hInfo->nfSeed);
                }
                sel_spec = igf_spec;
            }
            else
            {
                if (hPrivateData->n_noise_bands_off)
                {
                    IGF_replaceTCXNoise_2(src_spec,
                                          TCXNoise,
                                          hGrid->minSrcSubband,
                                          hGrid->startLine,
                                          hPrivateData->totalNoiseNrg_off,
                                          &hInfo->nfSeed);
                }
                sel_spec = src_spec;
            }
            for (sfb = hGrid->sfbWrap[tile_idx]; sfb < hGrid->sfbWrap[tile_idx + 1]; sfb++)
            {
                for (tb = swb_offset[sfb]; tb < swb_offset[sfb + 1]; tb++)
                {
                    igf_spec[tb] = sel_spec[strt_cpy];
                    strt_cpy++;
                }
            }
        }
    }
}

/**********************************************************************/ /*
calculates IGF energies
**************************************************************************/
static void IGF_calc(const IGF_DEC_PRIVATE_DATA_HANDLE                  hPrivateData,       /**< in:     | IGF private data handle                              */
                     const int                                          igfGridIdx,         /**< in:     | in case of CELP->TCX switching, use 1.25 framelength */
                     const float                                       *spectrum,           /**< in:     | MDCT spectrum                                        */
                     float                                             *igf_spec            /**< in:     | prepared IGF spectrum                                */
                    )
{
    H_IGF_GRID hGrid;
    float *igf_pN;
    float *igf_sN;
    float tmp[N_MAX_TCX] = {0.f};
    hGrid  = &hPrivateData->igfInfo.grid[igfGridIdx];
    igf_pN = hPrivateData->igf_pN;
    igf_sN = hPrivateData->igf_sN;

    IGF_getMDCTSquare(hGrid->startLine,  hGrid->stopLine, spectrum,          tmp);
    IGF_calcSfbEnergy(hGrid->startSfb,   hGrid->stopSfb,  hGrid->swb_offset, tmp, igf_sN);
    IGF_getMDCTSquare(hGrid->startLine,  hGrid->stopLine, igf_spec,          tmp);
    IGF_setLinesToZero(hGrid->startLine, hGrid->stopLine, spectrum,          tmp);
    IGF_calcSfbEnergy(hGrid->startSfb,   hGrid->stopSfb,  hGrid->swb_offset, tmp, igf_pN);
}

/**********************************************************************/ /*
apply IGF
**************************************************************************/
static void IGF_appl(IGF_DEC_PRIVATE_DATA_HANDLE                        hPrivateData,       /**< in:     | IGF private data handle                              */
                     const int                                          igfGridIdx,         /**< in:     | in case of CELP->TCX switching, use 1.25 framelength */
                     float                                             *pSpectralData,      /**< in: Q31 | MDCT spectrum                                        */
                     const float                                       *igf_spec,           /**< in: Q31 | prepared IGF spectrum                                */
                     float                                             *virtualSpec,        /**< out:Q31 | virtual IGF spectrum, used for temp flattening       */
                     int                                               *flag_sparse         /**< out: Q0 | temp flattening indicator                            */
                    )
{
    H_IGF_GRID hGrid;
    int tb;                                       /* target subband */
    int sfb;
    int s_sfb;
    int start_sfb;
    int stop_sfb;
    int *swb_offset;
    int hopsize;
    float tmp;
    float dE;
    float dN[IGF_MAX_SFB+1];
    float gain[IGF_MAX_SFB];
    float dS[IGF_MAX_SFB];
    float width;
    float sNlocal;
    float E;
    float sum;
    float val;
    float w0;
    float w1;
    float w2;
    float *sN;
    float *pN;
    float gFactor;                      /* general SCF adaption */
    float fFactor;                      /* first SCF adaption   */
    float lFactor;                      /* last SCF adaption    */

    /* initialize variables */
    w0          = 0.201f;
    w1          = 0.389f;
    w2          = 0.410f;
    dE          = 0.f;

    set_i(flag_sparse,   0.f, N_MAX_TCX-IGF_START_MN);
    set_f(virtualSpec, 0.f, N_MAX_TCX-IGF_START_MN);


    /* more inits */
    hGrid       = &hPrivateData->igfInfo.grid[igfGridIdx];
    sN          = hPrivateData->igf_sN;
    pN          = hPrivateData->igf_pN;
    start_sfb   = hGrid->startSfb;
    stop_sfb    = hGrid->stopSfb;
    gFactor     = hGrid->gFactor;
    fFactor     = hGrid->fFactor;
    lFactor     = hGrid->lFactor;
    swb_offset  = hGrid->swb_offset;

    /* collect energy below hGrid->startLine: */
    for (tb = hGrid->startLine-24; tb < hGrid->startLine; tb++)
    {
        dE += pSpectralData[tb] * pSpectralData[tb];
    }
    dE = (float)sqrt(dE/24.);

    hopsize = 2;
    hopsize = (hPrivateData->currWhiteningLevel[0] == IGF_WHITENING_OFF)    ? 4 : hopsize;
    hopsize = (hPrivateData->currWhiteningLevel[0] == IGF_WHITENING_MID)    ? 2 : hopsize;
    hopsize = (hPrivateData->currWhiteningLevel[0] == IGF_WHITENING_STRONG) ? 1 : hopsize;
    hopsize = min(hopsize, hPrivateData->igfInfo.maxHopsize);

    if(hopsize > 1)
    {
        for (sfb = start_sfb; sfb < stop_sfb; sfb += hopsize)
        {
            for (tb = sfb+1; tb < min((sfb+hopsize),stop_sfb); tb++)
            {
                sN[sfb] += sN[tb];
                pN[sfb] += pN[tb];
                sN[tb]   = 0.f;
                pN[tb]   = 0.f;
            }
        }
    }

    /* IGF_rescale_SCF */
    if (hGrid->infoIsRefined)
    {
        for (sfb = start_sfb; sfb < stop_sfb; sfb+=2)
        {
            width = (float)(swb_offset[sfb + 2] - swb_offset[sfb]);

            tmp = (float)hPrivateData->igf_curr[sfb >> 1];
            tmp = (float)pow(2.0, 0.25 * tmp - 4.0);
            tmp = tmp * tmp;

            sNlocal  = sN[sfb] + sN[sfb + 1];
            sNlocal /= width;

            tmp       = max(0.001 * sNlocal, tmp - sNlocal);
            dN[sfb]   = (float)sqrt(tmp);
            dN[sfb+1] = dN[sfb];
        }
    }
    else
    {
        for (sfb = start_sfb; sfb < stop_sfb; sfb++)
        {
            width = (float)(swb_offset[sfb + 1] - swb_offset[sfb]);

            tmp = (float)hPrivateData->igf_curr[sfb];
            tmp = (float)pow(2.0, 0.25 * tmp - 4.0);
            tmp = tmp * tmp;

            sNlocal  = sN[sfb];
            sNlocal /= width;

            tmp     = max(0.001 * sNlocal, tmp - sNlocal);
            dN[sfb] = (float)sqrt(tmp);
        }
    }

    dS[start_sfb] = dN[start_sfb];
    /* first value with adaption to core energy: */
    if (dE < dN[start_sfb])
    {
        dS[start_sfb] = dN[start_sfb] + fFactor * (dE - dN[start_sfb]);
    }
    /* last value with less energy: */
    dS[stop_sfb - 1] = lFactor * dN[stop_sfb - 1];

    if (hGrid->infoIsRefined && hopsize == 1)
    {
        /* apply filter to absolute energy values: */
        for (sfb = start_sfb+1; sfb < stop_sfb-1; sfb++)
        {
            dS[sfb] = w0 * dN[sfb - 1] + w1 * dN[sfb] + w2 * dN[sfb + 1];
        }
    }
    else
    {
        for (sfb = start_sfb + 1; sfb < stop_sfb - 1; sfb++)
        {
            dS[sfb] = dN[sfb];
        }
    }

    for (sfb = start_sfb; sfb < stop_sfb; sfb+=hopsize)
    {
        E   = 0.f;
        sum = 0;
        for (tb = 0; tb < hopsize; tb++)
        {
            width = (float)(swb_offset[min(sfb + tb + 1, stop_sfb)] - swb_offset[min(sfb + tb, stop_sfb)]);
            val   = dS[min(sfb + tb,stop_sfb - 1)];
            E    += val * val * width;
            sum  += width;
        }

        dS[sfb] = (float)sqrt((E * hopsize) / sum);
        dN[sfb] = gFactor * dS[sfb];

        width     = (float)(swb_offset[sfb + 1] - swb_offset[sfb]);
        dN[sfb]   = dN[sfb] * dN[sfb] * width;
        gain[sfb] = 0.f;

        if (pN[sfb] > 1.e-20f)
        {
            gain[sfb] = (float)sqrt(dN[sfb] / pN[sfb]);
        }

        for (s_sfb = sfb + 1; s_sfb < min(sfb + hopsize, stop_sfb); s_sfb++)
        {
            gain[s_sfb] = gain[sfb];
        }
    }

    /* tiling */
    for(sfb = start_sfb; sfb < stop_sfb; sfb++)
    {

        if (hPrivateData->frameLossCounter > 0)
        {
            gain[sfb] = min(gain[sfb], 12.f);

            if (hPrivateData->frameLossCounter < 5)
            {
                gain[sfb] -= gain[sfb] / 8 * hPrivateData->frameLossCounter;
            }
            else
            {
                gain[sfb] /= 2;
            }
        }

        for (tb = swb_offset[ sfb ]; tb < swb_offset[ sfb+1 ]; tb++)
        {
            if(pSpectralData[tb] == 0.f)
            {
                pSpectralData[tb] = igf_spec[tb] * gain[sfb];
                flag_sparse[tb-IGF_START_MN]   = 1;
            }
            else
            {
                virtualSpec[tb-IGF_START_MN]   = igf_spec[tb] * gain[sfb];
                flag_sparse[tb-IGF_START_MN]   = 2;

            }
        }
    }
}

/**********************************************************************/ /*
spectral whitening
**************************************************************************/
static void IGF_getWhiteSpectralData(const float                       *in,                 /**< in:     | MDCT spectrum              */
                                     float                             *out,                /**< out:    | whitened spectrum          */
                                     const int                          start,              /**< in:     | start MDCT subband index   */
                                     const int                          stop,               /**< in:     | stop MDCT subband index    */
                                     const int                          level               /**< in:     | whitening strength         */
                                    )
{
    int i;
    int n;
    int j;
    float div;
    float ak;

    /* inits */
    div = 0;

    for (i = start; i < stop - level; i++)
    {
        ak = 1e-3f;
        for (j = i - level; j < i + level + 1; j++)
        {
            ak += in[j] * in[j];
        }
        ak /= (float)(level * 2 + 1);


        n   = max(0.f, (int)(log(ak) * INV_LOG_2));           /* INV_LOG_2 = 1 / (float)log(2.0f)) */
        n >>= 1;                                              /* sqrt() */
        div = (float)(pow(2.0f, (float)(21 - n)));


        out[i] = in[i] * div;                                 /* same as shift */
    }

    for (; i < stop; i++)
    {
        ak = 1e-3f;

        for (j = i - level; j < stop; j++)
        {
            ak += in[j]*in[j];
        }
        ak /= (float)(stop - (i - level));


        n   = max(0.f,(int)(log(ak) * INV_LOG_2));            /* INV_LOG_2 = 1 / (float)log(2.0f)) */
        n >>= 1;                                              /* sqrt() */
        div = (float)(pow(2.0f, (float)(21 - n)));


        out[i] = in[i] * div;                                 /* same as shift */
    }
}

/**********************************************************************/ /*
refines the IGF grid
**************************************************************************/
static void IGF_RefineGrid(H_IGF_GRID                                   hGrid               /**< in/out: | IGF grid handle  */
                          )
{
    int a[IGF_MAX_SFB+1];           /* +1: because in for-loop one value too much will be extrapolated  */
    int sfb;

    set_i(a, 0, IGF_MAX_SFB+1);


    hGrid->infoIsRefined = 1;
    for (sfb = 0; sfb < hGrid->swb_offset_len; sfb++)
    {
        a[sfb*2 + 0] = hGrid->swb_offset[sfb];
        a[sfb*2 + 1] = round_f(hGrid->swb_offset[sfb] + 0.45f * (hGrid->swb_offset[sfb + 1] - hGrid->swb_offset[sfb]));
        if (a[sfb*2 + 1] & 1)
        {
            a[sfb*2 + 1]--;
        }
    }
    hGrid->stopSfb = hGrid->stopSfb * 2;
    for (sfb = 0; sfb <= hGrid->stopSfb; sfb++)
    {
        hGrid->swb_offset[sfb] = a[sfb];
    }

    for (sfb = 0; sfb <= hGrid->nTiles; sfb++)
    {
        hGrid->sfbWrap[sfb] *= 2;
    }
}

/**********************************************************************/ /*
reads whitening information from the bitstream
**************************************************************************/
void IGFDecReadData(const IGF_DEC_INSTANCE_HANDLE                    hInstance,          /**< in:     | instance handle of IGF Deccoder                      */
                    Decoder_State                                   *st,                 /**< in:     | decoder state                                        */
                    const int                                        igfGridIdx,         /**< in:     | in case of CELP->TCX switching, use 1.25 framelength */
                    const int                                        isIndepFrame        /**< in:     | if 1: arith dec force reset, if 0: no reset          */
                   )
{
    IGF_DEC_PRIVATE_DATA_HANDLE hPrivateData;
    H_IGF_GRID hGrid;
    int p;
    int nT;
    int tmp;

    if (hInstance != NULL)
    {
        hPrivateData = &hInstance->igfData;
        hGrid        = &hPrivateData->igfInfo.grid[igfGridIdx];
        nT           = hGrid->nTiles;
        tmp          = -1;

        for (p = 0; p < IGF_MAX_TILES; p++)
        {
            hPrivateData->currWhiteningLevel[p] = IGF_WHITENING_OFF;
        }

        if (isIndepFrame)
        {
            tmp = 0;
        }
        else
        {
            tmp = get_next_indice(st, 1);
        }
        if (tmp == 1)
        {
            for (p = 0; p < nT; p++)
            {
                hPrivateData->currWhiteningLevel[p] = hPrivateData->prevWhiteningLevel[p];
            }
        }
        else
        {
            IGF_decode_whitening_level(st, hPrivateData, 0);
            tmp = (int)get_next_indice(st, 1);
            if (tmp == 1)
            {
                for (p = 1; p < nT; p++)
                {
                    IGF_decode_whitening_level(st, hPrivateData, p);
                }
            }
            else
            {
                for (p = 1; p < nT; p++)
                {
                    hPrivateData->currWhiteningLevel[p] = hPrivateData->currWhiteningLevel[0];
                }
            }
        }
        for (p = 0; p < IGF_MAX_TILES; p++)
        {
            hPrivateData->prevWhiteningLevel[p] = hPrivateData->currWhiteningLevel[p];
        }
        IGF_decode_temp_flattening_trigger(st, hInstance);
    }
}

/**********************************************************************/ /*
read the IGF level information from the bitsream
**************************************************************************/
int IGFDecReadLevel(                                                                     /**< out:    | return igfAllZero flag indicating if no envelope is transmitted  */
    const IGF_DEC_INSTANCE_HANDLE                    hInstance,          /**< in:     | instance handle of IGF Deccoder                                  */
    Decoder_State                                   *st,                 /**< in:     | decoder state                                                    */
    const int                                        igfGridIdx,         /**< in:     | in case of CELP->TCX switching, use 1.25 framelength             */
    const int                                        isIndepFrame        /**< in:     | if 1: arith dec force reset, if 0: no reset                      */
)
{
    IGF_DEC_PRIVATE_DATA_HANDLE hPrivateData;
    H_IGF_GRID hGrid;
    int m_igfSfbStart;
    int IGFAllZero;

    IGFAllZero = 1;

    if (hInstance != NULL)
    {
        hPrivateData  = &hInstance->igfData;
        hGrid         = &hPrivateData->igfInfo.grid[igfGridIdx];
        m_igfSfbStart = hGrid->startSfb;
        IGFAllZero    = get_next_indice(st, 1);

        if (IGFAllZero == 0)
        {
            mvi2i(hPrivateData->igf_curr, hPrivateData->igf_prev, hGrid->stopSfb);
            IGFSCFDecoderDecode(&hPrivateData->hArithSCFdec,
                                st,
                                &hPrivateData->igf_curr[m_igfSfbStart],    /**< out: ptr to an array which will contain the decoded quantized coefficients   */
                                isIndepFrame);                             /**< in: if 1 on input reset will be forced, if 0 on input continue without reset */
        }
        else
        {
            IGFSCFDecoderReset(&hPrivateData->hArithSCFdec);
            set_i(&hPrivateData->igf_curr[m_igfSfbStart], 0, hGrid->stopSfb - m_igfSfbStart);
        }
    }

    hInstance->infoIGFAllZero = IGFAllZero;
    return IGFAllZero;
}

/**********************************************************************/ /*
apply the IGF decoder
**************************************************************************/
void IGFDecApplyMono(const IGF_DEC_INSTANCE_HANDLE                   hInstance,          /**< in:     | instance handle of IGF Decoder                       */
                     float                                          *spectrum,           /**< in/out: | MDCT spectrum                                        */
                     const int                                       igfGridIdx,         /**< in:     | in case of CELP->TCX switching, use 1.25 framelength */
                     int                                             bfi                 /**< in:     | frame loss == 1, frame good == 0                     */
                    )
{
    IGF_DEC_PRIVATE_DATA_HANDLE hPrivateData;
    H_IGF_GRID hGrid;
    int i;
    int whiteningLevel;
    float igf_spec[IGF_MAX_GRANULE_LEN] = {0.f};

    /* initialize variables */
    whiteningLevel                  = 7;
    hPrivateData                    = &hInstance->igfData;
    hGrid                           = &hPrivateData->igfInfo.grid[igfGridIdx];
    hPrivateData->totalNoiseNrg     = 0.f;
    hPrivateData->n_noise_bands     = 0;
    hPrivateData->totalNoiseNrg_off = 0.f;
    hPrivateData->n_noise_bands_off = 0;

    /* concealment counter */
    if (bfi)
    {
        hPrivateData->frameLossCounter++;
    }
    else
    {
        hPrivateData->frameLossCounter = 0;
    }

    /* skip IGF processing if all IGF levels are zero */
    if (!hInstance->infoIGFAllZero)
    {

        for (i = 0; i < hGrid->nTiles; i++)
        {
            if (hPrivateData->currWhiteningLevel[i] == IGF_WHITENING_MID)
            {
                IGF_getWhiteSpectralData(hPrivateData->pSpecFlat,
                                         igf_spec,
                                         hGrid->minSrcSubband,
                                         hGrid->startLine,
                                         whiteningLevel);

                hPrivateData->n_noise_bands = IGF_replaceTCXNoise_1(igf_spec,
                                              hInstance->infoTCXNoise,
                                              hGrid->minSrcSubband,
                                              hGrid->startLine,
                                              &hPrivateData->totalNoiseNrg);
                break;
            }
        }

        for (i = 0; i < hGrid->nTiles; i++)
        {
            if (hPrivateData->currWhiteningLevel[i] == IGF_WHITENING_OFF)
            {
                hPrivateData->n_noise_bands_off = IGF_replaceTCXNoise_1(hPrivateData->pSpecFlat,
                                                  hInstance->infoTCXNoise,
                                                  hGrid->minSrcSubband,
                                                  hGrid->startLine,
                                                  &hPrivateData->totalNoiseNrg_off);
                break;
            }
        }

        /* apply IGF in three steps: */
        IGF_prep(hPrivateData, igfGridIdx, hInstance->infoTCXNoise, igf_spec, hPrivateData->pSpecFlat);
        IGF_calc(hPrivateData, igfGridIdx, spectrum, igf_spec);
        IGF_appl(hPrivateData, igfGridIdx, spectrum, igf_spec, hInstance->virtualSpec, hInstance->flag_sparse);

    }

    /* reset TCX noise indicator vector */
    set_c((char*)(hInstance->infoTCXNoise), 0, IGF_START_MX);
}

/**********************************************************************/ /*
set mode is used to init the IGF dec with a new bitrate
**************************************************************************/
void IGFDecSetMode(const IGF_DEC_INSTANCE_HANDLE                     hInstance,          /**< in:     | instance handle of IGF Decoder */
                   const int                                         bitRate,            /**< in:     | bitrate                        */
                   const int                                         mode,               /**< in:     | bandwidth mode                 */
                   const int                                         defaultStartLine,   /**< in:     | default start subband index    */
                   const int                                         defaultStopLine     /**< in:     | default stop subband index     */
                   , const int                                         rf_mode             /**< in:     | flag to signal the RF mode */
                  )
{
    IGF_DEC_PRIVATE_DATA_HANDLE hPrivateData;

    hPrivateData           = &hInstance->igfData;
    hInstance->isIGFActive = 0;

    if (IGFCommonFuncsIGFConfiguration(bitRate, mode, &hPrivateData->igfInfo
                                       , rf_mode
                                      ))
    {
        IGFSCFDecoderOpen(&hPrivateData->hArithSCFdec,
                          hPrivateData->igfInfo.grid[0].stopSfb - hPrivateData->igfInfo.grid[0].startSfb,
                          bitRate,
                          mode
                          , rf_mode
                         );

        hInstance->infoIGFStopLine   = hPrivateData->igfInfo.grid[0].stopLine;
        hInstance->infoIGFStartLine  = hPrivateData->igfInfo.grid[0].startLine;
        hInstance->infoIGFStopFreq   = hPrivateData->igfInfo.grid[0].stopFrequency;
        hInstance->infoIGFStartFreq  = hPrivateData->igfInfo.grid[0].startFrequency;
        hInstance->infoIGFAllZero    = 0;
        hInstance->isIGFActive       = 1;

        if (hPrivateData->igfInfo.bitRateIndex <= IGF_BITRATE_SWB_48000 || hPrivateData->igfInfo.bitRateIndex <= IGF_BITRATE_FB_48000)
        {
            IGF_RefineGrid(&hPrivateData->igfInfo.grid[IGF_GRID_LB_NORM]);
            IGF_RefineGrid(&hPrivateData->igfInfo.grid[IGF_GRID_LB_TRAN]);
            IGF_RefineGrid(&hPrivateData->igfInfo.grid[IGF_GRID_LB_SHORT]);
        }
        /* IGFDecOutInformation(hInstance); */
    }
    else
    {
        hInstance->infoIGFStopLine   = defaultStopLine;
        hInstance->infoIGFStartLine  = defaultStartLine;
        hInstance->infoIGFStopFreq   = -1;
        hInstance->infoIGFStartFreq  = -1;
        fprintf(stderr,"IGFDecSetMode: initialization error!\n");
    }
}

/**********************************************************************/ /*
updates the start/stop frequency of IGF according to igfGridIdx
**************************************************************************/
void IGFDecUpdateInfo(const IGF_DEC_INSTANCE_HANDLE                  hInstance,          /**< in:     | instance handle of IGF Decoder */
                      const int                                      igfGridIdx          /**< in:     | IGF grid index                 */
                     )
{
    IGF_DEC_PRIVATE_DATA_HANDLE hPrivateData;
    H_IGF_GRID hGrid;

    hPrivateData = &hInstance->igfData;
    if (hInstance->isIGFActive)
    {
        hGrid                       = &hPrivateData->igfInfo.grid[igfGridIdx];
        hInstance->infoIGFStartFreq = hGrid->startFrequency;
        hInstance->infoIGFStopFreq  = hGrid->stopFrequency;
        hInstance->infoIGFStartLine = hGrid->startLine;
        hInstance->infoIGFStopLine  = hGrid->stopLine;
    }
}

/**********************************************************************/ /*
copy the LPC flat spectrum to IGF buffer
**************************************************************************/
void IGFDecCopyLPCFlatSpectrum(const IGF_DEC_INSTANCE_HANDLE         hInstance,          /**< in:     | instance handle of IGF Decoder     */
                               const float                          *pSpectrumFlat,      /**< in:     | LPC flattend spectrum from TCX dec */
                               const int                             igfGridIdx          /**< in:     | IGF grid index                     */
                              )
{
    IGF_DEC_PRIVATE_DATA_HANDLE hPrivateData;
    H_IGF_GRID hGrid;
    int i;

    if (hInstance)
    {
        hPrivateData = &hInstance->igfData;
        hGrid        = &hPrivateData->igfInfo.grid[igfGridIdx];


        for (i = hGrid->minSrcSubband - 7; i < hGrid->startLine; i++)
        {
            hPrivateData->pSpecFlat[i] = pSpectrumFlat[i] * 1024.f;
        }
    }
}

/**********************************************************************/ /*
store the IGF bitstream information for TCX10 subframes
**************************************************************************/
void IGFDecStoreTCX10SubFrameData(const IGF_DEC_INSTANCE_HANDLE      hInstance,          /**< in:     | instance handle of IGF Decoder */
                                  const int                          subFrameIdx         /**< in:     | index of subframe              */
                                 )
{
    IGF_DEC_PRIVATE_DATA_HANDLE hPrivateData;

    hPrivateData = &hInstance->igfData;

    /* store igf energies for subframe*/
    mvi2i(hPrivateData->igf_curr, hPrivateData->igf_curr_subframe[subFrameIdx][0], IGF_MAX_SFB);
    mvi2i(hPrivateData->igf_prev, hPrivateData->igf_prev_subframe[subFrameIdx], IGF_MAX_SFB);

    /* store spectral whitening information for current subframe */
    mvi2i(hPrivateData->currWhiteningLevel, hPrivateData->currWhiteningLevel_subframe[subFrameIdx], IGF_MAX_TILES);
    mvi2i(hPrivateData->prevWhiteningLevel, hPrivateData->prevWhiteningLevel_subframe[subFrameIdx], IGF_MAX_TILES);
    /* store flattening trigger for current subframe */
    hPrivateData->igf_flatteningTrigger_subframe[subFrameIdx] = hInstance->flatteningTrigger;
}

/**********************************************************************/ /*
restore the IGF bitstream information for TCX10 subframes
**************************************************************************/
void IGFDecRestoreTCX10SubFrameData(const IGF_DEC_INSTANCE_HANDLE    hInstance,          /**< in:     | instance handle of IGF Decoder */
                                    const int                        subFrameIdx         /**< in:     | index of subframe              */
                                   )
{
    IGF_DEC_PRIVATE_DATA_HANDLE hPrivateData;

    hPrivateData = &hInstance->igfData;

    /* store igf energies for subframe*/
    mvi2i(hPrivateData->igf_curr_subframe[subFrameIdx][0], hPrivateData->igf_curr, IGF_MAX_SFB);
    mvi2i(hPrivateData->igf_prev_subframe[subFrameIdx], hPrivateData->igf_prev, IGF_MAX_SFB);

    /* store spectral whitening information for current subframe */
    mvi2i(hPrivateData->currWhiteningLevel_subframe[subFrameIdx], hPrivateData->currWhiteningLevel, IGF_MAX_TILES);
    mvi2i(hPrivateData->prevWhiteningLevel_subframe[subFrameIdx], hPrivateData->prevWhiteningLevel, IGF_MAX_TILES);
    /* restore flattening trigger for current subframe */
    hInstance->flatteningTrigger = hPrivateData->igf_flatteningTrigger_subframe[subFrameIdx];
}

