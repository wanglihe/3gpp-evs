/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "rom_com.h"
#include "prot.h"
#include "stat_dec.h"

/*-------------------------------------------------------------------
 * Local constants
 *-------------------------------------------------------------------*/

#define LOBUF_NO_SMOOTHING_MODE             1

#define EPS                                 (1e-12f)
#define ENV_SCALE_OFFSET_1                  90.309f /* 10*log10(2^30) */
#define MAX_TEC_BW_LO                       (12)
#define MAX_NB_TEC_LOW_BAND                 (3)

const float TecSC[] = {0.3662f,    0.1078f,    0.1194f,    0.1289f,    0.1365f,    0.1412f};
const int TecLowBandTable[] = {0, 2, 4, 6};

#define TecSmoothingDeg                     5
#define NbTecLowBand                        3

#define ratioHiLoFac                        1.5894f
#define thRatio                             0.3649f
#define thRatio2                            2.0288f
#define thCorrCoef                          0.8795f
#define ratioHiLoFacDec                     (1.5894f * 0.75f)

/*-------------------------------------------------------------------
 * calcLoBufferDec()
 *
 *
 *-------------------------------------------------------------------*/

static void calcLoBufferDec(
    float **pCldfbReal,
    float **pCldfbImag,
    float* loBuffer,
    int startPos,
    int stopPos,
    int bandOffsetBottom,
    float scaleOffset
)
{
    int slot, k, lb, li, ui;
    float nrg, tmp;

    /* calc loTempEnv */
    for (slot = startPos; slot < stopPos; slot++)
    {
        tmp = 0;

        for (lb = 0; lb < NbTecLowBand; lb++)
        {
            li = TecLowBandTable[lb];
            ui = TecLowBandTable[lb+1];

            nrg = 0;
            for (k = li; k < ui; k++)
            {
                nrg += pCldfbReal[slot][k + bandOffsetBottom] * pCldfbReal[slot][k + bandOffsetBottom]
                       + pCldfbImag[slot][k + bandOffsetBottom] * pCldfbImag[slot][k + bandOffsetBottom];
            }
            tmp += (float) log10(nrg * normReciprocal[ui - li] + EPS);
        }

        loBuffer[slot] = 10.0f * tmp / NbTecLowBand + scaleOffset;
    }

    return;
}


/*-------------------------------------------------------------------
 * calcLoBufferEnc()
 *
 *
 *-------------------------------------------------------------------*/

static void calcLoBufferEnc(
    float** pCldfbPow,
    int startPos,
    int stopPos,
    int bandOffsetBottom,
    float scaleOffset,
    float* loBuffer
)
{
    int slot, k, lb, li, ui;
    float nrg, tmp;

    /* calc loTempEnv */
    for (slot = startPos; slot < stopPos; slot++)
    {
        tmp = 0;

        for (lb = 0; lb < NbTecLowBand; lb++)
        {
            li = TecLowBandTable[lb];
            ui = TecLowBandTable[lb+1];

            nrg = 0;
            for (k = li; k < ui; k++)
            {
                nrg += pCldfbPow[slot][k + bandOffsetBottom];
            }
            tmp += (float) log10(nrg * normReciprocal[ui - li] + EPS);
        }

        loBuffer[slot] = 10.0f * tmp / NbTecLowBand + scaleOffset;
    }

    return;
}


/*-------------------------------------------------------------------
 * calcLoTempEnv()
 *
 *
 *-------------------------------------------------------------------*/

static void calcLoTempEnv(
    float* loBuffer,
    int noCols,
    float* loTempEnv,
    float adjFac
)
{
    int slot, i;

    for (slot = 0; slot < noCols; slot++)
    {
        loTempEnv[slot] = TecSC[0] * loBuffer[slot];
        for (i = 1; i < TecSmoothingDeg + 1; i++)
        {
            loTempEnv[slot] += TecSC[i] * loBuffer[slot - i];
        }

        loTempEnv[slot] *= adjFac;
    }

    return;
}

/*-------------------------------------------------------------------
 * calcHiTempEnv()
 *
 *
 *-------------------------------------------------------------------*/

static void calcHiTempEnv(
    float** pCldfbPow,
    int startPos,
    int stopPos,
    int lowSubband,
    int highSubband,
    float scaleOffset,
    float* hiTempEnv
)
{
    int timeIndex, k;
    int bwHigh = highSubband - lowSubband;

    /* set hiTempEnv */
    for (timeIndex = startPos; timeIndex < stopPos; timeIndex++)
    {
        hiTempEnv[timeIndex] = 0;
        for (k = lowSubband; k < highSubband; k++)
        {
            hiTempEnv[timeIndex] += pCldfbPow[timeIndex][k];
        }
        hiTempEnv[timeIndex] = (float)(10 * log10(hiTempEnv[timeIndex] / bwHigh + EPS) + scaleOffset);
    }

    return;
}


/*-------------------------------------------------------------------
 * calcVar()
 *
 *
 *-------------------------------------------------------------------*/

static float calcVar(
    const float in[],
    int len,
    float* x
)
{
    float xx;
    int i;

    xx = 0;
    *x = 0;

    for(i=0; i<len; i++)
    {
        xx += in[i] * in[i];
        *x += in[i];
    }

    return (xx-(*x **x)/(float)len);
}


/*-------------------------------------------------------------------
 * calcCorrelationCoefficient2()
 *
 *
 *-------------------------------------------------------------------*/

static float calcCorrelationCoefficient2(
    const float in_vec1[],
    const float in_vec2[],
    int len,
    float var_x,
    float var_y,
    float x,
    float y
)
{

    int i;
    float xy;
    float ans;

    xy = 0;

    for(i = 0; i < len; i++)
    {
        xy += (in_vec1[i]*in_vec2[i]);
    }

    ans = (xy-x*y/(float)(len))/sqrt(var_x * var_y + EPS);

    return ans;
}


/*-------------------------------------------------------------------
 * resetTecDec()
 *
 *
 *-------------------------------------------------------------------*/

void resetTecDec(
    HANDLE_TEC_DEC hTecDec
)
{
    set_f(hTecDec->pGainTemp, 0,  CLDFB_NO_COL_MAX);

    set_f(hTecDec->loBuffer, 0.f, CLDFB_NO_COL_MAX + MAX_TEC_SMOOTHING_DEG);

    return;
}


/*-------------------------------------------------------------------
 * calcLoTempEnv_TBE()
 *
 *
 *-------------------------------------------------------------------*/

static void calcLoTempEnv_TBE(
    float* loBuffer,
    int noCols,
    float* loTempEnv,
    float adjFac
)
{
    int slot, i;
    int delay = 1;

    for (slot = 0; slot < noCols; slot++)
    {
        loTempEnv[slot] = TecSC[0] * loBuffer[slot - delay];

        for (i = 1; i < TecSmoothingDeg + 1; i++)
        {
            loTempEnv[slot] += TecSC[i] * loBuffer[slot - delay - i];
        }

        loTempEnv[slot] *= adjFac;
    }

    return;
}


/*-------------------------------------------------------------------
 * set_TEC_TFA_code()
 *
 *
 *-------------------------------------------------------------------*/

void set_TEC_TFA_code(
    const short corrFlag,
    short* tec_flag,
    short* tfa_flag
)
{
    *tec_flag = 0;
    if (*tfa_flag == 0)
    {
        if (corrFlag == 1)
        {
            *tec_flag = 1;
        }
        else if (corrFlag == 2)
        {
            *tec_flag = 1;
            *tfa_flag = 1;
        }
    }

    return;
}


/*-------------------------------------------------------------------
 * calcLoTempEnv_ns_TBE()
 *
 *
 *-------------------------------------------------------------------*/

static void calcLoTempEnv_ns_TBE(
    float* loBuffer,
    int noCols,
    float* loTempEnv
)
{
    int slot;
    int delay = 1;
    float fac = 1.4f;

    for (slot = 0; slot < noCols; slot++)
    {
        loTempEnv[slot] = fac * loBuffer[slot - delay];
    }

    return;
}


/*-------------------------------------------------------------------
 * calcLoTempEnv_ns()
 *
 *
 *-------------------------------------------------------------------*/

static void calcLoTempEnv_ns(
    float* loBuffer,
    int noCols,
    float* loTempEnv
)
{
    int slot;

    for (slot = 0; slot < noCols; slot++)
    {
        loTempEnv[slot] = loBuffer[slot];
    }

    return;
}

/*-------------------------------------------------------------------
 * calcGainLinear_TBE()
 *
 *
 *-------------------------------------------------------------------*/

static void calcGainLinear_TBE(
    const float* loTempEnv,
    const int startPos,
    const int stopPos,
    float* pGainTemp
)
{
    int timeIndex;
    float ftmp;

    for (timeIndex = startPos; timeIndex < stopPos; timeIndex++)
    {
        ftmp = loTempEnv[timeIndex];
        pGainTemp[timeIndex] = (float) pow(10.0, 0.1 * ftmp);
    }

    return;
}


/*-------------------------------------------------------------------
 * calcGainTemp_TBE()
 *
 *
 *-------------------------------------------------------------------*/

void calcGainTemp_TBE(
    float** pCldfbRealSrc,
    float** pCldfbImagSrc,
    float* loBuffer,
    int startPos,          /*!<  Start position of the current envelope. */
    int stopPos,           /*!<  Stop position of the current envelope. */
    int lowSubband,        /* lowSubband */
    float* pGainTemp,
    short code
)
{
    float loTempEnv[16];
    const int BW_LO = TecLowBandTable[NbTecLowBand];
    int slot;
    /*int lowSubband = pFreqBandTable[0];*/
    int noCols = stopPos - startPos;
    int bandOffset = 0;
    float scaleOffset = 0;

    assert(lowSubband >= BW_LO);


    bandOffset = lowSubband - BW_LO;

    calcLoBufferDec( pCldfbRealSrc, pCldfbImagSrc, loBuffer + MAX_TEC_SMOOTHING_DEG, startPos, stopPos, bandOffset, scaleOffset );

    if (code > 0)
    {
        if (code != 2)
        {
            calcLoTempEnv_TBE(loBuffer + MAX_TEC_SMOOTHING_DEG, noCols, loTempEnv, ratioHiLoFacDec);
        }
        else
        {
            calcLoTempEnv_ns_TBE(loBuffer + MAX_TEC_SMOOTHING_DEG, noCols, loTempEnv);
        }

        calcGainLinear_TBE(loTempEnv, startPos, stopPos, pGainTemp);
    }


    for (slot = 0; slot < MAX_TEC_SMOOTHING_DEG; slot++)
    {
        loBuffer[slot] = loBuffer[slot + stopPos];
    }

    return;
}


/*-------------------------------------------------------------------
 * setSubfrConfig()
 *
 *
 *-------------------------------------------------------------------*/

static void setSubfrConfig(
    int i_offset,
    int* k_offset,
    int* n_subfr,
    int l_subfr
)
{
    *n_subfr = N_TEC_TFA_SUBFR - i_offset;
    *k_offset = i_offset * l_subfr;

    return;
}


/*-------------------------------------------------------------------
 * calcSubfrNrg()
 *
 *
 *-------------------------------------------------------------------*/

static float calcSubfrNrg(
    const float* hb_synth,
    const int i_offset,
    float* enr,
    int k_offset,
    int l_subfr
)
{
    int i, j, k;
    float enr_all = 1e-12f;

    for (i=i_offset, k = k_offset; i<N_TEC_TFA_SUBFR; i++)
    {
        enr[i] = 1e-12f;
        for(j=0; j<l_subfr; j++)
        {
            enr[i] += hb_synth[k] * hb_synth[k];
            k++;
        }
        enr_all += enr[i];
    }

    return enr_all;
}


/*-------------------------------------------------------------------
 * procTfa()
 *
 *
 *-------------------------------------------------------------------*/

static void procTfa(
    float *hb_synth,
    int i_offset,
    int l_subfr
)
{
    int i,j,k;
    int k_offset, n_subfr;
    float enr[N_TEC_TFA_SUBFR], enr_ave;

    setSubfrConfig(i_offset, &k_offset, &n_subfr, l_subfr);

    enr_ave = calcSubfrNrg(hb_synth, i_offset, enr, k_offset, l_subfr) / n_subfr;

    for (i=i_offset, k=k_offset; i<N_TEC_TFA_SUBFR; i++)
    {
        float gain;

        gain = sqrt(enr_ave / enr[i]);

        for(j=0; j<l_subfr; j++)
        {
            hb_synth[k] *= gain;
            k++;
        }
    }

    return;
}


/*-------------------------------------------------------------------
 * procTec()
 *
 *
 *-------------------------------------------------------------------*/

static void procTec(
    float *hb_synth,
    int i_offset,
    int l_subfr,
    float* gain,
    short code
)
{
    int i,j,k;
    int k_offset, n_subfr;
    float enr[N_TEC_TFA_SUBFR], enr_ave, gain_ave;
    float inv_curr_enr[N_TEC_TFA_SUBFR];
    float max_inv_curr_enr, lower_limit_gain;
    float upper_limit_gain;
    float tmp;

    setSubfrConfig(i_offset, &k_offset, &n_subfr, l_subfr );

    enr_ave = calcSubfrNrg(hb_synth, i_offset, enr, k_offset, l_subfr ) / n_subfr;

    gain_ave = ( sum_f(&gain[i_offset], n_subfr) + 1.0e-12f) / n_subfr;

    max_inv_curr_enr = 10e-6f;;
    for (i=i_offset, k=k_offset; i<N_TEC_TFA_SUBFR; i++)
    {
        inv_curr_enr[i] = enr_ave / enr[i];
        gain[i] /= gain_ave;
        if (max_inv_curr_enr < inv_curr_enr[i])
        {
            max_inv_curr_enr = inv_curr_enr[i];
        }
    }

    lower_limit_gain = 0.1f;
    if (lower_limit_gain > (tmp = 1.0 / max_inv_curr_enr))
    {
        lower_limit_gain = tmp * 0.5;
    }

    upper_limit_gain = 1.2f;
    if (code == LOBUF_NO_SMOOTHING_MODE)
    {
        upper_limit_gain = 3.0;
    }

    for (i=i_offset, k=k_offset; i<N_TEC_TFA_SUBFR; i++)
    {
        if (lower_limit_gain > gain[i])
        {
            gain[i] = lower_limit_gain;
        }
        gain[i] *= inv_curr_enr[i];

        if (gain[i] > upper_limit_gain)
        {
            gain[i] = upper_limit_gain;
        }

        gain[i] = sqrt(gain[i]);
        for(j=0; j<l_subfr; j++)
        {
            hb_synth[k] *= gain[i];
            k++;
        }
    }

    return;
}


/*-------------------------------------------------------------------
 * procTecTfa_TBE()
 *
 *
 *-------------------------------------------------------------------*/

void procTecTfa_TBE(
    float *hb_synth,
    float *gain,
    short flat_flag,
    short last_core,
    int l_subfr,
    short code
)
{
    int i_offset = 0;

    if (flat_flag)
    {
        procTfa(hb_synth, i_offset, l_subfr);
    }
    else
    {
        if (last_core != ACELP_CORE)
        {
            i_offset = 1;
        }

        procTec(hb_synth, i_offset, l_subfr, gain, code);
    }

    return;
}



/*-------------------------------------------------------------------
 * resetTecEnc()
 *
 *
 *-------------------------------------------------------------------*/

void resetTecEnc(
    HANDLE_TEC_ENC hTecEnc,
    int flag
)
{
    if(flag == 0)
    {
        set_f(hTecEnc->loBuffer, 0.f, CLDFB_NO_COL_MAX + MAX_TEC_SMOOTHING_DEG + DELAY_TEMP_ENV_BUFF_TEC);
        set_f(hTecEnc->hiTempEnv, 0.f, CLDFB_NO_COL_MAX + DELAY_TEMP_ENV_BUFF_TEC + EXT_DELAY_HI_TEMP_ENV);
        set_f(hTecEnc->loTempEnv, 0.f, CLDFB_NO_COL_MAX);
        set_f(hTecEnc->loTempEnv_ns, 0.f, CLDFB_NO_COL_MAX);
    }
    else
    {
        set_f(hTecEnc->loBuffer, 0.f, MAX_TEC_SMOOTHING_DEG);
    }

    return;
}


/*-------------------------------------------------------------------
 * calcHiEnvLoBuff()
 *
 *
 *-------------------------------------------------------------------*/

void calcHiEnvLoBuff(
    int noCols,
    const int* pFreqBandTable,  /*!<  freqbandTable. */
    int nSfb,                   /*!<  Number of scalefactors. */
    float** pCldfbPow,
    float* loBuffer,
    float* hiTempEnvOrig
)
{
    const int BW_LO = TecLowBandTable[NbTecLowBand];
    const int lowSubband = pFreqBandTable[0];
    const int highSubband = pFreqBandTable[nSfb];

    int bandOffsetBottom = lowSubband - BW_LO;

    float scaleOffset = 0;
    float* hiTempEnv = hiTempEnvOrig + EXT_DELAY_HI_TEMP_ENV;

    assert(bandOffsetBottom > 0);

    /* calc hiTempEnv*/
    calcHiTempEnv( pCldfbPow, 0, noCols, lowSubband, highSubband, scaleOffset, hiTempEnv + DELAY_TEMP_ENV_BUFF_TEC );

    /* calc loBuffer */
    calcLoBufferEnc( pCldfbPow, 0, noCols, bandOffsetBottom, scaleOffset, loBuffer + MAX_TEC_SMOOTHING_DEG + DELAY_TEMP_ENV_BUFF_TEC );

    return;
}


/*-------------------------------------------------------------------
 * calcLoEnvCheckCorrHiLo()
 *
 *
 *-------------------------------------------------------------------*/

void calcLoEnvCheckCorrHiLo(
    int noCols,
    const int* pFreqBandTable,        /*!<  freqbandTable. */
    float* loBuffer,
    float* loTempEnv,
    float* loTempEnv_ns,
    float* hiTempEnvOrig,
    int* corrFlag
)
{
    const int BW_LO = TecLowBandTable[NbTecLowBand];
    const int lowSubband = pFreqBandTable[0];
    int i;
    float corrCoef;
    int bandOffsetBottom = lowSubband - BW_LO;
    float ratio;
    float hiVar, loVar;
    float hiSum, loSum;
    short code = 0; /* SET TENTATIVELY */
    float loVar_ns;
    float loSum_ns;
    float diff_hi_lo_sum;

    float* hiTempEnv = hiTempEnvOrig + EXT_DELAY_HI_TEMP_ENV;

    assert(bandOffsetBottom > 0);

    hiVar = calcVar(hiTempEnv, noCols, &hiSum);

    /* calc loTempEnv */
    calcLoTempEnv( loBuffer + MAX_TEC_SMOOTHING_DEG, noCols, loTempEnv, ratioHiLoFac );

    calcLoTempEnv_ns( loBuffer + MAX_TEC_SMOOTHING_DEG, noCols, loTempEnv_ns );

    loVar_ns = calcVar(loTempEnv_ns, noCols, &loSum_ns);
    diff_hi_lo_sum = loSum_ns - hiSum;

    if (hiVar > 800 && loVar_ns > 720 && diff_hi_lo_sum < 100)
    {
        code = 1;
    }

    *corrFlag = 0;

    assert(code == 0 || code == 1);

    if (code)
    {
        /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        /* ++++                     code == 1                      +++++*/
        /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        int maxPosHi, maxPosLo;
        float maxHi, maxLo;

        maxHi = hiTempEnv[0];
        maxLo = loTempEnv_ns[0];
        maxPosHi = maxPosLo = 0;

        for (i = 1; i < noCols; i++)
        {
            if (maxHi < hiTempEnv[i])
            {
                maxHi = hiTempEnv[i];
                maxPosHi = i;
            }

            if (maxLo < loTempEnv_ns[i])
            {
                maxLo = loTempEnv_ns[i];
                maxPosLo = i;
            }
        }

        if (abs(maxPosHi - maxPosLo) < 2)
        {
            *corrFlag = 2;
        }

        {
            float feature_max = 0;
            int pos_feature_max = 0;
            float feature[16];

            float min_local, max_local;
            int j;
            int len_window = EXT_DELAY_HI_TEMP_ENV + 1;
            float* curr_pos = hiTempEnv;

            feature_max = 0;
            pos_feature_max = 0;

            for (i = 0; i < 16; i++, curr_pos++)
            {
                max_local = min_local = curr_pos[0];

                for (j = 1; j < len_window; j++)
                {
                    if (max_local < curr_pos[-j])
                    {
                        max_local = curr_pos[-j];
                    }

                    if (min_local > curr_pos[-j])
                    {
                        min_local = curr_pos[-j];
                    }
                }

                feature[i] = max_local - min_local;

                if (feature_max < feature[i])
                {
                    feature_max = feature[i];
                    pos_feature_max = i;
                }
            }

            if (*corrFlag > 0)
            {
                if (!(feature_max > 20 && abs(pos_feature_max - maxPosHi) < 3 ))
                {
                    *corrFlag = 0;
                }
            }
        }
    }
    else
    {
        /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        /* ++++                     code == 0                       ++++*/
        /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

        /* calc the variance of loTempEnv */
        loVar = calcVar(loTempEnv, noCols, &loSum);

        /* calc correlation coefficient between loTempEnv and hiTempEnv */
        corrCoef = calcCorrelationCoefficient2( hiTempEnv, loTempEnv, noCols, hiVar, loVar, hiSum, loSum );

        ratio = hiVar / (loVar + EPS);

        if (corrCoef >=  thCorrCoef && ratio > thRatio && ratio < thRatio2)
        {
            *corrFlag = 1;
        }

    }

    for (i = 0; i < MAX_TEC_SMOOTHING_DEG + DELAY_TEMP_ENV_BUFF_TEC; i++)
    {
        loBuffer[i] = loBuffer[noCols + i];
    }

    for (i = 0; i < DELAY_TEMP_ENV_BUFF_TEC + EXT_DELAY_HI_TEMP_ENV; i++)
    {
        hiTempEnvOrig[i] = hiTempEnvOrig[noCols + i];
    }

    return;
}


/*-------------------------------------------------------------------
 * tecEnc_TBE()
 *
 *
 *-------------------------------------------------------------------*/

void tecEnc_TBE(
    int* corrFlag,
    const float *voicing,
    short coder_type
)
{
    float voice_sum;
    float voice_diff;

    /*-----------------------------------------------------------------*
     * TEC updates
     *-----------------------------------------------------------------*/

    voice_sum = voicing[0] + voicing[1];
    voice_diff = voicing[0] - voicing[1];

    if( voice_diff < 0 )
    {
        voice_diff *= -1.0f;
    }

    if( *corrFlag == 1 )
    {
        if( coder_type == INACTIVE || (((voice_sum > 0.35 * 2 && voice_sum < 0.55 * 2) && (voice_diff < 0.2))) )
        {
            *corrFlag = 0;
        }
    }
    if( voice_sum > 0.6 * 2 )
    {
        *corrFlag = 0;
    }

    return;
}
