/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"

/*-------------------------------------------------------------------*
* tcx_get_windows()
*
*
*-------------------------------------------------------------------*/

static void tcx_get_windows(
    TCX_config const * tcx_cfg,     /* i: TCX configuration                         */
    const short left_mode,          /* i: overlap mode of left window half          */
    const short right_mode,         /* i: overlap mode of right window half         */
    int *left_overlap,              /* o: left overlap length                       */
    float const **left_win,         /* o: left overlap window                       */
    int *right_overlap,             /* o: right overlap length                      */
    float const **right_win,        /* o: right overlap window                      */
    int fullband                    /* i: fullband flag                             */
)
{
    if (!fullband)
    {
        /* Left part */

        if (left_mode == TRANSITION_OVERLAP)
        {
            /* ACELP->TCX transition */
            *left_overlap = tcx_cfg->tcx_mdct_window_trans_length;
            *left_win = tcx_cfg->tcx_mdct_window_trans;
        }
        else if (left_mode == MIN_OVERLAP)
        {
            *left_overlap = tcx_cfg->tcx_mdct_window_min_length;
            *left_win = tcx_cfg->tcx_mdct_window_minimum;
        }
        else if (left_mode == HALF_OVERLAP)
        {
            *left_overlap = tcx_cfg->tcx_mdct_window_half_length;
            *left_win = tcx_cfg->tcx_mdct_window_half;
        }
        else if (left_mode == FULL_OVERLAP)
        {
            *left_overlap = tcx_cfg->tcx_mdct_window_length;
            *left_win = tcx_cfg->tcx_aldo_window_1_trunc;
        }
        else
        {
            assert(!"Not supported overlap");
        }

        /* Right part */

        if (right_mode == MIN_OVERLAP)
        {
            *right_overlap = tcx_cfg->tcx_mdct_window_min_length;
            *right_win = tcx_cfg->tcx_mdct_window_minimum;
        }
        else if (right_mode == HALF_OVERLAP)
        {
            *right_overlap = tcx_cfg->tcx_mdct_window_half_length;
            *right_win = tcx_cfg->tcx_mdct_window_half;
        }
        else if (right_mode == FULL_OVERLAP)
        {
            *right_overlap = tcx_cfg->tcx_mdct_window_delay;
            *right_win = tcx_cfg->tcx_aldo_window_2;
        }
        else
        {
            assert(!"Not supported overlap");
        }
    }
    else
    {
        /* Left part */

        if (left_mode == TRANSITION_OVERLAP)
        {
            /* ACELP->TCX transition */
            *left_overlap = tcx_cfg->tcx_mdct_window_trans_lengthFB;
            *left_win = tcx_cfg->tcx_mdct_window_transFB;
        }
        else if (left_mode == MIN_OVERLAP)
        {
            *left_overlap = tcx_cfg->tcx_mdct_window_min_lengthFB;
            *left_win = tcx_cfg->tcx_mdct_window_minimumFB;
        }
        else if (left_mode == HALF_OVERLAP)
        {
            *left_overlap = tcx_cfg->tcx_mdct_window_half_lengthFB;
            *left_win = tcx_cfg->tcx_mdct_window_halfFB;
        }
        else if (left_mode == RECTANGULAR_OVERLAP)
        {
            *left_overlap = 0;
            *left_win = NULL;
        }
        else if (left_mode == FULL_OVERLAP)
        {
            *left_overlap = tcx_cfg->tcx_mdct_window_lengthFB;
            *left_win = tcx_cfg->tcx_aldo_window_1_FB_trunc;
        }
        else
        {
            assert(!"Not supported overlap");
        }

        /* Right part */

        if (right_mode == MIN_OVERLAP)
        {
            *right_overlap = tcx_cfg->tcx_mdct_window_min_lengthFB;
            *right_win = tcx_cfg->tcx_mdct_window_minimumFB;
        }
        else if (right_mode == HALF_OVERLAP)
        {
            *right_overlap = tcx_cfg->tcx_mdct_window_half_lengthFB;
            *right_win = tcx_cfg->tcx_mdct_window_halfFB;
        }
        else if (right_mode == RECTANGULAR_OVERLAP)
        {
            *right_overlap = 0;
            *right_win = NULL;
        }
        else if (right_mode == FULL_OVERLAP)
        {
            *right_overlap = tcx_cfg->tcx_mdct_window_delayFB;
            *right_win = tcx_cfg->tcx_aldo_window_2_FB;
        }
        else
        {
            assert(!"Not supported overlap");
        }
    }

    return;
}

/*-------------------------------------------------------------------*
* tcx_windowing_analysis()
*
*
*-------------------------------------------------------------------*/

void tcx_windowing_analysis(
    float const *signal,        /* i: signal vector                              */
    int L_frame,                /* i: frame length                               */
    int left_overlap,           /* i: left overlap length                        */
    float const *left_win,      /* i: left overlap window                        */
    int right_overlap,          /* i: right overlap length                       */
    float const *right_win,     /* i: right overlap window                       */
    float *output               /* o: windowed signal vector                     */
)
{
    int w;

    /* Left overlap */
    for (w = 0; w < left_overlap; w++)
    {
        *output++ = *signal++ * left_win[w];
    }

    /* Non overlapping region */
    for (w = 0; w < L_frame-(left_overlap+right_overlap)/2; w++)
    {
        *output++ = *signal++;
    }

    /* Right overlap */
    for (w = 0; w < right_overlap; w++)
    {
        *output++ = *signal++ * right_win[right_overlap-1-w];
    }

    return;
}


/*-------------------------------------------------------------------*
* WindowSignal()
*
*
*-------------------------------------------------------------------*/

void WindowSignal(
    TCX_config const *tcx_cfg,                /* input: configuration of TCX              */
    int offset,                               /* input: left folding point offset relative to the input signal pointer */
    const short left_overlap_mode,            /* input: overlap mode of left window half  */
    const short right_overlap_mode,           /* input: overlap mode of right window half */
    int * left_overlap_length,                /* output: TCX window left overlap length   */
    int * right_overlap_length,               /* output: TCX window right overlap length  */
    float const in[],                         /* input: input signal                      */
    int * L_frame,                            /* input/output: frame length               */
    float out[],                              /* output: output windowed signal           */
    int fullband                              /* input: fullband flag                     */
)
{
    int l, r;
    float const * left_win;
    float const * right_win;

    /*-----------------------------------------------------------*
     * Init                                                      *
     *-----------------------------------------------------------*/

    tcx_get_windows(tcx_cfg, left_overlap_mode, right_overlap_mode, &l, &left_win, &r, &right_win, fullband );

    /* Init lengths */

    /* if past frame is ACELP */
    if (left_overlap_mode == TRANSITION_OVERLAP)
    {
        /* Increase frame size for 5ms */
        if (!fullband)
        {
            *L_frame += tcx_cfg->tcx5Size;
            offset = -tcx_cfg->tcx_mdct_window_trans_length/2;
        }
        else
        {
            *L_frame += tcx_cfg->tcx5SizeFB;
            offset = -tcx_cfg->tcx_mdct_window_trans_lengthFB/2;
        }
    }

    /*-----------------------------------------------------------*
     * Windowing                                                 *
     *-----------------------------------------------------------*/

    tcx_windowing_analysis(in-l/2+offset, *L_frame, l, left_win, r, right_win, out);

    if (left_overlap_mode == FULL_OVERLAP)
    {
        /* fade truncated ALDO window to avoid discontinuities */
        if (!fullband)
        {
            v_mult(out, tcx_cfg->tcx_mdct_window_minimum, out, tcx_cfg->tcx_mdct_window_min_length);
        }
        else
        {
            v_mult(out, tcx_cfg->tcx_mdct_window_minimumFB, out, tcx_cfg->tcx_mdct_window_min_lengthFB);
        }
    }

    *left_overlap_length = l;
    *right_overlap_length = r;

    return;
}


/*-------------------------------------------------------------------*
* tcx_windowing_synthesis_current_frame()
*
*
*-------------------------------------------------------------------*/

void tcx_windowing_synthesis_current_frame(
    float *signal,            /* i/o: signal vector                            */
    float *window,            /* i: TCX window vector                          */
    float *window_half,       /* i: TCX window vector for half-overlap window  */
    float *window_min,        /* i: TCX minimum overlap window                 */
    int window_length,        /* i: TCX window length                          */
    int window_half_length,   /* i: TCX half window length                     */
    int window_min_length,    /* i: TCX minimum overlap length                 */
    int left_rect,            /* i: left part is rectangular                   */
    int left_mode,            /* i: overlap mode of left window half           */
    float *acelp_zir,         /* i: acelp ZIR                                  */
    float *old_syn,           /* i: old synthesis                              */
    float *syn_overl,         /* i: overlap synthesis                          */
    float *A_zir,
    float *window_trans,
    int acelp_zir_len,
    int acelp_mem_len,
    int last_core_bfi,        /* i: last mode                                  */
    int last_is_cng,
    int fullbandScale
)
{
    int i, overlap;
    float tmp[L_FRAME_MAX/2];

    /* Init */
    overlap = window_length>>1;

    /* Past-frame is TCX concealed as CNG and current-frame is TCX */
    if ( last_is_cng==1 && left_rect==0 )
    {
        if (!fullbandScale)
        {
            set_zero(acelp_zir, acelp_zir_len);
            syn_filt(A_zir, M,acelp_zir, acelp_zir, acelp_zir_len, signal+overlap+acelp_mem_len-M, 0);
        }
        else
        {
            lerp(acelp_zir, tmp, acelp_zir_len, acelp_zir_len*FSCALE_DENOM/fullbandScale);
            acelp_zir = tmp;
        }

        for (i = 0; i < acelp_zir_len; i++)
        {
            signal[i] *= (float)(i)/(float)(acelp_zir_len);
            signal[i] += acelp_zir[i]*(float)(acelp_zir_len-i)/(float)(acelp_zir_len);
        }
    }
    else if ( left_rect==1 && last_core_bfi==ACELP_CORE ) /* Rectangular window (past-frame is ACELP) */
    {
        for (i=0; i<overlap-acelp_mem_len; i++)
        {
            signal[i] = 0;
        }

        if (fullbandScale == 0)
        {
            /*OLA with ACELP*/
            for (i = 0; i < 2*acelp_mem_len; i++)
            {
                /*window decoded TCX with aliasing*/
                signal[i+overlap-acelp_mem_len] *= window_trans[i];
                /*Time TDAC: 1)forward part of ACELP*/
                signal[i+overlap-acelp_mem_len] +=old_syn[acelp_zir_len-2*acelp_mem_len+i]*window_trans[2*acelp_mem_len-i-1]*window_trans[2*acelp_mem_len-i-1];

                /*Time TDAC: 1)reward part of ACELP*/
                signal[i+overlap-acelp_mem_len] +=old_syn[acelp_zir_len-i-1]*window_trans[i]*window_trans[2*acelp_mem_len-i-1];
            }

            for (i=0; i<M; i++)
            {
                signal[overlap+acelp_mem_len-M+i] -= old_syn[acelp_zir_len-M+i];
            }
        }

        /* ZIR at the end of the ACELP frame */
        acelp_zir_len=64;

        if (!fullbandScale)
        {
            set_zero(acelp_zir, acelp_zir_len);
            syn_filt(A_zir, M,acelp_zir, acelp_zir, acelp_zir_len, signal+overlap+acelp_mem_len-M, 0);
        }
        else
        {
            lerp(acelp_zir, tmp, acelp_zir_len * fullbandScale / FSCALE_DENOM, acelp_zir_len);
            acelp_zir_len = acelp_zir_len * fullbandScale / FSCALE_DENOM;
            acelp_zir = tmp;

            if(acelp_zir_len >= 2.0f * 64)
            {
                /* apply a simple low-pass to the ZIR, to avoid potentially unmasked HF content */
                for(i=2; i<acelp_zir_len; i++)
                {
                    acelp_zir[i] = 0.25f * acelp_zir[i-2] + 0.35f * acelp_zir[i-1] + 0.40f * acelp_zir[i];
                }
                acelp_zir[acelp_zir_len-1] = 0.40 * acelp_zir[acelp_zir_len-1] +  0.35f * acelp_zir[acelp_zir_len-1] + 0.25f * acelp_zir[acelp_zir_len-2];
                acelp_zir[acelp_zir_len-2] = 0.40 * acelp_zir[acelp_zir_len-2] +  0.35f * acelp_zir[acelp_zir_len-1] + 0.25f * acelp_zir[acelp_zir_len-1];
                for(i=acelp_zir_len-3; i>=0; i--)
                {
                    acelp_zir[i] = 0.40f * acelp_zir[i] + 0.35f * acelp_zir[i+1] + 0.25f * acelp_zir[i+2];
                }
            }
        }

        for (i = 0; i < acelp_zir_len; i++)
        {
            /*remove reconstructed ZIR and add ACELP ZIR*/
            signal[i+overlap+acelp_mem_len] -= acelp_zir[i]*(float)(acelp_zir_len-i)/(float)acelp_zir_len;
        }
    }
    else if ( left_rect==1 && last_core_bfi!=0 )          /* Rectangular window (past-frame is TCX) */
    {
        for (i=0; i<overlap+acelp_mem_len; i++)
        {
            signal[i] = 0;
        }
        for (i=0; i<window_length; i++)
        {
            signal[i+overlap+acelp_mem_len] *= window[i];
        }
    }
    else if (left_rect != 1 && last_core_bfi == 0)        /* Normal window (past-frame is ACELP) */
    {
        for (i=0; i<window_length; i++)
        {
            signal[i] *= window[i];
        }

        for (i=0; i<window_length; i++)
        {
            signal[i] += syn_overl[i];
        }
    }
    else                                                  /* Normal window (past-frame is TCX) */
    {
        if (left_mode == 2)
        {
            /* min. overlap */
            int w;
            for (i = 0; i < (window_length - window_min_length)/2; i++)
            {
                signal[i] = 0.0f;
            }
            for (w = 0; w < window_min_length; i++, w++)
            {
                signal[i] *= window_min[w];
            }
        }
        else if (left_mode == 3)
        {
            /* half OL */
            int w;
            for (i = 0; i < (window_length-window_half_length)/2; i++)
            {
                signal[i] = 0.0f;
            }
            for (w = 0; w < window_half_length; i++, w++)
            {
                signal[i] *= window_half[w];
            }
        }
        else
        {
            /* normal full/maximum overlap */
            for (i = 0; i < window_length; i++)
            {
                signal[i] *= window[i];
            }
        }
    }

    return;
}


/*-------------------------------------------------------------------*
* tcx_windowing_synthesis_past_frame()
*
*
*-------------------------------------------------------------------*/

void tcx_windowing_synthesis_past_frame(
    float *signal,           /* i/o: signal vector                            */
    float *window,           /* i: TCX window vector                          */
    float *window_half,      /* i: TCX window vector for half-overlap window  */
    float *window_min,       /* i: TCX minimum overlap window                 */
    int window_length,       /* i: TCX window length                          */
    int window_half_length,  /* i: TCX half window length                     */
    int window_min_length,   /* i: TCX minimum overlap length                 */
    int right_mode           /* i: overlap mode (left_mode of current frame)  */
)
{
    int i;

    if (right_mode == MIN_OVERLAP)
    {
        /* min. overlap */
        int w;
        for (i = (window_length - window_min_length)/2, w = 0; w < window_min_length; i++, w++)
        {
            signal[i] *= window_min[window_min_length-1-w];
        }
        for (; i < window_length; i++)
        {
            signal[i] = 0.0f;
        }
    }
    else if (right_mode == HALF_OVERLAP)
    {
        /* half OL */
        int w;
        for (i = (window_length-window_half_length)/2, w = 0; w < window_half_length; i++, w++)
        {
            signal[i] *= window_half[window_half_length-1-w];
        }
        for (; i < window_length; i++)
        {
            signal[i] = 0.0f;
        }
    }
    else if(right_mode == FULL_OVERLAP)
    {
        /* normal full/maximum overlap */
        for (i = 0; i < window_length; i++)
        {
            signal[i] *= window[window_length-1-i];
        }
    }

    return;
}


/*-------------------------------------------------------------------*
* lpc2mdct()
*
*
*-------------------------------------------------------------------*/

void lpc2mdct(
    float *lpcCoeffs,
    int lpcOrder,
    float *mdct_gains
)
{
    float RealData[2*FDNS_NPTS];
    float ImagData[2*FDNS_NPTS];
    float tmp;
    int i, sizeN;

    sizeN = 2*FDNS_NPTS;

    /*ODFT*/
    for(i=0; i<lpcOrder+1; i++)
    {
        tmp = (float)(((float)i)*EVS_PI/(float)(sizeN));
        RealData[i] = (float)( lpcCoeffs[i]*cos(tmp));
        ImagData[i] = (float)(-lpcCoeffs[i]*sin(tmp));
    }

    for(; i<sizeN; i++)
    {
        RealData[i] = 0.f;
        ImagData[i] = 0.f;
    }

    DoRTFTn(RealData, ImagData, sizeN);

    /*Get amplitude*/
    {
        for(i=0; i<FDNS_NPTS; i++)
        {
            mdct_gains[i] = (float)(1.0f/sqrt(RealData[i]*RealData[i] + ImagData[i]*ImagData[i]));
        }
    }

    return;
}


/*-------------------------------------------------------------------*
* mdct_noiseShaping()
*
*
*-------------------------------------------------------------------*/

void mdct_noiseShaping(
    float x[],
    int lg,
    const float gains[]
)
{
    int i, j, k, l;
    float g;
    int m, n, k1, k2;


    j = 0;
    k = lg/FDNS_NPTS;
    m = lg%FDNS_NPTS;

    if (m)
    {
        if ( m <= (FDNS_NPTS/2) )
        {
            n = FDNS_NPTS/m;
            k1 = k;
            k2 = k + 1;
        }
        else
        {
            n = FDNS_NPTS/(FDNS_NPTS-m);
            k1 = k + 1;
            k2 = k;
        }

        for (i=0; i<lg; )
        {
            if (j%n)
            {
                k = k1;
            }
            else
            {
                k = k2;
            }
            g = gains[j++];

            /* Limit number of loops, if end is reached */
            k = min(k, lg-i);

            for (l=0; l < k; l++)
            {
                x[i++] *= g;
            }
        }
    }
    else
    {
        for (i=0; i<lg; )
        {
            g = gains[j++];

            for (l=0; l < k; l++)
            {
                x[i++] *= g;
            }
        }
    }

    return;

}


/*-------------------------------------------------------------------*
* PsychAdaptLowFreqDeemph()
*
*
*-------------------------------------------------------------------*/

void PsychAdaptLowFreqDeemph(
    float x[],
    const float lpcGains[],
    float lf_deemph_factors[]
)
{
    int i;
    float max, fac, tmp;

    max = tmp = lpcGains[0];

    /* find minimum (tmp) and maximum (max) of LPC gains in low frequencies */
    for (i = 1; i < 9; i++)
    {
        if (tmp > lpcGains[i])
        {
            tmp = lpcGains[i];
        }
        if (max < lpcGains[i])
        {
            max = lpcGains[i];
        }
    }

    tmp *= 32.0f;

    if ((max < tmp) && (tmp > FLT_MIN))
    {
        fac = tmp = (float)pow(max / tmp, 0.0078125f);

        if (lf_deemph_factors)
        {
            /* gradual lowering of lowest 32 bins; DC is lowered by (max/tmp)^1/4 */
            for (i = 31; i >= 0; i--)
            {
                x[i] *= fac;
                lf_deemph_factors[i] *= fac;
                fac  *= tmp;
            }
        }
        else
        {
            /* gradual lowering of lowest 32 bins; DC is lowered by (max/tmp)^1/4 */
            for (i = 31; i >= 0; i--)
            {
                x[i] *= fac;
                fac  *= tmp;
            }
        }
    }

    return;
}


/*-------------------------------------------------------------------*
* AdaptLowFreqDeemph()
*
*
*-------------------------------------------------------------------*/

void AdaptLowFreqDeemph(
    float x[],
    short tcx_lpc_shaped_ari,
    const float lpcGains[],
    const int lg,
    float lf_deemph_factors[]
)
{

    int i, i_max_old,i_max;

    if(!tcx_lpc_shaped_ari)
    {
        /* 1. find first magnitude maximum in lower quarter of spectrum */
        i_max = -1;

        for (i = 0; i < lg/4; i++)
        {
            if ((x[i] <= -4.0f) || (x[i] >= 4.0f))
            {
                x[i] += (x[i] < 0.0f) ? 2.0f : -2.0f;
                i_max = i;
                break;
            }
        }
        /* 2. expand value range of all xi up to i_max: two extra steps */

        for (i = 0; i < i_max; i++)
        {
            x[i] *= 0.5f;
            lf_deemph_factors[i] *= 0.5f;
        }
        /* 3. find first magnitude maximum in lower quarter of spectrum */
        i_max_old = i_max;

        if (i_max_old > -1)
        {
            i_max = -1;

            for (i = 0; i < lg/4; i++)
            {
                if ((x[i] <= -4.0f) || (x[i] >= 4.0f))
                {
                    x[i] += (x[i] < 0.0f) ? 2.0f : -2.0f;
                    i_max = i;
                    break;
                }
            }
        }
        /* 4. expand value range of all xi up to i_max: two extra steps */
        for (i = 0; i < i_max; i++)
        {
            x[i] *= 0.5f;
            lf_deemph_factors[i] *= 0.5f;
        }
        /* 5. always expand two lines; lines could be at index 0 and 1! */
        if (i_max < i_max_old)
        {
            i_max = i_max_old;
        }
        i = i_max + 1;
        if (x[i] < 0.0f)
        {
            if (x[i] > -4.0f)
            {
                lf_deemph_factors[i] *= 0.5f;
            }
            x[i] = (x[i] <= -4.0f) ? x[i] + 2.0f : x[i] * 0.5f;
        }
        else
        {
            if (x[i] < 4.0f)
            {
                lf_deemph_factors[i] *= 0.5f;
            }
            x[i] = (x[i] >=  4.0f) ? x[i] - 2.0f : x[i] * 0.5f;
        }
        i++;
        if (x[i] < 0.0f)
        {
            if (x[i] > -4.0f)
            {
                lf_deemph_factors[i] *= 0.5f;
            }
            x[i] = (x[i] <= -4.0f) ? x[i] + 2.0f : x[i] * 0.5f;
        }
        else
        {
            if (x[i] < 4.0f)
            {
                lf_deemph_factors[i] *= 0.5f;
            }
            x[i] = (x[i] >=  4.0f) ? x[i] - 2.0f : x[i] * 0.5f;
        }
    }
    else
    {
        /*if(!tcx_lpc_shaped_ari)*/
        PsychAdaptLowFreqDeemph(x, lpcGains, lf_deemph_factors);
    }/*if(!tcx_lpc_shaped_ari)*/

    return;
}


/*-------------------------------------------------------------------*
* tcx_noise_filling()
*
*
*-------------------------------------------------------------------*/

void tcx_noise_filling(
    float *Q,
    const int noiseFillSeed,
    const int iFirstLine,
    const int lowpassLine,
    const int nTransWidth,
    const int L_frame,
    float tiltCompFactor,
    float fac_ns,
    unsigned char *infoTCXNoise
)
{
    int i, m, segmentOffset;
    int win; /* window coefficient */
    Word16 seed;
    float tilt_factor, nrg, tmp1, tmp2;

    /* 16-bit random number generator seed for generating filled lines */
    seed = (Word16)noiseFillSeed;

    tilt_factor = (float)pow(max(0.375f, tiltCompFactor), 1.0f/(float)L_frame);
    fac_ns /= (float)(nTransWidth * nTransWidth);

    /* find last nonzero line below iFirstLine, use it as start offset */
    for (i = iFirstLine; i > (iFirstLine >> 1); i--)
    {
        if (Q[i] != 0.0f)
        {
            break;
        }
    }
    fac_ns *= (float)pow(tilt_factor, (float)i);

    nrg = 1e-9f;
    win = 0;
    segmentOffset = ++i;

    for (; i < lowpassLine; i++)
    {
        fac_ns *= tilt_factor;

        if (Q[i] != 0.0f)
        {
            if (win > 0)
            {
                /* RMS-normalize current noise-filled segment */
                tmp1 = (float)sqrt((i - segmentOffset) / nrg);
                tmp2 = tmp1 * (float)nTransWidth;

                for (m = segmentOffset; m < i-win; m++)
                {
                    Q[m] *= tmp2;
                }

                for (; win > 0; win--)
                {
                    Q[m++] *= tmp1 * (float)win;
                }
                nrg = 1e-9f; /* start new segment: reset noise segment energy */
            }
            segmentOffset = i + 1;
        }
        else
        {
            /* line is zero, so fill line and update window and energy */
            if (win < nTransWidth)
            {
                win++;
            }
            tmp1 = (float)own_random(&seed);
            nrg += tmp1 * tmp1;   /* sum up energy of current noise segment */
            Q[i] = tmp1 * (float)win * fac_ns;

            if(infoTCXNoise)
            {
                /* set noiseflags for IGF */
                infoTCXNoise[i] = 1;
            }
        }
    }

    if (win > 0)
    {
        /* RMS-normalize uppermost noise-filled segment */
        tmp1 = (float)sqrt((lowpassLine - segmentOffset) / nrg);
        tmp2 = tmp1 * (float)nTransWidth;

        for (m = segmentOffset; m < lowpassLine; m++)
        {
            Q[m] *= tmp2;
        }
    }

    return;
}


/*-------------------------------------------------------------------*
* InitTnsConfigs()
*
*
*-------------------------------------------------------------------*/

void InitTnsConfigs(
    int nSampleRate,
    int L_frame,
    STnsConfig tnsConfig[2][2],
    int igfStopFreq,
    int bitrate
)
{
    if (bitrate > ACELP_32k)
    {
        InitTnsConfiguration(nSampleRate, L_frame/2, &tnsConfig[0][0], igfStopFreq, bitrate);
    }
    InitTnsConfiguration(nSampleRate, L_frame,   &tnsConfig[1][0], igfStopFreq, bitrate);
    InitTnsConfiguration(nSampleRate, L_frame  +L_frame/4, &tnsConfig[1][1], igfStopFreq, bitrate);

    return;
}



/*-------------------------------------------------------------------*
* SetTnsConfig()
*
*
*-------------------------------------------------------------------*/

void SetTnsConfig(
    TCX_config * tcx_cfg,
    int isTCX20,
    int isAfterACELP
)
{
    tcx_cfg->pCurrentTnsConfig = &tcx_cfg->tnsConfig[isTCX20][isAfterACELP];
    assert(tcx_cfg->pCurrentTnsConfig != NULL);

    return;
}

