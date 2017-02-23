/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "prot.h"
#include "options.h"

/*-------------------------------------------------------------------*
 * tcxGetNoiseFillingTilt()
 *
 *
 *-------------------------------------------------------------------*/

int tcxGetNoiseFillingTilt(
    float A[],
    int L_frame,
    int mode,
    float *noiseTiltFactor
)
{
    int firstLine;

    if (mode)
    {
        firstLine = L_frame / 6;
        *noiseTiltFactor = 0.5625f;
    }
    else
    {
        firstLine = L_frame / 8;
        *noiseTiltFactor = get_gain( A+1, A, M, NULL );
        *noiseTiltFactor = min(1.0f, (*noiseTiltFactor) + 0.09375f);
    }

    return firstLine;
}

/*-------------------------------------------------------------------*
 * tcxFormantEnhancement()
 *
 *
 *-------------------------------------------------------------------*/

void tcxFormantEnhancement(
    float xn_buf[],
    float gainlpc[],
    float spectrum[],
    int L_frame
)
{
    int k;
    int i, j, l = 0;
    float fac, step;

    k = L_frame / FDNS_NPTS;

    /* Formant enhancement via square root of the LPC gains */
    xn_buf[0] = (float)sqrt(gainlpc[0]);
    xn_buf[1] = (float)sqrt(gainlpc[1]);
    fac = 1.0f / min(xn_buf[0], xn_buf[1]);

    for (i = 1; i < FDNS_NPTS - 1; i++)
    {
        xn_buf[i+1] = (float)sqrt(gainlpc[i+1]);

        if ((xn_buf[i-1] <= xn_buf[i]) && (xn_buf[i+1] <= xn_buf[i]))
        {
            step = max(xn_buf[i-1], xn_buf[i+1]);
            step = (1.0f / step - fac) / (float)(i - l);
            xn_buf[l] = 1.0f;
            fac += step;
            for (j = l + 1; j < i; j++)
            {
                xn_buf[j] = min(1.0f, xn_buf[j] * fac);
                fac += step;
            }
            l = i;
        }
    }

    /* i = tcx_cfg->fdns_npts - 1; Completing changes to gains */
    step = min(xn_buf[i-1], xn_buf[i]);
    step = (1.0f / step - fac) / (float)(i - l);
    xn_buf[l] = 1.0f;
    fac += step;
    for (j = l + 1; j < i; j++)
    {
        xn_buf[j] = min(1.0f, xn_buf[j] * fac);
        fac += step;
    }
    xn_buf[i] = 1.0f;

    /* Application of changed gains onto decoded MDCT lines */
    for (i = j = 0; i < L_frame; j++)
    {
        for (l = 0; l < k; i++, l++)
        {
            spectrum[i] *= xn_buf[j];
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * tcxGetNoiseFillingTilt()
 *
 *
 *-------------------------------------------------------------------*/

void tcxInvertWindowGrouping(
    TCX_config *tcx_cfg,
    float xn_buf[],
    float spectrum[],
    int L_frame,
    int fUseTns,
    int last_core,
    int index,
    int frame_cnt,
    int bfi
)
{
    short i, w, t_integer;
    int L_win, L_spec;

    if (frame_cnt && !bfi && last_core!=0)
    {
        /* fix sub-window overlap */
        tcx_cfg->tcx_last_overlap_mode = tcx_cfg->tcx_curr_overlap_mode;
    }

    if (((!bfi) &&((tcx_cfg->tcx_last_overlap_mode != FULL_OVERLAP) ||
                   ((tcx_cfg->tcx_curr_overlap_mode == FULL_OVERLAP) && (frame_cnt == 0) && (index == 0))))
            ||
            ((bfi) &&((tcx_cfg->tcx_last_overlap_mode != FULL_OVERLAP) &&
                      !(tcx_cfg->tcx_curr_overlap_mode == FULL_OVERLAP))))
    {

        /* ungroup sub-windows: deinterleave MDCT bins into separate windows */
        for (t_integer = w = 0; w < 2; w++)
        {
            for (i = w; i < L_frame; i += 2)
            {
                xn_buf[t_integer++] = spectrum[i];
            }
        }

        mvr2r( xn_buf, spectrum, L_frame );

        if( tcx_cfg->fIsTNSAllowed && !bfi && fUseTns )
        {
            L_win = L_frame >> 1;
            L_spec = tcx_cfg->tnsConfig[0][0].iFilterBorders[0];

            /* rearrange LF sub-window lines prior to TNS synthesis filtering */
            if( L_spec < L_frame )
            {
                mvr2r( spectrum+8, spectrum+16, L_spec/2-8 );
                mvr2r( spectrum+L_frame/2, spectrum+8, 8 );
                mvr2r( spectrum+L_frame/2+8, spectrum+L_spec/2+8, L_spec/2-8 );
            }
            else
            {
                mvr2r( spectrum+8, xn_buf, L_win );
                mvr2r( xn_buf, spectrum+16, L_win-8 );
                mvr2r( xn_buf+L_win-8, spectrum+8, 8);
            }
        }
    }

    return;
}
