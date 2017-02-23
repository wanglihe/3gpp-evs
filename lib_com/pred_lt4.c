/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*-------------------------------------------------------------------*
 * pred_lt4()
 *
 * Compute the result of long term prediction with fractionnal
 * interpolation of resolution 1/4.
 *
 * On return, exc[0..L_subfr-1] contains the interpolated signal
 *   (adaptive codebook excitation)
 *-------------------------------------------------------------------*/

void pred_lt4(
    const float excI[],  /* i  : input excitation buffer  */
    float excO[],  /* o  : output excitation buffer */
    const short T0,      /* i  : integer pitch lag        */
    short frac,    /* i  : fraction of lag          */
    const short L_subfr, /* i  : subframe size            */
    const float *win,    /* i  : interpolation window     */
    const short nb_coef,  /* i  : nb of filter coef        */
    const short up_sample /* i  : up_sample        */
)
{
    short   i, j;
    float s;
    const float *x1, *x2, *x0, *c1, *c2;

    x0 = &excI[-T0];
    frac = -frac;

    if (frac < 0)
    {
        frac += up_sample;
        x0--;
    }

    for (j=0; j<L_subfr; j++)
    {
        x1 = x0++;
        x2 = x1+1;
        c1 = &win[frac];
        c2 = &win[up_sample-frac];

        s = 0.0f;
        for(i=0; i<nb_coef; i++, c1+=up_sample, c2+=up_sample)
        {
            s += (*x1--) * (*c1) + (*x2++) * (*c2);
        }
        excO[j] = s;
    }

    return;
}

/*-------------------------------------------------------------------*
 * pred_lt4_tc()
 *
 * adapt. search of the second impulse in the same subframe (when appears)
 * On return, exc[0..L_subfr-1] contains the interpolated signal
 *   (adaptive codebook excitation)
 *-------------------------------------------------------------------*/

void pred_lt4_tc(
    float exc[],   /* i/o: excitation buffer        */
    const short T0,      /* i  : integer pitch lag        */
    short frac,    /* i:   fraction of lag          */
    const float *win,    /* i  : interpolation window     */
    const short imp_pos, /* i  : glottal impulse position */
    const short i_subfr  /* i  : subframe index           */
)
{
    short i, j;
    float s;
    const float *x1, *x2, *x0, *c1, *c2;
    float excO[L_SUBFR+1];
    float excI[2*L_SUBFR];

    mvr2r( exc + i_subfr - L_SUBFR, excI, 2*L_SUBFR );

    if( ((T0+imp_pos-L_IMPULSE2) < L_SUBFR) && (T0 < L_SUBFR) )
    {
        set_f( excI + L_SUBFR - T0, 0 , T0 );
        set_f( excO, 0, L_SUBFR+1 );
        x0 = excI - T0 + L_SUBFR;
        x0 += T0;
        frac = -frac;

        if (frac < 0)
        {
            frac += PIT_UP_SAMP;
            x0--;
        }

        for (j=T0; j<L_SUBFR+1; j++)
        {
            x1 = x0++;
            x2 = x1+1;
            c1 = &win[frac];
            c2 = &win[PIT_UP_SAMP-frac];

            s = 0.0f;
            for(i=0; i<L_INTERPOL2; i++, c1+=PIT_UP_SAMP, c2+=PIT_UP_SAMP)
            {
                s += (*x1--) * (*c1) + (*x2++) * (*c2);
            }
            excO[j] = s;
        }

        for(i=T0; i<L_SUBFR; i++)
        {
            exc[i+i_subfr] += PIT_SHARP*excO[i];
        }

    }

    return;
}
