/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "prot.h"


/*-------------------------------------------------------------------*
 * interpolation()
 *
 * Fractional interpolation of signal at position (frac/resol)
 *-------------------------------------------------------------------*/

float interpolation(      /* o  : interpolated value   */
    const float *x,       /* i  : input vector         */
    const float *win,     /* i  : interpolation window */
    const short frac,     /* i  : fraction             */
    const short up_samp,  /* i  : upsampling factor    */
    const short nb_coef   /* i  : nb of filter coef    */
)
{
    short i;
    float s;
    const float *x1, *x2, *c1, *c2;

    x1 = &x[0];
    x2 = &x[1];
    c1 = &win[frac];
    c2 = &win[up_samp-frac];
    s = 0.0f;

    for (i=0; i< nb_coef; i++ )
    {
        s += (*x1--) * (*c1) + (*x2++) * (*c2);
        c1 += up_samp;
        c2 += up_samp;
    }

    return s;
}
