/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "options.h"
#include "jbm_pcmdsp_window.h"
#include "cnst.h"

/*-----------------------------------------------------------------------*
* hannWindow()
*
* Generates a Hann window (cos-shaped) of length n
*-----------------------------------------------------------------------*/

void hannWindow(
    uint16_t n,
    Float * w
)
{
    uint16_t i;
    Float arg;

    for (i = 0; i < n/2; i++)
    {
        arg = ((2.0f * EVS_PI) * i) / (Float) (n);
        w[i] = (Float) ((1.0f - cos (arg)) / 2.0f);
    }

    for ( ; i < n; i++)
    {
        w[i] = 1.0f - w[i-n/2];
    }

    return;
}


/*-----------------------------------------------------------------------*
* overlapAdd()
*
*  Overlap/Add of two signal with a given window
*-----------------------------------------------------------------------*/

void overlapAdd(
    const int16_t *fadeOut,
    const int16_t *fadeIn,
    int16_t *out,
    uint16_t n,
    uint16_t nChannels,
    const float *fadeOutWin,
    const float *fadeInWin
)
{
    float fdOutVal, fdInVal;
    int16_t i, j, hannIter;
    int32_t combinedVal;

    for(j = 0; j < nChannels; j++)
    {
        /* reset Hann window iterator to beginning (both channels use same window) */
        hannIter = 0;
        for(i = j; i < n; i += nChannels)
        {
            fdOutVal = fadeOut[i] * fadeOutWin[hannIter];
            fdInVal = fadeIn[i] * fadeInWin[hannIter];
            /* round combinedVal value (taking care of sign) */
            combinedVal = (int32_t)( (fdInVal + fdOutVal) + 0.5 );

            if( fdInVal + fdOutVal < 0.0 )
                combinedVal = (int32_t)( (fdInVal + fdOutVal) - 0.5 );
            /* saturate value */
            if (combinedVal > 32767)
            {
                combinedVal = 32767;
            }
            else if (combinedVal < -32768)
            {
                combinedVal = -32768;
            }
            out[i] = (int16_t) combinedVal;
            hannIter++;
        }
    }

    return;
}
