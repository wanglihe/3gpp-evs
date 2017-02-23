/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "prot.h"
#include "cnst.h"

/*-------------------------------------------------------------------*
 * corr_xh()
 *
 * Compute the correlation between the target signal and the impulse
 * response of the weighted synthesis filter.
 *
 * y[i] = sum(j=i,l-1) x[j]*h[j-i], i=0,l-1
 *-------------------------------------------------------------------*/

void corr_xh(
    const float *x,       /* i  : target signal                                   */
    float *y,       /* o  : correlation between x[] and h[]                 */
    const float *h,       /* i  : impulse response (of weighted synthesis filter) */
    const int   L_subfr   /* i  : length of the subframe                        */
)
{
    short i, j;
    float s;

    for (i = 0; i < L_subfr; i++)
    {
        s = 0.0f;
        for (j = i; j < L_subfr; j++)
        {
            s += x[j]*h[j-i];
        }

        y[i] = s;
    }

    return;
}
