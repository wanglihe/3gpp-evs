/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "cnst.h"
#include "prot.h"


/*-------------------------------------------------------------------*
 * TCX_MDCT()
 *
 *
 *-------------------------------------------------------------------*/

void TCX_MDCT(
    float const *x,
    float *y,
    int l,
    int m,
    int r
)
{
    short i;
    float dctInBuffer[N_MAX];

    /* Init */
    for(i=0; i<m/2; i++)
    {
        dctInBuffer[m/2 + r/2 + i]           = -1.0f * x [l + m/2 - 1 - i];
    }

    for(i=0; i<l/2; i++)
    {
        dctInBuffer[m/2 + r/2 + m/2 + i]     = x[i]  - x [l - 1 - i];
    }

    for(i=0; i<m/2; i++)
    {
        dctInBuffer[m/2 + r/2 - 1 - i]       = -1.0f * x[l + m/2 + i];
    }

    for(i=0; i<r/2; i++)
    {
        dctInBuffer[m/2 + r/2 - 1 - m/2 - i] = -1.0f * x[l + m   + i] - x[l + m + r - 1 - i];
    }

    edct( dctInBuffer, y, l/2 + m + r/2 );

    v_multc( y, (float)sqrt((float)NORM_MDCT_FACTOR / (l/2 + m + r/2)), y, l/2 + m + r/2 );

    return;
}


/*-------------------------------------------------------------------*
 * TCX_MDST()
 *
 *
 *-------------------------------------------------------------------*/

void TCX_MDST(
    float const *x,
    float *y,
    int l,
    int m,
    int r
)
{
    short i;
    float dctInBuffer[N_MAX];

    /* Init */
    for(i=0; i<m/2; i++)
    {
        dctInBuffer[m/2 + r/2 + i] = -1.0f * x[l + m/2 - 1 - i];
    }
    for(i=0; i<l/2; i++)
    {
        dctInBuffer[m/2 + r/2 + m/2 + i] = -1.0 * x[i] - x[l -1 - i];
    }
    for(i=0; i<m/2; i++)
    {
        dctInBuffer[m/2 + r/2 - 1 - i] = -1.0f * x[l + m/2 + i];
    }
    for(i=0; i<r/2; i++)
    {
        dctInBuffer[m/2 + r/2 - 1 - m/2 - i] = -1.0f * x[l+m+i] + x[l+m+r-1-i];
    }

    edst( dctInBuffer, y, l/2 + m + r/2 );

    v_multc( y, (float)sqrt((float)NORM_MDCT_FACTOR / (l/2 + m + r/2)), y, l/2 + m + r/2 );
}


/*-------------------------------------------------------------------*
 * TCX_MDCT_Inverse()
 *
 *
 *-------------------------------------------------------------------*/

void TCX_MDCT_Inverse(
    float *x,
    float *y,
    int l,
    int m,
    int r
)
{
    short i;
    const int L2 = l >> 1, R2 = r >> 1;
    float f;

    edct( x, y + L2, L2 + m + R2 );

    for (i = 0; i < R2; i++)
    {
        y[l + m + R2 + i] = -1.0f * y[L2 + i];   /* fold out right end of DCT */
    }

    mvr2r(y + L2 + m + R2, y, L2);  /* negate, fold out left end of DCT */

    for (i = 0; i < ((L2 + m + R2) >> 1); i++)
    {
        f = y[L2 + i];
        y[L2 + i] = -1.0f * y[l + m + R2 - 1 - i]; /* time-reverse mid of DCT */
        y[l + m + R2 - 1 - i] = -1.0f * f;
    }

    v_multc( y, (float)sqrt((float)(l/2 + m + r/2) / NORM_MDCT_FACTOR), y, l + m + r );

    return;
}
