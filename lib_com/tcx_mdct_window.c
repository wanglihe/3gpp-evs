/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <assert.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*-------------------------------------------------------------------
 * mdct_window_sine()
 *
 *
 *-------------------------------------------------------------------*/

void mdct_window_sine(
    float *window,
    int n
)
{
    int i;
    float c;

    c = EVS_PI / ( 2.0f * (float)n );

    for ( i = 0 ; i < n ; i++ )
    {
        window[i] = (float)sin( c * ( 0.5f + (float)i ) );
    }

    return;

}

/*-------------------------------------------------------------------
 * mdct_window_aldo()
 *
 *
 *-------------------------------------------------------------------*/

void mdct_window_aldo(
    float *window1,
    float *window2,
    int n
)
{
    int i, n1, n2, d;
    const float *p1, *p2;

    /* set table pointers and decimation factor */
    switch (n)
    {
    case 320/2:
        p1 = window_48kHz + 2;
        p2 = window_48kHz + 1110 - 3;
        d = 6;
        break;
    case 512/2:
        p1 = window_256kHz;
        p2 = window_256kHz + 592 - 1;
        d = 2;
        break;
    case 640/2:
        p1 = window_48kHz + 1;
        p2 = window_48kHz + 1110 - 2;
        d = 3;
        break;
    case 1024/2:
        p1 = window_256kHz;
        p2 = window_256kHz + 592 - 1;
        d = 1;
        break;
    case 1280/2:
        p1 = window_48kHz + 1;
        p2 = window_48kHz + 1110 - 2;
        d = 3;
        break;
    case 1920/2:
        p1 = window_48kHz;
        p2 = window_48kHz + 1110 - 1;
        d = 1;
        break;
    default:
        assert(0);
        return;
    }

    /* set lengths */
    n1 = n * 23 / 32; /* left slope length */
    n2 = n * 14 / 32; /* right slope length */

    /* first part (long slope) */
    if (n != 1280/2)
    {
        for (i = 0; i < n/2; i++)
        {
            *window1 = *p1;
            window1++;
            p1 += d;
        }

        if ((n == 512/2) || (n == 320 / 2)) p1++;

        for ( ; i < n1; i++)
        {
            *window1 = *p1;
            window1++;
            p1 += d;
        }
    }
    else
    {
        const float *pi = window_8_16_32kHz;

        for (i = 0; i < n/2; i+=2)
        {
            *window1 = *p1;
            window1++;
            p1 += d;

            *window1 = *pi;
            window1++;
            pi++;
        }
        for ( ; i < n1; i+=2)
        {
            *window1 = *pi;
            window1++;
            pi++;

            *window1 = *p1;
            window1++;
            p1 += d;
        }
    }

    /* second part (short slope) */

    if (n != 1280/2)
    {
        for (i = 0; i < n2/2; i++)
        {
            *window2 = *p2;
            window2++;
            p2 -= d;
        }

        if ((n == 512/2) || (n == 320 / 2)) p2--;

        for ( ; i < n2; i++)
        {
            *window2 = *p2;
            window2++;
            p2 -= d;
        }
    }
    else
    {
        const float *pi = window_8_16_32kHz + 370 - 1;

        for (i = 0; i < n2/2; i+=2)
        {
            *window2 = *p2;
            window2++;
            p2 -= d;

            *window2 = *pi;
            window2++;
            pi--;
        }

        for ( ; i < n2; i+=2)
        {
            *window2 = *pi;
            window2++;
            pi--;

            *window2 = *p2;
            window2++;
            p2 -= d;
        }
    }

    return;
}
