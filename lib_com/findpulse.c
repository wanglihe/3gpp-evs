/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "prot.h"
#include "cnst.h"


/*----------------------------------------------------------------------------------*
 *  findpulse()
 *
 *  Find first pitch pulse in a frame
 *----------------------------------------------------------------------------------*/

short findpulse(              /* o  : pulse position        */
    const short L_frame,      /* i  : length of the frame   */
    const float res[],        /* i  : residual signal       */
    const short T0,           /* i  : integer pitch         */
    const short enc_dec,      /* i  : flag enc/dec, 0 - enc, 1 - dec          */
    short *sign         /* i/o: sign of the maximum   */
)
{
    const float *ptr;
    float val, maxval;
    short i, maxi;
    float resf[L_FRAME16k];   /* Low pass filtered residual */

    if ( enc_dec == ENC )
    {
        /*-----------------------------------------------------------------*
         * Very simple LP filter
         *-----------------------------------------------------------------*/

        resf[0] = 0.50f * res[0] + 0.25f * res[1];
        for (i=1; i<L_frame-1; i++)
        {
            resf[i] = 0.25f * res[i-1] + 0.5f * res[i] + 0.25f * res[i+1];
        }
        resf[L_frame-1] = 0.25f * res[L_frame-2] + 0.50f * res[L_frame-1];

        /*-----------------------------------------------------------------*
         * Find "biggest" pulse in the last pitch section
         *-----------------------------------------------------------------*/

        ptr = resf + L_frame - 1;
        maxval = 0;
        maxi = 0;
        for (i=0; i<T0; i++)
        {
            val = (float)fabs(*ptr);
            if (val>maxval)
            {
                maxval = val;
                maxi = i;
                if(*ptr >= 0)
                {
                    *sign = 0;
                }
                else
                {
                    *sign = 1;
                }
            }
            ptr--;
        }
    }
    else
    {
        /*-----------------------------------------------------------------*
         * Find "biggest" pulse in the last pitch section according to the sign
         *-----------------------------------------------------------------*/

        ptr = res;
        maxval = 0;
        maxi = 0;

        if( *sign == 0 )
        {
            for (i=1; i<=T0; i++)
            {
                val = *ptr++;
                if ( val >= maxval )
                {
                    maxval = val;
                    maxi   = i;
                }
            }
        }
        else
        {
            for (i=1; i<=T0; i++)
            {
                val = *ptr++;
                if (val<=maxval)
                {
                    maxval = val;
                    maxi   = i;
                }
            }
        }
    }

    return(maxi);
}
