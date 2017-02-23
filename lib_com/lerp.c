/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "prot.h"


/*-------------------------------------------------------------*
 * procedure lerp()                                            *
 *                                                             *
 *                                                             *
 *-------------------------------------------------------------*/

static void lerp_proc(
    float *f,
    float *f_out,
    int bufferNewSize,
    int bufferOldSize
);


void lerp(
    float *f,
    float *f_out,
    int bufferNewSize,
    int bufferOldSize
)
{
    float maxFac;

    maxFac = 507.0/128.0;

    if( (float)bufferNewSize / bufferOldSize > maxFac )
    {
        int tmpNewSize = bufferOldSize*2;
        while(bufferNewSize > bufferOldSize)
        {
            if( (float)bufferNewSize / bufferOldSize <= maxFac )
            {
                tmpNewSize = bufferNewSize;
            }

            lerp_proc(f, f_out, tmpNewSize, bufferOldSize);

            f = f_out;
            bufferOldSize = tmpNewSize;
            tmpNewSize *= 2;
        }
    }
    else if( (float)bufferOldSize / bufferNewSize > maxFac )
    {
        int tmpNewSize = bufferOldSize/2;
        while(bufferNewSize < bufferOldSize)
        {
            if( (float)bufferOldSize / bufferNewSize <= maxFac )
            {
                tmpNewSize = bufferNewSize;
            }

            lerp_proc(f, f_out, tmpNewSize, bufferOldSize);

            f = f_out;
            bufferOldSize = tmpNewSize;
            tmpNewSize /= 2;
        }
    }
    else
    {
        lerp_proc(f, f_out, bufferNewSize, bufferOldSize);
    }
}

void lerp_proc(
    float *f,
    float *f_out,
    int bufferNewSize,
    int bufferOldSize
)
{
    int i, idx;
    float pos, shift, diff;
    float buf[2*L_FRAME_MAX];


    if( bufferNewSize == bufferOldSize )
    {
        mvr2r( f, buf, bufferNewSize );
        mvr2r( buf, f_out, bufferNewSize );
        return;
    }

    /* Using the basop code to avoid reading beyond end of input for bufferOldSize=320, bufferNewSize=640 */
    shift = (float)(L_shl(L_deposit_l(div_s( bufferOldSize, shl(bufferNewSize, 4))), 4-15+16))/65536.0f;
    pos = 0.5f * shift - 0.5f;

    if (shift < 0.3f)
    {
        pos = pos - 0.13f;
    }

    /* first point of interpolation */
    if( pos<0 )
    {
        buf[0]=f[0]+pos*(f[1]-f[0]);
    }
    else
    {
        idx=(int)pos;
        diff = pos - idx;
        buf[0] = f[idx] + diff * (f[idx+1]-f[idx]);
    }

    pos += shift;

    for ( i=1; i<bufferNewSize-1; i++ )
    {
        idx = (int)pos;
        diff = pos-idx;

        buf[i] = f[idx]+diff*(f[idx+1]-f[idx]);
        pos += shift;
    }


    /* last point */
    idx = (int)pos;

    if( pos > bufferOldSize-1 )
    {
        idx=bufferOldSize-2;
    }

    diff = pos - idx;

    buf[bufferNewSize-1] = f[idx]+diff*(f[idx+1]-f[idx]);

    mvr2r( buf, f_out, bufferNewSize );

    return;
}
