/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"

/*---------------------------------------------------------------------*
 * Local functions
 *---------------------------------------------------------------------*/

/*-------------------------------------------------------------------*
 * syn_output()
 *
 * Output synthesis signal with compensation for saturation
 * returns number of clipped samples
 *-------------------------------------------------------------------*/

unsigned int syn_output(      /* o  : number of clipped samples            */
    float *synth,       /* i/o: float synthesis signal               */
    const short output_frame, /* i  : output frame length                  */
    short *synth_out    /* o  : integer 16 bits synthesis signal     */
)
{
    /*-----------------------------------------------------------------*
     * float to integer conversion with saturation control
     *-----------------------------------------------------------------*/

    /* integer conversion */
    return mvr2s( synth, synth_out, output_frame );
}


/*-------------------------------------------------------------------*
 * AGC_dec()
 *
 * In-place saturation control (Automatic Gain Control)
 *-------------------------------------------------------------------*/

void AGC_dec(
    float x[],    /* i/o: input/output vector                   */
    float mem[],  /* i/o: mem[2] should be init to [0,0] */
    const short n       /* i  : vector size                    */
)
{
    short i;
    float fac, prev, tmp, frame_fac, max;

    /*-----------------------------------------------------------------*
     * calculate AGC factor to avoid saturation
     *-----------------------------------------------------------------*/

    max = 0.0f;

    for (i=0; i<n; i++)
    {
        tmp = (float)fabs(x[i]);
        if (tmp > max)
        {
            max = tmp;
        }
    }

    frame_fac = 0.0f;
    if ( max > 30000.0f )
    {
        frame_fac = 0.5f - (15000.0f/max);
    }

    fac = mem[0];
    prev = mem[1];

    /*-----------------------------------------------------------------*
     * AGC
     *-----------------------------------------------------------------*/

    for (i=0; i<n; i++)
    {
        /* update AGC factor (slowly) */
        fac = 0.99f * fac + 0.01f * frame_fac;

        /* convert float to integer with AGC */
        tmp = (1.0f - fac) *x[i] - fac * prev;
        prev = x[i];

        if (tmp > 32767.0f)
        {
            tmp = 32767.0f;
        }
        else if (tmp < -32768.0f)
        {
            tmp = -32768.0f;
        }

        x[i] = (short)floor(tmp + 0.5f);

    }

    mem[0] = fac;
    mem[1] = prev;

    return ;
}
