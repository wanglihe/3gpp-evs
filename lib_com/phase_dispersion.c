/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*-----------------------------------------------------------------------*
 * phase_dispersion()
 *
 * Post-processing to enhance noise at low bitrate
 *-----------------------------------------------------------------------*/

void phase_dispersion(
    const float gain_code,  /* i  : gain of code             */
    const float gain_pit,   /* i  : gain of pitch            */
    float code[],     /* i/o: code vector              */
    const short mode,       /* i  : level, 0=hi, 1=lo, 2=off */
    float disp_mem[]  /* i/o: static memory (size = 8) */
)
{
    short i, j, state;
    float *prev_gain_pit, *prev_gain_code, *prev_state;
    float code2[2*L_SUBFR];
    float h_disp[L_SUBFR], *code2_real, *code2_imag, *code_real, *code_imag, *h_real, *h_imag;

    prev_state = disp_mem;
    prev_gain_code = disp_mem+1;
    prev_gain_pit = disp_mem+2;

    state = 2;
    if (gain_pit < 0.6f)
    {
        state = 0;
    }
    else if (gain_pit < 0.9f)
    {
        state = 1;
    }

    for (i=5; i>0; i--)
    {
        prev_gain_pit[i] = prev_gain_pit[i-1];
    }
    prev_gain_pit[0] = gain_pit;

    if (gain_code - 3.0f * *prev_gain_code > 0.0f)
    {
        if (state < 2)
        {
            state++;
        }
    }
    else
    {
        j=0;
        for (i=0; i<6; i++)
        {
            if (prev_gain_pit[i] < 0.6f)
            {
                j++;
            }
        }

        if (j > 2)
        {
            state = 0;
        }

        if ((state - (short)*prev_state) > 1)
        {
            state--;
        }
    }

    *prev_gain_code = gain_code;
    *prev_state = (float)state;

    /*-----------------------------------------------------------------*
     * Circular convolution
     *-----------------------------------------------------------------*/

    state += mode;                        /* level of dispersion */
    if( state < 2 )
    {
        fft_rel( code, L_SUBFR, 6 );

        if (state == 0)
        {
            mvr2r( low_H, h_disp, L_SUBFR );
        }

        if (state == 1)
        {
            mvr2r( mid_H, h_disp, L_SUBFR );
        }

        code2_real = code2;
        code2_imag = code2 + L_SUBFR - 1;
        code_real = code;
        code_imag = code + L_SUBFR - 1;
        h_real = h_disp;
        h_imag = h_disp + L_SUBFR - 1;
        *code2_real++ = *code_real++ **h_real++;

        for (i=1; i<L_SUBFR/2; i++)
        {
            *code2_real++ = *code_real **h_real - *code_imag **h_imag;
            *code2_imag-- = *code_real++ **h_imag-- + *code_imag-- **h_real++;
        }

        *code2_real++ = *code_real++ **h_real++;
        ifft_rel( code2, L_SUBFR, 6 );

        mvr2r( code2, code, L_SUBFR );
    }

    return;
}
