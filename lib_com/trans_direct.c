/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*-------------------------------------------------------------------
 * direct_transform()
 *
 * Transformation of the signal to DCT domain
 *-------------------------------------------------------------------*/

void direct_transform(
    const float *in32,          /* i  : input signal                        */
    float *out32,         /* o  : transformation                      */
    const short is_transient,   /* i  : is transient                        */
    const short L               /* i  : length                              */
)
{
    short i;
    short seg;
    short segment_length;

    const float *wh;
    const float *wl;
    float *sh;
    float *sl;
    float *iseg;
    float *oseg;
    float dctin32[L_FRAME48k];
    float in32_r16[L_FRAME48k];

    const float *win;


    segment_length = L/2;

    if (is_transient)
    {
        if (L == L_FRAME48k)
        {
            win = wscw16q15;
        }
        else if (L == L_FRAME32k)
        {
            win = wscw16q15_32;
        }
        else if (L == L_FRAME8k)
        {
            win = wscw16q15_8;
        }
        else
        {
            win = wscw16q15_16;
        }

        for (i = 0; i < L/2; i++)
        {
            in32_r16[i]     = in32[L-1-i];
            in32_r16[L-1-i] = in32[i];
        }
        iseg = in32_r16 - segment_length/4;
        oseg = out32;

        wh = win + segment_length/4;
        wl = win + segment_length/4 - 1;
        sh = iseg + 3*segment_length/4;
        sl = iseg + 3*segment_length/4 - 1;
        for (i = 0; i < segment_length/4; i++)
        {
            dctin32[i] = ((*wl-- **sl--) - (*wh++ **sh++));
        }

        sl = iseg + segment_length/2 - 1;

        for (i = 0; i < segment_length/4; i++)
        {
            dctin32[segment_length/4 + i] = -(*sl--);
        }

        edct(dctin32, oseg, segment_length/2);

        iseg = iseg + segment_length/2;
        oseg = oseg + segment_length/2;

        for (seg = 1 ; seg <  NUM_TIME_SWITCHING_BLOCKS-1; seg++)
        {
            wh = win + segment_length/4;
            wl = win + segment_length/4 - 1;
            sh = iseg + 3*segment_length/4;
            sl = iseg + 3*segment_length/4 - 1;
            for (i = 0; i < segment_length/4; i++)
            {
                dctin32[i] = ((*wl-- **sl--) - (*wh++ **sh++));
            }

            sh = iseg;
            sl = iseg + segment_length/2 - 1;
            wh = win + segment_length/2 - 1;
            wl = win + 0;

            for (i = 0; i < segment_length/4; i++)
            {
                dctin32[segment_length/4 + i] = ((*wl++ **sl--) + (*wh-- **sh++));
            }

            edct(dctin32, oseg, segment_length/2);

            iseg = iseg + segment_length/2;
            oseg = oseg + segment_length/2;
        }

        sh = iseg + 3*segment_length/4 - 1;
        for (i = 0; i < segment_length /4; i++)
        {
            dctin32[i] = -(*sh--);
        }

        sh = iseg;
        sl = iseg + segment_length/2 - 1;
        wh = win + segment_length/2 - 1;
        wl = win + 0;

        for (i = 0; i < segment_length/4; i++)
        {
            dctin32[segment_length/4 + i] = ((*wh-- **sh++) + (*wl++ **sl--));
        }

        edct(dctin32, oseg, segment_length/2);
    }
    else
    {
        edct(in32, out32, L);
    }

    return;
}
