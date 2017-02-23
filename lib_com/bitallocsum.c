/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"


/*--------------------------------------------------------------------------
 *  bitallocsum()
 *
 *  Calculate the total number of bits allocated over frame
 *--------------------------------------------------------------------------*/

void bitallocsum(
    short *R,                      /* i  : bit-allocation vector                         */
    const short nb_sfm,                  /* i  : number of sub-vectors                         */
    short *sum,                    /* o  : total number of bits allocated                */
    short *Rsubband,               /* o  : rate per subband (Q3)                         */
    const short v,                       /* i  : bit rate                                      */
    const short length,                  /* i  : length of spectrum (32 or 48 kHz samplerate)  */
    const short *sfmsize                 /* i  : band length                                   */
)
{
    short i;
    short total, tmp;
    short diff;

    total = 0;
    for (i = 0; i < nb_sfm; i++)
    {
        tmp = R[i] * sfmsize[i];
        Rsubband[i] = tmp*8;
        total += tmp;
    }
    *sum = total;

    if ( length <= L_FRAME32k )
    {
        diff = v - *sum;
        i = 0;
        while ( diff > 0 )
        {
            if ( R[i] > 0 )
            {
                Rsubband[i] += 8;
                diff -= 1;
                *sum += 1;
            }
            i++;
            if ( i >= nb_sfm )
            {
                i = 0;
            }
        }
    }

    return;
}
