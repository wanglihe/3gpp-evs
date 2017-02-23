/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*--------------------------------------------------------------------------*
 * noise_adjust()
 *
 * Calculate attenuation
 *--------------------------------------------------------------------------*/

short noise_adjust(                 /* o  : index of noise attenuation */
    const float *coeffs_norm,       /* i  : normalized coefficients    */
    const short *bitalloc,          /* i  : bit allocation             */
    const short *sfm_start,
    const short *sfm_end,
    const short core_sfm            /* i  : index of the end band for core  */
)
{
    short nf_idx, sfm, bin, num_coeffs;
    int   E, diff, tmp;

    E = 0;
    num_coeffs = 0;

    for (sfm = 0; sfm <= core_sfm; sfm++)
    {
        if (bitalloc[sfm] == 0)
        {
            for (bin = sfm_start[sfm]; bin < sfm_end[sfm]; bin++)
            {
                if (coeffs_norm[bin] == 0)
                {
                    E = E - 1;
                }
                else
                {
                    tmp = (int)(floor(log10(fabs(coeffs_norm[bin] * L_FRAME))/0.301030f)+1);
                    if (tmp < 0)
                    {
                        tmp = 0;
                    }
                    E = E + tmp;
                }

                num_coeffs = num_coeffs + 1;
            }
        }
    }

    if (num_coeffs != 0)
    {
        E = E/num_coeffs;
    }
    else
    {
        E = 0;
    }

    diff = 7 - E;

    if ( diff >= 0 )
    {
        nf_idx = (short)diff;
        if ( diff > 3 )
        {
            nf_idx = 3;
        }
    }
    else
    {
        nf_idx = 0;
    }

    return nf_idx;
}
