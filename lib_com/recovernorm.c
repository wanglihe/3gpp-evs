/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "prot.h"
#include "cnst.h"
#include "rom_com.h"


/*--------------------------------------------------------------------------*
 * recovernorm()
 *
 * Recover reordered quantization indices and norms
 *--------------------------------------------------------------------------*/

void recovernorm(
    short *idxbuf,            /* i  : reordered quantization indices */
    short *ynrm,              /* o  : recovered quantization indices */
    short *normqlg2,          /* o  : recovered quantized norms      */
    short nb_sfm
)
{
    short i;
    const short *order = NULL;

    switch(nb_sfm)
    {
    case NB_SFM:
        order = norm_order_48;
        break;
    case SFM_N_SWB:
        order = norm_order_32;
        break;
    case SFM_N_WB:
        order = norm_order_16;
        break;
    default:
        order = norm_order_48;
        break;
    }

    for (i = 0; i < nb_sfm; i++)
    {
        ynrm[order[i]] = idxbuf[i];
        normqlg2[order[i]] = dicnlg2[idxbuf[i]];
    }

    return;
}
