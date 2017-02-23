/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "prot.h"
#include "cnst.h"
#include "rom_com.h"

/*--------------------------------------------------------------------------*/
/*  Function  reordernorm                                                   */
/*  ~~~~~~~~~~~~~~~~~~~~~                                                   */
/*                                                                          */
/*  Reorder quantization indices and norms                                  */
/*--------------------------------------------------------------------------*/
/*  short     *ynrm      (i)   quantization indices for norms               */
/*  short     *normqlg2  (i)   quantized norms                              */
/*  short     *idxbuf    (o)   reordered quantization indices               */
/*  short     *normbuf   (o)   reordered quantized norms                    */
/*  short     *nb_sfm    (i)   number of bands                              */
/*--------------------------------------------------------------------------*/

void reordernorm(
    const short *ynrm,
    const short *normqlg2,
    short *idxbuf,
    short *normbuf,
    const short nb_sfm
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
        idxbuf[i] = ynrm[order[i]];
        normbuf[i] = normqlg2[order[i]];
    }

    return;
}
