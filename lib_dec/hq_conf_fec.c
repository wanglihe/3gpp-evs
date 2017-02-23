/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*--------------------------------------------------------------------------*
 * hq_configure_bfi()
 *
 * Initialization of HQ bands and subframes
 *--------------------------------------------------------------------------*/

void hq_configure_bfi(
    short *nb_sfm,             /* o  : Number of sub bands               */
    short *num_Sb,             /* o  : Number of FEC sub bands ?         */
    short *num_bands_p,        /* o  : FEC sub bands                     */
    short const **sfmsize,     /* o  : Subband bandwidths                */
    short const **sfm_start,   /* o  : Subband start coefficients        */
    short const **sfm_end      /* o  : Subband end coefficients          */
)
{
    *num_Sb = MAX_SB_NB;
    *nb_sfm = SFM_N_NB;
    mvs2s( Num_bands_NB, num_bands_p, *num_Sb );
    *sfmsize = band_len_wb;
    *sfm_start = band_start_wb;
    *sfm_end = band_end_wb;

    return;
}
