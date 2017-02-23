/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prot.h"

/*------------------------------------------------------------------*
* vlpc_2st_dec()
*
*
*------------------------------------------------------------------*/

void vlpc_2st_dec(
    float *lsfq,    /* i/o:    i:1st stage   o:1st+2nd stage   */
    int *indx,      /* input:  index[] (4 bits per words)      */
    int mode,       /* input:  0=abs, >0=rel                   */
    float sr_core
)
{
    short i;
    float w[M];
    int xq[M];
    float scale = sr_core/INT_FS_12k8;

    /* weighting from the 1st stage */
    lsf_weight_2st( lsfq, w, mode, sr_core );

    /* quantize */
    AVQ_dec_lpc( indx, xq, 2 );

    /* quantized lsf */
    for( i=0; i<M; i++) lsfq[i] += scale*(w[i]*(float)xq[i] );

    /* reorder */
    v_sort( lsfq, 0, M-1 );
    reorder_lsf( lsfq, LSF_GAP, M, sr_core );

    return;
}
