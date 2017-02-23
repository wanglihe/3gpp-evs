/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <assert.h>
#include "prot.h"

extern float const dico_lsf_abs_8b[];

/*------------------------------------------------------------------*
* vlpc_1st_dec()
*
*
*------------------------------------------------------------------*/

void vlpc_1st_dec(
    int index,      /* input:  codebook index                  */
    float *lsfq,    /* i/o:    i:prediction   o:quantized lsf  */
    float sr_core
)
{
    short    i;
    const float *p_dico;
    float scale = sr_core/INT_FS_12k8;

    assert(index < 256);

    p_dico = &dico_lsf_abs_8b[index * M];
    for (i = 0; i < M; i++)
    {
        lsfq[i] += scale **p_dico++;
    }

    return;

}
